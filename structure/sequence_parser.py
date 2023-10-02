from django.conf import settings

from residue.models import Residue
from protein.models import Protein, ProteinSegment
from structure.models import Structure
from structure.functions import BlastSearch, BlastSearchOnline

from Bio import SeqIO, pairwise2
from Bio.PDB import PDBParser, PPBuilder
import Bio.PDB.Polypeptide as polypeptide

from collections import OrderedDict
import os, xlsxwriter

from datetime import datetime, date
startTime = datetime.now()

#Number of heavy atoms in each residue
atom_count = {
    "ALA": 5,
    "ARG": 11, 
    "ASN": 8,
    "ASP": 8, 
    "CYS": 6,
    "GLN": 9,
    "GLU": 9, 
    "GLY": 4,
    "HIS": 10,
    "ILE": 8,
    "LEU": 8,
    "LYS": 9,
    "MET": 8,
    "PHE": 11,
    "PRO": 7,
    "SER": 6,
    "THR": 7,
    "TRP": 14,
    "TYR": 12,
    "VAL": 7,
    }

class ParsedResidue(object):

    def __init__(self, res_name, res_num, gpcrdb=None, segment=None, coords='full'):
        
        self.resnum = None
        self.wt_num = res_num
        self.name = res_name
        self.mutation = None
        self.insertion = None
        self.deletion = False
        self.coords = coords
        self.gpcrdb = gpcrdb if gpcrdb else ''
        self.segment = segment
        self.seqres = False
        self.fusion = None


    def __repr__(self, **kwargs):
        return "Residue {}{} PDB {} ({})\tMutated: {} Insertion: {} Fusion: {} SEQRES: {}".format(self.name, self.wt_num, self.resnum if self.resnum else '-', self.gpcrdb, self.mutation if self.mutation else '-', self.insertion if self.insertion else '', self.fusion, self.seqres)


    def get_param_list(self):
        return [self.wt_num, self.name, self.resnum if self.resnum else '-', self.gpcrdb, self.mutation if self.mutation else '-', 'X' if self.seqres else '']


    def set_mutation(self, mutation):
        self.mutation = mutation


    def set_insertion(self, insertion):
        self.insertion = insertion

    def set_deletion(self, deletion=True):
        self.deletion = deletion        


    def set_fusion(self, fusion=True):
        self.fusion = fusion


    def set_seqres(self, seqres=True):
        self.seqres = seqres


    def set_coords_status(self, coords):
        self.coords = coords


    def set_gpcrdb(self, gpcrdb):
        self.gpcrdb = gpcrdb


    def set_wt_number(self, wt_num):
        self.wt_num = wt_num


    def set_pdb_res_num(self, res_num):
        self.resnum = res_num

class AuxProtein(object):
    """
    Class storing the mapping of the fusion/auxiliary protein.
    """

    residue_list = ["ARG","ASP","GLU","HIS","ASN","GLN","LYS","SER","THR", "HIS", "HID","PHE","LEU","ILE","TYR","TRP","VAL","MET","PRO","CYS","ALA","GLY"]

    def __init__(self, residues):
        
        self.residues = residues
        self.seq = self.get_peptide_sequence(self.residues)
        self.blast_online = BlastSearchOnline()

        self.id = ''
        self.mapping = {}
        self.map_aux()


    def get_peptide_sequence(self, residues):
        """
        Returns a sequence string of a given list of Bio.PDB.Residue objects.
        """
        return "".join([polypeptide.three_to_one(x.resname.replace('HID', 'HIS')) for x in residues if x.resname in self.residue_list])

    def map_aux(self):
        print('aux1',datetime.now() - startTime)
        alignments = self.blast_online.run(self.seq)
        print('aux2',datetime.now() - startTime)

        for alignment in alignments:
            self.id = alignment[0]
            for hsps in alignment[1].hsps:
                self.map_hsps(hsps)
    

    def map_hsps(self, hsps):
        """
        Analyzes the High Similarity Protein Segment.
        """
        offset = min([int(x.id[1]) for x in self.residues])
        q = hsps.query
        sbjct = hsps.sbjct
        sbjct_counter = hsps.sbjct_start	
        q_counter = hsps.query_start
        for s, q in zip(sbjct, q):
            if s == q:
                self.mapping[sbjct_counter] = offset - 1 + q_counter
                sbjct_counter += 1
                q_counter += 1
            elif s != '-' and q != '-':
                self.mapping[sbjct_counter] = offset - 1 + q_counter
                sbjct_counter += 1
                q_counter += 1
            elif s != '-' and q == '-':
                sbjct_counter += 1
            else:
                sbjct_counter += 1
                q_counter += 1

    def get_info(self):

        return  OrderedDict({
                "presence" : 'YES',
                "type" : "fusion",
                "uniprot" : self.id,
                "description" : "",
                "start" : min(self.mapping.values()),
                "end" : max(self.mapping.values()),
                "start (pdb)" : min(self.mapping.keys()),
                "end (pdb)" : max(self.mapping.keys()),
                })


class SequenceParser(object):
    """
    Class mapping the pdb, pdb_seqres, wildtype and any given sequence onto wt using blast with human sequences database. It produces a report with missing, mutated and inserted residues.
    """

    residue_list = ["ARG","ASP","GLU","HIS","ASN","GLN","LYS","SER","THR", "HIS", "HID","PHE","LEU","ILE","TYR","TRP","VAL","MET","PRO","CYS","ALA","GLY"]

    def __init__(self, pdb_file=None, sequence=None, wt_protein_id=None, db='protwis_blastdb'):

        # dictionary of 'ParsedResidue' object storing information about alignments and bw numbers
        self.mapping = {}
        self.residues = {}
        self.segments = {}
        self.blast = BlastSearch(blastdb=os.sep.join([settings.STATICFILES_DIRS[0], 'blast', db]))
        self.wt_protein_id = wt_protein_id

        if pdb_file is not None:
            self.pdb_struct = PDBParser(QUIET=True).get_structure('pdb', pdb_file)[0]
            # a list of SeqRecord objects retrived from the pdb SEQRES section
            try:
                self.seqres = list(SeqIO.parse(pdb_file, 'pdb-seqres'))
                self.struct_id = self.seqres[0].id.split(':')[0]
            except:
                self.seqres = None
                self.struct_id = None
            # SeqRecord id is a pdb_code:chain

        self.sequence = sequence
        if type(sequence) == "string":
            self.sequence = { x: y for x,y in enumerate(sequnece) }


        # If not specified, attempt to get wildtype from pdb.
        try:
            if not wt_protein_id and pdb_file is not None:
                self.wt = Structure.objects.get(pdb_code__index=self.struct_id).protein_conformation.protein.parent
            else:
                raise Exception()
        except:
            if not wt_protein_id:
                self.wt = None
                self.wt_seq = ''
            else:
                self.wt = Protein.objects.get(id=wt_protein_id)
                self.wt_seq = str(self.wt.sequence)
        self.fusions = []

        self.parse_pdb(self.pdb_struct)
        #if self.seqres:
        #    self.map_seqres()
        
        self.mark_deletions()


    def parse_pdb(self, pdb_struct):
        """
        extracting sequence and preparing dictionary of residues
        bio.pdb reads pdb in the following cascade: model->chain->residue->atom
        """

        for chain in pdb_struct:
            self.residues[chain.id] = []
            
            for res in chain:
            #in bio.pdb the residue's id is a tuple of (hetatm flag, residue number, insertion code)
                if res.resname.replace('HID', 'HIS') not in self.residue_list:
                    continue
                self.residues[chain.id].append(res)
            poly = self.get_chain_peptides(chain.id)
            for peptide in poly:
                #print("Start: {} Stop: {} Len: {}".format(peptide[0].id[1], peptide[-1].id[1], len(peptide)))
                self.map_to_wt_blast(chain.id, peptide, None, int(peptide[0].id[1]))


    def get_segments(self):

        #get the first chain
        c = list(self.mapping.keys())[0]

        for segment in ProteinSegment.objects.all():
            resi = []
            for r in Residue.objects.filter(protein_conformation__protein=self.wt.id, protein_segment=segment):
                if self.mapping[c][r.sequence_number].resnum is not None:
                    resi.append(self.mapping[c][r.sequence_number].resnum)
            if resi == []:
                continue
            self.segments[segment.slug] = [min(resi), max(resi)]
        return self.segments


    def get_chain_peptides(self, chain_id, gap_threshold=230):
        """
        Get peptides of sequential residue numbers (with respect to 230 aa gaps).
        The maximum length of ICL3 is 230 aa, and fusion proteins usualy have significantly different numbers, i.e. exceeding the 230 gap between TM5 and 6.

        The maximum allowed gap size can be evaluated automaticaly, but it is fairly costly:
        max([len(Residue.objects.filter(protein_segment=11, protein_conformation__protein=x)) for x in Protein.objects.filter(species=1)])
        """

        rnumbers = [int(x.id[1]) for x in self.residues[chain_id]]
        last_idx = len(rnumbers)-1
        peptides = []
        tmp = []
        for i, rnum in enumerate(rnumbers):
            if i == last_idx:
                #FIXME: Assuming that very last residue is actualy continuation of a chain
                tmp.append(self.residues[chain_id][i])
                peptides.append(tmp)
                break
            if rnumbers[i+1] != rnum+1 and abs(rnum+1 - rnumbers[i+1]) > gap_threshold:
                tmp.append(self.residues[chain_id][i])
                peptides.append(tmp)
                tmp = []
            else:
                tmp.append(self.residues[chain_id][i])
        return peptides


    def get_chain_sequence(self, chain):
        """
        Returns a sequence string of a given chain.
        """
        return "".join([polypeptide.three_to_one(x.resname.replace('HID', 'HIS')) for x in self.residues[chain] if x.resname in self.residue_list])

    def get_peptide_sequence(self, residues):
        """
        Returns a sequence string of a given list of Bio.PDB.Residue objects.
        """
        return "".join([polypeptide.three_to_one(x.resname.replace('HID', 'HIS')) for x in residues if x.resname in self.residue_list])
    
    def find_nonredundant_chains(self):
        """
        Returns a list of nonidentical chains.
        """
        nrc = []
        if len(self.mapping.keys()) == 1:
            return self.mapping.keys()

        for r_chain in self.mapping.keys():
            for chain in self.mapping.keys():
                if r_chain == chain:
                    continue
                if self.mapping[r_chain] != self.mapping[chain]:
                    nrc.append(r_chain)
        return nrc


    def map_to_wt_blast(self, chain_id, residues = None, sequence=None, starting_aa = 1, seqres = False):

        if residues:
            seq = self.get_peptide_sequence(residues)
        elif sequence:
            seq = sequence
        else:
            seq = self.get_chain_sequence(chain_id)
        alignments = self.blast.run(seq)
        
        if self.wt_protein_id!=None:
            self.wt = Protein.objects.get(id=self.wt_protein_id)
        else:
            self.wt = None
        for alignment in alignments:
            if self.wt==None:
                try:
                    self.wt = Protein.objects.get(entry_name=str(alignment[1].hit_def))
                    wt_resi = list(Residue.objects.filter(protein_conformation__protein=self.wt.id))
                    self.mapping[chain_id] = {x.sequence_number: ParsedResidue(x.amino_acid, x.sequence_number, str(x.display_generic_number) if x.display_generic_number else None, x.protein_segment) for x in wt_resi}
                except:
                    pass
            else:
                wt_resi = list(Residue.objects.filter(protein_conformation__protein=self.wt.id))
                self.mapping[chain_id] = {x.sequence_number: ParsedResidue(x.amino_acid, x.sequence_number, str(x.display_generic_number) if x.display_generic_number else None, x.protein_segment) for x in wt_resi}
            if alignment[1].hsps[0].expect > .5 and residues:
                # self.fusions.append(AuxProtein(residues))
                #The case when auxiliary protein is in a separate chain
                if self.get_chain_sequence(chain_id) == self.get_peptide_sequence(residues) and chain_id in self.mapping:
                    del self.mapping[chain_id]
                continue
            if self.wt.id != int(alignment[0]):
                continue
            for hsps in alignment[1].hsps:
                self.map_hsps(hsps, chain_id, starting_aa, seqres)
                # break
    

    def map_hsps(self, hsps, chain_id, offset = 1, seqres = False):
        """
        Analyzes the High Similarity Protein Segment.
        """
        q = hsps.query
        sbjct = hsps.sbjct
        sbjct_counter = hsps.sbjct_start	
        q_counter = hsps.query_start

        for s, q in zip(sbjct, q):
            if s == q:
                if seqres:
                    self.mapping[chain_id][sbjct_counter].set_seqres(True)
                else:
                    self.mapping[chain_id][sbjct_counter].set_pdb_res_num(offset - 1 + q_counter)
                sbjct_counter += 1
                q_counter += 1
            elif s != '-' and q != '-':
                self.mapping[chain_id][sbjct_counter].set_pdb_res_num(offset - 1 + q_counter)
                self.mapping[chain_id][sbjct_counter].set_mutation(q)
                sbjct_counter += 1
                q_counter += 1
            elif s == '-' and q != '-':
                self.mapping[chain_id][offset - 1 + q_counter].set_insertion(q)
                sbjct_counter += 1
                q_counter += 1
            elif s != '-' and q == '-':
                self.mapping[chain_id][sbjct_counter].set_deletion()
                sbjct_counter += 1
                q_counter += 1

    def map_to_wt_pw(self, chain_id, residues = None, sequence=None, starting_aa = 1):

        """
        @param sequence: a dictionary of residue number: residue one letter code pairs
        """

        if residues:
            seq = self.get_chain_sequence(residues)
        elif sequence:
            seq = sequence.values()
        else:
            return

        wt, chain_seq, score, start, end = pairwise2.align.localms(self.wt_seq, seq, 2, -4, -4, -.1, one_alignment_only=True)[0]

        offset = 0
        for w, c in zip(wt, chain_seq):
            if w == c:
                if seqres:
                    self.mapping[chain.id][starting_aa + offset].seqres=True
                r = Residue.objects.get(sequence_number=offset+self.wt_seq_start, protein_conformation__protein=self.wt.id)
                if r.display_generic_number is not None:
                    self.mapping[chain_id][starting_aa + offset].add_gpcrdb(r.display_generic_number)
                offset += 1
            elif c == '-' and w != '-':
                self.mapping[chain_id][starting_aa + offset].add_deletion()
            elif w != '-' and c != '-' and w != c:
                self.mapping[chain_id][starting_aa + offset].add_mutation(c)
                offset += 1
            elif w == '-' and c != '-':
                self.mapping[chain_id][starting_aa + offset].add_insertion(c)
                offset += 1


    def map_seqres(self):

        for sr in self.seqres:
            self.map_to_wt_blast(sr.annotations['chain'], sequence=sr.seq, seqres=True)

    def mark_deletions(self):
        for chain in self.mapping.keys():
            for num, res in self.mapping[chain].items():
                if res.resnum is None:
                    res.set_deletion()

    def get_mapping_dict(self, pdb_keys=False, seqres=False):

        if pdb_keys:
            return {x: {y: self.mapping[x][y].seqres if seqres else self.mapping[x][y].resnum for y in self.mapping[x].keys()} for x in self.mapping.keys()}
        else:
            if seqres:
                return {x: {y: self.mapping[x][y].resnum if self.mapping[x][y].seqres else '-' for y in self.mapping[x].keys()} for x in self.mapping.keys()}
            else:
                return {x: {y: self.mapping[x][y].resnum for y in self.mapping[x].keys()} for x in self.mapping.keys()}

    def get_fusions(self):

        if self.fusions == []:
            return {}
        fusion_dict = OrderedDict({"auxiliary": {}})
        count = 1
        for fusion in self.fusions:
            fusion_dict["auxiliary"]["aux{}".format(count)] = fusion.get_info()
        return fusion_dict

    def get_deletions(self):

        deletions_list = []

        for chain in self.find_nonredundant_chains():
            deletions = [x for x,y in self.mapping[chain].items() if y.deletion]
            deletion = deletions.reverse()
            tmp = []
            #for num, res in self.mapping[chain].items():
            #    if res.deletion:
            #        tmp.append(num)
            first = 0
            prev = 0
            while deletions != []:
                x = deletions.pop()
                #print("{}\t{}\t{}".format(x, first, prev))
                if first == 0:
                    tmp.append(x)
                    first = x
                    continue
                if prev == 0:
                    tmp.append(x)
                    prev = x
                    continue
                if abs(x - prev) == 1:
                    tmp.append(x)
                    prev = x
                else:
                    deletions_list.append(OrderedDict({
                        "start" : min(tmp),
                        "end" : max(tmp),
                        "type" : "single" if len(tmp) == 1 else "range",
                        "chain" : chain
                        }))
                    tmp = [x]
                    first = x
                    prev = x
            deletions_list.append(OrderedDict({
                        "start" : min(tmp),
                        "end" : max(tmp),
                        "type" : "single" if len(tmp) == 1 else "range",
                        "chain" : chain
                        }))

        return {"deletions" : deletions_list}

    def get_mutations(self):

        mutations_list = []
        for chain in self.find_nonredundant_chains():
            for num, res in self.mapping[chain].items():
                if res.mutation:
                    mutations_list.append(OrderedDict({
                        "wt" : res.name,
                        "mut" : res.mutation,
                        "pos (wt)" : num,
                        "pos (pdb)" : res.resnum,
                        "chain" : chain
                        }))
        return {"mutations" : mutations_list }


    def get_report(self):

        for chain in sorted(self.mapping.keys()):
            print("Chain {}".format(chain))
            for res in sorted(self.mapping[chain].keys()):
                print(self.mapping[chain][res])

    def save_excel_report(self, file_name):
        
        workbook = xlsxwriter.Workbook(file_name)
        
        for chain in sorted(self.mapping.keys()):
            worksheet = workbook.add_worksheet(chain)
            worksheet.write_row(0,0,["Protein number", "Residue name", "PDB number", "Generic number", "Mutation", "SEQRES"])

            row_id = 1
            for res in sorted(self.mapping[chain].keys()):
                tmp = self.mapping[chain][res]
                worksheet.write_row(row_id, 0, tmp.get_param_list())
                row_id += 1
        workbook.close()


class SequenceParserPW(object):
    """
    Class mapping the pdb, pdb_seqres, wildtype and any given sequence onto wt. It produces a report with missing, mutated and inserted residues.
    """

    residue_list = ["ARG","ASP","GLU","HIS","ASN","GLN","LYS","SER","THR","HID","PHE","LEU","ILE","TYR","TRP","VAL","MET","PRO","CYS","ALA","GLY"]

    def __init__(self, pdb_file, sequence=None):

        # dictionary of 'ParsedResidue' object storing information about alignments and bw numbers
        self.mapping = {}
        
        # a list of SeqRecord objects retrived from the pdb SEQRES section
        self.seqres = list(SeqIO.parse(pdb_file, 'pdb-seqres'))

        self.pdb_struct = PDBParser(QUIET=True).get_structure('pdb', pdb_file)[0]

        # SeqRecord id is a pdb_code:chain 
        self.struct_id = self.seqres[0].id.split(':')[0]
        self.wt = Structure.objects.get(pdb_code__index=self.struct_id).protein_conformation.protein.parent
        self.wt_seq = str(self.wt.sequence)
        self.wt_seq_start = Residue.objects.filter(protein_conformation__protein=self.wt.id).order_by("sequence_number")[0].sequence_number

        # a dictionary of per chain lists of peptides found in the pdb
        self.pdb_seq = {}
        for chain in self.pdb_struct:
            self.pdb_seq[chain.id] = self.get_res_list(chain)
            self.mapping[chain.id] = {}

            self.map_wildtype(chain=chain)

    def get_res_list(self, chain):

        #Both Polypeptide and and SeqIO suck at retrieving full aminoacid sequence. Have to do it the hard way.
        return [x  for x in chain if x.resname in self.residue_list]


    def get_chain_sequence(self, chain):
        return "".join([polypeptide.three_to_one(x.resname.replace('HID', 'HIS')) for x in chain if x.resname in self.residue_list])


    def align_to_wt(self, sequence):
        """
        Get the pairwise alignment between wildtype and a given sequence.
        """
        return pairwise2.align.localms(self.wt_seq, sequence, 2, -4, -4, -.1, one_alignment_only=True)[0]


    def map_wildtype(self, chain=None, seqres=None, sequence=None):

        if chain:
            query = self.get_chain_sequence(chain)
        elif seqres:
            query = seqres
        else:
            query = sequence

        wt, chain_seq, score, start, end = self.align_to_wt(query)
        offset = 0
        for w, c in zip(wt, chain_seq):
            if w == c:
                if seqres:
                    self.mapping[chain.id][offset+self.wt_seq_start].seqres=True
                r = Residue.objects.get(sequence_number=offset+self.wt_seq_start, protein_conformation__protein=self.wt.id)
                if r.display_generic_number is not None:
                    self.mapping[chain.id][offset+self.wt_seq_start].add_gpcrdb(r.display_generic_number)
                offset += 1
            elif c == '-' and w != '-':
                self.mapping[chain.id][offset+self.wt_seq_start].add_deletion()
                offset += 1
            elif w != '-' and c != '-' and w != c:
                self.mapping[chain.id][offset+self.wt_seq_start].add_mutation(c)
                offset += 1

