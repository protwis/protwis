from residue.models import Residue
from structure.models import Structure
from structure.functions import MappedResidue

from Bio import SeqIO, pairwise2
from Bio.PDB import PDBParser, PPBuilder

import enum

class ParsedResidue(object):

    def __init__(self, res_name, res_num):
        
        self.resnum = res_num
        self.name = res_name
        self.mutation = None
        self.deletion = False
        self.coords = enum.Enum(['full', 'partial', 'missing'])
        self.gpcrdb = ''
        self.seqres = False


    def add_mutation(self, mutation):
        self.mutation = mutation


    def add_deletion(self):
        self.deletion = True

    def coords(self, coords):
        self.coords = coords

    def add_gpcrdb(self, gpcrdb):
        self.gpcrdb = gpcrdb


class SequenceParser(object):
    """
    Class mapping the pdb, pdb_seqres, wildtype and any given sequence onto wt. It produces a report with missing, mutated and inserted residues.
    """

    residue_list = ["ARG","ASP","GLU","HIS","ASN","GLN","LYS","SER","THR","HID","PHE","LEU","ILE","TYR","TRP","VAL","MET","PRO","CYS","ALA","GLY"]

    def __init__(self, pdb_file, sequence=None):

        # dictionary of 'ParsedResidue' object storing information about alignments and bw numbers
        self.mapping = {}
        
        # a list of SeqRecord objects retrived from the pdb SEQRES section
        self.seqres = SeqIO.parse(pdb_file, 'pdb-seqres')

        self.pdb_struct = PDBParser(QUIET=True).get_structure('pdb', pdb_file)[0]

        # SeqRecord id is a pdb_code:chain 
        self.struct_id = self.seqres[0].id.split(':')[0]
        self.wt = Structure.get(pdb_code__index=struct_id).protein_conformation.protein.parent

        self.wt_seq = str(self.wt.sequence)
        self.wt_seq_start = Residue.objects.filter(protein_conformation__protein=self.wt.id).order_by("sequence_number")[0]

        # a dictionary of per chain lists of peptides found in the pdb
        self.pdb_seq = {}
        for chain in self.pdb_struct:
            self.pdb_seq[chain.id] = self.get_res_list(chain)
            self.mapping[chain.id] = {}


    def get_res_list(self, chain):

        #Both Polypeptide and and SeqIO suck at retrieving full aminoacid sequence. Have to do it the hard way.
        return [x  for x in chain if x.resname in self.residue_list]


    def get_chain_sequence(self, chain):
        return "".join([polypeptide.three_to_one(x.resname.replace('HID', 'HIS')) for x in chain if x.resname in self.residue_list])


    def align_to_wt(self, sequence):
        """
        Get the pairwise alignment between wildtype and a give sequence.
        """
        return pairwise2.align.localms(self.wt_seq, sequence, 2, -4, -4, -.1, one_alignment_only=True)[0]


    def map_wildtype(self, chain):
        
        wt, chain_seq, score, start, end = self.align_to_wt(self.get_chain_sequence(chain))

        for idx, w, c in enumerate(zip(wt.seq, chain_seq.seq)):
            if idx+self.wt_seq_start not in self.mapping[chain.id].keys():
                self.mapping[chain.id][idx+self.wt_seq_start] = ParsedResidue(w, idx+self.wt_seq_start)
            if w == c:
                r = Residue.objects.get(sequence_number=idx+self.wt_seq_start, protein_conformation__protein=self.wt.id)
                if r.display_generic_number is not None:
                    self.mapping[chain.id][idx+self.wt_seq_start].add_gpcrdb(r.display_generic_number)
            elif c == '-':
                self.mapping[chain.id][idx+self.wt_seq_start].add_deletion()
            elif w != '-' and c != '-':
                self.mapping[chain.id][idx+self.wt_seq_start].add_mutation(c)

