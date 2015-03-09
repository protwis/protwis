from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB import *
from Bio.PDB.PDBIO import Select

from common import BlastSearch, MappedResidue

import Bio.PDB.Polypeptide as polypeptide
import os


#==============================================================================

class generic_numbering(object):
    
    residue_list = ["ARG","ASP","GLU","HIS","ASN","GLN","LYS","SER","THR","HID","PHE","LEU","ILE","TYR","TRP","VAL","MET","PRO","CYS","ALA","GLY"]
  
    def __init__ (self, filename, local_gpcrdb = False):
    
        self.pdb_file = filename
        
        #dictionary of 'MappedResidue' object storing information about alignments and bw numbers
        self.residues = {}
        self.pdb_seq = {} #Seq('')
        #list of uniprot ids returned from blast
        self.up_id_list = []
        #setup for local blast search
        self.blast = BlastSearch()
        
        self.parse_pdb()
        

    def parse_pdb (self):
        #extracting sequence and preparing dictionary of residues
        #bio.pdb reads pdb in the following cascade: model->chain->residue->atom
        pdb_struct = PDBParser().get_structure('ref', self.pdb_file)[0]
    
        for chain in pdb_struct:
            self.residues[chain.id] = {}
            self.pdb_seq[chain.id] = Seq('')
            
            for res in chain:
            #in bio.pdb the residue's id is a tuple of (hetatm flag, residue number, insertion code)
                if res.resname == "HID":
                    resname = polypeptide.three_to_one('HIS')
                else:
                    if res.resname not in self.residue_list:
                        continue
                    self.residues[chain.id][res.id[1]] = MappedResidue(res.id[1], polypeptide.three_to_one(res.resname))
    
            self.pdb_seq[chain.id].seq = ''.join([self.residues[chain.id][x].name for x in sorted(self.residues[chain.id].keys())])
            
            for pos, res in enumerate(sorted(self.residues[chain.id].keys()), start=1):
                self.residues[chain.id][res].pos_in_aln = pos


    def locate_res_by_pos (self, chain, pos):
        for res in self.residues[chain].keys():
            if self.residues[chain][res].pos_in_aln == pos:
                return res
        return 0


    def map_blast_seq (self, upid, hsps, chain):
    
        #find uniprot residue numbers corresponding to those in pdb file
        q_seq = list(hsps.query)
        tmp_seq = list(hsps.sbjct)
        subj_counter = hsps.sbjct_start	
        q_counter = hsps.query_start
        
        #print "%s\n%s" %(hsps.query, hsps.sbjct)
        #print "%i\t%i" %(hsps.query_start, hsps.sbjct_start)

        while tmp_seq:
            #skipping position if there is a gap in either of sequences
            if q_seq[0] == '-' or q_seq[0] == 'X' or q_seq[0] == ' ':
                #print "Query puste: %s %s \nDone" %(q_seq[0], tmp_seq[0])
                subj_counter += 1
                #q_counter += 1
                tmp_seq.pop(0)
                q_seq.pop(0)
                continue
            if tmp_seq[0] == '-' or tmp_seq[0] == 'X' or tmp_seq[0] == ' ':
                #print "Tmp puste: %s %s \nDone" %(q_seq[0], tmp_seq[0])
                q_counter += 1
                tmp_seq.pop(0)
                q_seq.pop(0)
                continue
            if tmp_seq[0] == q_seq[0]:
                #print "Query and temp match %i:%s %i %s" %(q_counter,q_seq[0],subj_counter,tmp_seq[0])
                resn = self.locate_res_by_pos(chain, q_counter)
                #print "%i\t%i" %(resn, subj_counter)
                if resn != 0:
                    self.residues[chain][resn].add_mapping(upid, subj_counter)
                    
                    if upid not in self.up_id_list:
                        self.up_id_list.append(upid)
            q_counter += 1
            subj_counter += 1
            tmp_seq.pop(0)
            q_seq.pop(0)        
            
    
    def get_annotated_structure(self, gpcrdb=True, bw=False):
    
        pdb_struct = PDBParser().get_structure(os.path.splitext(os.path.basename(self.pdb_file))[0], self.pdb_file)
        
        for chain in pdb_struct[0]:
            for residue in chain:
                if self.residues[chain.id].has_key(residue.id[1]):
                    #print residue.id[1]
                    if self.residues[chain.id][residue.id[1]].gpcrdb != 0. and gpcrdb:
                        residue["CA"].set_bfactor(float(self.residues[chain.id][residue.id[1]].gpcrdb))
                        #print self.residues[chain.id][residue.id[1]].gpcrdb
                    if self.residues[chain.id][residue.id[1]].bw != 0. and bw:
                        residue["N"].set_bfactor(float(self.residues[chain.id][residue.id[1]].bw))
      
        return pdb_struct
  
  

    def save_gn_to_pdb(self, gpcrdb=True, bw=False):
    
        #replace bfactor field of CA atoms with b-w numbers and save structure to file
        #file name has '_GPCRDB' added before the extension
        pdb_struct = PDBParser().get_structure(os.path.splitext(os.path.basename(self.pdb_file))[0], self.pdb_file)
        for chain in pdb_struct[0]:
            for residue in chain:
                if self.residues[chain.id].has_key(residue.id[1]):
                    if self.residues[chain.id][residue.id[1]].gpcrdb != 0.:
                        residue["CA"].set_bfactor(float(self.residues[chain.id][residue.id[1]].gpcrdb))
                    if self.residues[chain.id][residue.id[1]].bw != 0.:
                        residue["N"].set_bfactor(float(self.residues[chain.id][residue.id[1]].bw))
                    r = self.residues[chain.id][residue.id[1]]
        #get the basename, extension and export the pdb structure with b-w numbers
        root, ext = os.path.splitext(self.pdb_file)
        io=PDBIO()
        io.set_structure(pdb_struct)
        io.save("%s_GPCRDB%s" %(root, ext))
        
    
    
    def assign_generic_numbers(self):
        
        alignments = {}
        #blast search goes first, looping through all the chains
        for chain in self.pdb_seq.keys():
            alignments[chain] = self.blast.run(self.pdb_seq[chain])
            
        #map the results onto pdb sequence for every sequence pair from blast
        for chain in self.pdb_seq.keys():
            for alignment in alignments[chain]:
                if alignment == []:
                    continue
                for hsps in alignment[1].hsps:
                    self.map_blast_seq(alignment[0], hsps, chain)
                    #now kiss
                    for res_num in self.residues[chain].keys():
                        for residue in residues:
                            if self.gpcrdb is not None:
                                if self.residues[chain][res_num].get_mapping(alignment[0]) == int(residue['residuenumber']) and residue["residuenumberfamilyalternate"] != 'None':
                                    self.residues[chain][res_num].add_gpcrdb_number(residue["residuenumberfamilyalternate"])
                                    #print residue["residuenumberfamilyalternate"]
                                if self.residues[chain][res_num].get_mapping(alignment[0]) == int(residue['residuenumber']) and residue["residuenumberfamilyalternate2"] != 'None':
                                    self.residues[chain][res_num].add_bw_number(residue["residuenumberfamilyalternate2"])

                            else:
                                if self.residues[chain][res_num].get_mapping(alignment[0]) == residue['residueNumber'] and hasattr(residue, "residueNumberFamilyAlternate"):
                                    self.residues[chain][res_num].add_gpcrdb_number(residue["residueNumberFamilyAlternate"])
                                if self.residues[chain][res_num].get_mapping(alignment[0]) == residue['residueNumber'] and hasattr(residue, "residueNumberFamilyAlternate2"):
                                    self.residues[chain][res_num].add_bw_number(residue["residueNumberFamilyAlternate2"])

        return self.get_annotated_structure()