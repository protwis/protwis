from subprocess import Popen, PIPE
from io import StringIO
from Bio.Blast import NCBIXML
import Bio.PDB.Polypeptide as polypeptide

from django.conf import settings
from common.selection import SimpleSelection
from common.alignment import Alignment
from protein.models import ProteinSegment
from residue.models import Residue
from structure.models import Structure

import os
import sys
import tempfile
import logging

logger = logging.getLogger("structure")

#==============================================================================
# I have put it into separate class for the sake of future uses
class BlastSearch(object):
    
    
    def __init__ (self, blast_path='blastp', blastdb=os.sep.join([settings.STATICFILES_DIRS[0], 'blast', 'protwis_blastdb']), top_results=1):
  
        self.blast_path = blast_path
        self.blastdb = blastdb
        #print(blastdb)
        #typicaly top scored result is enough, but for sequences with missing
        #residues it is better to use more results to avoid getting sequence of
        #e.g.  different species
        self.top_results = top_results
      
    #takes Bio.Seq sequence as an input and returns a list of tuples with the
    #alignments
    def run (self, input_seq):
    
        output = []
        #Windows has problems with Popen and PIPE
        if sys.platform == 'win32':
            tmp = tempfile.NamedTemporaryFile()
            logger.debug("Running Blast with sequence: {}".format(input_seq))
            tmp.write(bytes(str(input_seq) + '\n', 'latin1'))
            tmp.seek(0)
            blast = Popen('%s -db %s -outfmt 5' % (self.blast_path, self.blastdb), universal_newlines=True, shell=True, stdin=tmp, stdout=PIPE, stderr=PIPE)
            (blast_out, blast_err) = blast.communicate()
        else:
        #Rest of the world:
            blast = Popen('%s -db %s -outfmt 5' % (self.blast_path, self.blastdb), universal_newlines=True, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            (blast_out, blast_err) = blast.communicate(input=str(input_seq))
        if len(blast_err) != 0:
            logger.debug(blast_err)

        result = NCBIXML.read(StringIO(blast_out))
        for aln in result.alignments[:self.top_results]:         
            logger.debug("Looping over alignments, current hit: {}".format(aln.hit_id))
            output.append((aln.hit_id, aln))
        return output

#==============================================================================

#stores information about alignments and b-w numbers
class MappedResidue(object):
  
    def __init__ (self, res_num, res_name):
          
        self.number = res_num
        self.name = res_name
        self.pos_in_aln = 0
        self.mapping = {}
        self.bw = 0.
        self.gpcrdb = 0.       
  
    def add_bw_number (self, bw_number=''):
    
        self.bw = bw_number


    def add_gpcrdb_number (self, gpcrdb_number=''):

        #PDB format does not allow fractional part longer than 2 digits
        #so numbers x.xx1 are negative
        if len(gpcrdb_number) > 4:
          self.gpcrdb = '-' + gpcrdb_number[:4].replace('x', '.')
        else:
          self.gpcrdb = gpcrdb_number.replace('x', '.')

#==============================================================================

#turns selection into actual residues
class SelectionParser(object):

    def __init__ (self, selection):
    
        self.generic_numbers = []
        self.helices = []
        
        self.db_segments = [x.slug for x in ProteinSegment.objects.all()]
        
        for segment in selection.segments:
            if 'TM' in segment.item.slug:
                self.helices.append(int(segment.item.slug[-1]))
            elif segment not in self.db_segments:
                self.residues.append(segment)
        logger.debug("Helices selected: {}; Residues: {}".format(self.helices, self.generic_numbers))
    

#==============================================================================
class CASelector(object):

    def __init__ (self, parsed_selection, ref_pdbio_struct, alt_structs):
    
        self.ref_atoms = []
        self.alt_atoms = {}
    
        self.selection = parsed_selection
        try:
            self.ref_atoms.extend(self.select_generic_numbers(self.selection.generic_numbers, ref_pdbio_struct[0]))
            self.ref_atoms.extend(self.select_helices(self.selection.helices, ref_pdbio_struct[0]))
        except Exception as msg:
            logger.warning("Can't select atoms from the reference structure!\n{!s}".format(msg))
    
        for alt_struct in alt_structs:
            try:
                self.alt_atoms[alt_struct.id] = []
                self.alt_atoms[alt_struct.id].extend(self.select_generic_numbers(self.selection.generic_numbers, alt_struct[0]))
                self.alt_atoms[alt_struct.id].extend(self.select_helices(self.selection.helices, alt_struct[0]))
                
            except Exception as msg:
                logger.warning("Can't select atoms from structure {!s}\n{!s}".format(alt_struct.id, msg))

    def select_generic_numbers (self, gn_list, structure):
    
        if gn_list == []:
            return []
    
        atom_list = []
    
        for chain in structure:
            for res in chain:
                try:
                    if -8.1 < res['CA'].get_bfactor() < 8.1 and res["CA"].bfactor in gn_list:
                        atom_list.append(res['CA'])
                except:
                    continue

        if atom_list == []:
            logger.warning("No atoms with given generic numbers for  {!s}".format(structure.id))
        return atom_list
    
    
    def select_helices (self, helices_list, structure):
    
        if helices_list == []:
            return []
    
        atom_list = []
        for chain in structure:
            for res in chain:
                try:
                    if -8.1 < res['CA'].get_bfactor() < 8.1 and int(math.floor(abs(res['CA'].get_bfactor()))) in helices_list:
                        atom_list.append(res['CA'])
                except:
                    continue

        if atom_list == []:
            logger.warning("No atoms with given generic numbers for  {!s}".format(structure.id))

        return atom_list


    def get_consensus_sets (self, alt_id):
        
        tmp_ref = []
        tmp_alt = []
    
        for ref_at in self.ref_atoms:
            for alt_at in self.alt_atoms[alt_id]:
                if ref_at.get_bfactor() == alt_at.get_bfactor():
                    tmp_ref.append(ref_at)
                    tmp_alt.append(alt_at)
    
        if len(tmp_ref) != len(tmp_alt):
            return ([], [])
    
        return (tmp_ref, tmp_alt)
    
    
    def get_ref_atoms (self):
    
        return self.ref_atoms
  
  
    def get_alt_atoms (self, alt_id):
    
        try:
            return self.alt_atoms[alt_id]
        except:
            return []


    def get_alt_atoms_all (self):
    
        return self.alt_atoms


#==============================================================================
class BackboneSelector():
    """
    Selects backbone atoms from reference structure and rotamer for Superposition
    """

    similarity_dict = {
        "fragment_residue" : 0,
        "interaction_type" : 1,
        "target_residue" : 2
        }

    similarity_rules = [[['H', 'F', 'Y', 'W'], ['AEF', 'AFF'], ['H', 'F', 'Y', 'W']],
        [['Y'], ['AFE'], ['F']],
        [['S', 'T'], ['HBA', 'HBD'], ['S', 'T']],]

    def __init__ (self, ref_pdbio_struct, fragment):

        self.ref_atoms = []
        self.alt_atoms = []
        
        self.ref_atoms = self.select_ref_atoms(fragment, ref_pdbio_struct[0])
        self.alt_atoms = self.select_alt_atoms(PDBParser(PERMISSIVE=True).get_structure('ref', StringIO(fragment.rotamer.pdbdata)))
        
        
    def select_ref_atoms (self, fragment, ref_pdbio_struct, use_similar=False):

        for chain in ref_pdbio_struct:
            for res in chain:
                if self.get_generic_number(res) == fragment.rotamer.residue.generic_number:
                    if use_similar:
                        for rule in self.similarity_rules:
                            if polypeptide.three_to_one(res.resname) in rule[self.similarity_dict["target_residue"]] and fragment.residue.amino_acid in rule[self.similarity_dict["target_residue"]] and fragment.interaction_type.slug in rule[self.similarity_dict["interaction_type"]]:
                                return [res['CA'], res['N'], res['O']] 
                    else:
                        return [res['CA'], res['N'], res['O']] 

        return []                  


    def select_alt_atoms (self, rotamer_pdbio_struct):

        for chain in rotamer_pdbio_struct:
            for res in chain:
                try:
                    return [res['CA'], res['N'], res['O']]
                except:
                    continue
        return []
    

    def get_generic_number (self, res):

        if 0 < res['CA'].get_bfactor() < 8.1:
            return "{:2f}x{:2f}".format(res['N'].get_bfactor(), res['CA'].get_bfactor())
        if -8.1 < res['CA'].get_bfactor() < 0:
            return "{:2f}x{:2f}".format(res['N'].get_bfactor(), -res['CA'].get_bfactor() + 0.001)
        return 0.0


    def get_ref_atoms (self):
    
        return self.ref_atoms
  
  
    def get_alt_atoms (self):
    
        return self.alt_atoms



#==============================================================================
def check_gn (pdb_struct):
        
    gn_list = []
    for chain in pdb_struct:
        for residue in chain:
            try:
                if -8.1 < residue['CA'].get_bfactor() < 8.1:
                    gn_list.append(residue['CA'])
                    return True
            except:
                continue
    return False


#==============================================================================
def get_segment_template (protein, segments=['TM1', 'TM2', 'TM3', 'TM4','TM5','TM6', 'TM7']):

    a = Alignment()
    a.load_reference_protein(protein)
    #You are so gonna love it...
    a.load_proteins([x.protein_conformation.protein.parent for x in list(Structure.objects.order_by('protein_conformation__protein__parent','resolution').exclude(protein_conformation__protein=protein.id))])
    a.load_segments(ProteinSegment.objects.filter(slug__in=segments))
    a.build_alignment()
    a.calculate_similarity()

    return a.proteins[1]


#==============================================================================
def fetch_template_structure (template_protein):

    return Structure.objects.get(protein_conformation__protein__parent=template_protein.entry_name)