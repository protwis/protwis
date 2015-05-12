import os,sys,math,logging

from Bio.PDB import *
from structure_gpcr.common import SelectionParser
from structure_gpcr.assign_generic_numbers import GenericNumbering



#Basic functionality uses pdb residue numbers provided by input
#==============================================================================
class CASelector(object):

    logger = logging.getLogger("structure_gpcr")

    def __init__(self, parsed_selection, ref_struct, alt_structs):
    
        self.ref_atoms = []
        self.alt_atoms = {}
    
        self.selection = parsed_selection
        try:
            self.ref_atoms.extend(self.select_generic_numbers(self.selection.generic_numbers, ref_struct[0]))
            self.ref_atoms.extend(self.select_helices(self.selection.helices, ref_struct[0]))
            print(self.ref_atoms)
        except Exception as msg:
            self.logger.warning("Can't select atoms from the reference structure!\n{!s}".format(msg))
    
        for alt_struct in alt_structs:
            try:
                self.alt_atoms[alt_struct.id] = []
                self.alt_atoms[alt_struct.id].extend(self.select_generic_numbers(self.selection.generic_numbers, alt_struct[0]))
                self.alt_atoms[alt_struct.id].extend(self.select_helices(self.selection.helices, alt_struct[0]))
                
            except Exception as msg:
                self.logger.warning("Can't select atoms from structure {!s}\n{!s}".format(alt_struct.id, msg))

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
            self.logger.warning("No atoms with given generic numbers for  {!s}".format(structure.id))
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
            self.logger.warning("No atoms with given generic numbers for  {!s}".format(structure.id))

        return atom_list


    def get_consensus_sets(self, alt_id):
        
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
class ProteinSuperpose(object):
  
    logger = logging.getLogger("structure_gpcr")

    def __init__ (self, ref_file, alt_files, simple_selection):
    
        self.selection = SelectionParser(simple_selection)
    
        self.ref_struct = PDBParser().get_structure('ref', ref_file)
        assert self.ref_struct, self.logger.error("Can't parse the ref file %s".format(ref_file))
        if self.selection.generic_numbers != [] or self.selection.helices != []:
            if not self.check_gn(self.ref_struct[0]):
                gn_assigner = GenericNumbering(ref_file)
                self.ref_struct = gn_assigner.assign_generic_numbers()
      
        self.alt_structs = []
        for alt_id, alt_file in enumerate(alt_files):
            try:
                tmp_struct = PDBParser(PERMISSIVE=True).get_structure(alt_id, alt_file)
                if self.selection.generic_numbers != [] or self.selection.helices != []:
                    if not self.check_gn(tmp_struct[0]):
                        print("Assigning")
                        gn_assigner = GenericNumbering(alt_file)
                        self.alt_structs.append(gn_assigner.assign_generic_numbers())
                    else:
                        self.alt_structs.append(tmp_struct)
            except Exception as e:
                print(e)
                self.logger.warning("Can't parse the file {!s}\n{!s}".format(alt_id, e))
        self.selector = CASelector(self.selection, self.ref_struct, self.alt_structs)


    def run (self):
    
        if self.alt_structs == []:
            self.logger.error("No structures to align!")
            return []
    
        super_imposer = Superimposer()
        for alt_struct in self.alt_structs:
            try:
                ref, alt = self.selector.get_consensus_sets(alt_struct.id)
                super_imposer.set_atoms(ref, alt)
                super_imposer.apply(alt_struct[0].get_atoms())
                self.logger("RMS(first model, model {!s}) = {:d}".format(alt_struct.id, super_imposer.rms))
            except Exception as msg:
                self.logger.error("Failed to superpose structures {} and {}".format(self.ref_struct.id, alt_struct.id))

        return self.alt_structs


    def check_gn(self, pdb_struct):
        
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

class RotamerSuperpose(object):
    ''' Class to superimpose Atom objects on one-another. 

        @param original_rotamers: list of Atom objects of rotamers to be superposed on \n
        @param rotamers: list of Atom objects of rotamers to be superposed
    '''
    def __init__(self, reference_atoms, template_atoms):
        self.reference_atoms = reference_atoms
        self.template_atoms = template_atoms

    def run(self):
        ''' Run the superpositioning. 
        '''
        super_imposer = Superimposer()
        try:
            ref_backbone_atoms = [atom for atom in self.reference_atoms if atom.get_name() in ['N','CA','C','O']]
            temp_backbone_atoms = [atom for atom in self.template_atoms if atom.get_name() in ['N','CA','C','O']]
            super_imposer.set_atoms(ref_backbone_atoms, temp_backbone_atoms)
            super_imposer.apply(self.template_atoms)
            return self.template_atoms
        except Exception as msg:
            print("Failed to superpose atoms {} and {}".format(self.reference_atoms, self.template_atoms))
            return None

