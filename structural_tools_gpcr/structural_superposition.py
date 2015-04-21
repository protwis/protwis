import os,sys,math,logging

from Bio.PDB import *
from structural_tools_gpcr.common import SelectionParser
from structural_tools_gpcr.assign_generic_numbers import GenericNumbering



#Basic functionality uses pdb residue numbers provided by input
#==============================================================================
class CASelector(object):

    logger = logging.getLogger("structural_tools_gpcr")

    def __init__(self, parsed_selection, ref_struct, alt_structs):
    
        self.ref_atoms = []
        self.alt_atoms = {}
    
        self.selection = parsed_selection

        #Selecting the actual atoms from the selection string
        try:
              self.ref_atoms.extend(self.select_gn(self.selection.generic_numbers, ref_struct[0]))
              self.ref_atoms.extend(self.select_helices(self.selection.helices, ref_struct[0]))
        except:
              self.logger.warning("Can't select atoms from the reference structure!")
    
        for alt_struct in alt_structs:
            try:
                self.alt_atoms[alt_struct.id] = []
                self.alt_atoms[alt_struct.id].extend(self.select_gn(self.selection.generic_numbers, alt_struct[0]))
                self.alt_atoms[alt_struct.id].extend(self.select_helices(self.selection.helices, alt_struct[0]))
            except:
                self.logger.warning("Can't select atoms from structure {!s}".format(alt_struct.id))


    #GPCRdb numbers from annotated structure...
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
    
    
    #... and all residues from given helices
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


    #Just the tool functions for retrieving atom sets
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
  
    def __init__ (self, ref_file, alt_files, simple_selection):
    
        self.out_path = os.path.dirname(ref_file)
    
        self.selection = SelectionParser(simple_selection)
    
        #Setting up reference structure
        self.ref_struct = PDBParser().get_structure('ref', ref_file)
        #Kicking out if there is no reference
        assert self.ref_struct, "Can't parse the ref file %s" %ref_file
    
        #Checking if the reference structure is annotated with B-W numbers and running annotation if not
        if self.selection.generic_numbers != [] or self.selection.helices != []:
            if self.check_gn(self.ref_struct[0]) == False:
                gn_assigner = GenericNumbering(ref_file)
                self.ref_struct = gn_assigner.assign_generic_numbers()
      
        #Setting up structures to align
        self.alt_structs = []
        for alt_file in alt_files:
            self.out_path = os.path.dirname(alt_file)
      
        #Identifying structures by filenames
        alt_id = os.path.splitext(os.path.basename(alt_file))[0]
        #print alt_id
        try:
            tmp_struct = PDBParser().get_structure(alt_id, alt_file)
            #B-W check for alt structures as well
            if self.selection.bw_numbers != [] or self.selection.helices != []:
                if self.check_gn(tmp_struct[0]) == False:
                    gn_assigner = GenericNumbering(alt_file)
                    self.alt_structs.append(gn_assigner.assign_generic_numbers())
            else:
                self.alt_structs.append(tmp_struct)
    
        except:
            self.logger.warning("Can't parse the file {!s}".format(alt_file))
    
        self.selector = ca_selector(self.selection, self.ref_struct, self.alt_structs)


  def run (self):
    
    if self.alt_structs == []:
      print "No structures to align!"
      return None
    
    #The alignment itself
    super_imposer = Superimposer()
    for alt_struct in self.alt_structs:
      try:
        #print "Working on structure %s" %alt_struct.id
        #Selecting C-alphas based on selection string
        ref, alt = self.selector.get_consensus_sets(alt_struct.id)
        #super_imposer.set_atoms(self.selector.get_ref_atoms(), self.selector.get_alt_atoms(alt_struct.id))
        super_imposer.set_atoms(ref, alt)
        super_imposer.apply(alt_struct[0].get_atoms())
        print "RMS(first model, model %s) = %0.2f" % (alt_struct.id, super_imposer.rms)
        #print "Saving aligned structure as PDB file aligned-%s.pdb" % alt_struct.id
        #Saving aligned structures into files
        io=PDBIO()
        io.set_structure(alt_struct)
        io.save("%s/%s-superposed.pdb" % (self.out_path, alt_struct.id))
      except PDBExceptions.PDBException as e:
        print "Can't align the structures: %s" %e
        


  def check_gn(self, pdb_struct):
    gn_list = []
    
    #Checking for the first occurence of B-W notation and breaking the loop
    for chain in pdb_struct:
      for residue in chain:
        if residue.resname == "HID":
          resname = polypeptide.three_to_one('HIS')
        else:
          if residue.resname not in self.residue_list:
            continue
        if -8.1 < residue['CA'].get_bfactor() < 8.1:
          gn_list.append(residue['CA'])
          break
    
    if gn_list == []:
      return False
    return True

if __name__ == "__main__":
  (options, args) = parser.parse_args()
  
  aligner = prot_superpose(options.ref_file, args, options.sel_string)
  
  aligner.run()
  
