import os,sys,math,logging
from io import StringIO

import Bio.PDB.Polypeptide as polypeptide
from Bio.PDB import *
from Bio.Seq import Seq
from structure.functions import * #SelectionParser,BlastSearch,check_gn,get_segment_template
from structure_gpcr.assign_generic_numbers import GenericNumbering
from protein.models import Protein
from structure.models import Structure
from interaction.models import ResidueFragmentInteraction


#==============================================================================  
class ProteinSuperpose(object):
  
    logger = logging.getLogger("structure_gpcr")

    def __init__ (self, ref_file, alt_files, simple_selection):
    
        self.selection = SelectionParser(simple_selection)
    
        self.ref_struct = PDBParser().get_structure('ref', ref_file)
        assert self.ref_struct, self.logger.error("Can't parse the ref file %s".format(ref_file))
        if self.selection.generic_numbers != [] or self.selection.helices != []:
            if not check_gn(self.ref_struct[0]):
                gn_assigner = GenericNumbering(ref_file)
                self.ref_struct = gn_assigner.assign_generic_numbers()
      
        self.alt_structs = []
        for alt_id, alt_file in enumerate(alt_files):
            try:
                tmp_struct = PDBParser(PERMISSIVE=True).get_structure(alt_id, alt_file)
                if self.selection.generic_numbers != [] or self.selection.helices != []:
                    if not check_gn(tmp_struct[0]):
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
                self.logger.info("RMS(first model, model {!s}) = {:d}".format(alt_struct.id, super_imposer.rms))
            except Exception as msg:
                self.logger.error("Failed to superpose structures {} and {}".format(self.ref_struct.id, alt_struct.id))

        return self.alt_structs



#==============================================================================  
class FragmentSuperpose(object):

    logger = logging.getLogger("structure_gpcr")

    def __init__(self, pdb_file=None, pdb_filename=None):
        
        #pdb_file can be either a name/path or a handle to an open file
        self.pdb_file = pdb_file
        self.pdb_filename = pdb_filename
        self.pdb_seq = {}
        self.blast = BlastSearch()

        self.pdb_struct = self.parse_pdb()
        if not check_gn(self.pdb_struct):
            gn_assigner = GenericNumbering(self.pdb_file, self.pdb_filename)
            self.pdb_struct = gn_assigner.assign_generic_numbers()

        self.target = Protein.objects.get(pk=self.identify_receptor())


    def parse_pdb (self):

        pdb_struct = None
        #checking for file handle or file name to parse
        if self.pdb_file:
            pdb_struct = PDBParser(PERMISSIVE=True).get_structure('ref', self.pdb_file)[0]
        elif self.pdb_filename:
            pdb_struct = PDBParser(PERMISSIVE=True).get_structure('ref', self.pdb_file)[0]
        else:
            return None

        #extracting sequence and preparing dictionary of residues
        #bio.pdb reads pdb in the following cascade: model->chain->residue->atom
        for chain in pdb_struct:
            self.pdb_seq[chain.id] = Seq('')            
            for res in chain:
            #in bio.pdb the residue's id is a tuple of (hetatm flag, residue number, insertion code)
                if res.resname == "HID":
                    self.pdb_seq[chain.id] += polypeptide.three_to_one('HIS')
                else:
                    try:
                        self.pdb_seq[chain.id] += polypeptide.three_to_one(res.resname)
                    except Exception as msg:
                        continue
        return pdb_struct


    def identify_receptor(self):
        
        try:
            return self.blast.run(Seq(''.join([self.pdb_seq[x].seq for x in sorted(self.pdb_seq.keys())])))[0]        
        except Exception as msg:
            self.logger.error('Failed to identify protein for input file {!s}'.format(self.pdb_filename))
            return None


    def superpose_fragments(self, representative=False, use_similar=False):

        superposed_frags = [] #list of (fragment, superposed pdbdata) pairs
        if representative:
            fragments = self.get_representative_fragments()
        else:
            fragments = self.get_all_fragments()

        for fragment in fragments:
            atom_sel = BackboneSelector(self.pdb_struct, fragment, use_similar)
            if atom_sel.get_ref_atoms() == []:
                continue
            super_imposer = Superimposer()
            try:
                fragment_struct = PDBParser(PERMISSIVE=True).get_structure('alt', StringIO(fragment.get_pdbdata()))[0]
                super_imposer.set_atoms(atom_sel.get_ref_atoms(), atom_sel.get_alt_atoms())
                super_imposer.apply(fragment_struct)
                superposed_frags.append([fragment,fragment_struct])
            except Exception as msg:
                self.logger.error('Failed to superpose fragment {!s} with structure {!s}\nDebug message: {!s}'.format(fragment, self.pdb_filename, msg))
        return superposed_frags


    def get_representative_fragments(self):
        
        template = get_segment_template(self.target)
        return list(ResidueFragmentInteraction.filter(structure_ligand_pair__structure=template.id))


    def get_all_fragments(self):

        return list(ResidueFragmentInteraction.filter(structure_ligand_pair__structure__protein_conformation__protein__parent__not=self.target))



#==============================================================================  
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