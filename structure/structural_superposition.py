import os,sys,math,logging
from io import StringIO
from collections import OrderedDict
import numpy as np

import Bio.PDB.Polypeptide as polypeptide
from Bio.PDB import *
from Bio.Seq import Seq
from structure.functions import *
from structure.sequence_parser import *
from structure.assign_generic_numbers_gpcr import GenericNumbering
from protein.models import Protein
from structure.models import Structure
from interaction.models import ResidueFragmentInteraction

logger = logging.getLogger("protwis")

#==============================================================================
class ConvertSuperpose(object):

    def __init__ (self, ref_file, alt_files, simple_selection):

        self.selection = SelectionParser(simple_selection)
        print('Generic numbers beginning of ProteinSuperpose')
        self.ref_struct = PDBParser(PERMISSIVE=True).get_structure('ref', ref_file)[0]
        assert self.ref_struct, self.logger.error("Can't parse the ref file %s".format(ref_file))
        if self.selection.generic_numbers != [] or self.selection.helices != []:
            print('Assessing generic numbers')
            self.ref_struct = SequenceParser(self.ref_struct, db='g_protein_chimeras')

        self.alt_structs = []
        print('Parsing the targets')
        for alt_id, alt_file in enumerate(alt_files):
            try:
                tmp_struct = PDBParser(PERMISSIVE=True).get_structure(alt_id, alt_file)[0]
                if self.selection.generic_numbers != [] or self.selection.helices != []:
                    if not check_gn(tmp_struct):
                        print('Assessing generic numbers for target structure')
                        gn_assigner = GenericNumbering(structure=tmp_struct)
                        self.alt_structs.append(gn_assigner.assign_generic_numbers())
                        self.alt_structs[-1].id = alt_id
                    else:
                        self.alt_structs.append(tmp_struct)
            except Exception as e:
                print("Can't parse the file {!s}\n{!s}".format(alt_id, e))
                logger.warning("Can't parse the file {!s}\n{!s}".format(alt_id, e))
        print('Selection: {}. Ref Struct: {}. Alt Structs: {}'.format(self.selection, self.ref_struct, self.alt_structs))
        self.selector = CASelector(self.selection, self.ref_struct, self.alt_structs)

    def run (self):

        if self.alt_structs == []:
            logger.error("No structures to align!")
            return []

        super_imposer = Superimposer()
        for alt_struct in self.alt_structs:
            try:
                ref, alt = self.selector.get_consensus_atom_sets(alt_struct.id)
                super_imposer.set_atoms(ref, alt)
                super_imposer.apply(alt_struct.get_atoms())
                print("RMS(reference, model {!s}) = {:f}".format(alt_struct.id, super_imposer.rms))
                logger.info("RMS(reference, model {!s}) = {:f}".format(alt_struct.id, super_imposer.rms))
            except Exception as msg:
                print("Failed to superpose structures {} and {}\n{}".format(self.ref_struct.id, alt_struct.id, msg))
                logger.error("Failed to superpose structures {} and {}\n{}".format(self.ref_struct.id, alt_struct.id, msg))

        return self.alt_structs

#==============================================================================
class ProteinSuperpose(object):

    def __init__ (self, ref_file, alt_files, simple_selection):

        self.selection = SelectionParser(simple_selection)
        print('Generic numbers beginning of ProteinSuperpose')
        self.ref_struct = PDBParser(PERMISSIVE=True).get_structure('ref', ref_file)[0]
        assert self.ref_struct, self.logger.error("Can't parse the ref file %s".format(ref_file))
        if self.selection.generic_numbers != [] or self.selection.helices != []:
            if not check_gn(self.ref_struct):
                gn_assigner = GenericNumbering(structure=self.ref_struct)
                print('Assessing generic numbers')
                self.ref_struct = gn_assigner.assign_generic_numbers()

        self.alt_structs = []
        print('Parsing the targets')
        for alt_id, alt_file in enumerate(alt_files):
            try:
                tmp_struct = PDBParser(PERMISSIVE=True).get_structure(alt_id, alt_file)[0]
                if self.selection.generic_numbers != [] or self.selection.helices != []:
                    if not check_gn(tmp_struct):
                        print('Assessing generic numbers for target structure')
                        gn_assigner = GenericNumbering(structure=tmp_struct)
                        self.alt_structs.append(gn_assigner.assign_generic_numbers())
                        self.alt_structs[-1].id = alt_id
                    else:
                        self.alt_structs.append(tmp_struct)
            except Exception as e:
                print("Can't parse the file {!s}\n{!s}".format(alt_id, e))
                logger.warning("Can't parse the file {!s}\n{!s}".format(alt_id, e))
        print('Selection: {}. Ref Struct: {}. Alt Structs: {}'.format(self.selection, self.ref_struct, self.alt_structs))
        self.selector = CASelector(self.selection, self.ref_struct, self.alt_structs)

    def run (self):

        if self.alt_structs == []:
            logger.error("No structures to align!")
            return []

        super_imposer = Superimposer()
        for alt_struct in self.alt_structs:
            try:
                ref, alt = self.selector.get_consensus_atom_sets(alt_struct.id)
                super_imposer.set_atoms(ref, alt)
                super_imposer.apply(alt_struct.get_atoms())
                print("RMS(reference, model {!s}) = {:f}".format(alt_struct.id, super_imposer.rms))
                logger.info("RMS(reference, model {!s}) = {:f}".format(alt_struct.id, super_imposer.rms))
            except Exception as msg:
                print("Failed to superpose structures {} and {}\n{}".format(self.ref_struct.id, alt_struct.id, msg))
                logger.error("Failed to superpose structures {} and {}\n{}".format(self.ref_struct.id, alt_struct.id, msg))

        return self.alt_structs

#==============================================================================
class FragmentSuperpose(object):

    logger = logging.getLogger("structure")

    def __init__(self, pdb_file=None, pdb_filename=None):

        #pdb_file can be either a name/path or a handle to an open file
        self.pdb_file = pdb_file
        self.pdb_filename = pdb_filename
        self.pdb_seq = {}
        self.blast = BlastSearch()

        self.pdb_struct = self.parse_pdb()
        if not check_gn(self.pdb_struct):
            gn_assigner = GenericNumbering(structure=self.pdb_struct)
            self.pdb_struct = gn_assigner.assign_generic_numbers()
            self.target = Protein.objects.get(pk=gn_assigner.prot_id_list[0])
        else:
            self.target = Protein.objects.get(pk=self.identify_receptor())


    def parse_pdb (self):

        pdb_struct = None
        #checking for file handle or file name to parse
        if self.pdb_file:
            pdb_struct = PDBParser(PERMISSIVE=True, QUIET=True).get_structure('ref', self.pdb_file)[0]
        elif self.pdb_filename:
            pdb_struct = PDBParser(PERMISSIVE=True, QUIET=True).get_structure('ref', self.pdb_filename)[0]
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
            return self.blast.run(Seq(''.join([str(self.pdb_seq[x]) for x in sorted(self.pdb_seq.keys())])))[0][0]
        except Exception as msg:
            logger.error('Failed to identify protein for input file {!s}\nMessage: {!s}'.format(self.pdb_filename, msg))
            return None


    def superpose_fragments(self, representative=False, use_similar=False, state='inactive'):

        superposed_frags = [] #list of (fragment, superposed pdbdata) pairs
        if representative:
            fragments = self.get_representative_fragments(state)
        else:
            fragments = self.get_all_fragments()

        for fragment in fragments:
            atom_sel = BackboneSelector(self.pdb_struct, fragment, use_similar)
            if atom_sel.get_ref_atoms() == []:
                continue
            super_imposer = Superimposer()
            try:
                fragment_struct = PDBParser(PERMISSIVE=True, QUIET=True).get_structure('alt', StringIO(fragment.get_pdbdata()))[0]
                super_imposer.set_atoms(atom_sel.get_ref_atoms(), atom_sel.get_alt_atoms())
                super_imposer.apply(fragment_struct)
                superposed_frags.append([fragment,fragment_struct])
            except Exception as msg:
                logger.error('Failed to superpose fragment {!s} with structure {!s}\nDebug message: {!s}'.format(fragment, self.pdb_filename, msg))
        logger.info("Number of superimposed fragments: {}".format(len(superposed_frags)))
        return superposed_frags


    def get_representative_fragments(self, state):

        template = get_segment_template(self.target, state)
        return list(ResidueFragmentInteraction.objects.prefetch_related('rotamer__residue__display_generic_number', 'rotamer__residue', 'interaction_type').filter(structure_ligand_pair__structure__protein_conformation__protein=template.id))


    def get_all_fragments(self):

        return list(ResidueFragmentInteraction.objects.exclude(structure_ligand_pair__structure__protein_conformation__protein__parent=self.target).exclude(interaction_type__slug__in=['acc', 'hyd']).prefetch_related('rotamer__residue__display_generic_number', 'rotamer__residue', 'interaction_type'))

#==============================================================================
class RotamerSuperpose(object):
    ''' Class to superimpose Atom objects on one-another.

        @param reference_atoms: list of Atom objects of rotamers to be superposed on \n
        @param template_atoms: list of Atom objects of rotamers to be superposed
    '''
    def __init__(self, reference_atoms, template_atoms, TM_keys=None):
        self.reference_atoms = reference_atoms
        self.template_atoms = template_atoms
        self.backbone_rmsd, self.rmsd = None, None
        self.TM_keys = TM_keys
        self.num_atoms_used_for_superposition = 0

    def run(self):
        ''' Run the superpositioning.
        '''
        super_imposer = Superimposer()
        try:
            if not self.TM_keys:
                ref_backbone_atoms = [atom for atom in self.reference_atoms if atom.get_name() in ['N','CA','C','O']]
                temp_backbone_atoms = [atom for atom in self.template_atoms if atom.get_name() in ['N','CA','C','O']]
            else:
                ref_backbone_atoms = [atom for atom in self.reference_atoms if atom.get_name() in ['N','CA','C'] and atom.get_parent().get_full_id()[-1][1] in self.TM_keys]
                temp_backbone_atoms = [atom for atom in self.template_atoms if atom.get_name() in ['N','CA','C'] and atom.get_parent().get_full_id()[-1][1] in self.TM_keys]
            self.num_atoms_used_for_superposition = len(ref_backbone_atoms)
            super_imposer.set_atoms(ref_backbone_atoms, temp_backbone_atoms)
            super_imposer.apply(self.template_atoms)
            array1, array2 = np.array([0,0,0]), np.array([0,0,0])
            for atom1, atom2 in zip(ref_backbone_atoms, temp_backbone_atoms):
                array1 = np.vstack((array1, list(atom1.get_coord())))
                array2 = np.vstack((array2, list(atom2.get_coord())))
            diff = array1[1:]-array2[1:]
            self.backbone_rmsd = np.sqrt(sum(sum(diff**2))/array1[1:].shape[0])
            for atom1, atom2 in zip(self.reference_atoms, self.template_atoms):
                array1 = np.vstack((array1, list(atom1.get_coord())))
                array2 = np.vstack((array2, list(atom2.get_coord())))
            diff = array1[1:]-array2[1:]
            self.rmsd = np.sqrt(sum(sum(diff**2))/array1[1:].shape[0])
            return self.template_atoms
        except Exception as msg:
            if self.reference_atoms!='x':
                print("Failed rotamer superimposition:\n{}".format(msg))

#==============================================================================
class BulgeConstrictionSuperpose(object):
    ''' Class to superimpose bulge and constriction site.

        @param reference_dict: OrderedDict, dictionary of atoms to be superposed on, where keys are generic numbers
        and values are lists of atoms. \n
        @param template_dict: OrderedDict, dictionary of atoms to be superposed. Same format as reference_dict.
    '''
    def __init__(self, reference_dict, template_dict):
        self.reference_dict = reference_dict
        self.reference_gns = list(reference_dict.keys())
        self.template_dict = template_dict
        self.template_gns = list(template_dict.keys())
        self.starting_atom_type = template_dict[list(template_dict.keys())[0]][0].get_id()
        self.backbone_rmsd = None

    def run(self):
        ''' Run the superpositioning.
        '''
        super_imposer = Superimposer()
        ref_backbone_atoms = [atom for atom in self.reference_dict[self.reference_gns[0]] if atom.get_name() in
                                ['N','CA','C']] + [atom for atom in self.reference_dict[self.reference_gns[-1]] if
                                atom.get_name() in ['N','CA','C']]
        temp_backbone_atoms= [atom for atom in self.template_dict[self.template_gns[0]] if atom.get_name() in
                                ['N','CA','C']] + [atom for atom in self.template_dict[self.template_gns[-1]] if
                                atom.get_name() in ['N','CA','C']]
        all_template_atoms = []
        for gn, atoms in self.template_dict.items():
            all_template_atoms+=atoms
        super_imposer.set_atoms(ref_backbone_atoms, temp_backbone_atoms)
        super_imposer.apply(all_template_atoms)
        return self.rebuild_dictionary(all_template_atoms)

    def rebuild_dictionary(self, all_template_atoms):
        ''' Rebuild input ordered dictionary.
        '''
        residue = []
        temp_dict = OrderedDict()
        key_count = 0
        for atom in all_template_atoms:
            if atom.get_id()==self.starting_atom_type and residue!=[]:
                key_count+=1
                temp_dict[key_count] = residue
                residue = []
            residue.append(atom)
        temp_dict[key_count+1] = residue
        gn_count = 0
        for gn in self.template_gns:
            gn_count+=1
            self.template_dict[gn] = temp_dict[gn_count]
        return self.template_dict

    def calc_backbone_RMSD(self, ref_backbone_atoms, temp_backbone_atoms):
        ''' Calculate backbone RMSD.
        '''
        array1, array2 = np.array([0,0,0]), np.array([0,0,0])
        for atom1, atom2 in zip(ref_backbone_atoms, temp_backbone_atoms):
            array1 = np.vstack((array1, list(atom1.get_coord())))
            array2 = np.vstack((array2, list(atom2.get_coord())))
        diff = array1[1:]-array2[1:]
        return np.sqrt(sum(sum(diff**2))/array1[1:].shape[0])

#==============================================================================
class LoopSuperpose(BulgeConstrictionSuperpose):
    ''' Class to superpose loop regions on helix endings.
    '''
    def __init__(self, reference_dict, template_dict, ECL2=False, part=None):
        super(LoopSuperpose, self).__init__(reference_dict=reference_dict, template_dict=template_dict)
        self.ECL2 = ECL2
        self.part = part

    def run(self):
        ''' Run the superpositioning.
        '''
        super_imposer = Superimposer()
        ref_backbone_atoms, temp_backbone_atoms, all_template_atoms = [], [], []
        for gn, atoms in self.reference_dict.items():
            for atom in atoms:
                if atom.get_name() in ['N','CA','C']:
                    ref_backbone_atoms.append(atom)
        res_count=0
        array_length = len(self.template_dict.keys())
        edge1 = 4
        edge2 = 4
        if self.ECL2==True:
            if self.part==1:
                edge2 = 3
            elif self.part==2:
                edge1 = 3
        for gn, atoms in self.template_dict.items():
            res_count+=1
            for atom in atoms:
                if (res_count<=edge1 or array_length-edge2<res_count) and atom.get_name() in ['N','CA','C']:
                    temp_backbone_atoms.append(atom)
                all_template_atoms.append(atom)

        super_imposer.set_atoms(sorted(ref_backbone_atoms), sorted(temp_backbone_atoms))
        super_imposer.apply(all_template_atoms)
        self.backbone_rmsd = self.calc_backbone_RMSD(sorted(ref_backbone_atoms), sorted(temp_backbone_atoms))
        return self.rebuild_dictionary(all_template_atoms)

#==============================================================================
class OneSidedSuperpose(BulgeConstrictionSuperpose):
    ''' Class for one sided superposition. Used for helix ends and N- and C-terminus.
    '''
    def __init__(self, reference_dict, template_dict, num_frame, which_end):
        super(OneSidedSuperpose, self).__init__(reference_dict=reference_dict, template_dict=template_dict)
        self.num_frame = num_frame
        self.which_end = which_end

    def run(self):
        ''' Run the superpositioning.
        '''
        super_imposer = Superimposer()
        ref_backbone_atoms, temp_backbone_atoms, all_template_atoms = [], [], []
        for gn, atoms in self.reference_dict.items():
            for atom in atoms:
                if atom.get_name() in ['N','CA','C']:
                    ref_backbone_atoms.append(atom)
        res_count = 0
        if self.which_end==0:
            start = len(self.template_dict.keys())-self.num_frame
            end = start+self.num_frame
        elif self.which_end==1:
            start = 0
            end = self.num_frame-1
        for gn, atoms in self.template_dict.items():
            for atom in atoms:
                if start<=res_count<=end and atom.get_name() in ['N','CA','C']:
                    temp_backbone_atoms.append(atom)
                all_template_atoms.append(atom)
            res_count+=1
        super_imposer.set_atoms(ref_backbone_atoms, temp_backbone_atoms)
        super_imposer.apply(all_template_atoms)
        self.backbone_rmsd = self.calc_backbone_RMSD(ref_backbone_atoms, temp_backbone_atoms)
        return self.rebuild_dictionary(all_template_atoms)

#==============================================================================
class ECL2MidSuperpose(BulgeConstrictionSuperpose):
    ''' Class to superimpose 45x50-52 in ECL2 based on last residue of TM4, first residue of TM5 and 3x25 in TM3.
    '''

    def run(self):
        ''' Run the superpositioning.
        '''
        super_imposer = Superimposer()
        ref_backbone_atoms, temp_backbone_atoms, all_template_atoms = [], [], []
        for gn, atoms in self.reference_dict.items():
            for atom in atoms:
                if atom.get_name() in ['N','CA','C']:
                    ref_backbone_atoms.append(atom)
        res_count=0
        for gn, atoms in self.template_dict.items():
            res_count+=1
            for atom in atoms:
                if res_count<4 and atom.get_name() in ['N','CA','C']:
                    temp_backbone_atoms.append(atom)
                all_template_atoms.append(atom)
        self.backbone_rmsd = self.calc_backbone_RMSD(ref_backbone_atoms, temp_backbone_atoms)
        super_imposer.set_atoms(ref_backbone_atoms, temp_backbone_atoms)
        super_imposer.apply(all_template_atoms)
        return self.rebuild_dictionary(all_template_atoms)
