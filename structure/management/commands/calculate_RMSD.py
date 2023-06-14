# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 15:48:26 2015

@author: Gaspar Pandy
"""
from django.core.management.base import BaseCommand

import structure.structural_superposition as sp
import structure.assign_generic_numbers_gpcr as as_gn
from structure.functions import run_residue_flip, atoms_to_dict
from residue.models import Residue
from protein.models import Protein

import Bio.PDB as PDB
from Bio import SeqIO
from collections import OrderedDict
import numpy as np
from copy import deepcopy


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('files', help='Add any number of files as arguments. First one has to be the reference file.',
                            type=str, nargs='+')
        parser.add_argument('-p', help='Protein entry name.', type=str)
        parser.add_argument('-n', help='Specify residue sequence numbers to compare.', type=str, default=False, nargs='+')
        parser.add_argument('-r', help='Specify a range of residue sequence numbers to compare. Format: e.g. 1-300', type=str, default=False, nargs='+')
        parser.add_argument('-c', help='Specify chain ID. If not specified, the program will try to find one that matches.', type=str, default=False)
        parser.add_argument('--superpose_on', help='Specify superposition type. Supported options: base, nflag', type=str, default=False)
        parser.add_argument('--only_backbone', help='Calculate only the backbone atoms RMSD for the custom set.', action='store_true', default=False)

    def handle(self, *args, **options):
        v = Validation()
        if options['r']:
            seq_nums = []
            for r in options['r']:
                start, end = r.split('-')
                seq_nums+=[str(i) for i in list(range(int(start),int(end)+1))]
        else:
            seq_nums = options['n']
        if seq_nums==False:
            if options['c']==False:
                v.run_RMSD_list(options['files'], options['p'], superpose_on=options['superpose_on'], only_backbone=options['only_backbone'])
            else:
                v.run_RMSD_list(options['files'], options['p'], force_chain=options['c'], superpose_on=options['superpose_on'], only_backbone=options['only_backbone'])
        else:
            if options['c']==False:
                v.run_RMSD_list(options['files'], options['p'], seq_nums=seq_nums, superpose_on=options['superpose_on'], only_backbone=options['only_backbone'])
            else:
                v.run_RMSD_list(options['files'], options['p'], seq_nums=seq_nums, force_chain=options['c'], superpose_on=options['superpose_on'], only_backbone=options['only_backbone'])


class Validation():
    def __init__(self):
        pass

    def run_RMSD_list(self, files, protein, seq_nums=None, force_chain=None, superpose_on=False, only_backbone=False):
        """Calculates 3 RMSD values between a list of GPCR pdb files.

        It compares the files using sequence and generic numbers.
        First file in the list has to be the reference file.
        Params:
            @protein: UniProt entry name of GPCR or G protein, str
            @seq_nums: Specified list of sequence residue numbers for the Custom calculation, list
            @force_chain: Specify one letter chain name to use in the pdb files, str
            @superpose_on: Which backbone (N, CA, C) atoms to superimpose on (7TM, nflag), str
            @only_backbone: Calculate RMSD for only the backbone atoms, boolean
        """
        parser = PDB.PDBParser(QUIET=True)
        count = 0
        pdbs = []
        for f in files:
            count+=1
            pdb = parser.get_structure('struct{}'.format(count), f)[0]
            print(f)
            assign_gn = as_gn.GenericNumbering(pdb_file=f, sequence_parser=True)
            pdb = assign_gn.assign_generic_numbers_with_sequence_parser()
            pdbs.append(pdb)
        chains = []
        for p in pdbs:
            this = []
            for c in p.get_chains():
                this.append(c.get_id())
            chains.append(sorted(this))
        if force_chain:
            chains = [[force_chain] for i in chains]
        usable_chains = []
        for m in chains[1:]:
            for c in m:
                if c in chains[0]:
                    usable_chains.append(c)

        arrays = []
        model_counter = 0
        ### Creating full arrays
        for p in pdbs:
            try:
                if pdbs.index(p)==0 and len(usable_chains)==0:
                    chain = [c.get_id() for c in pdbs[0].get_chains()][0]
                else:
                    chain = p[usable_chains[0]].get_id()
            except KeyError:
                try:
                    chain = p[' '].get_id()
                except:
                    chain = p['A'].get_id()
            if force_chain and model_counter==0:
                chain = force_chain
            pdb_array1 = OrderedDict()
            for residue in p[chain]:
                if residue.get_full_id()[3][0]!=' ':
                    continue
                pdb_array1[int(residue.get_id()[1])] = residue
            arrays.append(pdb_array1)
            model_counter+=1

        ### Checking available residues in target and models
        all_deletes, all_keep = [], []
        for res in arrays[0]:
            for m in arrays[1:]:
                if res not in m:
                    all_deletes.append(res)
                else:
                    all_keep.append(res)
        ### Making unique list with residues present in all structures
        unique_nums = []
        for i in all_keep:
            if i not in unique_nums:
                unique_nums.append(i)
        all_keep = [i for i in unique_nums if i not in all_deletes]

        print('Residue sequence numbers present in all structures: {}'.format(len(all_keep)))
        print(all_keep)

        ### Checking available atoms in target and models
        atoms_to_keep, atoms_to_delete = OrderedDict(), OrderedDict()
        for target_resnum, target_res in arrays[0].items():
            atoms_to_keep[target_resnum] = []
            atoms_to_delete[target_resnum] = []
            for model in arrays[1:]:
                for model_resnum, model_res in model.items():
                    if target_resnum==model_resnum:
                        for atom in model_res:
                            if atom.id in target_res and atom.id not in atoms_to_keep[target_resnum]:
                                atoms_to_keep[target_resnum].append(atom.id)
                            elif atom.id not in target_res and atom.id not in atoms_to_delete[target_resnum]:
                                atoms_to_delete[target_resnum].append(atom.id)
                        for t_atom in target_res:
                            if t_atom.id in model_res and t_atom.id not in atoms_to_keep[target_resnum]:
                                atoms_to_keep[target_resnum].append(t_atom.id)
                            elif t_atom.id not in model_res and t_atom.id not in atoms_to_delete[target_resnum]:
                                atoms_to_delete[target_resnum].append(t_atom.id)
                        break

        ### Creating atom lists of structures
        atom_lists = []
        for m in arrays:
            atom_list = []
            for num, res in m.items():
                if num in all_keep and num not in all_deletes:
                    for atom in res:
                        if atom.id in atoms_to_keep[num] and atom.id not in atoms_to_delete[num]:
                            atom_list.append(atom)
            atom_lists.append(atom_list)

        ### Fetching data from GPCRdb
        protein = Protein.objects.get(entry_name=protein)
        # Base definition of receptors: 7TM
        if protein.family.slug.startswith('00'):
            base_nums = Residue.objects.filter(protein_conformation__protein=protein, protein_segment__slug__in=['TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7']).values_list('sequence_number', flat=True)
            base_type = '7TM'
        # Base definition of G proteins: HN and H5
        elif protein.family.slug.startswith('100'):
            base_nums = Residue.objects.filter(protein_conformation__protein=protein, protein_segment__slug__in=['H5']).values_list('sequence_number', flat=True)
            base_type = 'HN,H5'
        base_target_atom_list = [i for i in atom_lists[0] if i.get_parent().id[1] in base_nums]
        base_target_backbone_atom_list = [i for i in base_target_atom_list if i.id in ['N','CA','C']]
        base_atom_num = len(base_target_atom_list)
        print('{}_atom_num:'.format(base_type),base_atom_num)
        print('{}_backbone_atom_num:'.format(base_type),len(base_target_backbone_atom_list))
        print([i.get_parent() for i in base_target_backbone_atom_list])

        ### Running superposition and RMSD calculation
        target_dict = atoms_to_dict(atom_lists[0])
        c = 2
        for m in atom_lists[1:]:
            # Check for flipped residues
            atom_dict = deepcopy(atoms_to_dict(m))
            for i, j in zip(target_dict, atom_dict):
                if target_dict[i][0].get_parent().get_resname() in ['TYR','PHE','ARG','ASP','GLU']:
                    # print(i, target_dict[i][0].get_parent())
                    to_sup = deepcopy(atom_dict[j])
                    res_sup, _ = self.superpose(sorted(target_dict[i]), sorted(to_sup), [i])
                    rmsd = self.calc_RMSD(sorted(target_dict[i]), sorted(res_sup))
                    res_sup2 = run_residue_flip(res_sup)
                    rmsd2 = self.calc_RMSD(sorted(target_dict[i]), sorted(res_sup2))
                    if rmsd-rmsd2>=0.7:
                        atom_dict[j] = run_residue_flip(atom_dict[j])

            m = []
            for i, j in atom_dict.items():
                m+=j

            # print(len(m),len(m2))
            # for i,j in zip(m, m2):
            #     if i.get_coord()[0]!=j.get_coord()[0]:
            #         print(i.get_parent(), i,i.get_coord(),j.get_coord())
            print('########################################')
            print('Superposition base type for RMSD is: {}'.format(base_type))
            print('Model {}'.format(c-1))
            if seq_nums:
                seq_nums = [int(s) for s in seq_nums]
            else:
                seq_nums = all_keep
            ### Custom calculation
            if superpose_on:
                if superpose_on=='base':
                    custom_superposed, atoms_used_sp = self.superpose(sorted(atom_lists[0]), sorted(m), list(base_nums))
                elif superpose_on=='nflag':
                    custom_superposed, atoms_used_sp = self.superpose(sorted(atom_lists[0]), sorted(m), list(seq_nums))
                superposed = self.fetch_atoms_with_seqnum(custom_superposed, seq_nums, only_backbone)
                target_atoms = self.fetch_atoms_with_seqnum(atom_lists[0], seq_nums, only_backbone)
                rmsd = self.calc_RMSD(target_atoms, superposed)
            else:
                superposed, atoms_used_sp = self.superpose(sorted(atom_lists[0]), sorted(m))
                target_atoms = atom_lists[0]

            # for t, s in zip(sorted(target_atoms), sorted(superposed)):
            #     print(t, t.get_coord(), s, s.get_coord())
            rmsd = self.calc_RMSD(sorted(target_atoms), sorted(superposed))
            print('Num atoms sent for superposition: ', len(atom_lists[0]), len(m))
            print('Num atoms used for superposition: ', atoms_used_sp)
            print('Num atoms used for RMSD: ', len(target_atoms), len(superposed))
            print('Custom RMSD:', rmsd)

            # ### 7TM all atoms calculation
            TM_model_atom_list = [i for i in m if i.get_parent().id[1] in base_nums]
            superposed2, atoms_used_sp = self.superpose(sorted(base_target_atom_list), sorted(TM_model_atom_list), list(base_nums))
            rmsd = self.calc_RMSD(sorted(base_target_atom_list), sorted(superposed2))
            print('Num atoms sent for superposition: ', len(base_target_atom_list), len(TM_model_atom_list))
            print('Num atoms used for superposition: ', atoms_used_sp)
            print('Num atoms used for RMSD: ', len(base_target_atom_list), len(superposed2))
            print('{} all RMSD:'.format(base_type), rmsd)

            # ### 7TM only backbone (N, CA, C) calculation
            superposed, atoms_used_sp = self.superpose(sorted(base_target_atom_list), sorted(TM_model_atom_list), list(base_nums))
            superposed3 = self.fetch_atoms_with_seqnum(superposed, list(base_nums), True)
            rmsd = self.calc_RMSD(sorted(base_target_backbone_atom_list), sorted(superposed3))
            print('Num atoms sent for superposition: ', len(base_target_atom_list), len(TM_model_atom_list))
            print('Num atoms used for superposition: ', atoms_used_sp)
            print('Num atoms used for RMSD: ', len(base_target_backbone_atom_list), len(superposed3))
            print('{} backbone RMSD:'.format(base_type), rmsd)

            c+=1

    def fetch_atoms_with_seqnum(self, atom_list, seq_nums, only_backbone=False):
        """Gets atoms from list2 based on list1 resnums. Get only N, CA, C atoms of backbone when setting only_backbone to True."""
        out_list = []
        # print('Only backbone:', only_backbone)
        for i in atom_list:
            if i.get_parent().id[1] in seq_nums:
                if only_backbone and i.id in ['N','CA','C']:
                    out_list.append(i)
                elif not only_backbone:
                    out_list.append(i)
        return out_list

    def superpose(self, list1, list2, TM_keys=None):
        """Superposition, supply TM_keys (list of sequence residue numbers) when superimposing on only those atoms."""
        superpose = sp.RotamerSuperpose(list1, list2, TM_keys)
        return superpose.run(), superpose.num_atoms_used_for_superposition

    def calc_RMSD(self, list1, list2):
        """Calculates RMSD between two atoms lists. The two lists have to have the same length."""
        array1, array2 = np.array([0,0,0]), np.array([0,0,0])
        for a1, a2 in zip(list1, list2):
            array1 = np.vstack((array1, list(a1.get_coord())))
            array2 = np.vstack((array2, list(a2.get_coord())))
        rmsd = round(np.sqrt(sum(sum((array1[1:]-array2[1:])**2))/array1[1:].shape[0]),3)
        return rmsd

    ### Deprecated
    # def run_RMSD_list_archived(self, files, seq_nums=None, force_chain=None):
    #     """Calculates 4 RMSD values between a list of GPCR pdb files.

    #     It compares the files using sequence and generic numbers.
    #     First file in the list has to be the reference file.
    #         1. overall all atoms RMSD
    #         2. overall backbone atoms RMSD
    #         3. 7TM all atoms RMSD
    #         4. 7TM backbone atoms RMSD
    #     """
    #     c = 0
    #     for f in files:
    #         c+=1
    #         if c==1:
    #             self.number_of_residues_superposed['reference'] = OrderedDict()
    #             self.number_of_atoms_superposed['reference'] = OrderedDict()
    #             self.rmsds['reference'] = OrderedDict()
    #         else:
    #             self.number_of_residues_superposed['file{}'.format(str(c))] = OrderedDict()
    #             self.number_of_atoms_superposed['file{}'.format(str(c))] = OrderedDict()
    #             self.rmsds['file{}'.format(str(c))] = OrderedDict()
    #     parser = PDB.PDBParser(QUIET=True)
    #     count = 0
    #     pdbs = []
    #     for f in files:
    #         count+=1
    #         pdb = parser.get_structure('struct{}'.format(count), f)[0]
    #         assign_gn = as_gn.GenericNumbering(pdb_file=f, sequence_parser=True)
    #         pdb = assign_gn.assign_generic_numbers_with_sequence_parser()
    #         pdbs.append(pdb)
    #     chains = []
    #     for p in pdbs:
    #         this = []
    #         for c in p.get_chains():
    #             this.append(c.get_id())
    #         chains.append(this)
    #     usable_chains = []
    #     for m in chains[1:]:
    #         for c in m:
    #             if c in chains[0]:
    #                 usable_chains.append(c)
    #     if force_chain:
    #         chains[0] = [force_chain]

    #     arrays = []
    #     model_counter = 0
    #     for p in pdbs:
    #         try:
    #             if pdbs.index(p)==0 and len(usable_chains)==0:
    #                 chain = [c.get_id() for c in pdbs[0].get_chains()][0]
    #             else:
    #                 chain = p[usable_chains[0]].get_id()
    #         except:
    #             try:
    #                 chain = p[' '].get_id()
    #             except:
    #                 chain = p['A'].get_id()
    #         if force_chain and model_counter==0:
    #             chain = force_chain
    #         pdb_array1, pdb_array2 = OrderedDict(), OrderedDict()
    #         for residue in p[chain]:
    #             if residue.get_full_id()[3][0]!=' ':
    #                 continue
    #             if seq_nums!=None and str(residue.get_id()[1]) in seq_nums:
    #                 pdb_array1[int(residue.get_id()[1])] = residue
    #             elif seq_nums==None:
    #                 pdb_array1[int(residue.get_id()[1])] = residue
    #             try:
    #                 if -8.1 < residue['CA'].get_bfactor() < 8.1:
    #                     pdb_array2[int(residue.get_id()[1])] = residue
    #             except:
    #                 pass
    #         arrays.append([pdb_array1,pdb_array2])
    #         model_counter+=1

    #     all_deletes, TM_deletes = [], []
    #     all_keep, TM_keep = [], []
    #     for i in range(0,2):
    #         for res in arrays[0][i]:
    #             for m in arrays[1:]:
    #                 if res not in m[i]:
    #                     if i==0:
    #                         all_deletes.append(res)
    #                     else:
    #                         TM_deletes.append(res)
    #                 else:
    #                     if i==0:
    #                         all_keep.append(res)
    #                     else:
    #                         TM_keep.append(res)
    #     deletes = [all_deletes, TM_deletes]
    #     keeps = [all_keep, TM_keep]

    #     num_atoms1, num_atoms2 = OrderedDict(), OrderedDict()
    #     num_atoms = [num_atoms1, num_atoms2]
    #     mismatches = []
    #     resis_to_delete = []
    #     for m_i, m in enumerate(arrays):
    #         for i in range(0,2):
    #             for res in m[i]:
    #                 if res in deletes[i] or res not in keeps[i]:
    #                     resis_to_delete.append([m_i,i,res])
    #                 else:
    #                     try:
    #                         if m[i][res].get_resname()!=num_atoms[i][res][0].get_parent().get_resname():
    #                             del num_atoms[i][res]
    #                             mismatches.append(res)
    #                         else:
    #                             raise Exception()
    #                     except:
    #                         if res not in mismatches:
    #                             atoms = []
    #                             for atom in m[i][res]:
    #                                 atoms.append(atom)
    #                             if res not in num_atoms[i]:
    #                                 num_atoms[i][res] = atoms
    #                             else:
    #                                 if len(atoms)<len(num_atoms[i][res]):
    #                                     num_atoms[i][res] = atoms
    #     for i in resis_to_delete:
    #         del arrays[i[0]][i[1]][i[2]]
    #     atom_lists = []
    #     for m in arrays:
    #         this_model = []
    #         for i in range(0,2):
    #             this_list_all = []
    #             this_list_bb = []
    #             for res in m[i]:
    #                 if res in num_atoms[i]:
    #                     atoms = [a.get_id() for a in m[i][res].get_list()]
    #                     ref_atoms = [at.get_id() for at in num_atoms[i][res]]
    #                     for atom in sorted(atoms):
    #                         if atom in ref_atoms:
    #                             this_list_all.append(m[i][res][atom])
    #                             if atom in ['N','CA','C']:
    #                                 this_list_bb.append(m[i][res][atom])
    #             this_model.append(this_list_all)
    #             this_model.append(this_list_bb)
    #         atom_lists.append(this_model)
    #     _ = list(num_atoms[1].keys())
    #     c = 0
    #     for m in atom_lists:
    #         c+=1
    #         for i in range(0,4):
    #             if i<2:
    #                 j=0
    #             else:
    #                 j=1
    #             if c>1:
    #                 self.number_of_residues_superposed['file{}'.format(str(c))][self.four_scores[i]] = len(num_atoms[j])
    #                 self.number_of_atoms_superposed['file{}'.format(str(c))][self.four_scores[i]] = len(m[i])
    #                 rmsd = self.calc_RMSD(atom_lists[0][i], m[i])#, TM_keys)
    #                 self.rmsds['file{}'.format(str(c))][self.four_scores[i]] = rmsd
    #             else:
    #                 self.number_of_residues_superposed['reference'][self.four_scores[i]] = len(num_atoms[j])
    #                 self.number_of_atoms_superposed['reference'][self.four_scores[i]] = len(m[i])
    #                 self.rmsds['reference'][self.four_scores[i]] = None

    # def run_RMSD(self,file1,file2):
    #     """Calculates 4 RMSD values between two GPCR pdb files.

    #     It compares the two files using sequence numbers.
    #         1. overall all atoms RMSD
    #         2. overall backbone atoms RMSD
    #         3. 7TM all atoms RMSD
    #         4. 7TM backbone atoms RMSD
    #     """
    #     parser = PDB.PDBParser(QUIET=True)
    #     pdb1 = parser.get_structure('struct1', file1)[0]
    #     pdb2 = parser.get_structure('struct2', file2)[0]
    #     pdb_array1, pdb_array2, pdb_array3, pdb_array4 = OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict()

    #     assign_gn1 = as_gn.GenericNumbering(structure=pdb1)
    #     pdb1 = assign_gn1.assign_generic_numbers()
    #     assign_gn2 = as_gn.GenericNumbering(structure=pdb2)
    #     pdb2 = assign_gn2.assign_generic_numbers()

    #     for i in pdb1:
    #         for j in pdb2:
    #             if i.get_id()==j.get_id():
    #                 chain1 = i.get_id()
    #                 chain2 = i.get_id()
    #                 break

    #     if 'chain1' not in locals():
    #         for i in pdb1.get_chains():
    #             chain1 = i.get_id()
    #             break
    #     if 'chain2' not in locals():
    #         for i in pdb2.get_chains():
    #             chain2 = i.get_id()
    #             break

    #     for residue1 in pdb1[chain1]:
    #         if residue1.get_full_id()[3][0]!=' ':
    #             continue
    #         pdb_array1[int(residue1.get_id()[1])] = residue1
    #         try:
    #             if -8.1 < residue1['CA'].get_bfactor() < 8.1:
    #                 pdb_array3[int(residue1.get_id()[1])] = residue1
    #         except:
    #             pass
    #     for residue2 in pdb2[chain2]:
    #         if residue2.get_full_id()[3][0]!=' ':
    #             continue
    #         pdb_array2[int(residue2.get_id()[1])] = residue2
    #         try:
    #             if -8.1 < residue2['CA'].get_bfactor() < 8.1:
    #                 pdb_array4[int(residue2.get_id()[1])] = residue2
    #         except:
    #             pass
    #     overall_all1, overall_all2, overall_backbone1, overall_backbone2, o_a, o_b = self.create_lists(pdb_array1, pdb_array2)
    #     TM_all1, TM_all2, TM_backbone1, TM_backbone2, t_a, t_b = self.create_lists(pdb_array3, pdb_array4)

    #     rmsd1 = self.calc_RMSD(overall_all1, overall_all2,o_a)
    #     rmsd2 = self.calc_RMSD(overall_backbone1, overall_backbone2,o_b)
    #     rmsd3 = self.calc_RMSD(TM_all1, TM_all2,t_a)
    #     rmsd4 = self.calc_RMSD(TM_backbone1, TM_backbone2,t_b)
    #     return [rmsd1, rmsd2, rmsd3, rmsd4]

    # def create_lists(self, pdb_array1, pdb_array2):
    #     """Creates the 4 atom lists needed for run_RMSD()."""
    #     overall_all1, overall_all2, overall_backbone1, overall_backbone2 = [], [], [], []
    #     keys1, keys2 = [], []
    #     for num1, res1 in pdb_array1.items():
    #         for num2, res2 in pdb_array2.items():
    #             if num1==num2 and res1.get_resname()==res2.get_resname():
    #                 for atom1 in res1:
    #                     for atom2 in res2:
    #                         if atom1.get_id()==atom2.get_id():
    #                             overall_all1.append(atom1)
    #                             overall_all2.append(atom2)
    #                             keys1.append((num1,atom1.get_id()))
    #                             if atom1.get_id() in ['N','CA','C']:
    #                                 keys2.append((num1,atom1.get_id()))
    #                                 overall_backbone1.append(atom1)
    #                                 overall_backbone2.append(atom2)
    #                             break
    #                 break
    #     return overall_all1, overall_all2, overall_backbone1, overall_backbone2, keys1, keys2
