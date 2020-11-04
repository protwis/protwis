# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 15:48:26 2015

@author: Gaspar Pandy
"""
from django.core.management.base import BaseCommand

import structure.structural_superposition as sp
import structure.assign_generic_numbers_gpcr as as_gn

import Bio.PDB as PDB
from collections import OrderedDict
import numpy as np


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('files', help='Add any number of files as arguments. First one has to be the reference file.',
                            type=str, nargs='+')
        parser.add_argument('-n', help='Specify residue sequence numbers to compare.', type=str, default=False, nargs='+')
        parser.add_argument('-r', help='Specify a range of residue sequence numbers to compare. Format: e.g. 1-300', type=str, default=False, nargs='+')
        parser.add_argument('-c', help='Specify chain ID. If not specified, the program will try to find one that matches.', type=str, default=False)
        
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
                v.run_RMSD_list(options['files'])
            else:
                v.run_RMSD_list(options['files'], force_chain=options['c'])
        else:
            if options['c']==False:
                v.run_RMSD_list(options['files'], seq_nums=seq_nums)
            else:
                v.run_RMSD_list(options['files'], seq_nums=seq_nums, force_chain=options['c'])
        self.stdout.write('\nNumber of superposed residues:\n')
        for i,j in v.number_of_residues_superposed.items():
            self.stdout.write('{}: {}'.format(i,j))
        self.stdout.write('\nNumber of superposed atoms:\n')
        for i,j in v.number_of_atoms_superposed.items():
            self.stdout.write('{}: {}'.format(i,j))
        self.stdout.write('\nRMSD scores:\n')
        for i,j in v.rmsds.items():
            self.stdout.write('{}: {}'.format(i,j))
        self.stdout.write('\n',ending='')


class Validation():
    def __init__(self):
        self.number_of_residues_superposed = OrderedDict()
        self.number_of_atoms_superposed = OrderedDict()
        self.rmsds = OrderedDict()
        self.four_scores = ['overall_all','overall_backbone','TM_all','TM_backbone']
    
    def run_RMSD_list(self, files, seq_nums=None, force_chain=None):
        ''' Calculates 4 RMSD values between a list of GPCR pdb files. It compares the files using sequence and generic
        numbers. First file in the list has to be the reference file.
            1. overall all atoms RMSD
            2. overall backbone atoms RMSD
            3. 7TM all atoms RMSD
            4. 7TM backbone atoms RMSD
        '''
        c = 0
        for f in files:
            c+=1
            if c==1:
                self.number_of_residues_superposed['reference'] = OrderedDict()
                self.number_of_atoms_superposed['reference'] = OrderedDict()
                self.rmsds['reference'] = OrderedDict()
            else:
                self.number_of_residues_superposed['file{}'.format(str(c))] = OrderedDict()
                self.number_of_atoms_superposed['file{}'.format(str(c))] = OrderedDict()
                self.rmsds['file{}'.format(str(c))] = OrderedDict()
        parser = PDB.PDBParser(QUIET=True)
        count = 0
        pdbs = []
        for f in files:
            count+=1
            pdb = parser.get_structure('struct{}'.format(count), f)[0]
            assign_gn = as_gn.GenericNumbering(pdb_file=f, sequence_parser=True)
            pdb = assign_gn.assign_generic_numbers_with_sequence_parser()
            pdbs.append(pdb)
        chains = []
        for p in pdbs:
            this = []
            for c in p.get_chains():
                this.append(c.get_id())
            chains.append(this)
        usable_chains = []
        for m in chains[1:]:
            for c in m:
                if c in chains[0]:
                    usable_chains.append(c)
        if force_chain:
            chains[0] = [force_chain]
        
        arrays = []
        model_counter = 0
        for p in pdbs:
            try:
                if pdbs.index(p)==0 and len(usable_chains)==0:
                    chain = [c.get_id() for c in pdbs[0].get_chains()][0]
                else:
                    chain = p[usable_chains[0]].get_id()
            except:
                try:
                    chain = p[' '].get_id()
                except:
                    chain = p['A'].get_id()
            if force_chain and model_counter==0:
                chain = force_chain
            pdb_array1, pdb_array2 = OrderedDict(), OrderedDict()
            for residue in p[chain]:
                if residue.get_full_id()[3][0]!=' ':
                    continue
                if seq_nums!=None and str(residue.get_id()[1]) in seq_nums:
                    pdb_array1[int(residue.get_id()[1])] = residue
                elif seq_nums==None:
                    pdb_array1[int(residue.get_id()[1])] = residue
                try:
                    if -8.1 < residue['CA'].get_bfactor() < 8.1:
                        pdb_array2[int(residue.get_id()[1])] = residue
                except:
                    pass
            arrays.append([pdb_array1,pdb_array2])
            model_counter+=1

        all_deletes, TM_deletes = [], []
        all_keep, TM_keep = [], []
        for i in range(0,2):
            for res in arrays[0][i]:
                for m in arrays[1:]:
                    if res not in m[i]:
                        if i==0:
                            all_deletes.append(res)
                        else:
                            TM_deletes.append(res)
                    else:
                        if i==0:
                            all_keep.append(res)
                        else:
                            TM_keep.append(res)
        deletes = [all_deletes, TM_deletes]
        keeps = [all_keep, TM_keep]

        num_atoms1, num_atoms2 = OrderedDict(), OrderedDict()
        num_atoms = [num_atoms1, num_atoms2]
        mismatches = []
        resis_to_delete = []
        for m_i, m in enumerate(arrays):
            for i in range(0,2):
                for res in m[i]:
                    if res in deletes[i] or res not in keeps[i]:
                        resis_to_delete.append([m_i,i,res])
                    else:
                        try:
                            if m[i][res].get_resname()!=num_atoms[i][res][0].get_parent().get_resname():
                                del num_atoms[i][res]
                                mismatches.append(res)
                            else:
                                raise Exception()
                        except:
                            if res not in mismatches:
                                atoms = []
                                for atom in m[i][res]:
                                    atoms.append(atom)
                                if res not in num_atoms[i]:
                                    num_atoms[i][res] = atoms
                                else:
                                    if len(atoms)<len(num_atoms[i][res]):
                                        num_atoms[i][res] = atoms
        for i in resis_to_delete:
            del arrays[i[0]][i[1]][i[2]]
        atom_lists = []
        for m in arrays:
            this_model = []
            for i in range(0,2):
                this_list_all = []
                this_list_bb = []
                for res in m[i]:
                    if res in num_atoms[i]:
                        atoms = [a.get_id() for a in m[i][res].get_list()]
                        ref_atoms = [at.get_id() for at in num_atoms[i][res]]
                        for atom in sorted(atoms):
                            if atom in ref_atoms:
                                this_list_all.append(m[i][res][atom])
                                if atom in ['N','CA','C']:
                                    this_list_bb.append(m[i][res][atom])
                this_model.append(this_list_all)
                this_model.append(this_list_bb)
            atom_lists.append(this_model)
        TM_keys = list(num_atoms[1].keys())        
        c = 0
        for m in atom_lists:
            c+=1
            for i in range(0,4):
                if i<2:
                    j=0
                else:
                    j=1
                if c>1:
                    self.number_of_residues_superposed['file{}'.format(str(c))][self.four_scores[i]] = len(num_atoms[j])
                    self.number_of_atoms_superposed['file{}'.format(str(c))][self.four_scores[i]] = len(m[i])
                    rmsd = self.calc_RMSD(atom_lists[0][i], m[i])#, TM_keys)
                    self.rmsds['file{}'.format(str(c))][self.four_scores[i]] = rmsd
                else:
                    self.number_of_residues_superposed['reference'][self.four_scores[i]] = len(num_atoms[j])
                    self.number_of_atoms_superposed['reference'][self.four_scores[i]] = len(m[i])
                    self.rmsds['reference'][self.four_scores[i]] = None
    
    def run_RMSD(self,file1,file2):
        ''' Calculates 4 RMSD values between two GPCR pdb files. It compares the two files using sequence numbers.
            1. overall all atoms RMSD
            2. overall backbone atoms RMSD
            3. 7TM all atoms RMSD
            4. 7TM backbone atoms RMSD
        '''
        parser = PDB.PDBParser(QUIET=True)
        pdb1 = parser.get_structure('struct1', file1)[0]
        pdb2 = parser.get_structure('struct2', file2)[0]
        pdb_array1, pdb_array2, pdb_array3, pdb_array4 = OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict()
        
        assign_gn1 = as_gn.GenericNumbering(structure=pdb1)
        pdb1 = assign_gn1.assign_generic_numbers()
        assign_gn2 = as_gn.GenericNumbering(structure=pdb2)
        pdb2 = assign_gn2.assign_generic_numbers()

        for i in pdb1:
            for j in pdb2:
                if i.get_id()==j.get_id():
                    chain1 = i.get_id()                        
                    chain2 = i.get_id()
                    break
                
        if 'chain1' not in locals():
            for i in pdb1.get_chains():
                chain1 = i.get_id()
                break
        if 'chain2' not in locals():
            for i in pdb2.get_chains():
                chain2 = i.get_id()
                break

        for residue1 in pdb1[chain1]:
            if residue1.get_full_id()[3][0]!=' ':
                continue
            pdb_array1[int(residue1.get_id()[1])] = residue1
            try:
                if -8.1 < residue1['CA'].get_bfactor() < 8.1:
                    pdb_array3[int(residue1.get_id()[1])] = residue1
            except:
                pass
        for residue2 in pdb2[chain2]:
            if residue2.get_full_id()[3][0]!=' ':
                continue
            pdb_array2[int(residue2.get_id()[1])] = residue2
            try:
                if -8.1 < residue2['CA'].get_bfactor() < 8.1:
                    pdb_array4[int(residue2.get_id()[1])] = residue2
            except:
                pass
        overall_all1, overall_all2, overall_backbone1, overall_backbone2, o_a, o_b = self.create_lists(pdb_array1, pdb_array2)
        TM_all1, TM_all2, TM_backbone1, TM_backbone2, t_a, t_b = self.create_lists(pdb_array3, pdb_array4)

        rmsd1 = self.calc_RMSD(overall_all1, overall_all2,o_a)
        rmsd2 = self.calc_RMSD(overall_backbone1, overall_backbone2,o_b)
        rmsd3 = self.calc_RMSD(TM_all1, TM_all2,t_a)
        rmsd4 = self.calc_RMSD(TM_backbone1, TM_backbone2,t_b)  
        return [rmsd1, rmsd2, rmsd3, rmsd4]
        
    def create_lists(self, pdb_array1, pdb_array2):
        ''' Creates the 4 atom lists needed for run_RMSD().
        '''
        overall_all1, overall_all2, overall_backbone1, overall_backbone2 = [], [], [], []
        keys1, keys2 = [], []
        for num1, res1 in pdb_array1.items():
            for num2, res2 in pdb_array2.items():
                if num1==num2 and res1.get_resname()==res2.get_resname():
                    for atom1 in res1:
                        for atom2 in res2:
                            if atom1.get_id()==atom2.get_id():
                                overall_all1.append(atom1)
                                overall_all2.append(atom2)
                                keys1.append((num1,atom1.get_id()))
                                if atom1.get_id() in ['N','CA','C']:
                                    keys2.append((num1,atom1.get_id()))
                                    overall_backbone1.append(atom1)
                                    overall_backbone2.append(atom2)
                                break
                    break
        return overall_all1, overall_all2, overall_backbone1, overall_backbone2, keys1, keys2

    def calc_RMSD(self, list1, list2, TM_keys=None):
        ''' Calculates RMSD between two atoms lists. The two lists have to have the same length. 
        '''
        superpose = sp.RotamerSuperpose(list1, list2, TM_keys)
        print(list1)
        print(list2)
        print(TM_keys)
        list2 = superpose.run()
        array1, array2 = np.array([0,0,0]), np.array([0,0,0])
        for a1, a2 in zip(list1, list2):
            array1 = np.vstack((array1, list(a1.get_coord())))
            array2 = np.vstack((array2, list(a2.get_coord())))
        print(array1.shape,array2.shape)
        rmsd = np.sqrt(sum(sum((array1[1:]-array2[1:])**2))/array1[1:].shape[0])
        return rmsd
        
        