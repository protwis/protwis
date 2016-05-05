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
        parser.add_argument('file1')
        parser.add_argument('file2')
        
    def handle(self, *args, **options):
        v = Validation()
        self.stdout.write('1. Overall all, 2. Overall backbone, 3. 7TM all, 4. 7TM backbone', ending='\n')
        self.stdout.write(v.run_RMSD(options['file1'],options['file2']), ending='')
        self.stdout.write('\n',ending='')

class Validation():
    def __init__(self):
        pass
    
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
        
        for chain1 in pdb1:
            for residue1 in chain1:
                if residue1.get_full_id()[3][0]!=' ':
                    continue
                pdb_array1[int(residue1.get_id()[1])] = residue1
                try:
                    if -8.1 < residue1['CA'].get_bfactor() < 8.1:
                        pdb_array3[int(residue1.get_id()[1])] = residue1
                except:
                    pass
        for chain2 in pdb2:
            for residue2 in chain2:
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
        return rmsd1, rmsd2, rmsd3, rmsd4
        
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
                                if atom1.get_id() in ['N','CA','C','O']:
                                    keys2.append((num1,atom1.get_id()))
                                    overall_backbone1.append(atom1)
                                    overall_backbone2.append(atom2)
                                break
                    break
        return overall_all1, overall_all2, overall_backbone1, overall_backbone2, keys1, keys2

    def calc_RMSD(self, list1, list2, keys):
        ''' Calculates RMSD between two atoms lists. The two lists have to have the same length. 
        '''
        superpose = sp.RotamerSuperpose(list1, list2)
        list2 = superpose.run()
        array1, array2 = np.array([0,0,0]), np.array([0,0,0])
        for a1, a2, k in zip(list1, list2, keys):
            array1 = np.vstack((array1, list(a1.get_coord())))
            array2 = np.vstack((array2, list(a2.get_coord())))
        rmsd = np.sqrt(sum(sum((array1[1:]-array2[1:])**2))/array1[1:].shape[0])
        return rmsd
        
        