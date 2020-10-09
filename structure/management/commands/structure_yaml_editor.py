# -*- coding: utf-8 -*-
"""
Created on Mon Feb 06 16:50:57 2020

@author: Gaspar Pandy
"""
from build.management.commands.base_build import Command as BaseBuild
from django.db.models import Q

from protein.models import Protein, ProteinConformation, ProteinSequenceType, ProteinSource, ProteinState
from residue.models import Residue
from structure.models import Structure, PdbData, StructureType
from structure.sequence_parser import SequenceParser
from structure.functions import PdbChainSelector, PdbStateIdentifier
from construct.functions import *
from common.models import WebResource, WebLink, Publication

import Bio.PDB as PDB
from datetime import datetime
import urllib
import re
import os
import xmltodict
import yaml
import shlex
import subprocess
import pprint


class Command(BaseBuild):

    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('-y', help="Query specific yaml(s) with PDB code", default=False, type=str, nargs='+')
        parser.add_argument('-a', help="Check for auxiliary protein. Options: BRIL, T4-Lysozyme", default=False, type=str)
        parser.add_argument('--change', help="Change specific ligand value of specific field in yaml file. Arg1 field, arg2 old value, arg3 new value", default=False, type=str, nargs='+')
        parser.add_argument('-q', help="Add quotes to field", default=False, type=str)

    def handle(self, *args, **options):
        yaml_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'structures'])
        if options['y']:
            yamls = [i + '.yaml' for i in options['y']]
        else:
            yamls = os.listdir(yaml_dir)
        for y in yamls:
            sy = StructureYaml(y)
            if options['a']:
                aux_protein = options['a']
                out = sy.check_aux_protein(aux_protein)
                print(out)
            elif options['change']:
                if len(options['change'])!=3:
                    raise AssertionError('The --change flag takes exactly 3 arguments (field, old value, new value)')
                else:
                    field, old_value, new_value = options['change']
                    sy.change_ligand_value(y, field, old_value, new_value)
            elif options['q']:
                sy.add_quotes(y, options['q'])

class StructureYaml():
    def __init__(self, yaml):
        self.yaml = yaml
        self.yaml_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'structures'])
        self.pdb_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'pdbs'])
        self.aux_protein_table = {'BRIL':'P0ABE7', 'T4-Lysozyme':'P00720'}

    def check_aux_protein(self, aux_protein):
        ''' Return a two element list about the presence of the auxiliary protein. First element is true if the auxiliary protein is in
            the DBREF of the pdb file, second element is true if this info is in the structure yaml file.
            @param aux_protein: str, BRIL or T4-Lysozyme
        '''
        with open(os.sep.join([self.yaml_dir, self.yaml]), 'r') as yf:
            yaml_dict = yaml.load(yf, Loader=yaml.FullLoader)
        with open(os.sep.join([self.pdb_dir, self.yaml.split('.')[0]+'.pdb']), 'r') as pdb_file:
            lines = pdb_file.readlines()
            aux_pdb = False
            aux_line = ''
            for l in lines:
                if l.startswith('DBREF'):
                    l_split = l.split(' ')
                    if self.aux_protein_table[aux_protein] in l_split:
                        aux_pdb = True
                        aux_line = l
                        if yaml_dict and aux_protein not in yaml_dict['auxiliary_protein']:
                            return [True, False]
                        else:
                            return [True, True]
            return [False, False]

    def change_ligand_value(self, filename, field, old_value, new_value):
        with open(os.sep.join([self.yaml_dir, filename]), 'r') as yf:
            struct_yaml = yaml.load(yf, Loader=yaml.FullLoader)
        change = False
        if type(struct_yaml['ligand'])==type([]):
            for i, ligand in enumerate(struct_yaml['ligand']):
                if ligand[field]==old_value:
                    change = True
                    struct_yaml['ligand'][i][field] = new_value
            if change:
                with open(os.sep.join([self.yaml_dir, filename]), 'w') as yf:
                    yaml.dump(struct_yaml, yf, indent=4, default_flow_style=False)
        else:
            if struct_yaml['ligand'][field]==old_value:
                struct_yaml['ligand'][field] = new_value
                with open(os.sep.join([self.yaml_dir, filename]), 'w') as yf:
                    yaml.dump(struct_yaml, yf, indent=4, default_flow_style=False)

    def add_quotes(self, filename, field):
        change = False
        with open(os.sep.join([self.yaml_dir, filename]), 'r') as yf:
            struct_yaml = yaml.load(yf, Loader=yaml.FullLoader)
            pprint.pprint(struct_yaml)
            struct_yaml[field] = str(struct_yaml[field])
            pprint.pprint(struct_yaml)
            change = True
        if change:
            with open(os.sep.join([self.yaml_dir, filename]), 'w') as yf:
                yaml.dump(struct_yaml, yf, indent=4, default_flow_style=False)
