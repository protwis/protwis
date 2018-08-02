from build.management.commands.base_build import Command as BaseBuild
from build.management.commands.build_homology_models_zip import Command as UploadModel
from build_gpcr.management.commands.build_homology_models import CallHomologyModeling
from django.db.models import Q
from django.conf import settings

from protein.models import Protein, ProteinConformation, ProteinAnomaly, ProteinState, ProteinSegment
from residue.models import Residue
from residue.functions import dgn, ggn
from structure.models import *
from structure.functions import HSExposureCB, PdbStateIdentifier, HomologyModelingSupportFunctions, update_template_source, ImportHomologyModel
from common.alignment import AlignedReferenceTemplate, GProteinAlignment
from common.definitions import *
from common.models import WebLink
from signprot.models import SignprotComplex
import structure.structural_superposition as sp
import structure.assign_generic_numbers_gpcr as as_gn
import structure.homology_models_tests as tests
from structure.signprot_modeling import SignprotModeling, GPCRDBParsingPDB

import Bio.PDB as PDB
from modeller import *
from modeller.automodel import *
from collections import OrderedDict
import os
import shlex
import logging
import pprint
from io import StringIO, BytesIO
import sys
import re
import zipfile
import shutil
import math
from copy import deepcopy
from datetime import datetime, date
import yaml
import traceback
import subprocess
import pprint


startTime = datetime.now()
logger = logging.getLogger('homology_modeling')
hdlr = logging.FileHandler('./logs/homology_modeling.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
logger.setLevel(logging.INFO)
structure_path = './structure/'
pir_path = os.sep.join([structure_path, 'PIR'])

build_date = date.today()
atom_num_dict = {'E':9, 'S':6, 'Y':12, 'G':4, 'A':5, 'V':7, 'M':8, 'L':8, 'I':8, 'T':7, 'F':11, 'H':10, 'K':9, 
                         'D':8, 'C':6, 'R':11, 'P':7, 'Q':9, 'N':8, 'W':14, '-':0}
gprotein_segments = ProteinSegment.objects.filter(proteinfamily='Gprotein')
gprotein_segment_slugs = [i.slug for i in gprotein_segments]

import warnings
warnings.filterwarnings("ignore")

class Command(BaseBuild):  
    help = 'Build automated chimeric GPCR homology models'
    
    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--update', help='Upload model to GPCRdb, overwrites existing entry', default=False, 
                            action='store_true')
        parser.add_argument('-r', help='''Run program for specific receptor(s) by giving UniProt common name as 
                                          argument (e.g. 5ht2a_human) or build revised crystal by giving PDB code (e.g. 4K5Y)''', 
                            default=False, type=str, nargs='+')
        parser.add_argument('-z', help='Create zip file of model directory containing all built models', default=False,
                            action='store_true')
        parser.add_argument('-c', help='Select GPCR class (A, B1, B2, C, F)', default=False)
        parser.add_argument('-x', help='Select crystal structure refinement for all crystals in the db', default=False, action='store_true')
        parser.add_argument('--purge', help='Purge all existing records', default=False, action='store_true')
        parser.add_argument('-i', help='Number of MODELLER iterations for model building', default=1, type=int)
        parser.add_argument('--test_run', action='store_true', help='Build only a test set of homology models ', default=False)
        parser.add_argument('--debug', help='Debugging mode', default=False, action='store_true')
        parser.add_argument('--state', help='Specify state in debug mode', default=False, type=str, nargs='+')
        parser.add_argument('--complex', help='Build GPCR complex', default=False, action='store_true')
        parser.add_argument('--signprot', help='Specify signaling protein with UniProt name', default=False, type=str)
        parser.add_argument('--n_c_term', help='Model N- and C-termini', default=False, action='store_true')
        parser.add_argument('--force_main_temp', help='Build model using this xtal as main template', default=False, type=str)
        
    def handle(self, *args, **options):
        self.debug = options['debug']
        if not os.path.exists('./structure/homology_models/'):
            os.mkdir('./structure/homology_models')
        if not os.path.exists('./structure/PIR/'):
            os.mkdir('./structure/PIR')
        if not os.path.exists('./static/homology_models'):
            os.mkdir('./static/homology_models')
        open('./structure/homology_models/done_models.txt','w').close()
        if options['update']:
            self.update = True
        else:
            self.update = False
        if options['complex']:
            self.complex = True
        else:
            self.complex = False
        if not options['signprot']:
            self.signprot = False
        else:
            self.signprot = options['signprot']
        self.force_main_temp = options['force_main_temp']
        self.modeller_iterations = options['i']

        excludees = SignprotComplex.objects.all().values_list('structure__protein_conformation__protein__parent__entry_name', flat=True)
        classA_receptors = Protein.objects.filter(parent__isnull=True, accession__isnull=False, species__common_name='Human', family__parent__parent__parent__name='Class A (Rhodopsin)')
        gprotein_targets = {'Gs':['gnas2_human', 'gnal_human'], 'Gi/o':['gnai1_human', 'gnai2_human','gnai3_human','gnao_human','gnat1_human','gnat2_human','gnat3_human','gnaz_human'], 
                            'Gq/11':['gnaq_human','gna11_human','gna14_human','gna15_human'], 'G12/13':['gna12_human','gna13_human']}

        c=0
        classA_receptors = Protein.objects.filter(entry_name='ackr2_human')
        for receptor in classA_receptors:
            if receptor.entry_name not in excludees:
                c+=1
                print(receptor)
                first_in_subfam = True
                for gprotein_subfam, targets in gprotein_targets.items():
                    for target in targets:
                        if first_in_subfam:
                            mod = CallHomologyModeling(receptor.entry_name, 'Active', debug=True, update=True, complex_model=True, signprot='gnat1_human')
                            mod.run()
                        ihm = ImportHomologyModel(receptor.entry_name)
                        if ihm.find_files()!=None:
                            break
                    # mod = CallHomologyModeling('receptor.entry_name', 'Active', debug=True, update=True, complex_model=True, signprot='gnat1_human')
                    # mod.run(import_receptor=True)
        print(c)

        
        






