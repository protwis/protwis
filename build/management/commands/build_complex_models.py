from build.management.commands.base_build import Command as BaseBuild
from build.management.commands.build_homology_models_zip import Command as UploadModel
from build_gpcr.management.commands.build_homology_models import CallHomologyModeling
from django.db.models import Q
from django.conf import settings

from protein.models import Protein, ProteinConformation, ProteinAnomaly, ProteinState, ProteinSegment
from residue.models import Residue
from residue.functions import dgn, ggn
from structure.models import *
from structure.functions import HSExposureCB, PdbStateIdentifier, update_template_source
from common.alignment import AlignedReferenceTemplate, GProteinAlignment
from common.definitions import *
from common.models import WebLink
from signprot.models import SignprotComplex
import structure.structural_superposition as sp
import structure.assign_generic_numbers_gpcr as as_gn
import structure.homology_models_tests as tests
from structure.signprot_modeling import SignprotModeling 
from structure.homology_modeling_functions import SignprotFunctions, GPCRDBParsingPDB, ImportHomologyModel, Remodeling

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
gprotein_segments = ProteinSegment.objects.filter(proteinfamily='Alpha')
gprotein_segment_slugs = [i.slug for i in gprotein_segments]

import warnings
warnings.filterwarnings("ignore")

class Command(BaseBuild):  
    help = 'Build automated chimeric GPCR homology models'
    
    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--update', help='Upload models to GPCRdb, overwrites existing entry', default=False, 
                            action='store_true')
        parser.add_argument('-r', help='''Run program for specific receptor(s) by giving UniProt common name as 
                                          argument (e.g. 5ht2a_human) or build revised crystal by giving PDB code (e.g. 4K5Y)''', 
                            default=False, type=str, nargs='+')
        parser.add_argument('-z', help='Create zip file of model directory containing all built models', default=False,
                            action='store_true')
        # parser.add_argument('-c', help='Select GPCR class (A, B1, B2, C, F)', default=False)
        # parser.add_argument('-x', help='Select crystal structure refinement for all crystals in the db', default=False, action='store_true')
        parser.add_argument('--purge', help='Purge all existing records', default=False, action='store_true')
        parser.add_argument('-i', help='Number of MODELLER iterations for model building', default=1, type=int)
        # parser.add_argument('--test_run', action='store_true', help='Build only a test set of homology models ', default=False)
        parser.add_argument('--debug', help='Debugging mode', default=False, action='store_true')
        parser.add_argument('--signprot', help='Specify signaling protein with UniProt name', default=False, type=str)
        # parser.add_argument('--n_c_term', help='Model N- and C-termini', default=False, action='store_true')
        parser.add_argument('--force_main_temp', help='Build model using this xtal as main template', default=False, type=str)
        
    def handle(self, *args, **options):

        # rm = Remodeling('./structure/homology_models/ClassA_taar1_human-gnao_human_6G79_2018-12-04_GPCRDB.pdb', receptor=Protein.objects.get(entry_name='taar1_human'), signprot=Protein.objects.get(entry_name='gnao_human'))
        # rm.find_clashes()
        # rm.make_pirfile()
        # rm.run()
        # rm2 = Remodeling('./structure/homology_models/ClassA_taar1_human-gnao_human_6G79_2018-12-04_GPCRDB.pdb', receptor=Protein.objects.get(entry_name='taar1_human'), signprot=Protein.objects.get(entry_name='gnao_human'))
        # rm2.find_clashes()
        # print(datetime.now() - startTime)
        # return 0


        if options['purge']:
            self.purge_complex_entries()
        self.debug = options['debug']
        if not os.path.exists('./structure/homology_models/'):
            os.mkdir('./structure/homology_models')
        if not os.path.exists('./structure/PIR/'):
            os.mkdir('./structure/PIR')
        if not os.path.exists('./static/homology_models'):
            os.mkdir('./static/homology_models')
        if not os.path.exists('./structure/complex_models_zip/'):
            os.mkdir('./structure/complex_models_zip/')
        open('./structure/homology_models/done_models.txt','w').close()

        if options['update']:
            self.update = True
        else:
            self.update = False
        self.complex = True
        if not options['signprot']:
            self.signprot = False
        else:
            self.signprot = options['signprot']
        self.force_main_temp = options['force_main_temp']
        self.modeller_iterations = options['i']
        if options['debug']:
            self.debug = True
        else:
            self.debug = False

        sf = SignprotFunctions()

        # excludees = SignprotComplex.objects.all().values_list('structure__protein_conformation__protein__parent__entry_name', flat=True)
        gprots_with_structures = sf.get_subtypes_with_templates()
        if options['r']:
            self.receptor_list = Protein.objects.filter(entry_name__in=options['r'])
        else:
            self.receptor_list = Protein.objects.filter(parent__isnull=True, accession__isnull=False, species__common_name='Human', family__parent__parent__parent__name='Class A (Rhodopsin)')
        
        print(self.receptor_list)
        subfams = sf.get_subfamilies_with_templates()
        self.gprotein_targets = sf.get_subfam_subtype_dict(subfams)

        # del self.gprotein_targets['Gi/o']

        self.receptor_list = Protein.objects.filter(entry_name__in=['drd3_human'])

        self.processors = options['proc']
        self.prepare_input(self.processors, self.receptor_list)

        # delete unzipped folders in complex_models_zip
        for i in os.listdir('./structure/complex_models_zip/'):
            if i.startswith('Class') and not i.endswith('.zip'):
                shutil.rmtree('./structure/complex_models_zip/'+i)
        

    def main_func(self, positions, itearation, count, lock):
        processor_id = round(self.processors*positions[0]/len(self.receptor_list))+1
        i = 0
        while count.value<len(self.receptor_list):
            i += 1
            with lock:
                receptor = self.receptor_list[count.value]
                logger.info('Generating complex model for  \'{}\' ... ({} out of {}) (processor:{} count:{})'.format(receptor.entry_name, count.value+1, len(self.receptor_list),processor_id,i))
                count.value +=1 

            mod_startTime = datetime.now()
            self.build_all_complex_models_for_receptor(receptor, count, i, processor_id)
            logger.info('Complex model finished for  \'{}\' ... (processor:{} count:{}) (Time: {})'.format(receptor.entry_name, processor_id,i,datetime.now() - mod_startTime))

    def build_all_complex_models_for_receptor(self, receptor, count, i, processor_id):
        first_in_subfam = True
        for gprotein_subfam, targets in self.gprotein_targets.items():
            # print(gprotein_subfam, targets)
            for target in targets:
                # Only build gnat models with opsins
                if receptor.family.parent.name!='Opsins' and target in ['gnat1_human','gnat2_human','gnat3_human']:
                    continue
                # print(receptor, target)
                import_receptor = False
                if len(SignprotComplex.objects.filter(structure__protein_conformation__protein__parent__entry_name=receptor.entry_name, protein__entry_name=target))>0:
                    continue
                else:
                    if first_in_subfam:
                        mod = CallHomologyModeling(receptor.entry_name, 'Active', debug=self.debug, update=self.update, complex_model=True, signprot=target)
                        mod.run(fast_refinement=True)
                        first_in_subfam = False
                    else:
                        ihm = ImportHomologyModel(receptor.entry_name, target)
                        if ihm.find_files()!=None:
                            mod = CallHomologyModeling(receptor.entry_name, 'Active', debug=self.debug, update=self.update, complex_model=True, signprot=target)
                            mod.run(import_receptor=True, fast_refinement=True)
                            import_receptor = True
                        else:
                            mod = CallHomologyModeling(receptor.entry_name, 'Active', debug=self.debug, update=self.update, complex_model=True, signprot=target)
                            mod.run(fast_refinement=True)


    def purge_complex_entries(self):
        try:
            StructureComplexModelSeqSim.objects.all().delete()
        except:
            self.logger.warning('StructureComplexModelSeqSim data cannot be deleted')
        try:
            StructureComplexModelStatsRotamer.objects.all().delete()
        except:
            self.logger.warning('StructureComplexModelStatsRotamer data cannot be deleted')
        try:
            StructureComplexModel.objects.all().delete()
        except:
            self.logger.warning('StructureComplexModel data cannot be deleted')
        
        






