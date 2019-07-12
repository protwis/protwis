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

import warnings
warnings.filterwarnings("ignore")

class Command(BaseBuild):  
    help = 'Build automated chimeric GPCR homology models'
    
    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--update', help='Upload models to GPCRdb, overwrites existing entry', default=False, 
                            action='store_true')
        parser.add_argument('-r', help='''Run program for specific receptor(s) by giving UniProt common name as 
                                          argument (e.g. 5ht2a_human)''', 
                            default=False, type=str, nargs='+')
        parser.add_argument('--purge', help='Purge all existing db records', default=False, action='store_true')
        parser.add_argument('--purge_zips', help='Purge all existing local model zips', default=False, action='store_true')
        parser.add_argument('--test_run', action='store_true', help='Build only one complex model', default=False)
        parser.add_argument('--debug', help='Debugging mode', default=False, action='store_true')
        parser.add_argument('--signprot', help='Specify signaling protein with UniProt name', default=False, type=str, nargs='+')
        parser.add_argument('--force_main_temp', help='Build model using this xtal as main template', default=False, type=str)
        parser.add_argument('--skip_existing', help='Skip rebuilding models already in protwis/structure/complex_models_zip/', 
                            default=False, action='store_true')
        parser.add_argument('-z', help='Create zip file of complex model directory containing all built complex models', default=False,
                            action='store_true')
        
    def handle(self, *args, **options):
        if options['purge']:
            print("Delete existing db entries")
            self.purge_complex_entries()
        if options['purge_zips']:
            print("Delete existing local zips")
            complex_zip_path = './structure/complex_models_zip/'
            if os.path.exists(complex_zip_path):
                files = os.listdir(complex_zip_path)
                for f in files:
                    try:
                        os.unlink(complex_zip_path+f)
                    except:
                        shutil.rmtree(complex_zip_path+f)

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

        self.update = options['update']
        self.complex = True
        if not options['signprot']:
            self.signprot = False
        else:
            self.signprot = options['signprot']
        self.force_main_temp = options['force_main_temp']
        self.debug = options['debug']
        self.skip_existing = options['skip_existing']
        self.existing_list = []
        if self.skip_existing:
            existing = os.listdir('./structure/complex_models_zip/')
            for f in existing:
                if f.endswith('zip'):
                    split = f.split('_')
                    self.existing_list.append(split[1]+'_'+split[2]+'_'+split[3])

        sf = SignprotFunctions()

        gprots_with_structures = sf.get_subtypes_with_templates()
        if options['r']:
            self.receptor_list = Protein.objects.filter(entry_name__in=options['r'])
        else:
            self.receptor_list = Protein.objects.filter(parent__isnull=True, accession__isnull=False, species__common_name='Human', 
                                                        family__parent__parent__parent__name='Class A (Rhodopsin)')
        
        if options['signprot']:
            self.gprotein_targets = {'custom':options['signprot']}
        else:
            subfams = sf.get_subfamilies_with_templates()
            self.gprotein_targets = sf.get_subfam_subtype_dict(subfams)

        if options['test_run']:
            break_loop = False
            for receptor in self.receptor_list:
                for i,j in self.gprotein_targets.items():
                    for target in j:
                        if len(SignprotComplex.objects.filter(structure__protein_conformation__protein__parent__entry_name=receptor.entry_name, protein__entry_name=target))==0:
                            self.receptor_list = [receptor]
                            self.gprotein_targets = {i:[target]}
                            break_loop = True
                            break
                    if break_loop: break
                if break_loop: break

        ###
        # self.receptor_list = Protein.objects.filter(entry_name='gp139_human')
        # del self.gprotein_targets['Gi/o']
        # del self.gprotein_targets['Gs']
        # self.gprotein_targets['Gq/11'] = ['gna14_human', 'gna15_human']
        ###

        s_c = 0
        for i,j in self.gprotein_targets.items():
            for s in j:
                s_c+=1

        print('receptors to model: {}'.format(len(self.receptor_list)))
        print('signaling proteins per receptor: {}'.format(s_c))

        self.processors = options['proc']
        self.prepare_input(self.processors, self.receptor_list)

        # delete unzipped folders in complex_models_zip
        for i in os.listdir('./structure/complex_models_zip/'):
            if i.startswith('Class') and not i.endswith('.zip'):
                shutil.rmtree('./structure/complex_models_zip/'+i)
        
        #create master zip for archive
        if options['z']:
            os.chdir('./structure/')
            zipf = zipfile.ZipFile('../static/homology_models/GPCRdb_complex_homology_models_{}.zip'.format(str(build_date)),'w',zipfile.ZIP_DEFLATED)
            for root, dirs, files in os.walk('complex_models_zip'):
                for f in files:
                    zipf.write(os.path.join(root, f))
            zipf.close()

    def main_func(self, positions, itearation, count, lock):
        processor_id = round(self.processors*positions[0]/len(self.receptor_list))+1
        i = 0
        while count.value<len(self.receptor_list):
            i += 1
            with lock:
                receptor = self.receptor_list[count.value]
                logger.info('Generating complex model for  \'{}\' ... ({} out of {}) (processor:{} count:{})'.format(receptor.entry_name, 
                            count.value+1, len(self.receptor_list),processor_id,i))
                count.value +=1 

            mod_startTime = datetime.now()
            self.build_all_complex_models_for_receptor(receptor, count, i, processor_id)
            logger.info('Complex model finished for  \'{}\' ... (processor:{} count:{}) (Time: {})'.format(receptor.entry_name, 
                                                                                                    processor_id,i,datetime.now() - mod_startTime))

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
                # Skip models already in './structure/complex_models_zip/'
                if receptor.entry_name+'-'+target in self.existing_list:
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
        if os.path.exists('./structure/complex_models_zip/'):
            for i in os.listdir('./structure/complex_models_zip/'):
                shutil.rmtree('./structure/complex_models_zip/'+i)
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
        