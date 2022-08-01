from build.management.commands.base_build import Command as BaseBuild
from django.db import connection
from django.conf import settings


from protein.models import Protein, ProteinState
from structure.models import Structure, StructureModel, StructureComplexModel, StatsText, PdbData, StructureModelpLDDT
from residue.models import Residue
from common.definitions import *

import Bio.PDB as PDB
import os
import logging
import zipfile
import shutil
from datetime import datetime, date
import time
from io import StringIO


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
        parser.add_argument('-f', help='Specify file name to be uploaded to GPCRdb', default=False, type=str, nargs='+')
        parser.add_argument('-c', help='Upload only complex models to GPCRdb', default=False, action='store_true')
        parser.add_argument('--purge', help='Purge existing entries in GPCRdb', default=False, action='store_true')

    def get_structures(self,pdbname):
        if pdbname in self.cached_structures:
            return self.cached_structures[pdbname]
        else:
            if pdbname=='AF':
                return None
            else:
                s = Structure.objects.get(pdb_code__index=pdbname)
                self.cached_structures[pdbname] = s
                return s
        
    def handle(self, *args, **options):
        self.cached_structures = {}

        self.models_to_do = []

        if options['c'] and options['purge']:
            for s in StructureComplexModel.objects.all():
                s.main_template.refined = False
                s.main_template.save()
            StructureComplexModel.objects.all().delete()
        elif options['purge']:
            for s in StructureModel.objects.all():
                s.main_template.refined = False
                s.main_template.save()
            StructureModel.objects.all().delete()
        if options['c']:
            path = './structure/complex_models_zip/'
        else:
            path = './structure/homology_models_zip/'
        if not os.path.exists(path):
            os.mkdir(path)
        if options['f']:
            files = options['f']
        else:
            files = os.listdir(path)
        for f in files:
            if f.endswith('.zip'):
                modelname = f.split('.')[0]
                mod_dir = path+modelname
                if not os.path.exists(mod_dir):
                    os.mkdir(mod_dir)
                zip_mod = zipfile.ZipFile(path+f, 'r')
                zip_mod.extractall(mod_dir)
                zip_mod.close()
                self.models_to_do.append([modelname,path])

        self.processors = options['proc']
        self.prepare_input(options['proc'], self.models_to_do)

    def main_func(self, positions, iteration, count, lock):
        processor_id = round(self.processors*positions[0]/len(self.models_to_do))+1
        while count.value<len(self.models_to_do):
            with lock:
                if len(self.models_to_do)>count.value:
                    modelname,path = self.models_to_do[count.value]
                    count.value +=1
            start_import = time.time()
            start_connections = len(connection.queries)
            self.upload_to_db(modelname, path)
            mod_dir = path+modelname
            shutil.rmtree(mod_dir)
            # print("Done",modelname, time.time()-start_import,len(connection.queries)-start_connections)


    def upload_to_db(self, modelname, path):
        ''' Upload to model to StructureModel or StructureComplexModel
        '''
        name_list = modelname.split('_')
        if len(name_list)<3:
            return 0
        if name_list[3] in ['Inactive','Active','Intermediate']:
            self.complex = False
            self.revise_xtal = False
            gpcr_class = name_list[0][-1]
            gpcr_prot = '{}_{}'.format(name_list[1],name_list[2])
            state = name_list[3]
            main_structure = name_list[4]
            build_date = name_list[5]
        elif name_list[4]=='refined':
            self.complex = False
            self.revise_xtal = True
            gpcr_class = name_list[0][-1]
            gpcr_prot = name_list[3].lower()
            state = name_list[5]
            main_structure = name_list[3]
            build_date = name_list[6]
        elif name_list[5]=='refined':
            self.complex = True
            self.revise_xtal = True
            gpcr_class = name_list[0][-1]
            gpcr_prot = name_list[4].lower()
            sign_prot = '{}_{}'.format(name_list[2].split('-')[1], name_list[3])
            main_structure = name_list[4]
            build_date = name_list[6]
        else:
            self.complex = True
            self.revise_xtal = False
            gpcr_class = name_list[0][-1]
            gpcr_prot = '{}_{}'.format(name_list[1],name_list[2].split('-')[0])
            sign_prot = '{}_{}'.format(name_list[2].split('-')[1], name_list[3])
            main_structure = name_list[4]
            build_date = name_list[5]

        with open(os.sep.join([path, modelname, modelname+'.pdb']), 'r') as pdb_file:
            pdb_data = pdb_file.read()
        with open(os.sep.join([path, modelname, modelname+'.templates.csv']), 'r') as templates_file:
            templates = templates_file.readlines()
        # with open(os.sep.join([path, modelname, modelname+'.template_similarities.csv']), 'r') as sim_file:
        #     similarities = sim_file.readlines()

        if main_structure=='AF':
            stats_text = None
        else:
            stats_text = StatsText.objects.get_or_create(stats_text=''.join(templates))[0]
        pdb = PdbData.objects.get_or_create(pdb=pdb_data)[0]
        
        if self.complex:
            m_s = self.get_structures(main_structure)
            r_prot = Protein.objects.get(entry_name=gpcr_prot)
            s_prot = Protein.objects.get(entry_name=sign_prot)
            StructureComplexModel.objects.get_or_create(receptor_protein=r_prot, sign_protein=s_prot, main_template=m_s, pdb_data=pdb, version=build_date, stats_text=stats_text)
        else:
            s_state = ProteinState.objects.get(name=state)
            m_s = self.get_structures(main_structure)
            prot = Protein.objects.get(entry_name=gpcr_prot)
            sm, created = StructureModel.objects.get_or_create(protein=prot, state=s_state, main_template=m_s, pdb_data=pdb, version=build_date, stats_text=stats_text)
            if main_structure=='AF':
                try:
                    p = PDB.PDBParser().get_structure('model', os.sep.join([settings.DATA_DIR, 'structure_data', 'Alphafold', '{}_{}.pdb'.format(gpcr_prot, state.lower())]))[0]
                except:
                    print('ERROR: {}_{}.pdb is not in data folder'.format(gpcr_prot, state.lower()))
                resis = []
                for chain in p:
                    for res in chain:
                        plddt = res['CA'].get_bfactor()
                        r = StructureModelpLDDT(structure_model=sm, residue=Residue.objects.get(protein_conformation__protein=prot, sequence_number=res.get_id()[1]), pLDDT=plddt)
                        resis.append(r)
                StructureModelpLDDT.objects.bulk_create(resis)
        if self.revise_xtal:
            m_s.refined = True
            m_s.save()

