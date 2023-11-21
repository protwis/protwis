from build.management.commands.base_build import Command as BaseBuild
from django.db import connection
from django.conf import settings


from protein.models import Protein, ProteinState, ProteinConformation
from structure.models import Structure, StructureModel, StructureComplexModel, StatsText, PdbData, StructureModelpLDDT, StructureType, StructureExtraProteins
import structure.assign_generic_numbers_gpcr as as_gn
from residue.models import Residue
from common.models import WebResource, WebLink
from common.definitions import G_PROTEIN_DISPLAY_NAME as g_prot_dict
from signprot.models import SignprotComplex
from contactnetwork.cube import compute_interactions
from interaction.models import StructureLigandInteraction



import Bio.PDB as PDB
import os
import logging
import zipfile
import shutil
from datetime import datetime, date
import time


starttime = datetime.now()
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
    help = 'Upload GPCRdb structure models from zips to db'

    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('-f', help='Specify file name to be uploaded to GPCRdb', default=False, type=str, nargs='+')
        parser.add_argument('-c', help='Upload only complex models to GPCRdb', default=False, action='store_true')
        parser.add_argument('--purge', help='Purge existing entries in GPCRdb', default=False, action='store_true')
        parser.add_argument('--copy_from_data', help='Copy and rename AF models from data to upload folder', default=False, action='store_true')

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

        if options['copy_from_data']:
            af_folder = os.sep.join([settings.DATA_DIR, 'structure_data', 'AlphaFold'])
            data_files = os.listdir(af_folder)

            for f in data_files:
                try:
                    uni, species, state = f.split('.')[0].split('_')
                except ValueError:
                    uni, species, _, state = f.split('.')[0].split('_')
                if state.startswith('isoform'):
                    continue
                classname = Protein.objects.get(entry_name=uni+'_'+species).family.parent.parent.parent.name.split(' ')[1]
                new_name = 'Class{}_{}_{}_{}_AF_2022-08-16_GPCRdb'.format(classname, uni, species, state.capitalize())
                shutil.copyfile(os.sep.join([af_folder, f]), os.sep.join(['./structure/homology_models_zip/', new_name+'.pdb']))
                os.chdir('./structure/homology_models_zip/')
                assign_gn = as_gn.GenericNumbering(pdb_file=new_name+'.pdb', sequence_parser=True)
                pdb_struct = assign_gn.assign_generic_numbers_with_sequence_parser()
                ### Removing H-atoms from models
                for chain in pdb_struct:
                    for residue in chain:
                        for atom in residue.get_unpacked_list():
                            if atom.element=='H':
                                residue.detach_child(atom.get_id())
                io = PDB.PDBIO()
                io.set_structure(pdb_struct)
                io.save(new_name+'.pdb')
                with zipfile.ZipFile(new_name+'.zip', 'w') as zf:
                    zf.write(new_name+'.pdb')
                os.remove(new_name+'.pdb')
                os.chdir('../../')

        if options['c'] and options['purge']:
            for s in Structure.objects.filter(structure_type__slug='af-signprot-refined'):
                try:
                    PdbData.objects.filter(pdb=s.pdb_data.pdb).delete()
                    parent_struct = Structure.objects.get(pdb_code__index=s.pdb_code.index.split('_')[0])
                    parent_struct.refined = False
                    parent_struct.save()
                    s.stats_text.delete()
                except:
                    pass
            Structure.objects.filter(structure_type__slug='af-signprot-refined').delete()
        elif options['purge']:
            for s in StructureModel.objects.all():
                try:
                    s.pdb_data.delete()
                except:
                    pass
                try:
                    s.main_template.refined = False
                    s.main_template.save()
                except AttributeError:
                    pass
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
        _ = round(self.processors*positions[0]/len(self.models_to_do))+1
        while count.value<len(self.models_to_do):
            with lock:
                if len(self.models_to_do)>count.value:
                    modelname,path = self.models_to_do[count.value]
                    count.value +=1
            self.upload_to_db(modelname, path)
            mod_dir = path+modelname
            shutil.rmtree(mod_dir)

    def upload_to_db(self, modelname, path):
        ''' Upload to model to StructureModel or StructureComplexModel
        '''
        name_list = modelname.split('_')
        if len(name_list)<3:
            return 0
        if name_list[3] in ['Inactive','Active','Intermediate']:
            self.complex = False
            self.revise_xtal = False
            gpcr_prot = '{}_{}'.format(name_list[1],name_list[2])
            state = name_list[3]
            main_structure = name_list[4]
            build_date = name_list[5]
        elif name_list[4]=='refined':
            self.complex = False
            self.revise_xtal = True
            gpcr_prot = name_list[3].lower()
            state = name_list[5]
            main_structure = name_list[3]
            build_date = name_list[6]
        elif name_list[5]=='refined':
            self.complex = True
            self.revise_xtal = True
            gpcr_prot = name_list[4].lower()
            sign_prot = '{}_{}'.format(name_list[2].split('-')[1], name_list[3])
            main_structure = name_list[4]
            build_date = name_list[6]
        else:
            self.complex = True
            self.revise_xtal = False
            gpcr_prot = '{}_{}'.format(name_list[1],name_list[2].split('-')[0])
            sign_prot = '{}_{}'.format(name_list[2].split('-')[1], name_list[3])
            main_structure = name_list[4]
            build_date = name_list[5]

        try:
            ### Removing H atoms
            assign_gn = as_gn.GenericNumbering(pdb_file=os.sep.join([path, modelname, modelname+'.pdb']), sequence_parser=True)
            pdb_struct = assign_gn.assign_generic_numbers_with_sequence_parser()
            ### Removing H-atoms from models
            for chain in pdb_struct:
                for residue in chain:
                    for atom in residue.get_unpacked_list():
                        if atom.element=='H':
                            residue.detach_child(atom.get_id())
            io = PDB.PDBIO()
            io.set_structure(pdb_struct)
            io.save(os.sep.join([path, modelname, modelname+'.pdb']))
        except Exception as msg:
            print(modelname, msg)

        with open(os.sep.join([path, modelname, modelname+'.pdb']), 'r') as pdb_file:
            pdb_data = pdb_file.read()

        templates_file = os.sep.join([path, modelname, modelname+'.templates.csv'])
        if os.path.exists(templates_file):
            with open(templates_file, 'r') as templates_file:
                templates = templates_file.readlines()
        else:
            templates = ''

        if main_structure=='AF':
            stats_text = None
        else:
            stats_text = StatsText.objects.get_or_create(stats_text=''.join(templates))[0]
        pdb = PdbData.objects.get_or_create(pdb=pdb_data)[0]
        
        ### Alphafold refined structures
        if self.complex:
            m_s = self.get_structures(main_structure)
            r_prot = Protein.objects.get(entry_name=gpcr_prot)
            s_prot = Protein.objects.get(entry_name=sign_prot)
            # StructureComplexModel.objects.get_or_create(receptor_protein=r_prot, sign_protein=s_prot, main_template=m_s, pdb_data=pdb, version=build_date, stats_text=stats_text)
            parent_struct = Structure.objects.get(pdb_code__index=main_structure)
            protconf = ProteinConformation.objects.get(protein=parent_struct.protein_conformation.protein.parent)
            if parent_struct.structure_type.slug=='x-ray-diffraction':
                refined_type_slug = 'af-signprot-refined-xray'
                refined_type_name = 'Refined X-ray'
            elif parent_struct.structure_type.slug=='electron-microscopy':
                refined_type_slug = 'af-signprot-refined-cem'
                refined_type_name = 'Refined CEM'
            elif parent_struct.structure_type.slug=='electron-crystallography':
                refined_type_slug = 'af-signprot-refined-med'
                refined_type_name = 'Refined MED'
            signprotrefined, _ = StructureType.objects.get_or_create(slug=refined_type_slug, name=refined_type_name)
            webresource = WebResource.objects.get(slug='pdb')
            weblink, _ = WebLink.objects.get_or_create(index='{}_refined'.format(main_structure), web_resource=webresource)
            struct_obj, _ = Structure.objects.get_or_create(preferred_chain=parent_struct.preferred_chain, publication_date=build_date, pdb_data=pdb, pdb_code=weblink, build_check=True,
                                                            protein_conformation=protconf, state=parent_struct.state, structure_type=signprotrefined, author_state=parent_struct.author_state, stats_text=stats_text)
            signprot_complex, _ = SignprotComplex.objects.get_or_create(alpha=parent_struct.signprot_complex.alpha, protein=parent_struct.signprot_complex.protein, structure=struct_obj,
                                                                        beta_chain=parent_struct.signprot_complex.beta_chain, gamma_chain=parent_struct.signprot_complex.gamma_chain,
                                                                        beta_protein=parent_struct.signprot_complex.beta_protein, gamma_protein=parent_struct.signprot_complex.gamma_protein)
            struct_obj.signprot_complex = signprot_complex
            ligands = parent_struct.ligands.all()
            for l in ligands:
                try:
                    parent_sli = StructureLigandInteraction.objects.get(ligand=l, structure=parent_struct)
                except StructureLigandInteraction.MultipleObjectsReturned:
                    parent_sli = StructureLigandInteraction.objects.filter(ligand=l, structure=parent_struct)[0]
                sli, _ = StructureLigandInteraction.objects.get_or_create(pdb_reference=parent_sli.pdb_reference, annotated=True, ligand=parent_sli.ligand, ligand_role=parent_sli.ligand_role, structure=struct_obj)
                struct_obj.ligands.add(l)
            struct_obj.save()
            g_prot_dict[signprot_complex.protein.entry_name.split('_')[0].upper()]
            signprot_conf = ProteinConformation.objects.get(protein=signprot_complex.protein)
            sep, _ = StructureExtraProteins.objects.get_or_create(display_name=g_prot_dict[sign_prot.split('_')[0].upper()], note=None, chain=signprot_complex.alpha, category='G alpha', 
                                                                  wt_coverage=100, protein_conformation=signprot_conf, structure=struct_obj, wt_protein=signprot_complex.protein)
            parent_struct.refined = True
            parent_struct.save()
            try:
                compute_interactions(os.sep.join([path, modelname, modelname+'.pdb']), protein=struct_obj, signprot=signprot_complex.protein, do_complexes=True, save_to_db=True, file_input=True) # add do_complexes
            except Exception as msg:
                print('Error with interactions:', modelname, msg)
        else:
            s_state = ProteinState.objects.get(name=state)
            m_s = self.get_structures(main_structure)
            prot = Protein.objects.get(entry_name=gpcr_prot)
            sm, _ = StructureModel.objects.get_or_create(protein=prot, state=s_state, main_template=m_s, pdb_data=pdb, version=build_date, stats_text=stats_text)
            if main_structure=='AF':
                p = PDB.PDBParser().get_structure('model', os.sep.join([path, modelname, modelname+'.pdb']))[0]
                resis = []
                for chain in p:
                    for res in chain:
                        plddt = res['C'].get_bfactor()
                        try:
                            res_obj = Residue.objects.get(protein_conformation__protein=prot, sequence_number=res.get_id()[1])
                            r = StructureModelpLDDT(structure_model=sm, residue=res_obj, pLDDT=plddt)
                            resis.append(r)
                        except Residue.DoesNotExist:
                            if self.revise_xtal:
                                m_s.refined = False
                                m_s.save()
                            if sm.pdb_data:
                                sm.pdb_data.delete()
                            if sm.stats_text:
                                sm.stats_text.delete()
                            sm.delete()
                StructureModelpLDDT.objects.bulk_create(resis)
        if self.revise_xtal:
            m_s.refined = True
            m_s.save()
