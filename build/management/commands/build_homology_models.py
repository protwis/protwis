from build.management.commands.base_build import Command as BaseBuild
from build.management.commands.build_homology_models_zip import Command as UploadModel
from django.db.models import Q
from django.conf import settings

from protein.models import Protein, ProteinConformation, ProteinAnomaly, ProteinState, ProteinSegment, ProteinFamily
from residue.models import Residue, ResidueGenericNumberEquivalent, ResidueGenericNumber
from residue.functions import dgn, ggn
from structure.models import *
from structure.functions import HSExposureCB, PdbStateIdentifier, update_template_source, StructureSeqNumOverwrite, run_residue_flip
from common.alignment import AlignedReferenceTemplate, GProteinAlignment, Alignment
from common.definitions import *
from common.models import WebLink
from signprot.models import SignprotComplex
import structure.structural_superposition as sp
import structure.assign_generic_numbers_gpcr as as_gn
import structure.homology_models_tests as tests
from structure.signprot_modeling import SignprotModeling
from structure.homology_modeling_functions import GPCRDBParsingPDB, ImportHomologyModel, Remodeling, DummyStructure
from ligand.models import LigandPeptideStructure, AssayExperiment
from interaction.models import StructureLigandInteraction
from tools.management.commands.bioactivity_receptors import Command as Bioactivity

import Bio.PDB as PDB
from modeller import *
from modeller.automodel import *
from collections import OrderedDict
import glob
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
from sklearn.cluster import AffinityPropagation, MeanShift, KMeans
from sklearn import metrics
from sklearn.datasets import make_blobs
from sklearn.decomposition import PCA, TruncatedSVD
import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
import json



startTime = datetime.now()
logger = logging.getLogger('homology_modeling')
hdlr = logging.FileHandler('./logs/homology_modeling.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.INFO)
structure_path = './structure/'
pir_path = os.sep.join([structure_path, 'PIR'])
path = "./structure/homology_models/"
zip_path = "./structure/homology_models_zip/"

build_date = date.today()
atom_num_dict = {'E':9, 'S':6, 'Y':12, 'G':4, 'A':5, 'V':7, 'M':8, 'L':8, 'I':8, 'T':7, 'F':11, 'H':10, 'K':9,
                         'D':8, 'C':6, 'R':11, 'P':7, 'Q':9, 'N':8, 'W':14, '-':0}
gprotein_segments = ProteinSegment.objects.filter(proteinfamily='Alpha')
gprotein_segment_slugs = [i.slug for i in gprotein_segments]
structures_with_num_issues = [i.split('.')[0] for i in os.listdir(os.sep.join([settings.DATA_DIR, 'structure_data', 'wt_pdb_lookup']))]


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
        parser.add_argument('--purge_zips', help='Purge all zips in homology_models_zip dir', default=False, action='store_true')
        parser.add_argument('-i', help='Number of MODELLER iterations for model building', default=1, type=int)
        parser.add_argument('--test_run', action='store_true', help='Build only a test set of homology models ', default=False)
        parser.add_argument('--debug', help='Debugging mode', default=False, action='store_true')
        parser.add_argument('--state', help='Specify state in debug mode', default=False, type=str, nargs='+')
        parser.add_argument('--complex', help='Build GPCR complex', default=False, action='store_true')
        parser.add_argument('--signprot', help='Specify signaling protein with UniProt name', default=False, type=str)
        parser.add_argument('--n_c_term', help='Model N- and C-termini', default=False, action='store_true')
        parser.add_argument('--force_main_temp', help='Build model using this xtal as main template', default=False, type=str)
        parser.add_argument('--mutations', help='Add mutations based on array key', default=False, type=str, nargs='+')
        parser.add_argument('--fast_refinement', help='Chose fastest refinement option in MODELLER', default=False, action='store_true')
        parser.add_argument('--keep_hetatoms', help='Keep hetero atoms from main template, this includes ligands', default=False, action='store_true')
        parser.add_argument('--rerun', help='Skip models with matching zip archives and only run the missing models.', default=False, action='store_true')
        parser.add_argument('--alphafold', help='Run the alphafold pipeline.', default=False, action='store_true')
        parser.add_argument('--exclude_structures', help='Exclude structures as templates', default=False, type=str, nargs='+')
        parser.add_argument('--skip_existing', help='Skip model builds present in the homology_models_zip folder', default=False, action='store_true')


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
        self.fast_refinement = options['fast_refinement']
        self.keep_hetatoms = options['keep_hetatoms']
        self.rerun = options['rerun']
        self.added_mutations = options['mutations']
        self.alphafold = options['alphafold']
        self.exclude_structures = options['exclude_structures']

        if options['mutations']:
            mutations = options['mutations']
            if len(mutations)>0 and len(mutations)%2!=0:
                print('Error: Incorrect mutation input')
            else:
                self.added_mutations = OrderedDict()
                for i in range(0,len(mutations),2):
                    self.added_mutations[mutations[i]] = mutations[i+1]
                if self.debug:
                    print('Introducing mutations at: {}'.format(self.added_mutations))

        GPCR_class_codes = {'A':'001', 'B1':'002', 'B2':'003', 'C':'004', 'D1':'005', 'F':'006', 'T':'007'}
        self.modeller_iterations = options['i']
        self.build_all = False

        # Build all
        if options['purge']:
            # if updating all, then delete existing
            print("Delete existing homology model db entries")
            StructureModel.objects.all().delete()
        if options['purge_zips']:
            print("Delete existing local homology model zips")
            hommod_zip_path = './structure/homology_models_zip/'
            if os.path.exists(hommod_zip_path):
                files = os.listdir(hommod_zip_path)
                for f in files:
                    try:
                        os.unlink(hommod_zip_path+f)
                    except:
                        shutil.rmtree(hommod_zip_path+f)

        self.custom_selection = False
        if options['r']:
            all_receptors = Protein.objects.filter(entry_name__in=options['r'])
            self.custom_selection = True
        # Only refined structures
        elif options['x']:
            structs = Structure.objects.filter(annotated=True).exclude(structure_type__slug__startswith='af-').order_by('pdb_code__index')
            all_receptors = [i.protein_conformation.protein for i in structs]
        # Build all
        elif options['c']==False:
            self.build_all = True
            all_receptors = Protein.objects.filter(parent__isnull=True, accession__isnull=False, species__common_name='Human').filter(Q(family__slug__istartswith='001') |
                                                                                                                                      Q(family__slug__istartswith='002') |
                                                                                                                                      Q(family__slug__istartswith='003') |
                                                                                                                                      Q(family__slug__istartswith='004') |
                                                                                                                                      Q(family__slug__istartswith='005') |
                                                                                                                                      Q(family__slug__istartswith='006') |
                                                                                                                                      Q(family__slug__istartswith='007')).order_by('entry_name')
            structs = Structure.objects.filter(annotated=True).exclude(structure_type__slug__startswith='af-').order_by('pdb_code__index')
            all_receptors = list(all_receptors)+[i.protein_conformation.protein for i in structs]
        elif options['c'].upper() not in GPCR_class_codes:
            raise AssertionError('Error: Incorrect class name given. Use argument -c with class name A, B1, B2, D1, C, F or T')
        # Build one class
        else:
            all_receptors = Protein.objects.filter(parent__isnull=True, accession__isnull=False, species__common_name='Human',
                                                   family__slug__istartswith=GPCR_class_codes[options['c'].upper()])
        self.receptor_list, self.receptor_list_entry_names = [],[]

        # Find proteins and states for which there is no xtal yet
        for r in all_receptors:
            if r.accession==None:
                self.receptor_list.append([r, Structure.objects.get(pdb_code__index=r.entry_name.upper()).state.name])
                continue
            for st in self.get_states_to_model(r):
                self.receptor_list.append([r, st])

        if options['skip_existing']:
            existing_models = os.listdir(zip_path)
            models_to_skip = []
            for prot, state in self.receptor_list:
                if prot.parent:
                    for f in existing_models:
                        if prot.entry_name.upper() in f and 'refined' in f:
                            models_to_skip.append([prot,state])
                            break
                else:
                    for f in existing_models:
                        if prot.entry_name in f and state in f:
                            models_to_skip.append([prot,state])
                            break
            for i in models_to_skip:
                self.receptor_list.remove(i)

        # Only specified states
        if options['state']:
            self.receptor_list = [i for i in self.receptor_list if i[1] in options['state']]

        self.receptor_list_entry_names = [i[0].entry_name for i in self.receptor_list]

        # Test run, only 5 models
        if options['test_run']:
            self.receptor_list = self.receptor_list[:5]
            self.receptor_list_entry_names = self.receptor_list_entry_names[:5]

        if options['exclude_structures'] and options['r'] and options['state']:
            self.receptor_list = [[Protein.objects.get(entry_name=options['r'][0]), options['state'][0]]]

        # Model building
        print("receptors to do",len(self.receptor_list))
        self.processors = options['proc']
        self.prepare_input(options['proc'], self.receptor_list)

        # Cleanup
        missing_models = []
        with open('./structure/homology_models/done_models.txt') as f:
            for i in f.readlines():
                if i.split('\n')[0] not in self.receptor_list_entry_names:
                    missing_models.append(i.split('\n')[0])
        if len(missing_models)==0:
            print('All models were run')
        else:
            print('Missing models:')
            print(missing_models)
#            new_args = shlex.split('/env/bin/python3 manage.py build_homology_models -r {} -p {} -s {}'.format(' '.join(missing_models, options['proc'], options['s'])))
        os.remove('./structure/homology_models/done_models.txt')

        # Make zip file for archiving
        os.chdir('./structure/')
        if options['z']==True:
            if self.complex:
                zipf = zipfile.ZipFile('../static/homology_models/GPCRdb_complex_homology_models_{}.zip'.format(str(build_date)),'w',zipfile.ZIP_DEFLATED)
            else:
                zipf = zipfile.ZipFile('../static/homology_models/GPCRdb_homology_models_{}.zip'.format(str(build_date)),'w',zipfile.ZIP_DEFLATED)
            for root, dirs, files in os.walk('homology_models'):
                for f in files:
                    if 'post' not in f:
                        zipf.write(os.path.join(root, f))
            zipf.close()
        if not self.debug:
            shutil.rmtree('homology_models')
            shutil.rmtree('PIR')

    def main_func(self, positions, iteration, count, lock):
        processor_id = round(self.processors*positions[0]/len(self.receptor_list))+1
        i = 0
        while count.value<len(self.receptor_list):
            i += 1
            with lock:
                receptor = self.receptor_list[count.value]
                count.value +=1

            # RERUN: if a model zip file already exists, skip it and move to the next
            if self.rerun:
                # Init temporary model object for checks regarding signaling protein complexes etc.
                temp_model_check = HomologyModeling(receptor[0].entry_name, receptor[1], [receptor[1]], iterations=self.modeller_iterations, complex_model=self.complex, signprot=self.signprot, debug=self.debug,
                                                  force_main_temp=self.force_main_temp, fast_refinement=self.fast_refinement, keep_hetatoms=self.keep_hetatoms, mutations=self.added_mutations, exclude_structures=self.exclude_structures)

                path = './structure/complex_models_zip/' if temp_model_check.complex else './structure/homology_models_zip/'
                # Differentiate between structure refinement and homology modeling
                if temp_model_check.revise_xtal:
                    filepath = "{}*{}_refined_*.zip".format(path, receptor[0].entry_name.upper())
                else:
                    filepath = "{}*{}_{}_*.zip".format(path, receptor[0].entry_name, receptor[1])

                # Check if model zip file exists
                if len(glob.glob(filepath)) > 0:
                    continue

            mod_startTime = datetime.now()
            logger.info('Generating model for  \'{}\' ({})... ({} out of {}) (processor:{} count:{})'.format(receptor[0].entry_name, receptor[1],count.value, len(self.receptor_list),processor_id,i))
            chm = CallHomologyModeling(receptor[0].entry_name, receptor[1], iterations=self.modeller_iterations, debug=self.debug,
                                       update=self.update, complex_model=self.complex, signprot=self.signprot, force_main_temp=self.force_main_temp, keep_hetatoms=self.keep_hetatoms, mutations=self.added_mutations, alphafold=self.alphafold, exclude_structures=self.exclude_structures)
            chm.run(fast_refinement=self.fast_refinement)
            logger.info('Model finished for  \'{}\' ({})... (processor:{} count:{}) (Time: {})'.format(receptor[0].entry_name, receptor[1],processor_id,i,datetime.now() - mod_startTime))

    def get_states_to_model(self, receptor):
        if self.force_main_temp and self.custom_selection:
            return [Structure.objects.get(protein_conformation__protein__entry_name=self.force_main_temp.lower()).state.name]
        rec_class = ProteinFamily.objects.get(name=receptor.get_protein_class())
        if rec_class.name=='Class B2 (Adhesion)':
            rec_class = ProteinFamily.objects.filter(name__in=['Class B1 (Secretin)', 'Class B2 (Adhesion)'])
            structs_in_class = Structure.objects.filter(annotated=True).exclude(structure_type__slug__startswith='af-').filter(Q(protein_conformation__protein__parent__family__slug__startswith=rec_class[0].slug) |
                                                                               Q(protein_conformation__protein__parent__family__slug__startswith=rec_class[1].slug))
        else:
            structs_in_class = Structure.objects.filter(protein_conformation__protein__parent__family__slug__startswith=rec_class.slug, annotated=True).exclude(structure_type__slug__startswith='af-')
        possible_states = structs_in_class.exclude(protein_conformation__protein__parent=receptor).exclude(state__name='Other').values_list('state__name', flat=True).distinct()
        if len(possible_states)==0:
            if rec_class.name=='Class T (Taste 2)':
                rec_class = ProteinFamily.objects.get(name='Class A (Rhodopsin)')
                structs_in_class = Structure.objects.filter(protein_conformation__protein__parent__family__slug__startswith=rec_class.slug, annotated=True).exclude(structure_type__slug__startswith='af-')
            possible_states = structs_in_class.exclude(protein_conformation__protein__parent=receptor).exclude(state__name='Other').values_list('state__name', flat=True).distinct()
        structs = structs_in_class.filter(protein_conformation__protein__parent=receptor)
        li1 = list(possible_states)
        li2 = []
        for s in structs:
            if s.state.name in li1 and s.state.name not in li2:
                li2.append(s.state.name)
        li_dif = [i for i in li1 + li2 if i not in li1 or i not in li2]
        return li_dif


class CallHomologyModeling():
    def __init__(self, receptor, state, iterations=1, debug=False, update=False, complex_model=False, signprot=False, force_main_temp=False, keep_hetatoms=False, no_remodeling=False, mutations=False, alphafold=False, exclude_structures=False):
        self.receptor = receptor
        self.state = state
        self.modeller_iterations = iterations
        self.debug = debug
        self.update = update
        self.complex = complex_model
        self.signprot = signprot
        self.force_main_temp = force_main_temp
        self.keep_hetatoms = keep_hetatoms
        self.no_remodeling = no_remodeling
        self.mutations = mutations
        self.alphafold = alphafold
        self.alphafold_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'AlphaFold'])
        self.alphafold_refined_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'AlphaFold_refined'])
        self.alphafold_multimer_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'AlphaFold_multimer'])
        self.exclude_structures = exclude_structures


    def run(self, import_receptor=False, fast_refinement=False):
        try:
            # seq_nums_overwrite_cutoff_dict = {'4PHU':2000, '4LDL':1000, '4LDO':1000, '4QKX':1000, '5JQH':1000, '5TZY':2000, '6D26':2000, '6D27':2000, '6CSY':1000}

            ##### Ignore output from that can come from BioPDB! #####
            if not self.debug:
                _stdout = sys.stdout
                sys.stdout = open(os.devnull, 'w')

            Homology_model = HomologyModeling(self.receptor, self.state, [self.state], iterations=self.modeller_iterations, complex_model=self.complex, signprot=self.signprot, debug=self.debug,
                                              force_main_temp=self.force_main_temp, fast_refinement=fast_refinement, keep_hetatoms=self.keep_hetatoms, mutations=self.mutations, alphafold=self.alphafold)

            if import_receptor:
                ihm = ImportHomologyModel(self.receptor, self.signprot)
                model, templates, similarities = ihm.find_files()
                Homology_model.alignment.reference_dict, Homology_model.alignment.template_dict, Homology_model.alignment.alignment_dict, Homology_model.main_pdb_array = ihm.parse_model(model)
                Homology_model.template_source = ihm.parse_template_source(templates)
                Homology_model.similarity_table_all = ihm.parse_similarities(similarities)
                Homology_model.disulfide_pairs = ihm.find_disulfides()
                Homology_model.trimmed_residues = []
                Homology_model.main_template_preferred_chain = 'R'
                Homology_model.main_structure = Homology_model.template_source['TM1']['1x50'][0]

                cm = SignprotModeling(Homology_model.main_structure, self.signprot, Homology_model.template_source, Homology_model.trimmed_residues, Homology_model.alignment, Homology_model.main_pdb_array, Homology_model.debug)
                cm.run()
                Homology_model.template_source = cm.template_source
                Homology_model.trimmed_residues = cm.trimmed_residues
                Homology_model.alignment = cm.a
                Homology_model.main_pdb_array = cm.main_pdb_array
                Homology_model.target_signprot = cm.target_signprot
                Homology_model.signprot_protconf = cm.signprot_protconf
                Homology_model.signprot_complex = cm.signprot_complex
            elif self.alphafold and not Homology_model.revise_xtal:
                af_path = os.sep.join([self.alphafold_data_dir, Homology_model.reference_protein.entry_name+'_'+Homology_model.state.lower()+'.pdb'])
                ihm = ImportHomologyModel(self.receptor)
                ihm.path_to_pdb = af_path
                p = PDB.PDBParser()
                model = p.get_structure('receptor', af_path)[0]['A']
                Homology_model.alignment.reference_dict, Homology_model.alignment.template_dict, Homology_model.alignment.alignment_dict, Homology_model.main_pdb_array = deepcopy(ihm.parse_model(model))
                Homology_model.disulfide_pairs = ihm.find_disulfides()
                Homology_model.trimmed_residues = []
                Homology_model.main_template_preferred_chain = 'A'
                Homology_model.main_structure = DummyStructure('A')
                alignment = AlignedReferenceTemplate()
                alignment.run_hommod_alignment(Homology_model.reference_protein, ['TM1','TM2','TM3','TM4','TM5','TM6','TM7'], ['Active','Inactive'], 'similarity')
                Homology_model.similarity_table = alignment.similarity_table
                target_residues = Residue.objects.filter(protein_conformation__protein=Homology_model.reference_protein)
                rotamer_cluster_numbers = {}
                with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'rotamer_cluster_numbers.txt']), 'r') as rf:
                    lines = rf.readlines()
                    for l in lines:
                        s = l.split('\t')
                        rotamer_cluster_numbers[s[0]] = int(s[1])

                delete_list = []
                homologs = OrderedDict()
                protconfs = {}
                # Modality definition
                modality = {'Active':['Agonist', 'Full agonist', 'Agonist peptide', 'Agonist (partial)', 'Partial agonist'],
                            'Inactive':['Antagonist', 'Antagonist (neutral/silent)', 'Inverse agonist', 'Inverse agonist (partial)']}
                # ID "wrong" modality templates
                for i, j in Homology_model.similarity_table.items():
                    sli = StructureLigandInteraction.objects.filter(structure=i)
                    protconfs[i] = ProteinConformation.objects.get(protein=i.protein_conformation.protein.parent)
                    for s in sli:
                        if s.ligand_role.name not in modality[Homology_model.state]:
                            if len(sli)==1:
                                delete_list.append(i)
                            elif list(sli).index(s)==len(sli)-1 and i not in delete_list:
                                delete_list.append(i)
                        elif s.ligand_role.name in modality[Homology_model.state]:
                            if type(Homology_model.similarity_table[i])!=type([]):
                                Homology_model.similarity_table[i] = [Homology_model.similarity_table[i], modality[Homology_model.state].index(s.ligand_role.name)]
                            elif modality[Homology_model.state].index(s.ligand_role.name)<Homology_model.similarity_table[i][1]:
                                Homology_model.similarity_table[i][1] = modality[Homology_model.state].index(s.ligand_role.name)
                    if len(sli)==0:
                        delete_list.append(i)
                    if self.exclude_structures and i.pdb_code.index in self.exclude_structures:
                        delete_list.append(i)
                # Remove "wrong" templates from template list
                for i in delete_list:
                    del Homology_model.similarity_table[i]
                # Reorder based on seqsim, resolution modality
                resorted_keys = sorted(Homology_model.similarity_table.items(), key=lambda x: (-x[1][0],x[1][1],x[0].resolution))
                new_dict = OrderedDict()
                for r in resorted_keys:
                    new_dict[r[0]] = r[1][0]
                Homology_model.similarity_table = new_dict

                # ID most homologuous receptors
                for i, j in Homology_model.similarity_table.items():
                    if i.protein_conformation.protein.parent not in homologs:
                        homologs[i.protein_conformation.protein.parent] = [i]
                    elif i not in homologs[i.protein_conformation.protein.parent]:
                        homologs[i.protein_conformation.protein.parent].append(i)
                # ID interacting residues
                interacting_residues = OrderedDict()
                res_seqnums = []
                for seg_lab, seg in Homology_model.alignment.alignment_dict.items():
                    interacting_residues[seg_lab] = []
                    for gn, res in seg.items():
                        if res in ['G','A']:#['F','Y','W','R','E','K','N','H','Q','D']:
                            continue
                        if 'x' in gn:
                            if gn in ['45x50','3x25']:
                                continue
                            c = 0
                            found_in_receptor = False
                            receptor_count = 0
                            for receptor, structures in homologs.items():
                                for s in structures:
                                    for l in s.ligands.all():
                                        try:
                                            r = dgn(gn, protconfs[s])
                                            f = Fragment.objects.filter(structure=s, ligand=l, residue__display_generic_number__label=r)
                                        except Residue.DoesNotExist:
                                            continue
                                        if len(f)>0:
                                            c+=1
                                            found_in_receptor = True
                                            break
                                    if found_in_receptor:
                                        break
                                if c==3:
                                    interacting_residues[seg_lab].append(gn)
                                    res_seqnums.append(Residue.objects.get(display_generic_number__label=dgn(gn, Homology_model.prot_conf), protein_conformation=Homology_model.prot_conf).sequence_number)
                                    break
                                receptor_count+=1
                                if receptor_count==5:
                                    break

                Homology_model.alignment.reference_dict = deepcopy(Homology_model.alignment.reference_dict)
                Homology_model.alignment.template_dict = deepcopy(Homology_model.alignment.template_dict)

                # Check coverage per template
                coverage = OrderedDict()
                gn_per_struct = OrderedDict()
                struct_per_gn = OrderedDict()
                for seg,j in interacting_residues.items():
                    for gn in j:
                        struct_per_gn[gn] = []
                        for struct, sim in Homology_model.similarity_table.items():
                            if struct not in coverage:
                                coverage[struct] = [sim, 0]
                                gn_per_struct[struct] = []
                            try:
                                res = Fragment.objects.filter(structure=struct, residue__display_generic_number__label=dgn(gn, struct.protein_conformation), residue__amino_acid=Homology_model.alignment.reference_dict[seg][gn])
                                if len(res)>0:
                                    rot = Homology_model.right_rotamer_select(Rotamer.objects.filter(residue=res[0].residue))
                                    if not rot.missing_atoms:
                                        coverage[struct][1]+=1
                                        gn_per_struct[struct].append(gn)
                                        struct_per_gn[gn].append(struct)
                            except Residue.DoesNotExist:
                                continue

                bd = Bioactivity()
                bioactivity_data = bd.handle()
                bioactivity_structures = []
                if Homology_model.reference_protein.family.slug in bioactivity_data[Homology_model.state.lower()]:
                    if self.debug:
                        print(bioactivity_data[Homology_model.state.lower()][Homology_model.reference_protein.family.slug])
                    for i,j in bioactivity_data[Homology_model.state.lower()][Homology_model.reference_protein.family.slug].items():
                        bioactivity_structures+=[Structure.objects.get(pdb_code__index=j[0][0])]
                else:
                    if self.debug:
                        print('No Bioactivity data available')

                print(datetime.now() - startTime)
                resorted_keys = sorted(coverage.items(), key=lambda x: (-x[1][1],-x[1][0],x[0].resolution))
                new_dict = OrderedDict()
                for r in resorted_keys:
                    new_dict[r[0]] = r[1][0]
                Homology_model.similarity_table = new_dict
                prot_conf = ProteinConformation.objects.get(protein=Homology_model.reference_protein)

                ### Frequency benchmark TEMP MOD
                bioactivity_structures = []

                ### 7TM confidence score
                confidence_7TM = []
                for seg, res in Homology_model.main_pdb_array.items():
                    if seg.startswith('TM'):
                        for gn, atoms in res.items():
                            confidence_7TM.append(atoms[0].get_bfactor())

                alt_templates = OrderedDict()

                ### Pocket alignment
                pocket_alignment = Alignment()
                pocket_alignment.load_reference_protein(Homology_model.reference_protein)
                pocket_alignment.load_proteins([i.protein_conformation.protein for i in list(Homology_model.similarity_table.keys())])
                pocket_segments = []
                new_segment_dict = OrderedDict([('Custom', [])])
                new_generic_numbers = OrderedDict([('Custom', OrderedDict())])
                for gn in struct_per_gn:
                    new_segment_dict['Custom'].append(gn)
                    new_generic_numbers['Custom'][gn] = []
                pocket_alignment.segments = new_segment_dict
                pocket_alignment.generic_numbers['gpcrdba'] = new_generic_numbers
                pocket_alignment.use_residue_groups = True
                pocket_alignment.build_alignment()
                pocket_alignment.calculate_similarity(normalized=True)
                for p in pocket_alignment.proteins[1:]:
                    s = Structure.objects.get(protein_conformation__protein=p.protein)
                    Homology_model.similarity_table[s] = int(p.identity)

                resorted_keys = sorted(Homology_model.similarity_table.items(), key=lambda x: (-x[1],x[0].resolution))
                new_dict = OrderedDict()
                for r in resorted_keys:
                    new_dict[r[0]] = r[1]
                Homology_model.similarity_table = new_dict
                for gn, structs in struct_per_gn.items():
                    res = target_residues.get(display_generic_number__label=dgn(gn, Homology_model.prot_conf))
                    for s, i in Homology_model.similarity_table.items():
                        if s in structs:
                            alt_templates[gn] = [res.amino_acid, res.sequence_number, s]
                            break

                ### Clustering
                # for gn, structs in struct_per_gn.items():
                #     fix_st = False
                #     if len(bioactivity_structures)>0:
                #         fix_st = bioactivity_structures[0]
                #     # fix_st = Structure.objects.get(pdb_code__index='4IAR')
                #     if fix_st in struct_per_gn[gn]:
                #         structs = [fix_st, fix_st]
                #         fix_st = True
                #     ### Bioactivity benchmark TEMP MOD
                #     # else:
                #     #     continue
                #     ###
                #     if len(structs)==0:
                #         continue
                #     data = np.array([0,0])
                #     target_residue = Residue.objects.get(protein_conformation=prot_conf, display_generic_number__label=dgn(gn, prot_conf))
                #     print(gn, structs)
                #     first = True
                #     ref_atoms = []
                #     if len(structs)>1:
                #         for s in structs:
                #             this_struct_data = np.array([0,0,0])
                #             chi_angles = np.array([])
                #             rot_obj = Rotamer.objects.filter(structure=s, residue__display_generic_number__label=dgn(gn, s.protein_conformation))
                #             rot_obj = Homology_model.right_rotamer_select(rot_obj)
                #             if rot_obj.missing_atoms:
                #                 continue
                #             rot = PDB.PDBParser().get_structure('rot', StringIO(rot_obj.pdbdata.pdb))
                #             rot.atom_to_internal_coordinates()

                #             for res in rot.get_residues():
                #                 for i in range(1,6):
                #                     chi = res.internal_coord.get_angle('chi{}'.format(i))
                #                     # print(gn,s,res,chi,i)
                #                     if chi:
                #                         chi_angles = np.append(chi_angles, chi)
                #                 if first and len(chi_angles)>0:
                #                     ref_atoms = [a for a in res]
                #                     sorted_atoms = ref_atoms
                #                     first = False
                #                     if len(chi_angles)==1:
                #                         data = np.array([0]*2)
                #                     else:
                #                         data = np.array([0]*len(chi_angles))
                #                     print(ref_atoms[0].get_parent())
                #                     print(ref_atoms)
                #                     print(chi_angles)
                #                 else:
                #                     sup = sp.RotamerSuperpose(sorted(ref_atoms), sorted([a for a in res]))
                #                     sorted_atoms = sup.run()
                #                 for a in sorted_atoms:
                #                     if a.get_id() in ['N','CA','C','O']:
                #                         continue
                #                     this_struct_data = np.vstack((this_struct_data, list(a.get_coord())))
                #             # pprint.pprint(this_struct_data[1:])

                #             if len(chi_angles)>2:
                #                 if len(chi_angles)==5:
                #                     array_ext = np.array([0,0,0,0,0])
                #                 elif len(chi_angles)==4:
                #                     array_ext = np.array([0,0,0,0])
                #                 elif len(chi_angles)==3:
                #                     array_ext = np.array([0,0,0])
                #                 # tsvd = TruncatedSVD(n_components=2, random_state=1)
                #                 # chi_angles = np.vstack((chi_angles, array_ext))
                #                 # tsvd_out = tsvd.fit_transform(chi_angles)
                #                 # chi_angles = tsvd_out[0]
                #             elif len(chi_angles)==1:
                #                 chi_angles = np.append(chi_angles, 0)
                #             if len(chi_angles)>0:
                #                 data = np.vstack((data, chi_angles))

                #         data = data[1:]
                #         pprint.pprint(data)

                #         # Affinity Propagation
                #         af = AffinityPropagation(preference=-50, random_state=0).fit(data)
                #         cluster_centers_indices = af.cluster_centers_indices_
                #         labels = af.labels_

                #         if len(cluster_centers_indices)>rotamer_cluster_numbers[PDB.Polypeptide.three_to_one(ref_atoms[0].get_parent().get_resname())]:
                #             # KMeans
                #             kmeans = KMeans(n_clusters=rotamer_cluster_numbers[PDB.Polypeptide.three_to_one(ref_atoms[0].get_parent().get_resname())], random_state=0).fit(data)
                #             labels = kmeans.labels_
                #             cluster_centers = kmeans.cluster_centers_
                #             cluster_centers_indices = []
                #             print('KMeans')
                #             print(cluster_centers)
                #             cluster_centers_indices = []
                #             for c in cluster_centers:
                #                 closest = None
                #                 closest_i = None
                #                 min_val = 10000
                #                 for i,j in enumerate(data):
                #                     dist = np.linalg.norm(c-j)
                #                     if dist<min_val:
                #                         closest = j
                #                         min_val = dist
                #                         closest_i = i
                #                 cluster_centers_indices.append(closest_i)

                #         n_clusters_ = len(cluster_centers_indices)
                #         if n_clusters_==0:
                #             ms = MeanShift(bin_seeding=True).fit(data)
                #             labels = ms.labels_
                #             cluster_centers = ms.cluster_centers_
                #             print('MeanShift', cluster_centers)
                #             cluster_centers_indices = []
                #             for c in cluster_centers:
                #                 closest = None
                #                 closest_i = None
                #                 min_val = 10000
                #                 for i,j in enumerate(data):
                #                     dist = np.linalg.norm(c-j)
                #                     if dist<min_val:
                #                         closest = j
                #                         min_val = dist
                #                         closest_i = i
                #                 cluster_centers_indices.append(closest_i)

                #         print(cluster_centers_indices)
                #         print(labels)
                #         print(n_clusters_)
                #         cluster_repr_struct = structs[cluster_centers_indices[max(set(labels), key = list(labels).count)]]
                #         clusters_by_receptor = OrderedDict()
                #         for i, l in enumerate(labels):
                #             if l not in clusters_by_receptor:
                #                 clusters_by_receptor[l] = [structs[i].protein_conformation.protein.parent]
                #             else:
                #                 rec = structs[i].protein_conformation.protein.parent
                #                 if rec not in clusters_by_receptor[l]:
                #                     clusters_by_receptor[l].append(rec)
                #         print(clusters_by_receptor)

                #         ### If largest cluster only has 1 receptor, skip
                #         largest_cluster_label = max(clusters_by_receptor, key=lambda k: len(clusters_by_receptor[k]))
                #         if len(clusters_by_receptor[largest_cluster_label])==1:
                #             if not fix_st:
                #                 continue

                #         largest_clusters = [largest_cluster_label]
                #         for l, r in clusters_by_receptor.items():
                #             if l!=largest_cluster_label and len(clusters_by_receptor[l])==len(clusters_by_receptor[largest_cluster_label]):
                #                 largest_clusters.append(l)
                #         print('LARGEST clusters', largest_clusters)
                #         if len(largest_clusters)>1:
                #             ## When receptor number is equal, select cluster with most structures
                #             # struct_counts = OrderedDict()
                #             # for l in largest_clusters:
                #             #     struct_counts[l] = list(labels).count(l)
                #             #     print(l, list(labels).count(l))
                #             # counts = [l for l in struct_counts if struct_counts[l]==max(struct_counts.values())]
                #             # print(counts, max(struct_counts.values()))
                #             # if len(counts)==1:
                #             #     largest_cluster_label = counts[0]
                #             # else:
                #             #     print('same number of receptors and structures in class!')
                #             # cluster_repr_struct = structs[cluster_centers_indices[largest_cluster_label]]

                #             ## When receptor number is equal, select cluster representative with higher similarity
                #             sim_struct = {}
                #             for l in largest_clusters:
                #                 s = structs[cluster_centers_indices[l]]
                #                 sim_struct[s] = Homology_model.similarity_table[s]
                #             cluster_repr_struct = sorted(sim_struct.items(), key=lambda k: -k[1])[0][0]
                #         else:
                #             cluster_repr_struct = structs[cluster_centers_indices[largest_cluster_label]]

                #         print('CLUSTER Repr:', cluster_repr_struct)
                #         alt_templates[gn] = [target_residue.amino_acid, target_residue.sequence_number, cluster_repr_struct]
                #         # plt.close("all")
                #         # plt.figure(1)
                #         # plt.clf()

                #         # colors = cycle("bgrcmykbgrcmykbgrcmykbgrcmyk")
                #         # for k, col in zip(range(n_clusters_), colors):
                #         #     class_members = labels == k
                #         #     cluster_center = data[cluster_centers_indices[k]]
                #         #     plt.plot(data[class_members, 0], data[class_members, 1], col + ".")
                #         #     plt.plot(
                #         #         cluster_center[0],
                #         #         cluster_center[1],
                #         #         "o",
                #         #         markerfacecolor=col,
                #         #         markeredgecolor="k",
                #         #         markersize=14,
                #         #     )
                #         #     for x in data[class_members]:
                #         #         plt.plot([cluster_center[0], x[0]], [cluster_center[1], x[1]], col)

                #         # plt.title("Estimated number of clusters: %d" % n_clusters_)
                #         # plt.show()
                #         # plt.savefig('plot.png')
                if self.debug:
                    print('ALT TEMPLATES')
                    pprint.pprint(alt_templates)

                used_rotamer_templates = OrderedDict()
                confidence_scores = []
                for gn, struct in alt_templates.items():
                    struct = struct[2]
                    rot_obj = Rotamer.objects.filter(structure=struct, residue__display_generic_number__label=dgn(gn, struct.protein_conformation))
                    rot_obj = Homology_model.right_rotamer_select(rot_obj)
                    rot = PDB.PDBParser().get_structure('rot', StringIO(rot_obj.pdbdata.pdb))
                    seg = rot_obj.residue.protein_segment.slug
                    struct_atoms = Homology_model.main_pdb_array[seg][gn.replace('x','.')]
                    confidence_scores.append(struct_atoms[0].get_bfactor())
                    alt_atoms = deepcopy([a for a in rot.get_atoms()])
                    sup = sp.RotamerSuperpose(sorted(struct_atoms), sorted(alt_atoms))
                    new_atoms = sup.run()

                    alt_rmsd = deepcopy(sup.rmsd)
                    if alt_atoms[0].get_parent().get_resname() in ['TYR','PHE','ARG','ASP','GLU']:
                        new_alt_atoms = run_residue_flip(alt_atoms)
                        superpose2 = sp.RotamerSuperpose(sorted(struct_atoms), sorted(new_alt_atoms))
                        new_atoms2 = superpose2.run()
                        if alt_rmsd>superpose2.rmsd:
                            alt_rmsd = superpose2.rmsd
                            new_atoms = new_atoms2
                    # if alt_rmsd>0.2:
                    Homology_model.main_pdb_array[seg][gn.replace('x','.')] = new_atoms
                    Homology_model.template_source = update_template_source(Homology_model.template_source,[gn],struct,seg,just_rot=True)
                    used_rotamer_templates[gn] = [struct, struct_atoms[0].get_bfactor(), Homology_model.similarity_table[struct]]

                Homology_model.statistics.add_info('possible_rotamer_templates', alt_templates)
                Homology_model.statistics.add_info('used_rotamer_templates', used_rotamer_templates)

                # Label positions to change
                # c = 0
                # for seg_lab, seg in interacting_residues.items():
                #     for gn in seg:
                #         Homology_model.alignment.alignment_dict[seg_lab][gn] = '.'
                #         Homology_model.alignment.template_dict[seg_lab][gn] = 'G'
                #         c+=1
                # print('Number of positions to change: ',c)
                # Homology_model.statistics.info_dict['pdb_db_inconsistencies'] = []

                # h = Homology_model.run_non_conserved_switcher(Homology_model.main_pdb_array, Homology_model.alignment.reference_dict, Homology_model.alignment.template_dict, Homology_model.alignment.alignment_dict)
                if self.debug:
                    pprint.pprint(Homology_model.statistics)
            ### Refined alphafold structure model
            elif self.alphafold and Homology_model.revise_xtal:
                ### structure
                Homology_model.main_structure = Structure.objects.get(protein_conformation__protein=Homology_model.reference_protein)
                if Homology_model.complex:
                    target_residues = Residue.objects.filter(protein_conformation__protein__entry_name__in=[Homology_model.reference_protein.parent.entry_name, Homology_model.signprot]).order_by('protein_conformation__protein__family__slug', 'sequence_number')
                    alpha_protobj = Protein.objects.get(entry_name=Homology_model.reference_protein.entry_name+'_a')
                    Homology_model.target_signprot = alpha_protobj.parent
                    Homology_model.signprot_protconf = ProteinConformation.objects.get(protein=Homology_model.target_signprot)
                    Homology_model.signprot_complex = SignprotComplex.objects.get(structure=Homology_model.main_structure)
                    template_residues = Residue.objects.filter(protein_conformation__protein__in=[Homology_model.reference_protein, alpha_protobj], protein_segment__isnull=False).order_by('protein_conformation__protein__family__slug', 'sequence_number')
                else:
                    target_residues = Residue.objects.filter(protein_conformation__protein=Homology_model.reference_protein.parent)
                    template_residues = Residue.objects.filter(protein_conformation__protein=Homology_model.reference_protein, protein_segment__isnull=False)

                wt_lookup_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'wt_pdb_lookup', Homology_model.reference_entry_name+'.json'])
                wt_lookup_needed = False
                Homology_model.main_structure = Structure.objects.get(pdb_code__index=Homology_model.reference_entry_name.upper())
                Homology_model.main_template_preferred_chain = Homology_model.main_structure.preferred_chain
                alignment = AlignedReferenceTemplate()
                ### AF
                if Homology_model.complex:
                    complex_name = Homology_model.reference_protein.parent.entry_name+'-'+Homology_model.signprot
                    af_path =  os.sep.join([self.alphafold_multimer_data_dir, complex_name, complex_name+'.pdb'])
                    if not os.path.exists(af_path):
                        af_path = os.sep.join([self.alphafold_multimer_data_dir, complex_name+'_gbb1_human_gbg2_human', complex_name+'_gbb1_human_gbg2_human.pdb'])
                        if not os.path.exists(af_path):
                            return 0
                else:
                    af_path = os.sep.join([self.alphafold_refined_data_dir, Homology_model.main_structure.pdb_code.index+'_refined.pdb'])
                
                ihm = ImportHomologyModel(Homology_model.reference_protein.parent)
                ihm.path_to_pdb = af_path
                p = PDB.PDBParser()
                model_receptor = p.get_structure('receptor', af_path)[0]['A']
                af_reference_dict, af_template_dict, af_alignment_dict, af_main_pdb_array = deepcopy(ihm.parse_model(model_receptor))

                if Homology_model.complex:
                    signprot = Protein.objects.get(entry_name=Homology_model.signprot)
                    model_signprot = p.get_structure('signprot',  af_path)[0]['B']
                    ihm2 = ImportHomologyModel(signprot, 'Alpha')
                    spaf_reference_dict, spaf_template_dict, spaf_alignment_dict, spaf_main_pdb_array = deepcopy(ihm2.parse_model(model_signprot))
                    
                    for i, j, k, l in zip(spaf_reference_dict, spaf_template_dict, spaf_alignment_dict, spaf_main_pdb_array):
                        af_reference_dict[i] = spaf_reference_dict[i]
                        af_template_dict[j] = spaf_template_dict[j]
                        af_alignment_dict[k] = spaf_alignment_dict[k]
                        af_main_pdb_array[l] = spaf_main_pdb_array[l]

                print('pLDDT')
                for seg, resis in af_main_pdb_array.items():
                    for gn, atoms in resis.items():
                        print(gn, atoms[0].get_parent(), atoms[0].get_bfactor())
                # Homology_model.disulfide_pairs = ihm.find_disulfides()
                Homology_model.trimmed_residues = []
                parse = GPCRDBParsingPDB()

                ### get wt_pdb_lookup data when there are incorrect sequence numbers in structure
                if os.path.exists(wt_lookup_file):
                    ssno = StructureSeqNumOverwrite(Homology_model.main_structure)
                    wt_lookup_needed = True

                ### build reference_dict
                reference_dict = OrderedDict()
                for r in target_residues:
                    if r.protein_segment.slug not in reference_dict:
                        reference_dict[r.protein_segment.slug] = OrderedDict()
                    if r.display_generic_number:
                        if r.protein_segment.proteinfamily=='GPCR':
                            reference_dict[r.protein_segment.slug][ggn(r.display_generic_number.label)] = r
                        else:
                            reference_dict[r.protein_segment.slug][r.display_generic_number.label] = r
                    else:
                        reference_dict[r.protein_segment.slug][str(r.sequence_number)] = r

                segment_labels = [i for i in list(reference_dict.keys()) if len(reference_dict[i])>0]

                ### build template_dict
                template_dict = OrderedDict()
                distorted_residues = OrderedDict()
                for seglab, seg in reference_dict.items():
                    template_dict[seglab] = OrderedDict()
                    distorted_residues[seglab] = []
                    for gn, res in seg.items():
                        gn = str(gn)
                        if wt_lookup_needed and res.protein_segment.proteinfamily=='GPCR':
                            if res.sequence_number in ssno.wt_pdb_table:
                                seqnum = ssno.wt_pdb_table[res.sequence_number]
                            elif res.sequence_number in ssno.pdb_wt_table:
                                seqnum = 0
                            else:
                                seqnum = res.sequence_number
                        elif Homology_model.complex and res.protein_segment.proteinfamily=='Alpha':
                            try:
                                seqnum = template_residues.get(display_generic_number=res.display_generic_number).sequence_number
                            except Residue.DoesNotExist:
                                seqnum = 0
                        else:
                            seqnum = res.sequence_number
                        try:
                            temp_res = template_residues.get(sequence_number=seqnum, protein_segment__proteinfamily=res.protein_segment.proteinfamily)
                            ### Checking annotation
                            if res.display_generic_number and not temp_res.display_generic_number:
                                distorted_residues[res.protein_segment.slug].append(ggn(res.display_generic_number.label))
                        except Residue.DoesNotExist:
                            temp_res = '-'
                        if temp_res!='-':
                            template_dict[seglab][gn] = temp_res
                        else:
                            template_dict[seglab][gn] = '-'

                ### Custom fixes
                if Homology_model.main_structure.pdb_code.index in ['2YCW','2YCZ']:
                    distorted_residues['TM6']+=['6x22','6x23']

                ### build alignment_dict
                alignment_dict = OrderedDict()
                for seglab, seg in template_dict.items():
                    alignment_dict[seglab] = OrderedDict()
                    for gn, res in seg.items():
                        gn = str(gn)
                        if template_dict[seglab][gn]=='-':
                            alignment_dict[seglab][gn] = '-'
                        elif reference_dict[seglab][gn].amino_acid!=template_dict[seglab][gn].amino_acid:
                            alignment_dict[seglab][gn] = '.'
                        else:
                            alignment_dict[seglab][gn] = res

                ### build main_pdb_array
                main_pdb_array = OrderedDict()
                for seglab, seg in template_dict.items():
                    main_pdb_array[seglab] = OrderedDict()
                    for gn, res in seg.items():
                        mpa_gn = str(gn).replace('x','.')
                        if res=='-':
                            main_pdb_array[seglab][mpa_gn] = 'x'
                        else:
                            rot_obj = Homology_model.right_rotamer_select(Rotamer.objects.filter(residue=res))
                            rot = PDB.PDBParser().get_structure('rot', StringIO(rot_obj.pdbdata.pdb))[0]
                            if rot_obj.missing_atoms:
                                alignment_dict[seglab][gn] = '.'
                            for chain in rot:
                                for res in chain:
                                    atoms = [a for a in res if a.get_name()[0]!='H']
                            main_pdb_array[seglab][mpa_gn] = atoms

                ### reformat dictionaries
                for ref_seg, temp_seg, ali_seg in zip(reference_dict, template_dict, alignment_dict):
                    for ref_res, temp_res, ali_res in zip(reference_dict[ref_seg], template_dict[temp_seg], alignment_dict[ali_seg]):
                        if type(reference_dict[ref_seg][ref_res])!=type(''):
                            reference_dict[ref_seg][ref_res] = reference_dict[ref_seg][ref_res].amino_acid
                        if type(template_dict[temp_seg][temp_res])!=type(''):
                            template_dict[temp_seg][temp_res] = template_dict[temp_seg][temp_res].amino_acid
                        if type(alignment_dict[ali_seg][ali_res])!=type(''):
                            alignment_dict[ali_seg][ali_res] = alignment_dict[ali_seg][ali_res].amino_acid

                ### build template source
                template_source = OrderedDict()
                for seglab, seg in alignment_dict.items():
                    template_source[seglab] = OrderedDict()
                    for gn, res in seg.items():
                        if res!='-':
                            template_source[seglab][gn] = [Homology_model.main_structure, Homology_model.main_structure]
                        else:
                            template_source[seglab][gn] = [None, None]

                ### Identify missing sections
                missing_sections = []
                new_start, new_end = None, None
                section = []
                prev_seg = None
                for seg, resis in template_source.items():
                    if seg=='G.HN':
                        break
                    for gn, res in resis.items():
                        if res[0]!=Homology_model.main_structure:
                            if not new_start:
                                new_start = [seg, gn]
                                new_end = [seg, gn]
                                section = [[seg, gn]]
                            else:
                                new_end = [seg, gn]
                                section.append([seg, gn])
                        else:
                            new_start = None
                            if section not in missing_sections:
                                missing_sections.append(section)
                if section not in missing_sections:
                    missing_sections.append(section)

                if Homology_model.complex:
                    new_start, new_end = None, None
                    section = []
                    prev_seg = None
                    for seg, resis in template_source.items():
                        if '.' not in seg:
                            continue
                        for gn, res in resis.items():
                            if res[0]!=Homology_model.main_structure:
                                if not new_start:
                                    new_start = [seg, gn]
                                    new_end = [seg, gn]
                                    section = [[seg, gn]]
                                else:
                                    new_end = [seg, gn]
                                    section.append([seg, gn])
                            else:
                                new_start = None
                                if section not in missing_sections:
                                    missing_sections.append(section)
                    if section not in missing_sections:
                        missing_sections.append(section)

                    ### Chimera check
                    for section in missing_sections:
                        H5_resis = [s[0] for s in section if s[0]=='G.H5']
                        if len(H5_resis)>20:
                            if self.debug:
                                print('ERROR: incorrect chimera mapping for {}'.format(Homology_model.reference_protein))
                                print(H5_resis, len(H5_resis))
                            logger.error('Incorrect chimera mapping for {}'.format(Homology_model.reference_protein))
                            raise AssertionError

                    ### Shorten ICL3
                    for i in reference_dict:
                        if i.startswith('ICL3'):
                            label = i
                            break

                    if len(reference_dict[label])>14:
                        delete_ts, delete_r, delete_t, delete_a, delete_m = set(),set(),set(),set(),set()
                        chain_break = False
                        icl3_c = 0
                        keys = list(template_source['ICL3'].keys())
                        length = len(template_dict[label])
                        ref_prot = Homology_model.reference_protein
                        for r_s,t_s,a_s,ar_s in zip(reference_dict[label],template_dict[label],
                                                    alignment_dict[label],main_pdb_array[label]):
                            icl3_c+=1
                            if 7<icl3_c<length-6:
                                if chain_break==False:
                                    reference_dict[label][r_s] = '/'
                                    template_dict[label][t_s] = '/'
                                    alignment_dict[label][a_s] = '/'
                                    main_pdb_array[label][ar_s] = '/'
                                    delete_ts.add(('ICL3',keys[icl3_c-1]))
                                    chain_break = True
                                else:
                                    delete_r.add((label,r_s))
                                    delete_t.add((label,t_s))
                                    delete_a.add((label,a_s))
                                    delete_m.add((label,ar_s))
                                    delete_ts.add(('ICL3',keys[icl3_c-1]))
                        Homology_model.icl3_delete['A'] = []
                        for i,ii in delete_ts:
                            Homology_model.icl3_delete['A'].append(int(ii))
                            del template_source[i][ii]
                        for i,ii in delete_r:
                            del reference_dict[i][ii]
                        for i,ii in delete_t:
                            del template_dict[i][ii]
                        for i,ii in delete_a:
                            del alignment_dict[i][ii]
                        for i,ii in delete_m:
                            del main_pdb_array[i][ii]

                if self.debug:
                    print('Missing sections')
                    pprint.pprint(missing_sections)

                ### Swap sections
                seg_dict = {}
                for i in missing_sections:
                    if len(i)==0:
                        continue
                    if (i[0][0]=='N-term' and i[0][1]=='1') or ('N-term' not in reference_dict and i[0][0]=='TM1' and int(i[0][1].split('x')[1])<45) or (i[0][0]=='G.HN' and i[0][1]=='G.HN.01'):
                        i, plus_gns, last_seg = parse.find_after_4_gns(i, reference_dict, distorted_residues, segment_labels)
                        orig_residues = parse.fetch_residues_from_array(main_pdb_array[last_seg], plus_gns)
                        alt_residues_end = parse.fetch_residues_from_array(af_main_pdb_array[last_seg], plus_gns)

                        swapped_dict = OrderedDict()
                        for seg, pos in i:
                            if seg not in swapped_dict:
                                swapped_dict[seg] = [pos]
                            else:
                                swapped_dict[seg].append(pos)
                            template_source[seg][pos] = [None, None]
                            seg_dict[pos] = seg

                        alt_residues = OrderedDict()
                        for seg, posis in swapped_dict.items():
                            alt_mid = parse.fetch_residues_from_array(af_main_pdb_array[seg], posis)
                            alt_residues = parse.add_two_ordereddict(alt_residues, alt_mid)
                        alt_residues = parse.add_two_ordereddict(alt_residues, alt_residues_end)

                        superpose = sp.OneSidedSuperpose(orig_residues,alt_residues,4,0)
                        new_residues = superpose.run()

                        key_list = list(new_residues.keys())[:-4]
                        Homology_model.trimmed_residues.append(key_list[-1])
                        Homology_model.trimmed_residues.append(list(new_residues.keys())[-4])

                    elif (i[-1][0]=='C-term' and i[-1][1]==list(template_source['C-term'].keys())[-1]) or (i[-1][0]=='G.H5' and i[-1][1]==list(template_source['G.H5'].keys())[-1]):
                        i, minus_gns, first_seg = parse.find_before_4_gns(i, reference_dict, distorted_residues, segment_labels)
                        orig_residues = parse.fetch_residues_from_array(main_pdb_array[first_seg], minus_gns)
                        alt_residues_start = parse.fetch_residues_from_array(af_main_pdb_array[first_seg], minus_gns)

                        swapped_dict = OrderedDict()
                        for seg, pos in i:
                            if seg not in swapped_dict:
                                swapped_dict[seg] = [pos]
                            else:
                                swapped_dict[seg].append(pos)
                            seg_dict[pos] = seg

                        alt_residues = OrderedDict()
                        for seg, posis in swapped_dict.items():
                            alt_mid = parse.fetch_residues_from_array(af_main_pdb_array[seg], posis)
                            alt_residues = parse.add_two_ordereddict(alt_residues, alt_mid)
                        alt_residues = parse.add_two_ordereddict(alt_residues_start, alt_residues)

                        superpose = sp.OneSidedSuperpose(orig_residues,alt_residues,4,1)
                        new_residues = superpose.run()

                        key_list = list(new_residues.keys())[4:]
                        Homology_model.trimmed_residues.append(key_list[0])
                        Homology_model.trimmed_residues.append(list(new_residues.keys())[3])

                    else:
                        ### find right before and after 4 gns
                        i, minus_gns, first_seg = parse.find_before_4_gns(i, reference_dict, distorted_residues, segment_labels)
                        i, plus_gns, last_seg = parse.find_after_4_gns(i, reference_dict, distorted_residues, segment_labels)

                        ### Adjusting swappable section for truncated ICL3
                        if len(Homology_model.icl3_delete)>0:
                            icl3_positions = []
                            for e, pos in enumerate(i):
                                if pos[0]=='ICL3':
                                    icl3_positions.append(e)
                            if len(icl3_positions)>0:
                                start = icl3_positions[7]
                                end = icl3_positions[-7]
                                i = i[:start]+i[end:]

                        orig_residues1 = parse.fetch_residues_from_array(main_pdb_array[first_seg], minus_gns)
                        orig_residues2 = parse.fetch_residues_from_array(main_pdb_array[last_seg], plus_gns)
                        orig_residues = parse.add_two_ordereddict(orig_residues1, orig_residues2)

                        alt_residues1 = parse.fetch_residues_from_array(af_main_pdb_array[first_seg], minus_gns)
                        alt_residues2 = parse.fetch_residues_from_array(af_main_pdb_array[last_seg], plus_gns)

                        swapped_dict = OrderedDict()
                        for seg, pos in i:
                            if seg not in swapped_dict:
                                swapped_dict[seg] = [pos]
                            else:
                                swapped_dict[seg].append(pos)
                            seg_dict[pos] = seg

                        for seg, posis in swapped_dict.items():
                            alt_mid = parse.fetch_residues_from_array(af_main_pdb_array[seg], posis)
                            alt_residues1 = parse.add_two_ordereddict(alt_residues1, alt_mid)
                        alt_residues = parse.add_two_ordereddict(alt_residues1, alt_residues2)

                        superpose = sp.LoopSuperpose(orig_residues, alt_residues)
                        new_residues = superpose.run()

                        key_list = list(new_residues.keys())[4:-4]
                        Homology_model.trimmed_residues.append(key_list[0])
                        Homology_model.trimmed_residues.append(key_list[-1])
                        Homology_model.trimmed_residues.append(list(new_residues.keys())[3])
                        Homology_model.trimmed_residues.append(list(new_residues.keys())[-4])

                    for key in key_list:
                        if len(key.split('.'))>2:
                            main_pdb_array[seg_dict[key]][key] = new_residues[key]
                            template_dict[seg_dict[key]][key] = PDB.Polypeptide.three_to_one(new_residues[key][0].get_parent().get_resname())
                            alignment_dict[seg_dict[key]][key] = PDB.Polypeptide.three_to_one(new_residues[key][0].get_parent().get_resname())
                        else:
                            main_pdb_array[seg_dict[key.replace('.','x')]][key] = new_residues[key]
                            template_dict[seg_dict[key.replace('.','x')]][key.replace('.','x')] = PDB.Polypeptide.three_to_one(new_residues[key][0].get_parent().get_resname())
                            alignment_dict[seg_dict[key.replace('.','x')]][key.replace('.','x')] = PDB.Polypeptide.three_to_one(new_residues[key][0].get_parent().get_resname())

                ### Rotamer swap (missing atoms and non-wt)
                for seglab, seg in alignment_dict.items():
                    for gn, res in seg.items():
                        if res=='.' and gn in af_template_dict[seglab]:
                            mpa_gn = gn.replace('x','.')
                            orig_res = main_pdb_array[seglab][mpa_gn]
                            alt_res = af_main_pdb_array[seglab][mpa_gn]
                            superpose = sp.RotamerSuperpose(orig_res, alt_res)
                            new_atoms = superpose.run()
                            main_pdb_array[seglab][mpa_gn] = new_atoms
                            new_aa = PDB.Polypeptide.three_to_one(new_atoms[0].get_parent().get_resname())
                            template_dict[seglab][gn] = new_aa
                            alignment_dict[seglab][gn] = new_aa
                            template_source[seglab][gn][1] = None

                if Homology_model.complex:
                    # Add Beta and Gamma chains
                    p = PDB.PDBParser(QUIET=True).get_structure('structure', StringIO(Homology_model.main_structure.pdb_data.pdb))[0]
                    beta = p[Homology_model.signprot_complex.beta_chain]
                    gamma = p[Homology_model.signprot_complex.gamma_chain]
                    reference_dict['Beta'] = OrderedDict()
                    template_dict['Beta'] = OrderedDict()
                    alignment_dict['Beta'] = OrderedDict()
                    main_pdb_array['Beta'] = OrderedDict()
                    template_source['Beta'] = OrderedDict()
                    reference_dict['Gamma'] = OrderedDict()
                    template_dict['Gamma'] = OrderedDict()
                    alignment_dict['Gamma'] = OrderedDict()
                    main_pdb_array['Gamma'] = OrderedDict()
                    template_source['Gamma'] = OrderedDict()
                    for b_res in beta:
                        key = str(b_res.get_id()[1])
                        # remove waters
                        if b_res.get_resname()=='HOH':
                            continue
                        reference_dict['Beta'][key] = PDB.Polypeptide.three_to_one(b_res.get_resname())
                        template_dict['Beta'][key] = PDB.Polypeptide.three_to_one(b_res.get_resname())
                        alignment_dict['Beta'][key] = PDB.Polypeptide.three_to_one(b_res.get_resname())
                        atoms = [atom for atom in b_res]
                        main_pdb_array['Beta'][key] = atoms
                        template_source['Beta'][key] = [Homology_model.main_structure, Homology_model.main_structure]
                    for g_res in gamma:
                        key = str(g_res.get_id()[1])
                        # remove waters
                        if g_res.get_resname()=='HOH':
                            continue
                        reference_dict['Gamma'][key] = PDB.Polypeptide.three_to_one(g_res.get_resname())
                        template_dict['Gamma'][key] = PDB.Polypeptide.three_to_one(g_res.get_resname())
                        alignment_dict['Gamma'][key] = PDB.Polypeptide.three_to_one(g_res.get_resname())
                        atoms = [atom for atom in g_res]
                        main_pdb_array['Gamma'][key] = atoms
                        template_source['Gamma'][key] = [Homology_model.main_structure, Homology_model.main_structure]

                # C-term truncation
                if 'C-term' in reference_dict and len(reference_dict['C-term'])>10:
                    keys = list(reference_dict['C-term'].keys())[::-1]
                    for i, key in enumerate(keys):
                        del reference_dict['C-term'][key]
                        del template_dict['C-term'][key]
                        del alignment_dict['C-term'][key]
                        del main_pdb_array['C-term'][key]
                        del template_source['C-term'][key]
                        if i==len(keys)-11:
                            break

                Homology_model.alignment.reference_dict = reference_dict
                Homology_model.alignment.template_dict = template_dict
                Homology_model.alignment.alignment_dict = alignment_dict
                Homology_model.main_pdb_array = main_pdb_array
                Homology_model.template_source = template_source
            else:
                alignment = Homology_model.run_alignment([self.state])
                Homology_model.build_homology_model(alignment)
                print('icl3', Homology_model.icl3_delete)

            Homology_model.build_homology_model_second_part()
            formatted_model = Homology_model.format_final_model()
            if formatted_model==None:
                raise ValueError('Incorrect assignment of generic numbers in {} {}'.format(self.receptor, self.state))

            # # if Homology_model.changes_on_db:
            ssno = StructureSeqNumOverwrite(Homology_model.main_structure)
            if len(ssno.pdb_wt_table)>0 and not self.alphafold:
                ssno.seq_num_overwrite('wt')

            # Run clash and break test
            p = PDB.PDBParser()
            post_model = p.get_structure('model','./structure/homology_models/{}.pdb'.format(Homology_model.modelname))
            if Homology_model.signprot:
                hse = HSExposureCB(post_model, radius=11, check_chain_breaks=True, check_knots=True, receptor=self.receptor, signprot=Homology_model.signprot, check_hetatoms=True)
                # Run remodeling
                if len(hse.remodel_resis)>0 and not self.no_remodeling and not Homology_model.revise_xtal:
                    rm = Remodeling('./structure/homology_models/{}.pdb'.format(Homology_model.modelname), gaps=hse.remodel_resis, receptor=self.receptor, signprot=Homology_model.signprot,
                                    icl3_delete=Homology_model.icl3_delete)
                    rm.make_pirfile()
                    rm.run()
                    logger.info('Remodeled {} {} at {}'.format(self.receptor, Homology_model.signprot, hse.remodel_resis))
                    with open('./structure/homology_models/{}.pdb'.format(Homology_model.modelname), 'r') as remodeled_pdb:
                        formatted_model = remodeled_pdb.read()
                    formatted_model = Homology_model.format_final_model()
            else:
                hse = HSExposureCB(post_model, radius=11, check_chain_breaks=True, receptor=self.receptor, check_hetatoms=True)

            # Remove not interacting HETATM residues
            if self.debug:
                print('HETRESIS to remove: {}'.format(hse.hetresis_to_remove))
            if len(hse.hetresis_to_remove)>0:
                post_model2 = p.get_structure('model','./structure/homology_models/{}.pdb'.format(Homology_model.modelname))
                for chain in post_model2[0]:
                    for hetres in hse.hetresis_to_remove:
                        try:
                            chain.detach_child(hetres.id)
                            logger.info('{} detached from model due to no interaction'.format(hetres))
                        except KeyError:
                            print('Cannot detach {} from model'.format(hetres))
                            logger.info('Cannot detach {} from model'.format(hetres))
                io = PDB.PDBIO()
                io.set_structure(post_model2)
                io.save('./structure/homology_models/{}.pdb'.format(Homology_model.modelname))
                logger.info('HETRESIS removed: {}'.format(hse.hetresis_to_remove))
                formatted_model = Homology_model.format_final_model()

            # Check for residue shifts in model
            residue_shift = False
            db_res = ''
            for chain in post_model[0]:
                if Homology_model.complex and chain.id!='R':
                    continue
                for res in chain:
                    try:
                        if Homology_model.revise_xtal==False:
                            db_res = Residue.objects.get(protein_conformation__protein=Homology_model.reference_protein, sequence_number=int(res.get_id()[1]))
                        else:
                            db_res = Residue.objects.get(protein_conformation__protein=Homology_model.reference_protein.parent, sequence_number=int(res.get_id()[1]))
                        # Skip amino acid HETRESIS
                        if res.get_id()[0].startswith('H_'):
                            continue
                        if PDB.Polypeptide.three_to_one(res.get_resname())!=db_res.amino_acid:
                            residue_shift = True
                            break
                    except:
                        pass

            # if residue_shift==True:
            #     if self.debug:
            #         print('Residue shift in model {} at {}'.format(Homology_model.reference_entry_name, db_res))
            #     logger.warning('Residue shift in model {} at {}'.format(Homology_model.reference_entry_name, db_res))
            #     raise ValueError('Error: Residue shift in model {} at {}'.format(Homology_model.reference_entry_name, db_res))

            # Check for clashes in model
            if len(hse.clash_pairs)>0:
                if self.debug:
                    print('Remaining clashes in {}:'.format(Homology_model.reference_entry_name))
                for i in hse.clash_pairs:
                    if self.debug:
                        print(i)
                    if i[0][1]==i[1][1]-1 or i[0][1]==i[1][1]:
                        hse.clash_pairs.remove(i)
                logger.warning('Remaining clashes in {}\n{}'.format(Homology_model.reference_entry_name,hse.clash_pairs))

            # Check for chain breaks in model
            if len(hse.chain_breaks)>0:
                if self.debug:
                    print('Chain breaks in {}:'.format(Homology_model.reference_entry_name))
                    for j in hse.chain_breaks:
                        print(j)
                logger.warning('Chain breaks in {}\n{}'.format(Homology_model.reference_entry_name,hse.chain_breaks))


            ##### Resume output #####
            if not self.debug:
                sys.stdout = _stdout
                sys.stdout.close()

            # Output as zip
            if Homology_model.complex:
                path = './structure/complex_models_zip/'
            else:
                path = './structure/homology_models_zip/'
            if not os.path.exists(path):
                os.mkdir(path)
            out_stream = BytesIO()
            zipf = zipfile.ZipFile('{}{}.zip'.format(path, Homology_model.modelname),'w',zipfile.ZIP_DEFLATED)
            zipf.writestr(Homology_model.modelname+'.pdb', formatted_model)
            with open('./structure/homology_models/'+Homology_model.modelname+'.templates.csv','r') as templates_csv:
                tc = templates_csv.read()
            zipf.writestr(Homology_model.modelname+'.templates.csv', tc)
            with open('./structure/homology_models/'+Homology_model.modelname+'.template_similarities.csv','r') as temp_sim_csv:
                tsc = temp_sim_csv.read()
            zipf.writestr(Homology_model.modelname+'.template_similarities.csv', tsc)
            zipf.close()


            # Upload to db
            if self.update and not residue_shift:
                # Handling local vs live builds
                manage_call = "./manage.py"
                if os.path.exists("/env/bin/python3"):
                    manage_call = "/env/bin/python3 manage.py"

                if Homology_model.complex:
                    new_args = shlex.split((manage_call + ' build_homology_models_zip -c -f {}.zip').format(Homology_model.modelname))
                else:
                    new_args = shlex.split((manage_call + ' build_homology_models_zip -f {}.zip').format(Homology_model.modelname))

                subprocess.call(new_args)

                # logger.info('{} ({}) homology model uploaded to db'.format(Homology_model.reference_entry_name,self.state))
                if self.debug:
                    print('{} homology model uploaded to db'.format(Homology_model.reference_entry_name))

            with open('./structure/homology_models/done_models.txt','a') as f:
                f.write(self.receptor+'\n')

        except Exception as msg:
            try:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                if self.debug:
                    print('Error on line {}: Failed to build model {} (main structure: {})\n{}'.format(exc_tb.tb_lineno, self.receptor,
                                                                                            Homology_model.main_structure,msg))
                    print(''.join(traceback.format_tb(exc_tb)))
                logger.error('Failed to build model {} {}\n    {}'.format(self.receptor, self.state, msg))
                # Return main template structure residue table sequence numbering to original
                ssno = StructureSeqNumOverwrite(Homology_model.main_structure)
                if len(ssno.pdb_wt_table)>0 and not self.alphafold:
                # # if Homology_model.changes_on_db:
                    ssno = StructureSeqNumOverwrite(Homology_model.main_structure)
                    ssno.seq_num_overwrite('wt')
                    self.logger.info('Structure {} residue table sequence number overwrite wt to pdb'.format(structure))
                t = tests.HomologyModelsTests()
                if 'Number of residues in the alignment and  pdb files are different' in str(msg):
                    t.pdb_alignment_mismatch(Homology_model.alignment, Homology_model.main_pdb_array,
                                             Homology_model.main_structure)
                    if Homology_model.complex:
                        t.pdb_pir_mismatch(os.sep.join([structure_path, 'homology_models', '{}_{}_post.pdb'.format(Homology_model.reference_entry_name, Homology_model.target_signprot.entry_name)]),
                                           os.sep.join([pir_path, '{}_{}.pir'.format(Homology_model.reference_protein.accession, Homology_model.target_signprot.entry_name)]))
                    else:
                        t.pdb_pir_mismatch(os.sep.join([structure_path, 'homology_models', '{}_{}_post.pdb'.format(Homology_model.reference_entry_name, Homology_model.state)]),
                                           os.sep.join([pir_path, '{}_{}.pir'.format(Homology_model.reference_protein.accession, Homology_model.state)]))
                elif 'No such residue:' in str(msg):
                    if self.debug:
                        t.pdb_pir_mismatch(Homology_model.main_pdb_array, Homology_model.model_sequence)
                with open('./structure/homology_models/done_models.txt','a') as f:
                    f.write(self.receptor+'\n')
            except:
                try:
                    Protein.objects.get(entry_name=self.receptor)
                except:
                    logger.error('Invalid receptor name: {}'.format(self.receptor))
                    print('Invalid receptor name: {}'.format(self.receptor))


class HomologyModeling(object):
    ''' Class to build homology models for GPCRs.

        @param reference_entry_name: str, protein entry name \nco
        @param state: str, endogenous ligand state of reference \n
        @param query_states: list, list of endogenous ligand states to be applied for template search \n
        @param iterations: int, number of MODELLER iterations
    '''
    segment_coding = {1:'TM1',2:'TM2',3:'TM3',4:'TM4',5:'TM5',6:'TM6',7:'TM7',8:'H8', 12:'ICL1', 23:'ECL1', 34:'ICL2',
                      45:'ECL2'}

    def __init__(self, reference_entry_name, state, query_states, iterations=1, complex_model=False, signprot=False, debug=False,
                 force_main_temp=False, fast_refinement=False, keep_hetatoms=False, mutations=False, alphafold=False):
        self.debug = debug
        self.complex = complex_model
        self.modelname = ''
        self.signprot = signprot
        self.target_signprot = None
        self.force_main_temp = force_main_temp
        self.keep_hetatoms = keep_hetatoms
        self.fast_refinement = fast_refinement
        self.version = build_date
        self.reference_entry_name = reference_entry_name.lower()
        self.state = state
        self.query_states = query_states
        self.modeller_iterations = iterations
        self.statistics = CreateStatistics(self.reference_entry_name)
        self.reference_protein = Protein.objects.get(entry_name=self.reference_entry_name)
        self.reference_class = self.reference_protein.family.parent.parent.parent
        self.segments = []
        self.similarity_table = OrderedDict()
        self.similarity_table_all = OrderedDict()
        self.similarity_table_other_states = OrderedDict()
        self.main_structure = None
        self.signprot_complex = None
        self.main_template_preferred_chain = ''
        self.loop_template_table = OrderedDict()
        self.loops = OrderedDict()
        self.icl3_delete = OrderedDict()
        self.changes_on_db = False
        self.alphafold = alphafold
        if len(self.reference_entry_name)==4:
            self.prot_conf = ProteinConformation.objects.get(protein=self.reference_protein.parent)
            self.uniprot_id = self.reference_protein.parent.accession
            self.revise_xtal = True
            if self.reference_entry_name.upper()=='5LWE':
                raise AssertionError('Skipping 5LWE refined structure model')
            try:
                # Added beta protein check to avoid main tempalte switch for structures that only have the alpha subunit (5G53, 6FUF) - FIX later
                sc = SignprotComplex.objects.get(structure__pdb_code__index=self.reference_entry_name.upper(), beta_protein__isnull=False)
                self.complex = True
                self.signprot = sc.protein.entry_name
            except SignprotComplex.DoesNotExist:
                if self.debug:
                    print('Warning: {} not in SignprotComplex, creating non-complex refined model'.format(self.reference_entry_name))
        else:
            self.prot_conf = ProteinConformation.objects.get(protein=self.reference_protein)
            self.uniprot_id = self.reference_protein.accession
            self.revise_xtal = False
        class_tree = {'001':'A', '002':'B1', '003':'B2', '004':'C', '005':'D1', '006':'F', '007':'T'}
        self.class_name = 'Class'+class_tree[Protein.objects.get(entry_name=self.reference_entry_name).family.parent.slug[:3]]
        self.statistics.add_info('uniprot_id',self.uniprot_id)
        self.statistics.add_info('state',self.state)
        self.template_source = OrderedDict()
        self.helix_end_mods = None
        self.alignment = OrderedDict()
        self.main_pdb_array = OrderedDict()
        self.disulfide_pairs = []
        self.trimmed_residues = []
        self.peptide_ligands = []
        for r in Residue.objects.filter(protein_conformation=self.prot_conf):
            if r.protein_segment.slug not in self.template_source:
                self.template_source[r.protein_segment.slug] = OrderedDict()
            try:
                self.template_source[r.protein_segment.slug][ggn(r.display_generic_number.label)] = [None,None]
            except:
                self.template_source[r.protein_segment.slug][str(r.sequence_number)] = [None,None]
        self.mutations = mutations

    def __repr__(self):
        return "<Hommod: {}, {}>".format(self.reference_entry_name, self.state)

    def right_rotamer_select(self, rotamer):
        ''' Filter out compound rotamers.
        '''
        if len(rotamer)>1:
            for i in rotamer:
                if i.pdbdata.pdb.startswith('COMPND')==False:
                    rotamer = i
                    break
        else:
            rotamer=rotamer[0]
        return rotamer

    def pdb_line_whitespace(self, pos_list, i, group3, whitespace):
        if len(str(pos_list[i]))-len(group3)==0:
            whitespace = whitespace*' '
        elif len(str(pos_list[i]))-len(group3)==1:
            whitespace = (whitespace-1)*' '
        elif len(str(pos_list[i]))-len(group3)==2:
            whitespace = (whitespace-2)*' '
        elif len(str(pos_list[i]))-len(group3)==-1:
            whitespace = (whitespace+1)*' '
        elif len(str(pos_list[i]))-len(group3)==-2:
            whitespace = (whitespace+2)*' '
        elif len(str(pos_list[i]))-len(group3)==-3:
            whitespace = (whitespace+3)*' '
        else:
            whitespace = (whitespace-3)*' '
        return whitespace

    def format_final_model(self):
        ''' Do final formatting on homology model pdb file. Adds REMARK line, correct residue numbering and
            class-specific generic numbers. Returns the pdb in string format.
        '''
        pos_list = []
        if self.complex:
            first_signprot_res = False
        gpcr_segments = [i.slug for i in ProteinSegment.objects.filter(proteinfamily='GPCR') if not i.name.startswith("ECD")] #### FIX ME, for class B structures there will be issues because of same name of ECD H1 and G alpha H1
        for seg in self.template_source:
            for num in self.template_source[seg]:
                try:
                    if seg in gpcr_segments:
                        num = str(Residue.objects.get(protein_conformation=self.prot_conf,
                                                      display_generic_number__label=dgn(num,self.prot_conf)).sequence_number)
                    else:
                        raise AssertionError
                except:
                    if self.complex and seg not in gpcr_segments:
                        try:
                            num = str(Residue.objects.get(protein_conformation=self.signprot_protconf,
                                                          display_generic_number__label=num).sequence_number)
                            if not first_signprot_res:
                                first_signprot_res = num
                        except Exception as msg:
                            pass
                pos_list.append(num)
        if self.complex:
            first_beta_res = list(self.template_source['Beta'].keys())[0]
            first_gamma_res = list(self.template_source['Gamma'].keys())[0]

        if self.complex:
            # first_signprot_i = [(i,j) for i,j in enumerate(pos_list) if j==first_signprot_res][-1][0]
            prev_p = int(pos_list[0])
            sp_first_indeces = []
            for i,p in enumerate(pos_list):
                if int(p)<prev_p:
                    sp_first_indeces.append(i)
                prev_p = int(p)
        i, pep_i = 0, 0
        path = './structure/homology_models/'

        with open (path+self.modelname+'.pdb', 'r+') as f:
            pdblines = f.readlines()
            out_list = []
            prev_num, prev_chain = None, None
            first_hetatm = False
            water_count = 0
            first_signprot_res_found = False
            first_ligand_pep_res = False
            atom_num = 1
            atom_num_offset = []
            pref_chain = str(self.main_structure.preferred_chain)

            if self.debug:
                if self.complex:
                    print('sp_first_indeces:', sp_first_indeces)
                print('pos_list:\n', pos_list)
            for line in pdblines:
                try:
                    if prev_num==None:
                        pdb_re = re.search('(ATOM[A-Z\s\d]{13}\S{3})([\sA-Z]+)(\d+)([A-Z\s\d.-]{49,53})',line)
                        prev_num = int(pdb_re.group(3))
                        prev_chain = pdb_re.group(2).strip()
                    pdb_re = re.search('(ATOM[A-Z\s\d]{13}\S{3})([\sA-Z]+)(\d+)([A-Z\s\d.-]{49,53})',line)
                    if int(pdb_re.group(3))>prev_num or (int(pdb_re.group(3))<prev_num and prev_chain!=pdb_re.group(2).strip()):
                        i+=1
                        prev_num = int(pdb_re.group(3))
                        prev_chain = pdb_re.group(2).strip()
                    whitespace = self.pdb_line_whitespace(pos_list, i, pdb_re.group(3), len(pdb_re.group(2)))
                    whitespace = ' '+pref_chain+whitespace[2:]
                    group1 = pdb_re.group(1)
                    if 'OXT' in group1:
                        if not self.complex:
                            atom_num_offset.append(int(pos_list[i]))
                        else:
                            if i<sp_first_indeces[0]:
                                atom_num_offset.append(int(pos_list[i])) ### only add OXT induced atom num offset from Receptor, not from sign proteins
                    if self.complex:
                        if self.revise_xtal:
                            gpcr_chain, alpha_chain, beta_chain, gamma_chain = self.signprot_complex.structure.preferred_chain, self.signprot_complex.alpha, self.signprot_complex.beta_chain, self.signprot_complex.gamma_chain
                        else:
                            gpcr_chain, alpha_chain, beta_chain, gamma_chain = 'R','A','B','G'
                        if i<sp_first_indeces[0]:
                            if len(whitespace)==2:
                                whitespace = whitespace[0]+gpcr_chain
                            else:
                                whitespace = whitespace[0]+gpcr_chain+whitespace[2:]
                        elif sp_first_indeces[0]<=i<sp_first_indeces[1]:
                            if len(whitespace)==2:
                                whitespace = whitespace[0]+alpha_chain
                            else:
                                whitespace = whitespace[0]+alpha_chain+whitespace[2:]
                        elif sp_first_indeces[1]<=i<sp_first_indeces[2]:
                            if len(whitespace)==2:
                                whitespace = whitespace[0]+beta_chain
                            else:
                                whitespace = whitespace[0]+beta_chain+whitespace[2:]
                        elif sp_first_indeces[2]<=i:
                            if len(whitespace)==2:
                                whitespace = whitespace[0]+gamma_chain
                            else:
                                whitespace = whitespace[0]+gamma_chain+whitespace[2:]
                    out_line = group1+whitespace+pos_list[i]+pdb_re.group(4)
                    ### Needed when complex remodel is run
                    if out_line[-1]!='\n':
                        out_line+='\n'
                    out_list.append(out_line)
                except Exception as msg:
                    if line.startswith('TER'):
                        # pdb_re = re.search('(TER\s+\d+\s+\S{3})([\sA-Z]+)(\d+)',line)
                        # out_list.append(pdb_re.group(1)+len(pdb_re.group(2))*' '+pos_list[i]+"\n")
                        out_list.append('TER\n')
                        atom_num+=1
                    else:
                        try:
                            # Peptide ligands
                            if len(self.peptide_ligands)>0 and line.startswith('ATOM'):
                                print(line)
                                for l in self.peptide_ligands:
                                    if l[0].chain==line[21]:
                                        break
                                pdb_re = re.search('(ATOM[A-Z\s\d]{13}\S{3})([\sA-Z]+)(\d+)([A-Z\s\d.-]{49,53})',line)
                                if not first_ligand_pep_res:
                                    prev_num = pdb_re.group(3)
                                    prev_chain = pdb_re.group(2).strip()
                                    first_ligand_pep_res = True
                                if pdb_re.group(3)>prev_num:
                                    pep_i+=1
                                print(pep_i, l[-2][pep_i])
                                prev_num = pdb_re.group(3)
                                whitespace = self.pdb_line_whitespace(l[-2], pep_i, pdb_re.group(3), len(pdb_re.group(2)))
                                whitespace = ' '+l[0].chain+whitespace[2:]
                                out_list.append(pdb_re.group(1)+whitespace+str(l[-2][pep_i])+pdb_re.group(4))
                            ###
                            if len(pref_chain)>1:
                                pref_chain = pref_chain[0]
                            pdb_re = re.search('(HETATM[0-9\sA-Z{apo}]{{11}})([A-Z0-9\s]{{3}})([\sA-Z]+)(\d+)([\s0-9.A-Z-]+)'.format(apo="'"),line)
                            if not pdb_re:
                                continue
                            alternate_water = False
                            whitespace3 = (len(pdb_re.group(3))-2)*' '
                            if first_hetatm==False:
                                prev_hetnum = int(pdb_re.group(4))
                                first_hetatm = True
                                atom_num = int(pdb_re.group(1)[7:11])+1
                                if self.complex:
                                    num = int(pos_list[sp_first_indeces[0]-1])+100
                                else:
                                    num = int(pos_list[i])+100
                                if 'HOH' in pdb_re.group(2):
                                    water_count+=1
                                    if water_count in self.alternate_water_positions:
                                        if len(str(num))==3:
                                            whitespace1 = ' '
                                            whitespace2 = 5*' '
                                        else:
                                            whitespace1 = ''
                                            whitespace2 = 4*' '
                                        bwater = 'HETATM {}  O  BHOH {}{}{}{}'.format(str(atom_num+1), pref_chain, whitespace1[:-1], num+1, whitespace2)+self.alternate_water_positions[water_count][31:]
                                        alternate_water = True
                                if alternate_water==True:
                                    out_list.append(pdb_re.group(1)[:7]+str(atom_num)+pdb_re.group(1)[11:-1]+'A'+pdb_re.group(2)+' '+pref_chain+whitespace3+str(num+1)+pdb_re.group(5))
                                    out_list.append(bwater)
                                    atom_num+=2
                                else:
                                    out_list.append(pdb_re.group(1)[:7]+str(atom_num)+pdb_re.group(1)[11:]+pdb_re.group(2)+' '+pref_chain+whitespace3+str(num)+pdb_re.group(5))
                                    atom_num+=1
                            else:
                                if int(pdb_re.group(4))!=prev_hetnum:
                                    if 'HOH' in pdb_re.group(2):
                                        water_count+=1
                                        if water_count in self.alternate_water_positions:
                                            if len(str(num))==3:
                                                whitespace1 = ' '
                                                whitespace2 = 5*' '
                                            else:
                                                whitespace1 = ''
                                                whitespace2 = 4*' '
                                            bwater = 'HETATM {}  O  BHOH {}{}{}{}'.format(str(atom_num+1), pref_chain, whitespace1[:-1], num+1, whitespace2)+self.alternate_water_positions[water_count][31:]
                                            alternate_water = True
                                    if alternate_water==True:
                                        out_list.append(pdb_re.group(1)[:7]+str(atom_num)+pdb_re.group(1)[11:-1]+'A'+pdb_re.group(2)+' '+pref_chain+whitespace3+str(num+1)+pdb_re.group(5))
                                        out_list.append(bwater)
                                        atom_num+=2
                                    else:
                                        out_list.append(pdb_re.group(1)[:7]+str(atom_num)+pdb_re.group(1)[11:]+pdb_re.group(2)+' '+pref_chain+whitespace3+str(num+1)+pdb_re.group(5))
                                        atom_num+=1
                                    prev_hetnum+=1
                                    num+=1

                                else:
                                    out_list.append(pdb_re.group(1)+pdb_re.group(2)+' '+pref_chain+whitespace3+str(num)+pdb_re.group(5))
                                    atom_num+=1
                        except:
                            out_list.append(line)

        with open (path+self.modelname+'.pdb', 'w') as f:
            f.write(''.join(out_list))

        with open (path+self.modelname+'.pdb', 'r') as f:
            pdb_struct = f.read()

        pdb_struct = StringIO(pdb_struct)
        assign_gn = as_gn.GenericNumbering(pdb_file=pdb_struct, sequence_parser=True)
        pdb_struct = assign_gn.assign_generic_numbers_with_sequence_parser()

        # Check if gn assigning is ok
        if self.revise_xtal==None:
            helices = ProteinSegment.objects.filter(category='helix', proteinfamily='GPCR')
            misassigned = []
            for chain in pdb_struct:
                if self.complex:
                    if chain.get_id()!='R':
                        continue
                for res in chain:
                    residue = Residue.objects.get(protein_conformation=self.prot_conf, sequence_number=int(res.get_id()[1]))
                    if residue.protein_segment not in helices:
                        continue
                    gn = str(res['CA'].get_bfactor()).replace('.','x')
                    if len(gn.split('x')[1])==1:
                        gn = gn+'0'
                    if gn.startswith('-'):
                        gn = gn[1:]+'1'
                    display_gn = dgn(gn, self.prot_conf)
                    if residue.display_generic_number.label!=display_gn:
                        misassigned.append([residue.display_generic_number.label, display_gn])
            if len(misassigned)>0:
                print("MISASSIGNED GNS: ")
                print(misassigned)
                logger.error("Misassigned gns in {}:\n{}".format(self.reference_entry_name, misassigned))
            if len(misassigned)!=0:
                return None

        # add SSBOND and CONECT
        ssbond, conect = '',''
        count = 0
        if self.debug:
            pprint.pprint(self.disulfide_pairs)
        for i, j in enumerate(self.disulfide_pairs):
            if j[0]!=0:
                count+=1
                if j[0].startswith('ECL3'):
                    c1 = Residue.objects.filter(protein_conformation=self.prot_conf, protein_segment__slug='ECL3')[int(j[0][-1])-1]
                    c2 = Residue.objects.filter(protein_conformation=self.prot_conf, protein_segment__slug='ECL3')[int(j[1][-1])-1]
                elif j[0].startswith('6x'):
                    c1 = Residue.objects.get(protein_conformation=self.prot_conf, display_generic_number__label=dgn(j[0], self.prot_conf))
                    c2 = Residue.objects.filter(protein_conformation=self.prot_conf, protein_segment__slug='ECL3')[int(j[1][-1])-1]
                elif j[0].startswith('3x2'):
                    c1 = Residue.objects.get(protein_conformation=self.prot_conf, display_generic_number__label=dgn(j[0], self.prot_conf))
                    c2 = Residue.objects.get(protein_conformation=self.prot_conf, display_generic_number__label=dgn(j[1], self.prot_conf))
                else:
                    c1 = Residue.objects.get(protein_conformation=self.prot_conf, sequence_number=int(j[0]))
                    c2 = Residue.objects.get(protein_conformation=self.prot_conf, sequence_number=int(j[1]))
                for c in pdb_struct:
                    if self.debug:
                        for res in c:
                            print(res)
                    chain = c[int(c1.sequence_number)].get_parent().get_id()
                    break
                ws1 = ' '*(5-len(str(c1.sequence_number)))
                ws2 = ' '*(5-len(str(c2.sequence_number)))
                ssbond+='SSBOND    {} CYS {}{}{}    CYS {}{}{}\n'.format(count, chain, ws1, str(c1.sequence_number), chain, ws2, str(c2.sequence_number))
                # for res in pdb_struct[chain]:
                    # print(res)
                for atom in pdb_struct[chain][c1.sequence_number]:
                    print('c1',atom, atom.get_serial_number(), c1.sequence_number)
                for atom in pdb_struct[chain][c2.sequence_number]:
                    print('c2',atom, atom.get_serial_number(), c2.sequence_number)
                sg1 = pdb_struct[chain][c1.sequence_number]['SG'].get_serial_number()
                sg2 = pdb_struct[chain][c2.sequence_number]['SG'].get_serial_number()
                offset = 0
                for a in atom_num_offset:
                    if a<c2.sequence_number:
                        offset+=1
                print(sg1,sg2,atom_num_offset, offset)
                sg1-=offset
                sg2-=offset
                print(sg1,sg2)
                ws3 = ' '*(5-len(str(sg1)))
                ws4 = ' '*(5-len(str(sg2)))
                conect+='CONECT{}{}{}{}\n'.format(ws3, str(sg1), ws4, str(sg2))

        # Put in peptide ligands
        if len(self.peptide_ligands)>0:
            for l in self.peptide_ligands:
                pdb_struct.detach_child(l[0].chain)
                pdb_struct.add(l[2][l[0].chain])

        ### FIXME later - fetch UNK unnatural residues from original file and keep the modeled natural ones
        # if len(self.peptide_ligands)>0:
        #     for l in self.peptide_ligands:
        #         model_peptide_resis = [i for i in pdb_struct[l[0].chain]]
        #         print(model_peptide_resis)
        #         orig_peptide_resis = [i for i in l[2][l[0].chain]]
        #         print(orig_peptide_resis)
        #         for i, res in enumerate(pdb_struct[l[0].chain]):
        #             if res.get_resname()=='UNK':
        #                 print(res, model_peptide_resis[i].get_id()[1], pdb_struct[l[0].chain][model_peptide_resis[i].get_id()[1]])
        #                 pdb_struct[l[0].chain].detach_child(model_peptide_resis[i].get_id())
        #                 pdb_struct[l[0].chain].add(l[2][l[0].chain][orig_peptide_resis[i].get_id()])
        # pprint.pprint(self.peptide_ligands)

        io = PDB.PDBIO()
        io.set_structure(pdb_struct)
        io.save(path+self.modelname+'.pdb')

        with open (path+self.modelname+'.pdb', 'r+') as f:
            content = f.read()
            first_line  = 'REMARK    1 MODEL FOR {} CREATED WITH GPCRDB HOMOLOGY MODELING PIPELINE, VERSION {}\n'.format(self.reference_entry_name, build_date)
            second_line = 'REMARK    2 MAIN TEMPLATE: {}\n'.format(self.main_structure)
            f.seek(0,0)
            f.write(first_line+second_line+ssbond+content[:-7]+conect+'END')

        return first_line+second_line+ssbond+content[:-7]+conect+'END'

    def run_alignment(self, query_states, core_alignment=True,
                      segments=['TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','TM5','TM6','TM7','H8'],
                      order_by='similarity'):
        ''' Creates pairwise alignment between reference and target receptor(s).
            Returns Alignment object.

            @param core_alignment: boolean, False if only create core alignment (no loops)
            @param query_states: list, list of endogenous ligand states to be applied for template search
            @param segments: list, list of segments to use, e.g.: ['TM1','ICL1','TM2','ECL1'] \n
            @param order_by: str, order results by identity, similarity or simscore
        '''
        alignment = AlignedReferenceTemplate()
        alignment.run_hommod_alignment(self.reference_protein, segments, query_states, order_by, complex_model=self.complex, signprot=self.signprot, force_main_temp=self.force_main_temp,
                                       core_alignment=core_alignment)
        parse = GPCRDBParsingPDB()
        main_pdb_array = OrderedDict()
        if core_alignment==True:
            if self.debug:
                print('Alignment: ',datetime.now() - startTime)
            self.changes_on_db = alignment.changes_on_db
            alignment.enhance_alignment(alignment.reference_protein, alignment.main_template_protein)
            if self.debug:
                print('Enhanced alignment: ',datetime.now() - startTime)
            self.segments = segments
            self.main_structure = alignment.main_template_structure
            if self.debug:
                print('Main structure: {}'.format(self.main_structure))
            self.similarity_table = alignment.similarity_table
            self.similarity_table_all = self.run_alignment(["Inactive","Intermediate","Active"], core_alignment=False)[0].similarity_table
            for i,j in self.similarity_table_all.items():
                if i not in self.similarity_table:
                    self.similarity_table_other_states[i] = j
            self.main_template_preferred_chain = str(self.main_structure.preferred_chain)[0]
            self.statistics.add_info("main_template", self.main_structure)
            self.statistics.add_info("preferred_chain", self.main_template_preferred_chain)


            main_pdb_array = parse.pdb_array_creator(structure=self.main_structure)

            ### main_pdb_array sanity check
            if self.debug:
                print('main_pdb_array sanity check:')
                for i, j in main_pdb_array.items():
                    for k,l in j.items():
                        try:
                            print(k, l[0].get_parent())
                        except:
                            print(k, l)
            if self.main_structure.pdb_code.index=='4OR2':
                main_pdb_array['H8'] = OrderedDict()

            try:
                if len(alignment.reference_dict['H8'])==0:
                    del alignment.reference_dict['H8']
                    del alignment.template_dict['H8']
                    del alignment.alignment_dict['H8']
                    del main_pdb_array['H8']
                elif 'H8' not in list(main_pdb_array.keys()):
                    main_pdb_array['H8'] = OrderedDict()
            except:
                pass

            for seg_l, seg in main_pdb_array.items():
                for gn, res in seg.items():
                    self.template_source = update_template_source(self.template_source, [gn.replace('.','x')],self.main_structure,seg_l)

            helixends = HelixEndsModeling(self.similarity_table_all, self.template_source, self.main_structure, self.revise_xtal, self.debug)

            try:
                if (len(main_pdb_array['H8'])==0 and len(list(Residue.objects.filter(protein_conformation=self.prot_conf, protein_segment__slug='H8')))>0 or
                   (self.reference_protein.family.slug.startswith('004') and self.main_structure.pdb_code.index not in ['4OO9']) or
                   (self.main_structure.pdb_code.index in ['5UNF','5UNG','5UNH','5O9H'] and self.revise_xtal==False)):
                    if self.debug:
                        print('Main template {} H8 is not usable, skipping'.format(self.main_structure))
                    helixends.correct_helix_ends(self.main_structure, main_pdb_array, alignment,
                                                 self.template_source, separate_H8=True)
                    main_pdb_array = helixends.main_pdb_array
                    alignment = helixends.alignment

                    self.template_source = helixends.template_source
                    self.helix_end_mods = helixends.helix_end_mods
                    if self.reference_protein.family.slug.startswith('004'):
                        struct = Structure.objects.get(pdb_code__index='4OO9')
                        alt_simtable = self.similarity_table_all
                        alt_simtable[struct] = 0
                        gn_list = list(Residue.objects.filter(protein_conformation=struct.protein_conformation,
                                                              protein_segment__slug='H8'))
                    else:
                        alt_simtable = self.similarity_table_all
                        for struct in alt_simtable:
                            if struct.pdb_code.index in ['5UNF','5UNG','5UNH','5O9H']:
                            ## NMUR2
                            # if struct.pdb_code.index in ['5UNF','5UNG','5UNH','5O9H','7F9Z', '6KO5', '7F9Y', '4BUO', '4XES', '4GRV', '4XEE', '5T04', '7L0P', '7L0R', '7L0Q', '7L0S', '6OSA','6PWC','6OS9','6UP7','6B73','6VI4','4DJH','4N6H','6PT2','5XSZ','4RWD','6PT3','4EJ4','5X93','6TO7','5WS3','7DDZ','6ME2','5DHG','7L1V','7DB6','5WQC','6IGK','6TOD','6TOS','6TOT','6TQ4','6TP6','5GLI','4S0V','6TQ6','6TPN','6TQ7','6TQ9','6IGL','6K1Q','6TPG','6TPJ','4ZJ8','5GLH','4ZJC','6ME3','6LRY','5DHH','4EA3','6TP4','6TP3','7L1U','6ME4','6ME5','6PS8','6V9S','5XPR','5C1M','7F8Y','6BQH','7CMV','7JVR','4ZUD','4DKL','6CM4','3PBL','7EZM','6JOD','7MBX','7MBY','7BVQ','7BU7','7BU6','6OS2','6OS1','5UNF','5UNG','7F8U','6OS0','6DO1','4YAY','5UNH','6BQG','7CMU','7DFP','6LUQ','7EZK','7BTS','5XJM']:
                            # if struct.pdb_code.index in ['5UNF','5UNG','5UNH','5O9H','7F9Z']:
                            ## Kappa opioid
                            # if struct.pdb_code.index in ['5UNF','5UNG','5UNH','5O9H','6B73', '6VI4', '4DJH', '5C1M', '4DKL', '6DDF']:
                            ## NPY1R
                            # if struct.pdb_code.index in ['5ZBH','5ZBQ','7DDZ','6TO7','5WS3','6ME2','7L1V','7F8W','7DB6','7L1U','6HLP','7F8Y','7EZM']:#,'6B73','6VI4','5C1M','4DKL','6DDE','4N6H','6PT2','7E32','5DHG','7EXD','1GZM','7JVR','6CM4','5JQH','3SN6','7F9Z','6KO5','4ZWJ','6OFJ','7F9Y','5X93','5ZKC','6OL9','6BQH','6KUX','7CMV','4IAQ','4ZUD','3PBL','7E33']:
                                continue
                            try:
                                gn_list = list(Residue.objects.filter(protein_conformation=struct.protein_conformation,
                                                                      protein_segment__slug='H8'))
                                if len(gn_list)>0:
                                    break
                            except:
                                pass
                    for i in alignment.ordered_proteins:
                        if i.protein.entry_name==struct.protein_conformation.protein.parent.entry_name:
                            break
                    H8_alignment = AlignedReferenceTemplate()
                    H8_alignment.enhance_alignment(alignment.ordered_proteins[0],i)

        ######### temporary
                    reference_dict, template_dict, alignment_dict = OrderedDict(),OrderedDict(),OrderedDict()
                    for i,j,k in zip(H8_alignment.reference_dict['H8'],H8_alignment.template_dict['H8'],H8_alignment.alignment_dict['H8']):
                        if i in self.template_source['H8']:
                            reference_dict[i] = H8_alignment.reference_dict['H8'][i]
                            template_dict[i] = H8_alignment.template_dict['H8'][i]
                            alignment_dict[i] = H8_alignment.alignment_dict['H8'][i]
        ###################
        ######### change values
                    alignment.reference_dict['H8'] = reference_dict
                    alignment.template_dict['H8'] = template_dict
                    alignment.alignment_dict['H8'] = alignment_dict
        #######################

                    gn_num_list = [ggn(i.display_generic_number.label) for i in gn_list if i.display_generic_number!=None]
                    found_match = False
                    c1 = -4
                    c2 = None
                    while found_match==False:
                        refs = list(main_pdb_array['TM7'].keys())[c1:c2]
                        try:
                            for gn in refs:
                                Residue.objects.get(protein_conformation=struct.protein_conformation,
                                                    display_generic_number__label=dgn(gn.replace('.','x'),struct.protein_conformation))
                            found_match=True
                        except:
                            c1-=1
                            if c2==None:
                                c2 = -1
                            else:
                                c2-=1
                            if c1<-10:
                                break

                    refs = [i.replace('.','x') for i in refs]
                    H8_reference = parse.fetch_residues_from_array(main_pdb_array['TM7'], refs)
                    H8_template = parse.fetch_residues_from_pdb(struct, refs+gn_num_list)
                    superpose = sp.OneSidedSuperpose(H8_reference,H8_template,4,1)
                    sup_residues = superpose.run()
                    H8_array = OrderedDict()
                    for i,j in alignment.template_dict['H8'].items():
                        if j not in ['-','x']:
                            try:
                                H8_array[i.replace('x','.')] = sup_residues[i.replace('x','.')]
                            except:
                                H8_array[i.replace('x','.')] = 'x'
                        else:
                            H8_array[i.replace('x','.')] = 'x'
                    main_pdb_array['H8'] = H8_array
                    for gn, res in main_pdb_array['H8'].items():
                        try:
                            if gn.replace('.','x') in gn_num_list:
                                self.template_source = update_template_source(self.template_source,[gn.replace('.','x')],struct,'H8')
                        except:
                            pass
                    if self.debug:
                        print('Correcting helix ends of only H8')
                    pprint.pprint(self.template_source['H8'])
                    helixends.correct_helix_ends(self.main_structure, main_pdb_array, alignment,
                                                 self.template_source, separate_H8=False)
                    self.helix_end_mods['added']['H8'] = helixends.helix_end_mods['added']['H8']
                    self.helix_end_mods['removed']['H8'] = helixends.helix_end_mods['removed']['H8']
                    self.template_source = helixends.template_source
                else:
                    raise Exception()
            except:
                if len(list(Residue.objects.filter(protein_conformation=self.prot_conf, protein_segment__slug='H8')))==0:
                    sep_H8 = True
                else:
                    sep_H8 = None
                if self.debug:
                    if sep_H8:
                        print('H8 build exception: skipping H8 ')
                    else:
                        print('H8 build exception: correcting all helices')
                helixends.correct_helix_ends(self.main_structure, main_pdb_array, alignment,
                                             self.template_source, separate_H8=sep_H8)
                self.helix_end_mods = helixends.helix_end_mods
                self.template_source = helixends.template_source

            self.statistics.add_info('helix_end_mods',self.helix_end_mods)

            if self.debug:
                print('Corrected helix ends: ',datetime.now() - startTime)

            main_pdb_array = helixends.main_pdb_array
            alignment = helixends.alignment

            loops_in_ref = [i for i in list(self.template_source) if i[0] not in ['N','C','T','H']]
            for loop in loops_in_ref:
                loop_alignment = AlignedReferenceTemplate()
                loop_alignment.run_hommod_alignment(self.reference_protein, [loop], ['Inactive','Intermediate','Active'],
                                                    order_by='similarity',
                                                    provide_main_template_structure=self.main_structure,
                                                    provide_similarity_table=self.similarity_table_all,
                                                    main_pdb_array=main_pdb_array, provide_alignment=alignment, complex_model=self.complex)
                self.loop_template_table[loop] = loop_alignment.loop_table
                try:
                    if loop in list(alignment.alignment_dict.keys()) and self.main_structure in loop_alignment.loop_table:
                        temp_loop_table = OrderedDict([('aligned',100)])
                        try:
                            for lab, val in loop_alignment.loop_table.items():
                                temp_loop_table[lab] = val
                            self.loop_template_table[loop] = temp_loop_table
                        except:
                            pass
                except:
                    pass
            self.statistics.add_info('similarity_table', self.similarity_table)
            self.statistics.add_info('loops',self.loop_template_table)
            if self.debug:
                print('Loop alignment: ',datetime.now() - startTime)

        return alignment, main_pdb_array


    def build_homology_model(self, ref_temp_alignment, switch_bulges=True, switch_constrictions=True, loops=True,
                             switch_rotamers=True, N_and_C_termini=True):
        ''' Function to identify and switch non conserved residues in the alignment. Optionally,
            it can identify and switch bulge and constriction sites too.

            @param ref_temp_alignment: AlignedReferenceAndTemplate, alignment of reference and main template with
            alignment string. \n
            @param switch_bulges: boolean, identify and switch bulge sites. Default = True.
            @param switch_constrictions: boolean, identify and switch constriction sites. Default = True.
            @param loops: boolean, set it to True if you want loop modeling. Default = True.
            @param switch_rotamers: boolean, set it to True if you want alternative rotamer tempaltes. Default = True.
            @param N_and_C_termini: boolean, set it to True if you want to model N/C-termini. Only applies for refining
            crystals. Default = True.
        '''
        a = ref_temp_alignment[0]
        main_pdb_array = ref_temp_alignment[1]
        ref_bulge_list, temp_bulge_list, ref_const_list, temp_const_list = [],[],[],[]
        parse = GPCRDBParsingPDB()
        # for i,j,k,l in zip(a.reference_dict, a.template_dict, a.alignment_dict, main_pdb_array):
        #     print(a.reference_dict[i])
        #     print(a.template_dict[j])
        #     print(a.alignment_dict[k])
        #     for o in main_pdb_array[l]:
        #         try:
        #             print(main_pdb_array[l][o][0].get_parent())
        #         except:
        #             print(main_pdb_array[l][o])
        # raise AssertionError
        # Delete H8 from dictionaries if it's not present in reference (e.g. gnrhr_human)
        if self.revise_xtal:
            del_H8_prot = self.reference_protein.parent
        else:
            del_H8_prot = self.reference_protein
        if len(Residue.objects.filter(protein_conformation__protein=del_H8_prot, protein_segment__slug='H8'))==0:
            try:
                del a.reference_dict['H8']
                del a.template_dict['H8']
                del a.alignment_dict['H8']
                del main_pdb_array['H8']
            except:
                pass
        trimmed_residues=[]

        # loops
        if loops==True:
            c3x25 = {'001':'3x25','002':'3x29','003':'3x29','004':'3x29','006':'3x25'}
            model_loops = []
            loop_stat = OrderedDict()
            for label, structures in self.loop_template_table.items():
                if label in ['ICL1','ECL1','ICL2']:
                    x50_present = False
                    l_gns = list(Residue.objects.filter(protein_conformation=self.prot_conf, protein_segment__slug=label))
                    for i in l_gns:
                        try:
                            if 'x50' in i.display_generic_number.label:
                                structures = self.similarity_table_all
                                x50_present = True
                                break
                        except:
                            pass

                loop = Loops(self.reference_protein, label, structures, self.main_structure, self.helix_end_mods,
                             list(self.template_source), self.revise_xtal)

                loop_template = loop.fetch_loop_residues(main_pdb_array)

                if (loop.loop_output_structure not in [self.main_structure,None] and label in ['ICL1','ECL1','ICL2'] and
                    x50_present==True):
                    al = AlignedReferenceTemplate()
                    t = al.run_hommod_alignment(self.reference_protein, [label], ['Inactive','Intermediate','Active'],
                                                order_by='similarity',
                                                only_output_alignment=loop.loop_output_structure.protein_conformation.protein.parent)
                    al.enhance_alignment(t[0],t[1])
                    a.reference_dict[label] = al.reference_dict[label]
                    a.template_dict[label] = al.template_dict[label]
                    a.alignment_dict[label] = al.alignment_dict[label]

                if label=='ECL2' and (loop.partialECL2_1==True or loop.partialECL2_2==True):
                    al = AlignedReferenceTemplate()
                    t = al.run_hommod_alignment(self.reference_protein, [label], ['Inactive','Intermediate','Active'],
                                                order_by='similarity',
                                                only_output_alignment=loop.loop_output_structure[1].protein_conformation.protein.parent)
                    al.enhance_alignment(t[0],t[1],keep_all=True)
                    a.reference_dict[label] = al.reference_dict[label]
                    a.template_dict[label] = al.template_dict[label]
                    a.alignment_dict[label] = al.alignment_dict[label]

                if type(loop.loop_output_structure)!=type([]):
                    loop_insertion = loop.insert_loop_to_arrays(loop.loop_output_structure, main_pdb_array, loop_template,
                                                                a.reference_dict, a.template_dict, a.alignment_dict)
                else:
                    print(loop_template, loop.loop_output_structure, loop.partialECL2_1, loop.partialECL2_2)
                    loop_insertion = loop.insert_ECL2_to_arrays(loop.loop_output_structure, main_pdb_array, loop_template,
                                                                a.reference_dict, a.template_dict, a.alignment_dict,
                                                                loop.partialECL2_1, loop.partialECL2_2)
                if loop.model_loop==True and loop.new_label!=None:
                    model_loops.append(loop.new_label)
                main_pdb_array = loop_insertion.main_pdb_array
                a.reference_dict = loop_insertion.reference_dict
                a.template_dict = loop_insertion.template_dict
                a.alignment_dict = loop_insertion.alignment_dict

                # update template_source with backbone template of loop
                if loop.new_label!=None:
                    change_i, change_template_list = [], []
                    loop_stat[loop.new_label] = loop.loop_output_structure
                    ref_loop_ids = a.reference_dict[loop.new_label]
                    c = 0
                    for i, v in ref_loop_ids.items():
                        if v!='x':
                            if '?' not in i:
                                change_i.append(c)
                            c+=1
                    for i, v in enumerate(list(self.template_source[label])):
                        if i in change_i:
                            change_template_list.append(v)
                    self.template_source = update_template_source(self.template_source, change_template_list,loop.loop_output_structure,label)
                else:
                    loop_stat[label] = loop.loop_output_structure
                    add_ECL2_disulfide = False
                    if label=='ECL2' and loop.loop_output_structure!=None:
                        change_i1, change_i2, change_i3, change_templates1, change_templates2, change_templates3 = [],[],[],[],[],[]
                        x50 = list(self.template_source[label].keys()).index('45x50')
                        c = 0
                        for i, v in a.reference_dict[label].items():
                            if v!='x':
                                if '?' not in i:
                                    if c<x50:
                                        change_i1.append(c)
                                    elif x50<=c<x50+3:
                                        change_i2.append(c)
                                    elif x50+2<c:
                                        change_i3.append(c)
                                c+=1
                        for i, v in enumerate(list(self.template_source[label])):
                            if i in change_i1:
                                change_templates1.append(v)
                            elif i in change_i2:
                                change_templates2.append(v)
                            elif i in change_i3:
                                change_templates3.append(v)
                        self.template_source = update_template_source(self.template_source,change_templates1,loop.loop_output_structure[0],label)
                        self.template_source = update_template_source(self.template_source,change_templates2,loop.loop_output_structure[1],label)
                        self.template_source = update_template_source(self.template_source,change_templates3,loop.loop_output_structure[2],label)
                        add_ECL2_disulfide = True
                    # if chain break in ECL2_2 is predicted, let MODELLER optimize conserved middle part, while keeping disulfide bridge
                    if label=='ECL2' and loop.evade_chain_break:
                        trimmed_residues+=['45.50','45.51','45.52']
                        add_ECL2_disulfide = True
                    # add 3x25-45x50 disulfide bond
                    if label=='ECL2' and add_ECL2_disulfide:
                        # if loop.loop_output_structure[1]!=self.main_structure:
                            self.reference_class
                            self.disulfide_pairs.append([c3x25[self.reference_class.slug[:3]],'45x50'])
                        # else:
                        #     self.disulfide_pairs.append([0,0])

            self.statistics.add_info('loops', loop_stat)
            self.loops = loop_stat

        if self.debug:
            print(loop_stat)
            print('Integrate loops: ',datetime.now() - startTime)

        # for i,j,k,l in zip(a.reference_dict, a.template_dict, a.alignment_dict, main_pdb_array):
        #     print(a.reference_dict[i])
        #     print(a.template_dict[j])
        #     print(a.alignment_dict[k])
        #     for o in main_pdb_array[l]:
        #         try:
        #             print(main_pdb_array[l][o][0].get_parent())
        #         except:
        #             print(main_pdb_array[l][o])
        # raise AssertionError

        # bulges and constrictions
        if switch_bulges==True or switch_constrictions==True:
            delete_r = set()
            delete_t = set()
            delete_a = set()

            ref_seg_keys = list(a.reference_dict.keys())
            temp_seg_keys = list(a.template_dict.keys())
            aligned_seg_keys = list(a.alignment_dict.keys())

            for ref_seg, temp_seg, aligned_seg in zip(ref_seg_keys, temp_seg_keys, aligned_seg_keys):
                if ref_seg[0]=='T':
                    ref_res_keys = list(a.reference_dict[ref_seg].keys())
                    temp_res_keys = list(a.template_dict[temp_seg].keys())
                    aligned_res_keys = list(a.alignment_dict[aligned_seg].keys())
                    for ref_res, temp_res, aligned_res in zip(ref_res_keys, temp_res_keys,
                                                              aligned_res_keys):
                        gn = ref_res
                        gn_num = parse.gn_num_extract(gn, 'x')[1]

                        if a.alignment_dict[aligned_seg][aligned_res]=='-':
                            if (a.reference_dict[ref_seg][ref_res]=='-' and
                                a.reference_dict[ref_seg][parse.gn_indecer(gn,'x',-1)] not in
                                ['-','/'] and a.reference_dict[ref_seg][parse.gn_indecer(gn,'x',+1)] not in ['-','/']):

                                # bulge in template
                                if len(str(gn_num))==3:
                                    if switch_bulges==True:
                                        try:
                                            Bulge = Bulges(gn)
                                            bulge_template = Bulge.find_bulge_template(self.similarity_table_all,
                                                                                       bulge_in_reference=False)
                                            l = list(main_pdb_array[temp_seg].keys())
                                            this = l.index(gn.replace('x','.'))
                                            bulge_site = OrderedDict([(l[this-2],main_pdb_array[ref_seg][l[this-2]]),
                                                                      (l[this-1],main_pdb_array[ref_seg][l[this-1]]),
                                                                      (l[this],main_pdb_array[ref_seg][l[this]]),
                                                                      (l[this+1],main_pdb_array[ref_seg][l[this+1]]),
                                                                      (l[this+2],main_pdb_array[ref_seg][l[this+2]])])
                                            superpose = sp.BulgeConstrictionSuperpose(bulge_site, bulge_template)
                                            new_residues = superpose.run()
                                            switch_res = 0
                                            for gen_num, atoms in bulge_template.items():
                                                if switch_res!=0 and switch_res!=3:
                                                    gn__ = gen_num.replace('.','x')
                                                    self.template_source = update_template_source(self.template_source,[gn__],Bulge.template,ref_seg)
                                                    main_pdb_array[ref_seg][gen_num] = new_residues[gen_num]
                                                    a.template_dict[temp_seg][gn__] = PDB.Polypeptide.three_to_one(
                                                                                       atoms[0].get_parent().get_resname())
                                                    if a.template_dict[temp_seg][gn__]==a.reference_dict[ref_seg][gn__]:
                                                        a.alignment_dict[aligned_seg][gn__]=a.template_dict[temp_seg][gn__]
                                                    else:
                                                        a.alignment_dict[aligned_seg][gn__]='.'
                                                switch_res+=1
                                            del main_pdb_array[ref_seg][gn.replace('x','.')]
                                            delete_r.add((ref_seg,gn))
                                            delete_t.add((temp_seg,gn))
                                            delete_a.add((aligned_seg,gn))

                                            temp_bulge_list.append({gn:Bulge.template})
                                        except:
                                            temp_bulge_list.append({gn:None})

                                # constriction in reference
                                else:
                                    if switch_constrictions==True:
                                        try:
                                            Const = Constrictions(gn)
                                            constriction_template = Const.find_constriction_template(
                                                                                            self.similarity_table_all,
                                                                                            constriction_in_reference=True)
                                            constriction_site = OrderedDict([
                                                (parse.gn_indecer(gn,'x',-2).replace('x','.'),
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',-2).replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',-1).replace('x','.'),
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',-1).replace('x','.')]),
                                                (gn.replace('x','.'),
                                                 main_pdb_array[ref_seg][gn.replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',+1).replace('x','.'),
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',+1).replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',+2).replace('x','.'),
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',+2).replace('x','.')])])
                                            superpose = sp.BulgeConstrictionSuperpose(constriction_site,
                                                                                      constriction_template)
                                            new_residues = superpose.run()
                                            switch_res = 0
                                            for gen_num, atoms in constriction_template.items():
                                                if switch_res!=0 and switch_res!=3:
                                                    gn__ = gen_num.replace('.','x')
                                                    self.template_source = update_template_source(self.template_source,[gn__],Const.template,ref_seg)
                                                    main_pdb_array[ref_seg][gen_num] = new_residues[gen_num]
                                                    a.template_dict[gn__] = PDB.Polypeptide.three_to_one(
                                                                                       atoms[0].get_parent().get_resname())
                                                    if a.template_dict[temp_seg][gn__]==a.reference_dict[ref_seg][gn__]:
                                                        a.alignment_dict[aligned_seg][gn__]=a.template_dict[temp_seg][gn__]
                                                switch_res+=1
                                            ref_const_list.append({gn:Const.template})
                                            del main_pdb_array[ref_seg][gn.replace('x','.')]
                                            delete_r.add((ref_seg,gn))
                                            delete_t.add((temp_seg,gn))
                                            delete_a.add((aligned_seg,gn))
                                        except:
                                            ref_const_list.append({gn:None})
                            elif (a.template_dict[ref_seg][temp_res]=='-' and
                                  a.template_dict[temp_seg][parse.gn_indecer(gn,'x',-1)] not in
                                  ['-','/'] and a.template_dict[temp_seg][parse.gn_indecer(gn,'x',+1)] not in ['-','/']):

                                # bulge in reference
                                if len(str(gn_num))==3:
                                    if switch_bulges==True:
                                        try:
                                            Bulge = Bulges(gn)
                                            bulge_template = Bulge.find_bulge_template(self.similarity_table_all,
                                                                                       bulge_in_reference=True)
                                            bulge_site = OrderedDict([
                                                (parse.gn_indecer(gn,'x',-2).replace('x','.'),
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',-2).replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',-1).replace('x','.'),
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',-1).replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',+1).replace('x','.'),
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',+1).replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',+2).replace('x','.'),
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',+2).replace('x','.')])])
                                            superpose = sp.BulgeConstrictionSuperpose(bulge_site, bulge_template)
                                            new_residues = superpose.run()
                                            switch_res = 0
                                            for gen_num, atoms in bulge_template.items():
                                                if switch_res!=0 and switch_res!=4:
                                                    gn__ = gen_num.replace('.','x')
                                                    self.template_source = update_template_source(self.template_source,[gn__],Bulge.template,ref_seg)
                                                    main_pdb_array[ref_seg][gen_num] = new_residues[gen_num]
                                                    a.template_dict[temp_seg][gn__] = PDB.Polypeptide.three_to_one(
                                                                                       atoms[0].get_parent().get_resname())
                                                    if a.template_dict[temp_seg][gn__]==a.reference_dict[ref_seg][gn__]:
                                                        a.alignment_dict[aligned_seg][gn__]=a.template_dict[temp_seg][gn__]
                                                switch_res+=1
                                            ref_bulge_list.append({gn:Bulge.template})
                                            if a.reference_dict[ref_seg][gn] == a.template_dict[temp_seg][gn]:
                                                a.alignment_dict[ref_seg][gn] = a.reference_dict[ref_seg][gn]
                                            else:
                                                a.alignment_dict[ref_seg][gn] = '.'
                                        except:
                                            ref_bulge_list.append({gn:None})

                                # constriction in template
                                else:
                                    if switch_constrictions==True:
                                        try:
                                            Const = Constrictions(gn)
                                            constriction_template = Const.find_constriction_template(
                                                                                           self.similarity_table_all,
                                                                                           constriction_in_reference=False)
                                            constriction_site = OrderedDict([
                                                (parse.gn_indecer(gn,'x',-2).replace('x','.'),
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',-2).replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',-1).replace('x','.'),
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',-1).replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',+1).replace('x','.'),
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',+1).replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',+2).replace('x','.'),
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',+2).replace('x','.')])])
                                            superpose = sp.BulgeConstrictionSuperpose(constriction_site,
                                                                                      constriction_template)
                                            new_residues = superpose.run()
                                            switch_res = 0
                                            for gen_num, atoms in constriction_template.items():
                                                if switch_res!=0 and switch_res!=4:
                                                    gn__ = gen_num.replace('.','x')
                                                    self.template_source = update_template_source(self.template_source,[gn__],Const.template,ref_seg)
                                                    main_pdb_array[ref_seg][gen_num] = new_residues[gen_num]
                                                    a.template_dict[temp_seg][gn__] = PDB.Polypeptide.three_to_one(
                                                                                       atoms[0].get_parent().get_resname())
                                                    if a.template_dict[temp_seg][gn__]==a.reference_dict[ref_seg][gn__]:
                                                        a.alignment_dict[aligned_seg][gn__]=a.template_dict[temp_seg][gn__]
                                                switch_res+=1
                                            temp_const_list.append({gn:Const.template})
                                            if a.reference_dict[ref_seg][gn] == a.template_dict[temp_seg][gn]:
                                                a.alignment_dict[ref_seg][gn] = a.reference_dict[ref_seg][gn]
                                            else:
                                                a.alignment_dict[ref_seg][gn] = '.'
                                        except:
                                            temp_const_list.append({gn:None})

            for i,ii in delete_r:
                del a.reference_dict[i][ii]
            for i,ii in delete_t:
                del a.template_dict[i][ii]
            for i,ii in delete_a:
                del a.alignment_dict[i][ii]

            self.statistics.add_info('reference_bulges', ref_bulge_list)
            self.statistics.add_info('template_bulges', temp_bulge_list)
            self.statistics.add_info('reference_constrictions', ref_const_list)
            self.statistics.add_info('template_constrictions', temp_const_list)

            # insert bulge to array in the right place
            if ref_bulge_list!=[]:
                out_pdb_array = OrderedDict()
                bulge_gns = []
                for bulge in ref_bulge_list:
                    if list(bulge.values())[0]!=None:
                        gn = list(bulge.keys())[0].replace('x','.')
                        bulge_gns.append(gn)
                for seg_id, residues in main_pdb_array.items():
                    seg = OrderedDict()
                    for key, value in residues.items():
                        seg[key] = value
                        if str(key)+'1' in bulge_gns:
                            seg[str(key)+'1'] = main_pdb_array[seg_id][str(key)+'1']
                    out_pdb_array[seg_id] = seg
                main_pdb_array = out_pdb_array
            if temp_const_list!=[]:
                out_pdb_array = OrderedDict()
                const_gns = []
                for const in temp_const_list:
                    gn_ = list(const.keys())[0].replace('x','.')
                    const_gns.append(gn_)
                for seg_id, residues in main_pdb_array.items():
                    seg = OrderedDict()
                    for key, value in residues.items():
                        seg[key] = value
                        if parse.gn_indecer(key, '.', +1) in const_gns:
                            gn = parse.gn_indecer(key, '.', +1)
                            try:
                                seg[gn] = main_pdb_array[seg_id][gn]
                            except:
                                seg[gn.replace('.','?')] = '-'
                                a.alignment_dict[seg_id][gn.replace('.','?')] = '-'
                                a.reference_dict[seg_id] = OrderedDict([(g.replace('x','?'), v) if g==gn.replace('.','x') else (g, v) for g, v in a.reference_dict[seg_id].items()])
                                a.template_dict[seg_id] = OrderedDict([(g.replace('x','?'), v) if g==gn.replace('.','x') else (g, v) for g, v in a.template_dict[seg_id].items()])
                    out_pdb_array[seg_id] = seg
                main_pdb_array = out_pdb_array
        if self.debug:
            print('Integrate bulges/constrictions: ',datetime.now() - startTime)

        # check for inconsitencies with db
        pdb_db_inconsistencies = []
        for seg_label, segment in a.template_dict.items():
            try:
                for gn, res in segment.items():
                    try:
                        if (res==PDB.Polypeptide.three_to_one(
                                            main_pdb_array[seg_label][gn.replace('x','.')][0].get_parent().get_resname()) and
                                            main_pdb_array[seg_label][gn.replace('x','.')][0].get_parent().get_resname()!='YCM'):
                            pass
                        elif 'x' in gn:
                            try:
                                Residue.objects.get(
                                        protein_conformation__protein=self.main_structure.protein_conformation.protein,
                                        display_generic_number__label=dgn(gn,self.main_structure.protein_conformation))
                                pdb_db_inconsistencies.append({gn:a.template_dict[seg_label][gn]})
                            except:
                                pass
                        else:
                            pass
                    except:
                        pass
            except:
                pass

        if pdb_db_inconsistencies!=[]:
            for incons in pdb_db_inconsistencies:
                seg = self.segment_coding[int(list(incons.keys())[0][0])]
                seq_num = Residue.objects.get(
                                        protein_conformation__protein=self.main_structure.protein_conformation.protein,
                                        display_generic_number__label=dgn(list(incons.keys())[0],self.main_structure.protein_conformation))
                temp_segment, temp_array = OrderedDict(), OrderedDict()
                for key, value in main_pdb_array[seg].items():
                    if key==str(seq_num.sequence_number):
                        temp_segment[list(incons.keys())[0].replace('x','.')] = value
                    else:
                        temp_segment[key] = value
                for seg_id, segment in main_pdb_array.items():
                    if seg_id==seg:
                        temp_array[seg_id] = temp_segment
                    else:
                        temp_array[seg_id] = segment
                main_pdb_array = temp_array
                a.template_dict[seg][list(incons.keys())[0]] = PDB.Polypeptide.three_to_one(
                            main_pdb_array[seg][list(incons.keys())[0].replace('x','.')][0].get_parent().get_resname())
                if a.reference_dict[seg][list(incons.keys())[0]]==a.template_dict[seg][list(incons.keys())[0]]:
                    a.alignment_dict[seg][list(incons.keys())[0]] = a.reference_dict[seg][list(incons.keys())[0]]
                else:
                    a.alignment_dict[seg][list(incons.keys())[0]] = '.'

        for seg in main_pdb_array:
            for gn, atoms in main_pdb_array[seg].items():
                try:
                    if atoms[0].get_parent().get_resname() in ['YCM','CSD','SEP','TYS']:
                        non_ess_res = atoms[0].get_parent().get_resname()
                        if self.debug:
                            print(gn, atoms[0].get_parent().get_resname(), atoms[0].get_parent().get_id())
                        if non_ess_res=='SEP':
                            a.template_dict[seg][gn.replace('.','x')] = 'S'
                        elif non_ess_res=='TYS':
                            a.template_dict[seg][gn.replace('.','x')] = 'Y'
                        else:
                            a.template_dict[seg][gn.replace('.','x')] = 'C'
                            a.alignment_dict[seg][gn.replace('.','x')] = '.'
                except:
                    pass

        self.statistics.add_info('pdb_db_inconsistencies', pdb_db_inconsistencies)

        if not os.path.exists(path):
            os.mkdir(path)

        if self.debug:
            print('Check inconsistencies: {}'.format(pdb_db_inconsistencies),datetime.now() - startTime)

        # inserting loops for free modeling
        for label, template in loop_stat.items():
            if template==None:
                modeling_loops = Loops(self.reference_protein, label, self.similarity_table_all, self.main_structure,
                                       self.helix_end_mods, list(self.template_source), self.revise_xtal)
                modeling_loops.insert_gaps_for_loops_to_arrays(main_pdb_array, a.reference_dict, a.template_dict,
                                                               a.alignment_dict)
                main_pdb_array = modeling_loops.main_pdb_array
                a.reference_dict = modeling_loops.reference_dict
                a.template_dict = modeling_loops.template_dict
                a.alignment_dict = modeling_loops.alignment_dict
        if self.debug:
            print('Free loops: ',datetime.now() - startTime)

        # Adjust H8 if needed
        if 'H8' in main_pdb_array and 'ICL4' not in main_pdb_array and len(self.helix_end_mods['removed']['TM7'][1])>0:
            unwind_num = math.ceil(len(self.helix_end_mods['removed']['TM7'][1])/2)
            trimmed_residues+=list(main_pdb_array['TM7'].keys())[(unwind_num*-1):]+list(main_pdb_array['H8'].keys())[:unwind_num]

        # N- and C-termini
        if N_and_C_termini==True: #and self.prot_conf.protein==self.main_structure.protein_conformation.protein.parent:
            N_struct = self.template_source['TM1'][list(self.template_source['TM1'])[0]][0]
            N_term = Residue.objects.filter(protein_conformation=self.prot_conf, protein_segment__slug='N-term')
            if N_struct!=None:
                N_term_temp = Residue.objects.filter(protein_conformation=N_struct.protein_conformation,
                                                     protein_segment__slug='N-term')
                last_five = [i.sequence_number for i in list(N_term_temp) if i.sequence_number<1000]
                if self.prot_conf.protein!=self.main_structure.protein_conformation.protein.parent:
                    last_five = last_five[-5:]
            else:
                N_term_temp = []
                last_five = []
            if self.main_structure==N_struct and len(self.helix_end_mods['removed']['TM1'][0])==0:
                try:
                    temp_coo = list(parse.fetch_residues_from_pdb(N_struct,last_five).values())
                except:
                    temp_coo = None
            elif len(last_five)==5 or len(self.helix_end_mods['removed']['TM1'][0])>0:
                try:
                    temp_nums = last_five + [i for i in range(last_five[-1]+1,last_five[-1]+5)]
                    template = parse.fetch_residues_from_pdb(N_struct,temp_nums)
                    ref_nums = list(main_pdb_array['TM1'])[:4]
                    reference = OrderedDict()
                    for i in ref_nums:
                        reference[i] = main_pdb_array['TM1'][i]
                    superpose = sp.OneSidedSuperpose(reference,template,4,0)
                    sup_residues = superpose.run()
                    n_count2 = 0
                    temp_coo = []
                    for num, atoms in sup_residues.items():
                        if n_count2<5:
                            temp_coo.append(atoms)
                        n_count2+=1
                except:
                    temp_coo = None
            else:
                temp_coo = None

            prev_num = 0

            # check for gaps in N-term, usual in class B structures
            N_gap = False
            for n in N_term_temp:
                if prev_num!=0 and n.sequence_number-prev_num!=1:
                    N_gap = True
                prev_num = n.sequence_number

            r_i, t_i, a_i, arr_i = OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict()
            N_r, N_t, N_a, N_arr = OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict()
            n_count = 0

            # Gapped N-term residues aren't currently modeled, they are skipped ### FIXME
            if N_gap:
                if self.main_structure==N_struct:
                    for n in N_term:
                        n_count+=1
                        N_r[str(n.sequence_number)] = n.amino_acid
                        n_t = N_term_temp.filter(sequence_number=n.sequence_number)
                        if len(n_t)==1:
                            N_t[str(n.sequence_number)] = n_t[0].amino_acid
                            self.template_source['N-term'][str(n.sequence_number)][0] = N_struct
                            self.template_source['N-term'][str(n.sequence_number)][1] = N_struct
                            N_arr[str(n.sequence_number)] = list(parse.fetch_residues_from_pdb(N_struct,[n.sequence_number]).values())[0]
                            N_a[str(n.sequence_number)] = '-'
                        else:
                            N_t[(str(n.sequence_number))] = '-'
                            self.template_source['N-term'][str(n.sequence_number)][0] = None
                            self.template_source['N-term'][str(n.sequence_number)][1] = None
                            N_arr[str(n.sequence_number)] = '-'
                            N_a[str(n.sequence_number)] = '-'
                            trimmed_residues.append(n.sequence_number)
                else:
                    logger.warning('Gaps in N-term in template {} for {} {} model, skipping N-term'.format(N_struct, self.reference_protein, self.state))
                    del self.template_source['N-term']
                    if self.debug:
                        print('Gaps in N-term in template {} for {} {} model, skipping N-term'.format(N_struct, self.reference_protein, self.state))
            else:
                for n in N_term:
                    n_count+=1
                    N_r[str(n.sequence_number)] = n.amino_acid
                    N_a[str(n.sequence_number)] = '-'

                    try:
                        N_arr[str(n.sequence_number)] = temp_coo[-1*(len(N_term)-n_count+1)]
                        N_t[str(n.sequence_number)] = list(N_term_temp)[-1*(len(N_term)-n_count+1)].amino_acid
                        self.template_source['N-term'][str(n.sequence_number)][0] = N_struct
                        self.template_source['N-term'][str(n.sequence_number)][1] = N_struct
                    except:
                        N_t[str(n.sequence_number)] = '-'
                        N_arr[str(n.sequence_number)] = '-'

            r_i['N-term'] = N_r
            t_i['N-term'] = N_t
            a_i['N-term'] = N_a
            arr_i['N-term'] = N_arr
            for r,t,al,arr in zip(a.reference_dict,a.template_dict,a.alignment_dict,main_pdb_array):
                r_i[r]=a.reference_dict[r]
                t_i[t]=a.template_dict[t]
                a_i[al]=a.alignment_dict[al]
                arr_i[arr]=main_pdb_array[arr]
            a.reference_dict = r_i
            a.template_dict = t_i
            a.alignment_dict = a_i
            main_pdb_array = arr_i

            # N-term alignment sanity check
            # for i,j,k,l in zip(a.reference_dict['N-term'], a.template_dict['N-term'], a.alignment_dict['N-term'],main_pdb_array['N-term']):
            #     print(i,j,k,l, a.reference_dict['N-term'][i], a.template_dict['N-term'][j],a.alignment_dict['N-term'][k],main_pdb_array['N-term'][l])

            try:
                index = -1
                while self.template_source['H8'][list(self.template_source['H8'])[index]][0]==None:
                    index-=1
                C_struct = self.template_source['H8'][list(self.template_source['H8'])[index]][0]
                last_seg = 'H8'
            except:
                C_struct = self.template_source['TM7'][list(self.template_source['TM7'])[-1]][0]
                last_seg = 'TM7'
            C_term = Residue.objects.filter(protein_conformation=self.prot_conf, protein_segment__slug='C-term')
            C_term_temp = Residue.objects.filter(protein_conformation=C_struct.protein_conformation,
                                                 protein_segment__slug='C-term')

            first_five = [i.sequence_number for i in list(C_term_temp) if i.sequence_number<1000]
            if self.main_structure==C_struct:
                try:
                    temp_coo2 = list(parse.fetch_residues_from_pdb(C_struct,first_five).values())
                except:
                    temp_coo2 = None
            elif len(first_five)==5:
                try:
                    temp_nums2 = [i for i in range(first_five[0]-4,first_five[0])] + first_five
                    template2 = parse.fetch_residues_from_array(C_struct,temp_nums2)
                    ref_nums2 = list(main_pdb_array[last_seg])[-4:]
                    reference2 = OrderedDict()
                    for i in ref_nums2:
                        reference2[i] = main_pdb_array[last_seg][i]
                    superpose2 = sp.OneSidedSuperpose(reference2,template2,4,1)
                    sup_residues2 = superpose2.run()
                    c_count2 = 0
                    temp_coo2 = []
                    for num, atoms in sup_residues2.items():
                        if c_count2<5:
                            temp_coo2.append(atoms)
                        c_count2+=1
                except:
                    temp_coo2 = None
            else:
                temp_coo2 = None

            a.reference_dict['C-term'],a.template_dict['C-term'] = OrderedDict(),OrderedDict()
            a.alignment_dict['C-term'],main_pdb_array['C-term'] = OrderedDict(),OrderedDict()
            c_count = -1
            for c in C_term:
                c_count+=1
                a.reference_dict['C-term'][str(c.sequence_number)] = c.amino_acid
                a.alignment_dict['C-term'][str(c.sequence_number)] = '-'
                try:
                    if self.revise_xtal and self.main_structure.pdb_code.index=='5F8U' and c.sequence_number==359:
                        raise Exception()
                    main_pdb_array['C-term'][str(c.sequence_number)] = temp_coo2[c_count]
                    a.template_dict['C-term'][str(c.sequence_number)] = list(C_term_temp)[c_count].amino_acid
                    self.template_source['C-term'][str(c.sequence_number)][0] = C_struct
                    self.template_source['C-term'][str(c.sequence_number)][1] = C_struct
                except:
                    a.template_dict['C-term'][str(c.sequence_number)] = '-'
                    main_pdb_array['C-term'][str(c.sequence_number)] = '-'

            # Shorten N- and C-termini
            n_count=1
            delete_termini = set()
            # n_term_no_coord = [i for i in a.template_dict['N-term'] if a.template_dict['N-term'][i]=='-']
            # n_term_coord = [i for i in a.template_dict['N-term'] if a.template_dict['N-term'][i]!='-']
            for num in list(a.template_dict['N-term'])[:-5]:
                if a.template_dict['N-term'][num]=='-':
                    delete_termini.add(('N-term', num))
                else:
                    break
                n_count+=1

            c_count=1
            for num in list(a.template_dict['C-term'])[5:]:
                if a.template_dict['C-term'][num]=='-':
                    delete_termini.add(('C-term', num))
                c_count+=1

            for i, ii in delete_termini:
                del a.reference_dict[i][ii]
                del a.template_dict[i][ii]
                del a.alignment_dict[i][ii]
                del main_pdb_array[i][ii]
                del self.template_source[i][ii]

            if len(a.reference_dict['N-term'])==0:
                del a.reference_dict['N-term']
                del a.template_dict['N-term']
                del a.alignment_dict['N-term']
                del main_pdb_array['N-term']

        # Shorten ICL3
        for i in a.reference_dict:
            if i.startswith('ICL3'):
                label = i
                break
        try:
            if '_dis' not in label and len(a.reference_dict[label])>10:
                delete_ts, delete_r, delete_t, delete_a, delete_m = set(),set(),set(),set(),set()
                chain_break = False
                icl3_c = 0
                keys = list(self.template_source['ICL3'].keys())
                length = len(a.template_dict[label])
                if self.revise_xtal==True:
                    ref_prot = self.reference_protein.parent
                else:
                    ref_prot = self.reference_protein
                for r_s,t_s,a_s,ar_s in zip(a.reference_dict[label],a.template_dict[label],
                                            a.alignment_dict[label],main_pdb_array[label]):
                    icl3_c+=1
                    if 5<icl3_c<length-4:
                        if chain_break==False and self.main_structure.protein_conformation.protein.parent==ref_prot:
                            a.reference_dict[label][r_s] = '/'
                            a.template_dict[label][t_s] = '/'
                            a.alignment_dict[label][a_s] = '/'
                            main_pdb_array[label][ar_s] = '/'
                            delete_ts.add(('ICL3',keys[icl3_c-1]))
                            chain_break = True
                        else:
                            delete_r.add((label,r_s))
                            delete_t.add((label,t_s))
                            delete_a.add((label,a_s))
                            delete_m.add((label,ar_s))
                            delete_ts.add(('ICL3',keys[icl3_c-1]))
                self.icl3_delete['R'] = []
                for i,ii in delete_ts:
                    self.icl3_delete['R'].append(int(ii))
                    del self.template_source[i][ii]
                for i,ii in delete_r:
                    del a.reference_dict[i][ii]
                for i,ii in delete_t:
                    del a.template_dict[i][ii]
                for i,ii in delete_a:
                    del a.alignment_dict[i][ii]
                for i,ii in delete_m:
                    del main_pdb_array[i][ii]
        except:
            pass

        # introducing mutations
        if self.mutations:
            for i, j in self.mutations.items():
                for seg, resis in a.reference_dict.items():
                    for gn, res in resis.items():
                        if i==gn:
                            a.reference_dict[seg][gn] = j
                            break

        # non-conserved residue switching
        if switch_rotamers==True:
            non_cons_switch = self.run_non_conserved_switcher(main_pdb_array,a.reference_dict,a.template_dict,
                                                              a.alignment_dict)
            main_pdb_array = non_cons_switch[0]
            a.reference_dict = non_cons_switch[1]
            a.template_dict = non_cons_switch[2]
            a.alignment_dict = non_cons_switch[3]
            trimmed_residues+=non_cons_switch[4]
        else:
            for seg_id, seg in main_pdb_array.items():
                for key in seg:
                    if a.reference_dict[seg_id][str(key).replace('.','x')]!='-':
                        trimmed_residues.append(key)

        if 'ICL4_free' in main_pdb_array:
            freeICL4=True
        else:
            freeICL4=False
        if freeICL4==True:
            for i in list(main_pdb_array['H8']):
                if i not in trimmed_residues:
                    trimmed_residues.append(i)
        array_keys = list(main_pdb_array.keys())

        for i,j in self.helix_end_mods['added'].items():
            try:
                if j[0][-1].replace('x','.') not in trimmed_residues:
                    trimmed_residues.append(j[0][-1].replace('x','.'))
            except:
                pass
            try:
                if j[1][0].replace('x','.') not in trimmed_residues:
                    trimmed_residues.append(j[1][0].replace('x','.'))
            except:
                pass
            try:
                if j[0][0].replace('x','.') not in trimmed_residues and '_cont' in array_keys[array_keys.index(i)-1]:
                    trimmed_residues.append(j[0][0].replace('x','.'))
            except:
                pass
            try:
                if j[1][-1].replace('x','.') not in trimmed_residues and array_keys[array_keys.index(i)+1]+'_cont' in array_keys:
                    trimmed_residues.append(j[1][-1].replace('x','.'))
            except:
                pass
            try:
                trimmed_residues.append(parse.gn_indecer(j[0][-1],'x',1).replace('x','.'))
            except:
                pass
            try:
                trimmed_residues.append(parse.gn_indecer(j[1][0],'x',-1).replace('x','.'))
            except:
                pass
        for i in ref_bulge_list+temp_bulge_list+ref_const_list+temp_const_list:
            i = list(i.keys())[0].replace('x','.')
            if parse.gn_indecer(i,'.',-2) not in trimmed_residues:
                trimmed_residues.append(parse.gn_indecer(i,'.',-2))
            if parse.gn_indecer(i,'.',-1) not in trimmed_residues:
                trimmed_residues.append(parse.gn_indecer(i,'.',-1))
            if parse.gn_indecer(i,'.',1) not in trimmed_residues:
                trimmed_residues.append(parse.gn_indecer(i,'.',1))
            if parse.gn_indecer(i,'.',2) not in trimmed_residues:
                trimmed_residues.append(parse.gn_indecer(i,'.',2))
        for s in a.reference_dict:
            if 'dis' in s:
                key_list=list(a.reference_dict[s].keys())
                i, j = 0, -1
                found_first, found_second = False, False
                while found_first==False:
                    if a.reference_dict[s][key_list[i]]!='x':

                        found_first=True
                    i+=1
                while found_second==False:
                    if a.reference_dict[s][key_list[j]]!='x':
                        trimmed_residues.append(key_list[j].replace('x','.'))
                        found_second=True
                    j-=1
                segs = list(a.reference_dict.keys())
                if i>0:
                    trimmed_residues.append(list(a.reference_dict[segs[segs.index(s)-1]].keys())[-1].replace('x','.'))
                if j<-1:
                    trimmed_residues.append(list(a.reference_dict[segs[segs.index(s)+1]].keys())[0].replace('x','.'))
        if self.reference_entry_name.startswith('taar') and str(self.main_structure)=='4IAR':
            trimmed_residues.append('5.36')

        if self.debug:
            print('Rotamer switching: ',datetime.now() - startTime)

        for i in model_loops:
            for j in a.reference_dict[i]:
                trimmed_residues.append(j.replace('x','.'))
        if self.reference_protein.family.slug.startswith('004'):
            for i in a.template_dict['H8']:
                trimmed_residues.append(i.replace('x','.'))
        # Add last res of TM7 to model when start of H8 was added from alt template
        if 'H8' in self.template_source and self.template_source['H8'][list(self.template_source['H8'])[0]][0]!=self.main_structure:
            trimmed_residues.append(list(main_pdb_array['TM7'])[-1])

        self.statistics.add_info('trimmed_residues', trimmed_residues)

        # check if ECL3 might have a disulfide bridge
        self.disulfide_pairs.append(self.ECL3_disulfide(a.reference_dict))
        self.statistics.add_info('disulfide_pairs', self.disulfide_pairs)

        # building complex model, adding signaling protein and structure data to dictionaries
        if self.complex:
            cm = SignprotModeling(self.main_structure, self.signprot, self.template_source, trimmed_residues, a, main_pdb_array, self.debug)
            cm.run()
            self.template_source = cm.template_source
            trimmed_residues = cm.trimmed_residues
            a = cm.a
            main_pdb_array = cm.main_pdb_array
            self.target_signprot = cm.target_signprot
            self.signprot_protconf = cm.signprot_protconf
            self.signprot_complex = cm.signprot_complex
        ######## end of complex modeling
        self.alignment = a
        self.main_pdb_array = main_pdb_array
        self.trimmed_residues = trimmed_residues

        return self

    def build_homology_model_second_part(self):
        # write to file
        if self.complex:
            post_file = '{}{}_{}_post.pdb'.format(path, self.reference_entry_name, self.target_signprot)
        else:
            post_file = path+self.reference_entry_name+'_'+self.state+"_post.pdb"

        ## NMUR2
        # for i in ['ECL2|2','ECL2|3','ECL2|4','ECL2|5','ECL2|6','ECL2|7','ECL2|8','ECL2|9','ECL2|10','ECL2|11','ECL2|12','ECL2|13','ECL2|14','ECL2|15','ECL2|16','ECL2|17','ECL2|23','ECL2|24','ECL2|25']:
        #     self.trimmed_residues.remove(i)
        ## NPY1R
        # for i in [ 'ECL2|2','ECL2|3','ECL2|4','ECL2|5','ECL2|6','ECL2|7','ECL2|8','ECL2|9','ECL2|10','ECL2|11','ECL2|12','ECL2|13','ECL2|14','ECL2|15','ECL2|17','ECL2|18','ECL2|19','ECL2|25','ECL3|2','ECL3|3','ECL3|4']:
        #     self.trimmed_residues.remove(i)

        trimmed_res_nums, helix_restraints, icl3_mid, disulfide_nums, complex_start, beta_start, gamma_start = self.write_homology_model_pdb(post_file,
                                                                                                                    self.main_pdb_array, self.alignment,
                                                                                                                    trimmed_residues=self.trimmed_residues,
                                                                                                                    disulfide_pairs=self.disulfide_pairs, complex=self.complex)
        self.statistics.add_info('template_source',self.template_source)

        # Adding HETATMs when revising xtal or keeping hetatms
        hetatm_count = 0
        water_count = 0
        if self.revise_xtal==True:
            ref_prot = self.reference_protein.parent
        else:
            ref_prot = self.reference_protein
        if (self.revise_xtal and ref_prot==self.main_structure.protein_conformation.protein.parent) or self.keep_hetatoms:
            pdb = PDB.PDBList()
            self.alternate_water_positions = OrderedDict()
            last_seqnum = list(Residue.objects.filter(protein_conformation__protein=self.reference_protein.parent))[-1].sequence_number
            with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'pdbs', str(self.main_structure)+'.pdb']), 'r') as f:
                lines = f.readlines()
            with open(post_file, 'a') as model:
                hetatm = 1
                for line in lines:
                    if line.startswith('HETATM'):
                        if 'YCM' in line or 'CSD' in line or 'SEP' in line or 'TYS' in line or 'ACE' in line:
                            continue
                        if self.main_structure.pdb_code.index=='4BUO' and 'GLY' in line:
                            continue
                        pref_chain = str(self.main_structure.preferred_chain)
                        if len(pref_chain)>1:
                            pref_chain = pref_chain[0]
                        pdb_re = re.search('(HETATM[0-9\sA-Z{apo}]{{11}})([A-Z0-9\s]{{3}})\s({pref})([0-9\s]{{4}})'.format(apo="'",pref=pref_chain), line)
                        if not pdb_re:
                            continue
                        if pdb_re.group(2)!='HOH':
                            if hetatm!=pdb_re.group(4):
                                hetatm_count+=1
                                hetatm = pdb_re.group(4)
                                ### Hetatm seqnums are too low, same as normal residues - FIXME
                                # print("HETATM STUFF")
                                # print(int(hetatm), last_seqnum)
                                # if int(hetatm)<=last_seqnum:
                                #     num_hetatm = last_seqnum+hetatm_count
                                #     if len(hetatm)>len(str(num_hetatm)):
                                #         num_hetatm = ' '*(len(hetatm)-len(str(num_hetatm)))+str(num_hetatm)
                                #     line.replace(hetatm, num_hetatm)
                                #     print(num_hetatm)
                                #     print(line)
                        else:
                            if pdb_re.group(1)[-1]==' ' or pdb_re.group(1)[-1]==pref_chain:
                                water_count+=1
                            elif pdb_re.group(1)[-1] in ['B','C','D']:
                                self.alternate_water_positions[water_count] = line
                        if pdb_re!=None:
                            model.write(line)
                try:
                    ligand_pep_structs = LigandPeptideStructure.objects.filter(structure=self.main_structure)
                    for l in ligand_pep_structs:
                        model.write('TER\n')
                        ### Skip peptide ligands that are not present in structure
                        if l.chain==' ':
                            continue
                        if self.debug:
                            print(l)
                        pdbdata = ''
                        # pdbdata = '\n'.join([i for i in self.main_structure.pdb_data.pdb.split('\n') if i.startswith('ATOM') or i.startswith('HETATM') and len(i)>20 and i[21]==l.chain])
                        for i in self.main_structure.pdb_data.pdb.split('\n'):
                            if 'HOH' in i:
                                continue
                            if (i.startswith('ATOM') or i.startswith('HETATM')) and len(i)>20 and i[21]==l.chain:
                                pdbdata+=i+'\n'
                        model.write(pdbdata+'\n')
                        pep = PDB.PDBParser().get_structure('pep', StringIO(pdbdata))[0]
                        pep_resnums = []
                        pep_seq = ''
                        ppb = PDB.PPBuilder()
                        for pp in ppb.build_peptides(pep, aa_only=False):
                            pep_seq = pp.get_sequence()
                        for r in pep[l.chain]:
                            pep_resnums.append(r.get_id()[1])
                        if len(pep_seq)!=len(pep_resnums):
                            pep_seq = ''
                            for r in pep[l.chain]:
                                try:
                                    pep_seq+=PDB.Polypeptide.three_to_one(r.get_resname())
                                except KeyError:
                                    pep_seq+='X'
                        self.peptide_ligands.append([l, pdbdata, pep, pep_resnums, pep_seq])
                except LigandPeptideStructure.DoesNotExist:
                    pass
                model.write('TER\n')
                model.write('END')
        print(self.peptide_ligands)
        # correcting for side chain clashes
        p = PDB.PDBParser()
        post_model = p.get_structure('post', post_file)[0]
        hse = HSExposureCB(post_model, radius=11)
        clash_pairs = hse.clash_pairs
        trimmed_res_nums['clashes'] = OrderedDict()
        for i in clash_pairs:
            trimmed_res_nums['clashes'][i[0][1]] = i[0][1]
            trimmed_res_nums['clashes'][i[1][1]] = i[1][1]

        self.statistics.add_info('clashing_residues', clash_pairs)

        # Check improved sequence identity
        self.identicals = 0
        counter = 0
        for r_s, t_s in zip(self.alignment.reference_dict, self.alignment.template_dict):
            for r, t in zip(self.alignment.reference_dict[r_s], self.alignment.template_dict[t_s]):
                if self.alignment.reference_dict[r_s][r]==self.alignment.template_dict[t_s][t]:
                    self.identicals+=1
                counter+=1
        # print(self.state, counter, self.identicals)

        ### Alignment sanity check
        # if self.debug:
        #     for i,j,k in zip(self.alignment.reference_dict, self.alignment.template_dict, self.main_pdb_array):
        #         print(i)
        #         pprint.pprint(self.alignment.reference_dict[i])
        #         pprint.pprint(self.alignment.template_dict[j])
        #         pprint.pprint(self.alignment.alignment_dict[i])
        #         pprint.pprint(self.main_pdb_array[k])

        # Model with MODELLER
        self.create_PIR_file(self.alignment.reference_dict, self.alignment.template_dict, post_file, hetatm_count, water_count)

        if self.revise_xtal==False:
            if self.complex:
                self.modelname = "{}_{}-{}_{}_{}_GPCRdb".format(self.class_name, self.reference_entry_name, self.target_signprot.entry_name, self.main_structure, build_date)
            else:
                self.modelname = "{}_{}_{}_{}_{}_GPCRdb".format(self.class_name, self.reference_entry_name, self.state, self.main_structure, build_date)
        else:
            if self.complex:
                self.modelname = "{}_{}-{}_{}_refined_{}_GPCRdb".format(self.class_name, self.reference_protein.parent.entry_name, self.target_signprot.entry_name, self.main_structure, build_date)
            else:
                self.modelname = "{}_{}_{}_refined_{}_{}_GPCRdb".format(self.class_name, self.reference_protein.parent.entry_name, self.main_structure, self.main_structure.state.name, build_date)
            try:
                os.remove('./pdb{}.ent'.format(self.reference_entry_name))
            except:
                pass

        # Ignore output from modeller!
        if not self.debug:
            _stdout = sys.stdout
            sys.stdout = open(os.devnull, 'w')

        if self.complex:
            pir_file = "./structure/PIR/{}_{}.pir".format(self.uniprot_id, self.target_signprot.entry_name)
        else:
            pir_file = "./structure/PIR/"+self.uniprot_id+"_"+self.state+".pir"

        if self.alphafold and len([trimmed_res_nums[i] for i in trimmed_res_nums if len(trimmed_res_nums[i])>0])==0:
            trimmed_res_nums['N-term']['1'] = 1

        if self.debug:
            print('Atom selection for MODELLER:')
            print(trimmed_res_nums)

        self.run_MODELLER(pir_file, post_file,
                          self.uniprot_id, self.modeller_iterations, path+self.modelname+'.pdb',
                          atom_dict=trimmed_res_nums, helix_restraints=helix_restraints, icl3_mid=icl3_mid, disulfide_nums=disulfide_nums, complex_start=complex_start, beta_start=beta_start,
                          gamma_start=gamma_start, peptide_ligands=self.peptide_ligands)
        # Resume output
        if not self.debug:
            sys.stdout.close()
            sys.stdout = _stdout

        if not self.debug:
            os.remove(post_file)

        # stat file
        with open(path+self.modelname+'.templates.csv','w') as s_file:
            rot_table = []
            s_file.write('Segment,Sequence_number,Generic_number,Backbone_template,Rotamer_template\n')
            for seg, resis in self.template_source.items():
                list_keys = list(resis)
                if len(list_keys)==0:
                    continue
                first_gn = list_keys[0]
                first_temp = self.template_source[seg][first_gn][0]
                if 'x' in first_gn:
                    try:
                        first_seqnum = Residue.objects.get(protein_conformation=self.prot_conf,display_generic_number__label=dgn(list_keys[0],self.prot_conf)).sequence_number
                    except:
                        try:
                            first_seqnum = int(list_keys[0])
                        except:
                            continue
                elif self.complex and first_gn!=None and len(first_gn.split('.'))==3:
                    first_seqnum = Residue.objects.get(protein_conformation=self.signprot_protconf,display_generic_number__label=list_keys[0]).sequence_number
                else:
                    first_seqnum = int(first_gn)
                    first_gn = None
                for gn, res in resis.items():
                    key = gn
                    if 'x' in gn:
                        seq_num = Residue.objects.get(protein_conformation=self.prot_conf,display_generic_number__label=dgn(gn,self.prot_conf)).sequence_number
                        curr_seqnum = seq_num
                    elif self.complex and first_gn!=None and len(first_gn.split('.'))==3:
                        seq_num = Residue.objects.get(protein_conformation=self.signprot_protconf,display_generic_number__label=gn).sequence_number
                        curr_seqnum = seq_num
                    else:
                        seq_num = int(gn)
                        curr_seqnum = seq_num
                        gn = None

                    if seg in gprotein_segment_slugs:
                        rot_table.append([seg,seq_num,gn,res[0],res[1]])
                    elif seg=='Beta':
                        rot_table.append([seg,seq_num,gn,res[0],res[1]])
                    elif seg=='Gamma':
                        rot_table.append([seg,seq_num,gn,res[0],res[1]])
                    else:
                        rot_table.append([seg,seq_num,gn,res[0],res[1]])

            for rot in rot_table:
                try:
                    bb = rot[3].pdb_code.index
                except:
                    bb = rot[3]
                try:
                    rt = rot[4].pdb_code.index
                except:
                    rt = rot[4]
                l = "{},{},{},{},{}\n".format(rot[0],rot[1],rot[2],bb,rt)
                s_file.write(l)

        self.template_stats = rot_table


        # template seq sim file
        with open(path+self.modelname+'.template_similarities.csv','w') as s_file:
            template_list, self.template_list = [], OrderedDict()
            for r in self.template_stats:
                if r[3] not in template_list and r[3]!=None:
                    template_list.append(r[4])
                if r[4] not in template_list and r[4]!=None:
                    template_list.append(r[4])
            s_file.write('Template,Sequence_similarity,Resolution,Representative,State\n')
            for temp, sim in self.similarity_table_all.items():
                if temp in template_list:
                    self.template_list[temp] = sim
            for t, s in self.template_list.items():
                s_file.write('{},{},{},{},{}\n'.format(t.pdb_code.index, s, t.resolution, t.representative, t.state.slug))


        if self.debug:
            print('MODELLER build: ',datetime.now() - startTime)
            pprint.pprint(self.statistics)
            print('################################')
        return self

    def run_non_conserved_switcher(self, main_pdb_array, reference_dict, template_dict, alignment_dict):
        ''' Switches non-conserved residues with best possible template. Returns refreshed main_pdb_array
            (atom coordinates), reference_dict (reference generic numbers and residue ids), template_dict (template
            generic numbers and residue ids) and alignment_dict (aligned reference and template dictionary).

            @param main_pdb_array: nested OrderedDict(), output of GPCRDBParsingPDB().pdb_array_creator()
            @param reference_dict: reference dictionary of AlignedReferenceTemplate.
            @param template_dict: template dictionary of AlignedReferenceTemplate.o2
            @param alignment_dict: alignment dictionary of AlignedReferenceTemplate.
        '''
        atom_num_dict = {'E':9, 'S':6, 'Y':12, 'G':4, 'A':5, 'V':7, 'M':8, 'L':8, 'I':8, 'T':7, 'F':11, 'H':10, 'K':9,
                         'D':8, 'C':6, 'R':11, 'P':7, 'Q':9, 'N':8, 'W':14, '-':0}
        parse = GPCRDBParsingPDB()
        ref_length = 0
        conserved_count = 0
        non_cons_count = 0
        trimmed_res_num = 0
        switched_count = 0
        non_cons_res_templates, conserved_residues = OrderedDict(), OrderedDict()
        trimmed_residues = []
        inconsistencies = []

        ###
        ## Default Class A set
        self.same_state_residues = ['1x49','1x50','1x53','12x50','12x51','2x37','2x39','2x40','2x42','2x43','2x47','2x50','3x40','3x41','3x43','3x46','3x49','3x50','4x49','5x50','5x51','5x54','5x58','6x34','6x37','6x38','6x40','6x41','6x44','6x47','6x48','7x40','7x46','7x49','7x50','7x52','7x53','7x54','8x49','8x50','8x51',]
        ## Extended Class A set with TM5e, TM6b, TM7e, H8b
        self.same_state_residues = ['1x49','1x50','1x53','12x50','12x51','2x37','2x39','2x40','2x42','2x43','2x47','2x50','3x40','3x41','3x43','3x46','3x49','3x50','4x49','5x50','5x51','5x54','5x57','5x58','5x59','5x60','5x61','5x62','5x63','5x64','5x65','5x66','5x67','5x68','5x69','5x70','5x71','5x72','5x73','6x23','6x24','6x25','6x26','6x27','6x28','6x29','6x30','6x31','6x32','6x33','6x34','6x35','6x36','6x37','6x38','6x39','6x40','6x41','6x44','6x47','6x48','7x40','7x46','7x49','7x50','7x52','7x53','7x54','7x56','8x47','8x48','8x49','8x50','8x51',]
        ###

        if self.revise_xtal==True:
            ref_prot = self.reference_protein.parent
        else:
            ref_prot = self.reference_protein
        for incons in self.statistics.info_dict['pdb_db_inconsistencies']:
            inconsistencies.append(list(incons.keys())[0])
        for ref_seg, temp_seg, aligned_seg in zip(reference_dict, template_dict, alignment_dict):
            if len(ref_seg)>4:
                segment = ref_seg[:4]
            else:
                segment = ref_seg
            for ref_res, temp_res, aligned_res in zip(reference_dict[ref_seg], template_dict[temp_seg],
                                                      alignment_dict[aligned_seg]):
                ### State change
                # if temp_res not in self.same_state_residues and 'x' in temp_res:
                #     alignment_dict[aligned_seg][aligned_res] = '.'
                #     # template_dict[temp_seg][temp_res] = '.'
                #     # inconsistencies.append(temp_res)
                #     main_pdb_array[ref_seg][ref_res.replace('x','.')] = main_pdb_array[ref_seg][ref_res.replace('x','.')][:5]
                ###

                if self.revise_xtal==True and reference_dict[ref_seg][ref_res]!=template_dict[temp_seg][temp_res]:
                    alignment_dict[aligned_seg][aligned_res]='.'
                if template_dict[temp_seg][temp_res]=='/':
                    continue
                if reference_dict[ref_seg][ref_res]!='-':
                    ref_length+=1
                else:
                    trimmed_residues.append(ref_res.replace('x','.'))
                if '?' in temp_res:
                    trimmed_residues.append(ref_res.replace('x','.'))
                    trimmed_res_num+=1
                    non_cons_count+=1
                    continue
                if '-term' in ref_seg and (template_dict[temp_seg][temp_res]=='-' or
                                           reference_dict[ref_seg][ref_res]!=template_dict[temp_seg][temp_res] or
                                           len(main_pdb_array[ref_seg][ref_res])<atom_num_dict[template_dict[temp_seg][temp_res]]):
                    trimmed_residues.append(ref_res.replace('x','.'))
                    trimmed_res_num+=1
                    non_cons_count+=1
                    continue
                if (ref_res not in inconsistencies and
                    alignment_dict[aligned_seg][aligned_res]!='.' and
                    alignment_dict[aligned_seg][aligned_res]!='x' and
                    alignment_dict[aligned_seg][aligned_res]!='-' and
                    alignment_dict[aligned_seg][aligned_res]!='/' and
                    len(main_pdb_array[ref_seg][ref_res.replace('x','.')])>=atom_num_dict[template_dict[temp_seg][temp_res]]):
                    try:
                        rot_test = Rotamer.objects.filter(structure=self.main_structure,
                                                          residue__display_generic_number__label=dgn(ref_res,
                                                          self.main_structure.protein_conformation))
                        rot_test = self.right_rotamer_select(rot_test)
                        if rot_test.missing_atoms==True:
                            alignment_dict[aligned_seg][aligned_res]='.'
                            template_dict[temp_seg][temp_res]='G'
                        else:
                            raise Exception()
                    except:
                        conserved_residues[ref_res] = alignment_dict[aligned_seg][aligned_res]
                        conserved_count+=1
                        if 'x' not in ref_res:
                            num_in_loop = parse.gn_num_extract(ref_res,'|')[1]
                            try:
                                this_res = list(Residue.objects.filter(protein_conformation=self.prot_conf,
                                                                       protein_segment__slug=segment))[num_in_loop-1]
                            except:
                                trimmed_residues.append(ref_res.replace('x','.'))
                                continue
                            seq_num = str(this_res.sequence_number)
                            try:
                                self.template_source = update_template_source(self.template_source,[seq_num],self.template_source[segment][seq_num][0],segment,
                                                            just_rot=True)
                                key_in_template_source = seq_num
                            except:
                                self.template_source = update_template_source(self.template_source,[ggn(this_res.display_generic_number.label)],
                                                            self.template_source[segment][ggn(this_res.display_generic_number.label)][0],
                                                            segment,just_rot=True)
                                key_in_template_source = ggn(this_res.display_generic_number.label)
                        else:
                            try:
                                self.template_source = update_template_source(self.template_source,[ref_res],self.template_source[segment][ref_res][0],segment,
                                                            just_rot=True)
                                key_in_template_source = ref_res
                            except:
                                missing_i = list(reference_dict[ref_seg].keys()).index(ref_res)
                                gaps_before = [x for x in list(reference_dict[ref_seg].keys())[:missing_i] if reference_dict[ref_seg][x]=='-']
                                this_loop = Residue.objects.filter(protein_conformation__protein=self.reference_protein, protein_segment__slug=ref_seg[:4])
                                right_res = str(this_loop[missing_i-len(gaps_before)].sequence_number)
                                self.template_source = update_template_source(self.template_source,[right_res],self.template_source[segment][right_res][0],segment,
                                                            just_rot=True)
                                key_in_template_source = right_res
                        if '_dis' in ref_seg or (ref_seg=='ECL2' and self.template_source['ECL2'][key_in_template_source][0]!=self.main_structure
                                                 and '|' in ref_res):
                            trimmed_residues.append(ref_res.replace('x','.'))
                gn = ref_res

                if ((gn in inconsistencies or alignment_dict[aligned_seg][aligned_res]=='.' and
                    reference_dict[ref_seg][gn]!=template_dict[temp_seg][gn]) or (template_dict[temp_seg][temp_res]!='x' and
                    len(main_pdb_array[ref_seg][ref_res.replace('x','.')])<atom_num_dict[template_dict[temp_seg][temp_res]])):
                    print(gn, alignment_dict[aligned_seg][aligned_res], template_dict[temp_seg][gn])
                ###########
                # if (gn in inconsistencies) or (template_dict[temp_seg][temp_res]!='x' and
                #     len(main_pdb_array[ref_seg][ref_res.replace('x','.')])<atom_num_dict[template_dict[temp_seg][temp_res]]) or (reference_dict[ref_seg][gn]==template_dict[temp_seg][gn] or
                #     gn in inconsistencies or alignment_dict[aligned_seg][aligned_res]=='.'):
                ###########
                    non_cons_count+=1
                    gn_ = str(ref_res).replace('x','.')
                    no_match = True
                    if '|' in gn_:
                        try:
                            list_num = int(gn.split('|')[1])-1
                            gn = ggn(list(Residue.objects.filter(protein_conformation__protein=ref_prot,
                                         protein_segment__slug=ref_seg.split('_')[0]))[list_num].display_generic_number.label)
                            gn_ = gn.replace('x','.')
                        except:
                            pass
                    this_state_structs_with_resi, other_state_structs_with_resi = [],[]

                    # if gn in self.same_state_residues:
                    #     main_pdb_array, template_dict, non_cons_res_templates, switched_count, no_match = self.find_and_switch_rotamer(self.similarity_table, gn, gn_,
                    #                 reference_dict, ref_seg, ref_res, main_pdb_array, atom_num_dict, template_dict, temp_seg, temp_res, non_cons_res_templates, switched_count, no_match, segment)
                    # else:
                    main_pdb_array, template_dict, non_cons_res_templates, switched_count, no_match = self.find_and_switch_rotamer(self.similarity_table, gn, gn_,
                                reference_dict, ref_seg, ref_res, main_pdb_array, atom_num_dict, template_dict, temp_seg, temp_res, non_cons_res_templates, switched_count, no_match, segment)
                    if no_match==True:
                        main_pdb_array, template_dict, non_cons_res_templates, switched_count, no_match = self.find_and_switch_rotamer(self.similarity_table_other_states, gn, gn_,
                                reference_dict, ref_seg, ref_res, main_pdb_array, atom_num_dict, template_dict, temp_seg, temp_res, non_cons_res_templates, switched_count, no_match, segment)
                        if no_match==True:
                            print('No rotamer match: ', gn, reference_dict[ref_seg][gn])
                            try:
                                if 'free' not in ref_seg:
                                    residue = main_pdb_array[ref_seg][str(ref_res).replace('x','.')]
                                    main_pdb_array[ref_seg][str(ref_res).replace('x','.')] = residue[0:5]
                                    trimmed_residues.append(gn_)
                                    trimmed_res_num+=1
                                    if '_' in ref_seg:
                                        orig_ref_seg = ref_seg.split('_')[0]
                                    else:
                                        orig_ref_seg = ref_seg
                                    if ref_res in self.template_source[orig_ref_seg]:
                                        self.template_source = update_template_source(self.template_source, [ref_res], None, orig_ref_seg, just_rot=True)
                                    else:
                                        fetched_key = list(self.template_source[orig_ref_seg].keys())[int(ref_res.split('|')[1])-1]
                                        self.template_source = update_template_source(self.template_source,[fetched_key], None, orig_ref_seg, just_rot=True)
                                elif 'free' in ref_seg:
                                    trimmed_residues.append(gn_)
                                    trimmed_res_num+=1
                            except:
                                print("Missing atoms in {} at {}".format(self.main_structure,gn))
                                logging.warning("Missing atoms in {} at {}".format(self.main_structure,gn))
                elif alignment_dict[aligned_seg][aligned_res]=='x':
                    trimmed_residues.append(gn.replace('x','.'))
                    trimmed_res_num+=1
        self.statistics.add_info('ref_seq_length', ref_length)
        self.statistics.add_info('conserved_num', conserved_count)
        self.statistics.add_info('non_conserved_num', non_cons_count)
        self.statistics.add_info('trimmed_residues_num', trimmed_res_num)
        self.statistics.add_info('non_conserved_switched_num', switched_count)
        self.statistics.add_info('conserved_residues', conserved_residues)
        self.statistics.add_info('non_conserved_residue_templates', non_cons_res_templates)

        return [main_pdb_array, reference_dict, template_dict, alignment_dict, trimmed_residues]

    def find_and_switch_rotamer(self, similarity_table, gn, gn_, reference_dict, ref_seg, ref_res, main_pdb_array, atom_num_dict,
                                template_dict, temp_seg, temp_res, non_cons_res_templates, switched_count, no_match, segment):
        parse = GPCRDBParsingPDB()
        ######################
        # if self.revise_xtal and self.force_main_temp:
        #     s = Structure.objects.get(pdb_code__index='6CM4')
        #     temp_st = OrderedDict([(s, 100)])
        #     for st, sim in similarity_table.items():
        #         temp_st[st] = sim
        #     similarity_table = temp_st
        ######################
        for struct in similarity_table:
            try:
                alt_temp = parse.fetch_residues_from_pdb(struct, [gn])
                if reference_dict[ref_seg][ref_res]==PDB.Polypeptide.three_to_one(
                                                        alt_temp[gn_][0].get_parent().get_resname()):
                    orig_res = main_pdb_array[ref_seg][str(ref_res).replace('x','.')]
                    alt_res = alt_temp[gn_]
                    if len(alt_res)!=atom_num_dict[reference_dict[ref_seg][ref_res]]:
                        continue
                    superpose = sp.RotamerSuperpose(orig_res, alt_res)
                    new_atoms = superpose.run()
                    if self.debug:
                        print(struct, gn_, superpose.backbone_rmsd, superpose.rmsd)
                    if superpose.backbone_rmsd>0.45:
                        continue
                    main_pdb_array[ref_seg][str(ref_res).replace('x','.')] = new_atoms
                    template_dict[temp_seg][temp_res] = reference_dict[ref_seg][ref_res]
                    non_cons_res_templates[gn] = struct
                    switched_count+=1
                    no_match = False
                    if 'x' not in ref_res:
                        num_in_loop = parse.gn_num_extract(ref_res,'|')[1]
                        seq_num = str(list(Residue.objects.filter(protein_conformation=self.prot_conf,
                                                                  protein_segment__slug=segment))[num_in_loop-1].sequence_number)
                        self.template_source = update_template_source(self.template_source,[seq_num],struct,segment,just_rot=True)
                    else:
                        self.template_source = update_template_source(self.template_source,[ref_res],struct,segment,just_rot=True)
                    break
            except:
                pass
        return main_pdb_array, template_dict, non_cons_res_templates, switched_count, no_match

    def ECL3_disulfide(self, reference_dict):
        c61, c62 = False, False
        try:
            if reference_dict['TM6']['6x61']=='C':
                c61 = True
        except:
            pass
        try:
            if reference_dict['TM6']['6x62']=='C':
                c62 = True
        except:
            pass
        ecl3_lab = [i for i in reference_dict if i.startswith('ECL3')][0]
        ecl3_c = []
        for gn, res in reference_dict[ecl3_lab].items():
            if res=='C':
                ecl3_c.append(gn)
        if c61==True and len(ecl3_c)>0:
            return ['6x61', ecl3_c[0]]
        elif c62==True and len(ecl3_c)>0:
            return ['6x62', ecl3_c[0]]
        elif len(ecl3_c)>=2:
            return [ecl3_c[0], ecl3_c[1]]
        else:
            return [0,0]

    def write_homology_model_pdb(self, filename, main_pdb_array, alignment, trimmed_residues=[], disulfide_pairs=[], complex=False):
        ''' Write PDB file from pdb array to file.

            @param filename: str, filename of output file \n
            @param main_pdb_array: OrderedDict(), of atoms of pdb, where keys are generic numbers/residue numbers and
            values are list of atoms. Output of GPCRDBParsingPDB.pdb_array_creator().
            @param alignment: AlignedReferenceTemplate class, alignment of reference and template.
            @trimmed_residues: list, list of generic numbers that are trimmed/to be modeled by MODELLER.
        '''
        key = ''
        res_num = 0
        counter_num = 0
        atom_num = 0
        trimmed_resi_nums = OrderedDict()
        helix_restraints = []
        prev_seg = '0'
        icl3_mid = None
        disulfide_nums = [[0,0],[0,0]]
        complex_start, beta_start, gamma_start = None, None, None
        with open(filename,'w+') as f:
            for seg_id, segment in main_pdb_array.items():
                if seg_id!='TM1' and prev_seg!='0' and seg_id.startswith('T') and prev_seg.startswith('T'):
                    atom_num+=1
                    # f.write("\nTER{}      {} {}{}".format(str(atom_num).rjust(8),atom.get_parent().get_resname(),str(self.main_template_preferred_chain)[0],str(res_num).rjust(4)))
                trimmed_segment = OrderedDict()
                for key in segment:
                    res_num+=1
                    counter_num+=1
                    for i, d_p in enumerate(disulfide_pairs):
                        for j, d in enumerate(d_p):
                            try:
                                if key==d.replace('x','.'):
                                    disulfide_nums[i][j] = res_num
                                    break
                            except:
                                pass
                    try:
                        if alignment.reference_dict[seg_id][key.replace('.','x')] in ['-','x']:
                            counter_num-=1
                            res_num-=1
                            continue
                    except:
                        try:
                            if alignment.reference_dict[seg_id][key] in ['-','x']:
                                counter_num-=1
                                res_num-=1
                                continue
                        except:
                            pass
                    if seg_id=='G.HN':
                        if complex_start==None:
                            complex_start = counter_num
                    if seg_id=='Beta':
                        if beta_start==None:
                            beta_start = counter_num
                    if seg_id=='Gamma':
                        if gamma_start==None:
                            gamma_start = counter_num
                    if segment[key]=='/':
                        atom_num+=1
                        icl3_mid = counter_num
                        res_num-=1
                        counter_num-=1
                        continue
                    if key in trimmed_residues:
                        trimmed_segment[key] = counter_num
                        if 'x' in segment[key]:
                            if '?' in key:
                                atom_num+=1
                                continue
                            else:
                                helix_restraints.append(counter_num)
                                continue
                    if 'x' in segment[key]:
                        atom_num+=1
                        continue
                    if '?' in key and '-' in segment[key]:
                        atom_num+=1
                        continue
                    if '-term' in seg_id and segment[key]=='-':
                        continue
                    for atom in main_pdb_array[seg_id][key]:
                        atom_num+=1
                        coord = list(atom.get_coord())
                        coord1 = "%8.3f"% (coord[0])
                        coord2 = "%8.3f"% (coord[1])
                        coord3 = "%8.3f"% (coord[2])
                        if str(atom.get_id())=='CA':
                            if len(key)==4:
                                bfact = "%6.2f"% (float(key))
                            elif '.' not in key:
                                bfact = "%6.2f"% (float(atom.get_bfactor()))
                            else:
                                key_split = key.split('.')
                                if '.' in key and len(key_split[1])==3:
                                    bfact = " -%4.2f"% (float(key))
                                elif len(key_split)==3:
                                    bfact = "%6.2f"% (float(atom.get_bfactor()))
                                else:
                                    bfact = " %5.2f"% (float(key))
                        else:
                            bfact = "%6.2f"% (float(atom.get_bfactor()))
                        occupancy = "%6.2f"% (atom.get_occupancy())
                        template="""
ATOM{atom_num}  {atom}{res} {chain}{res_num}{coord1}{coord2}{coord3}{occupancy}{bfactor}{atom_s}  """
                        context={"atom_num":str(atom_num).rjust(7), "atom":str(atom.get_id()).ljust(4),
                                 "res":atom.get_parent().get_resname(),
                                 "chain":str(self.main_template_preferred_chain)[0],
                                 "res_num":str(res_num).rjust(4), "coord1":coord1.rjust(12),
                                 "coord2":coord2.rjust(8), "coord3":coord3.rjust(8),
                                 "occupancy":str(occupancy).rjust(3),
                                 "bfactor":str(bfact).rjust(4), "atom_s":str(str(atom.get_id())[0]).rjust(12)}
                        f.write(template.format(**context))
                trimmed_resi_nums[seg_id] = trimmed_segment
                prev_seg = seg_id[:4]
            f.write("\nTER\n")
            print(self.main_structure)
            if self.reference_entry_name!=self.main_structure.protein_conformation.protein.parent.entry_name:
                atom_num+=1
                # f.write("\nTER{}      {} {}{}".format(str(atom_num).rjust(8),atom.get_parent().get_resname(),str(self.main_template_preferred_chain)[0],str(res_num).rjust(4)))

        return trimmed_resi_nums, helix_restraints, icl3_mid, disulfide_nums, complex_start, beta_start, gamma_start

    def create_PIR_file(self, reference_dict, template_dict, template_file, hetatm_count, water_count):
        ''' Create PIR file from reference and template alignment (AlignedReferenceAndTemplate).

            @param reference_dict: AlignedReferenceAndTemplate.reference_dict
            @param template_dict: AlignedReferenceAndTempalte.template_dict
            @template_file: str, name of template file with path
            @param hetatm_count: int, number of hetero atoms
            @param water_count: int, number of water atoms

        '''
        ref_sequence, temp_sequence = '',''
        res_num = 1
        with open(template_file,'r') as f:
            lines = f.readlines()
            for line in lines:
                try:
                    pdb_re = re.search('(ATOM[A-Z\s\d]{13}\S{3}\s\S\s+)(\d+)([A-Z\s\d.-]{49,53})',line)
                    start_num = pdb_re.group(2)
                    break
                except:
                    try:
                        pdb_re = re.search('(ATOM[A-Z\s\d]{13}\S{3}\s+)(\d+)([A-Z\s\d.-]{49,53})',line)
                        start_num = pdb_re.group(2)
                        break
                    except:
                        pass
        for ref_seg, temp_seg in zip(reference_dict, template_dict):
            if ref_seg in ['G.HN','Beta','Gamma']:
                ref_sequence+='/'
                temp_sequence+='/'
            for ref_res, temp_res in zip(reference_dict[ref_seg], template_dict[temp_seg]):
                if reference_dict[ref_seg][ref_res] in ['-','x']:
                    continue
                else:
                    ref_sequence+=reference_dict[ref_seg][ref_res]
                if template_dict[temp_seg][temp_res] in ['-','x']:
                    temp_sequence+='-'
                else:
                    temp_sequence+=template_dict[temp_seg][temp_res]
                res_num+=1
        for i in range(hetatm_count):
            ref_sequence+='.'
            temp_sequence+='.'
        for i in range(water_count):
            ref_sequence+='w'
            temp_sequence+='w'
        if len(self.peptide_ligands)>0:
            for l in self.peptide_ligands:
                ref_sequence+='/'+l[-1]
                temp_sequence+='/'+l[-1]

        self.model_sequence = temp_sequence

        if self.complex:
            pir_file = "./structure/PIR/{}_{}.pir".format(self.uniprot_id, self.target_signprot.entry_name)
        else:
            pir_file = "./structure/PIR/"+self.uniprot_id+"_"+self.state+".pir"
        with open(pir_file, 'w+') as output_file:
            template="""
>P1;{temp_file}
structure:{temp_file}:{start}:{chain}:{res_num}:{chain}::::
{temp_sequence}*

>P1;{uniprot}
sequence:{uniprot}::::::::
{ref_sequence}*
            """
            context={"temp_file":template_file,
                     "start":start_num,
                     "chain":self.main_template_preferred_chain,
                     "res_num":res_num,
                     "temp_sequence":temp_sequence,
                     "uniprot":self.uniprot_id,
                     "ref_sequence":ref_sequence}
            output_file.write(template.format(**context))

    def run_MODELLER(self, pir_file, template, reference, number_of_models, output_file_name, atom_dict=None,
                     helix_restraints=[], icl3_mid=None, disulfide_nums=[], complex_start=None, beta_start=None, gamma_start=None, peptide_ligands=[]):
        ''' Build homology model with MODELLER.

            @param pir_file: str, file name of PIR file with path \n
            @param template: str, file name of template with path \n
            @param reference: str, Uniprot code of reference sequence \n
            @param number_of_models: int, number of models to be built \n
            @param output_file_name: str, name of output file
            @param atom_dict: nested OrderedDict(), atoms to model with MODELLER organized by segments and generic
            numbers, default=None
            @param helix_restraints: list, list of generic numbers that should be modelled as helical regions by
            MODELLER, default=[]
            @param icl3_mid: int, position of the break in the middle of ICL3, default=None
        '''
        log.none()
        env = environ(rand_seed=1028) #!!random number generator
        if self.revise_xtal==True:
            ref_prot = self.reference_protein.parent
        else:
            ref_prot = self.reference_protein
        if ref_prot==self.main_structure.protein_conformation.protein.parent or self.keep_hetatoms:
            env.io.hetatm = True
            env.io.water = True
        if atom_dict==None:
            a = automodel(env, alnfile = pir_file, knowns = template, sequence = reference,
                          assess_methods=(assess.DOPE))
        else:
            a = HomologyMODELLER(env, alnfile = pir_file, knowns = template, sequence = reference,
                                 assess_methods=(assess.DOPE), atom_selection=atom_dict,
                                 helix_restraints=helix_restraints, icl3_mid=icl3_mid, disulfide_nums=disulfide_nums, complex_start=complex_start, beta_start=beta_start, gamma_start=gamma_start, peptide_ligands=peptide_ligands)

        a.starting_model = 1
        a.ending_model = number_of_models
        if self.fast_refinement:
            a.md_level = refine.very_fast
        else:
            a.md_level = refine.slow
        path = "./structure/homology_models/"
        if not os.path.exists(path):
            os.mkdir(path)
        a.make()

        # Get a list of all successfully built models from a.outputs
        ok_models = [x for x in a.outputs if x['failure'] is None]
        if len(ok_models)==0:
            os.rename("./"+template, output_file_name)
            return 0

        # Rank the models by DOPE score
        key = 'DOPE score'
        if sys.version_info[:2] == (2,3):
            # Python 2.3's sort doesn't have a 'key' argument
            ok_models.sort(lambda a,b: cmp(a[key], b[key]))
        else:
            ok_models.sort(key=lambda a: a[key])

        # Get top model
        m = ok_models[0]
#        print("Top model: %s (DOPE score %.3f)" % (m['name'], m[key]))

        for file in os.listdir("./"):
            if file==m['name']:
                os.rename("./"+file, output_file_name)
            elif file.startswith(self.uniprot_id):
                os.remove("./"+file)


class SilentModeller(object):
    ''' No text to console.
    '''
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, *args):
        sys.stdout.close()
        sys.stdout = self._stdout


class HomologyMODELLER(automodel):
    def __init__(self, env, alnfile, knowns, sequence, assess_methods, atom_selection, helix_restraints=[], icl3_mid=None, disulfide_nums=[], complex_start=None, beta_start=None, gamma_start=None, peptide_ligands=[]):
        super(HomologyMODELLER, self).__init__(env, alnfile=alnfile, knowns=knowns, sequence=sequence,
                                               assess_methods=assess_methods)
        self.alnfile = alnfile
        self.atom_dict = atom_selection
        self.helix_restraints = helix_restraints
        self.icl3_mid = icl3_mid
        print('ICL3', self.icl3_mid)
        self.disulfide_nums = disulfide_nums
        print('DISULFIDE')
        print(self.disulfide_nums)
        self.complex_start = complex_start
        self.beta_start = beta_start
        self.gamma_start = gamma_start
        print(self.complex_start, self.beta_start, self.gamma_start)
        self.peptide_ligands = peptide_ligands

    def identify_chain(self, seq_num):
        if len(self.peptide_ligands)>0:
            i = len(self.peptide_ligands)
        else:
            i = 0
        if self.complex_start!=None:
            if len(self.chains)==2+i:
                if seq_num<self.complex_start:
                    return 'A'
                elif seq_num>=self.complex_start:
                    return 'B'
                else:
                    return ''
            elif len(self.chains)==3+i:
                ValueError('DEBUG me: 3 chains during modeller run. Complex model with long ICL3')
            elif len(self.chains)==4+i:
                # Normal complex model with one receptor chain and 3 g protein subunits
                if seq_num<self.complex_start:
                    return 'A'
                elif self.complex_start<=seq_num<self.beta_start:
                    return 'B'
                elif self.beta_start<=seq_num<self.gamma_start:
                    return 'C'
                elif seq_num>=self.gamma_start:
                    return 'D'
            elif len(self.chains)==5+i:
                # Long ICL3 complex model with two receptor chains and 3 g protein subunits
                if seq_num<self.icl3_mid:
                    return 'A'
                elif self.icl3_mid<=seq_num<self.complex_start:
                    return 'B'
                elif self.complex_start<=seq_num<self.beta_start:
                    return 'C'
                elif self.beta_start<=seq_num<self.gamma_start:
                    return 'D'
                elif seq_num>=self.gamma_start:
                    return 'E'
            else:
                return ''

        else:
            if len(self.chains)==2+i:
                if self.icl3_mid==None:
                    return 'A'
                elif seq_num<self.icl3_mid:
                    return 'A'
                elif seq_num>=self.icl3_mid:
                    return 'B'
                else:
                    return ''
            elif len(self.chains)==1+i and i!=0:
                return 'A'
            else:
                return ''

    def find_helix_restraints(self):
        print(self.helix_restraints)
        start = 0
        prev = 0
        out = []
        for i in self.helix_restraints:
            if prev==0:
                start = i
            if i==prev+1:
                pass
            elif prev!=0:
                out.append([start,prev])
                start = i
            prev=i
        if len(self.helix_restraints)>0:
            out.append([start,prev])
        print('helix restraints', out)
        return out

    def select_atoms(self):
        for r in self.residues:
            print(r)
        selection_out = []
        for seg_id, segment in self.atom_dict.items():
            for gn, atom in segment.items():
                chain = self.identify_chain(atom)
                print(seg_id, gn, atom, chain)
                selection_out.append(self.residues[str(atom)+':{}'.format(chain)])
        return selection(selection_out)

    def special_restraints(self, aln):
        rsr = self.restraints
        for i in self.find_helix_restraints():
            chain = self.identify_chain(i[0])
            for j, k in self.atom_dict.items():
                segment = sorted(list(k.items()),key=lambda e: e[0])
                if segment==[]:
                    continue
                # No helix extension for Gprotein helices
                if j in gprotein_segment_slugs:
                    continue
                if j=='H8' and i[0]==segment[0][1] and i[1]==segment[-1][1]:
                    rsr.add(secondary_structure.alpha(self.residue_range('{}:{}'.format(i[0]-4,chain),'{}:{}'.format(i[1],chain))))
                    break
                elif i[0]==segment[0][1]:
                    rsr.add(secondary_structure.alpha(self.residue_range('{}:{}'.format(i[0],chain),'{}:{}'.format(i[1]+4,chain))))
                    break
                elif i[1]==segment[-1][1]:
                    rsr.add(secondary_structure.alpha(self.residue_range('{}:{}'.format(i[0]-4,chain),'{}:{}'.format(i[1],chain))))
                    break
        # Custom special distance based restraints
        with open(self.alnfile, 'r') as alnfile:
            template_line = alnfile.readlines()[1]
        # To fix ECL1 loop knot with N-term in 5NX2
        if '5nx2_Intermediate_post.pdb' in template_line:
            at = self.atoms
            rsr.add(forms.gaussian(group=physical.xy_distance,
                                   feature=features.distance(at['CA:63'],
                                                             at['CA:186']),
                                   mean=15.0, stdev=0.1))
            rsr.add(forms.gaussian(group=physical.xy_distance,
                                   feature=features.distance(at['CA:63'],
                                                             at['CA:187']),
                                   mean=15.0, stdev=0.1))
            rsr.add(forms.gaussian(group=physical.xy_distance,
                                   feature=features.distance(at['CA:63'],
                                                             at['CA:188']),
                                   mean=15.0, stdev=0.1))
            rsr.add(forms.gaussian(group=physical.xy_distance,
                                   feature=features.distance(at['CA:63'],
                                                             at['CA:189']),
                                   mean=15.0, stdev=0.1))
            rsr.add(forms.gaussian(group=physical.xy_distance,
                                   feature=features.distance(at['CA:63'],
                                                             at['CA:190']),
                                   mean=15.0, stdev=0.1))
            rsr.add(forms.gaussian(group=physical.xy_distance,
                                   feature=features.distance(at['CA:63'],
                                                             at['CA:191']),
                                   mean=15.0, stdev=0.1))

    def special_patches(self, aln):
        for d in self.disulfide_nums:
            if d[0]==0:
                continue
            try:
                if self.identify_chain(d[0]) in ['A','B','C','D']:
                    formatted_res_id1 = '{}:{}'.format(d[0], self.identify_chain(d[0]))
                else:
                    formatted_res_id1 = str(d[0])
                if self.identify_chain(d[1]) in ['A','B','C','D']:
                    formatted_res_id2 = '{}:{}'.format(d[1], self.identify_chain(d[1]))
                else:
                    formatted_res_id2 = str(d[1])
                self.patch(residue_type='DISU', residues=(self.residues[formatted_res_id1],
                                                          self.residues[formatted_res_id2]))
            except Exception as msg:
                print("DISULFIDE ISSUE")
                print(msg)
                logger.error('Disulfide issue: {} {}'.format(d[0], d[1]))
                pass

    # def make(self):
    #     with SilentModeller():
    #         super(HomologyMODELLER, self).make()


class SegmentEnds(object):
    def __init__(self):
        self.start = None
        self.end = None
        self.protein_segment = None

    def __repr__(self):
        return "<{},{},{}>".format(self.start,self.end,self.protein_segment)


class HelixEndsModeling(HomologyModeling):
    ''' Class for modeling the helix ends of GPCRs.
    '''
    def __init__(self, similarity_table, template_source, main_structure, revise_xtal=False, debug=False):
        self.helix_ends = OrderedDict()
        self.helix_end_mods = OrderedDict()
        self.main_pdb_array = OrderedDict()
        self.alignment = OrderedDict()
        self.similarity_table = similarity_table
        self.template_source = template_source
        self.main_structure = main_structure
        self.revise_xtal =  revise_xtal
        self.templates_to_skip = OrderedDict([('TM1',[['6BQG','6BQH'],[]]),('TM2',[[],[]]),('TM3',[[],[]]),('TM4',[[],[]]),
                                              ('TM5',[[],[]]),('TM6',[[],[]]),('TM7',[[],['5UNF','5UNG','5UNH']]),('H8',[[],[]])])
        ## NMUR2
        # self.templates_to_skip = OrderedDict([('TM1',[['6BQG','6BQH'],[]]),('TM2',[[],[]]),('TM3',[[],[]]),('TM4',[[],[]]),
        #                                       ('TM5',[[],[]]),('TM6',[[],[]]),('TM7',[[],['5UNF','5UNG','5UNH','7F9Z', '6KO5', '7F9Y', '4BUO', '4XES', '4GRV', '4XEE', '5T04', '7L0P', '7L0R', '7L0Q', '7L0S', '6OSA','6PWC','6OS9','6UP7','6B73','6VI4','4DJH','4N6H','6PT2','5XSZ','4RWD','6PT3','4EJ4','5X93','6TO7','5WS3','7DDZ','6ME2','5DHG','7L1V','7DB6','5WQC','6IGK','6TOD','6TOS','6TOT','6TQ4','6TP6','5GLI','4S0V','6TQ6','6TPN','6TQ7','6TQ9','6IGL','6K1Q','6TPG','6TPJ','4ZJ8','5GLH','4ZJC','6ME3','6LRY','5DHH','4EA3','6TP4','6TP3','7L1U','6ME4','6ME5','6PS8','6V9S','5XPR','5C1M','7F8Y','6BQH','7CMV','7JVR','4ZUD','4DKL','6CM4','3PBL','7EZM','6JOD','7MBX','7MBY','7BVQ','7BU7','7BU6','6OS2','6OS1','5UNF','5UNG','7F8U','6OS0','6DO1','4YAY','5UNH','6BQG','7CMU','7DFP','6LUQ','7EZK','7BTS','5XJM']]),('H8',[[],[]])])
        ## Kappa opioid
        # self.templates_to_skip = OrderedDict([('TM1',[['6BQG','6BQH'],['6B73', '6VI4', '4DJH', '5C1M', '4DKL', '6DDF']]),('TM2',[['6B73', '6VI4', '4DJH', '5C1M', '4DKL', '6DDF'],[]]),('TM3',[[],['6B73', '6VI4', '4DJH', '5C1M', '4DKL', '6DDF']]),('TM4',[['6B73', '6VI4', '4DJH', '5C1M', '4DKL', '6DDF'],[]]),
        #                                       ('TM5',[[],['6B73', '6VI4', '4DJH', '5C1M', '4DKL', '6DDF']]),('TM6',[['6B73', '6VI4', '4DJH', '5C1M', '4DKL', '6DDF'],[]]),('TM7',[[],['5UNF','5UNG','5UNH','6B73', '6VI4', '4DJH', '5C1M', '4DKL', '6DDF']]),('H8',[[],[]])])
        ## NPY1R - OX2R based
        # self.templates_to_skip = OrderedDict([('TM1',[['6BQG','6BQH'],[]]),('TM2',[[],[]]),('TM3',[[],[]]),('TM4',[[],[]]),
        #                                       ('TM5',[[],['5ZBH','5ZBQ','7DDZ','6TO7','5WS3','6ME2','7L1V','7F8W','7DB6','7L1U','6HLP','7F8Y','7EZM','7EZH','6B73','6VI4','5C1M','4DKL','6DDE','4N6H','6PT2','7E32','5DHG','7EXD','1GZM']]),('TM6',[['5ZBH','5ZBQ','7DDZ','6TO7','5WS3','6ME2','7L1V','7F8W','7DB6','7L1U','6HLP','7F8Y','7EZM'],[]]),('TM7',[[],['5UNF','5UNG','5UNH','5ZBH','5ZBQ','7DDZ','6TO7','5WS3','6ME2','7L1V','7F8W','7DB6','7L1U','6HLP','7F8Y','7EZM']]),('H8',[[],[]])])
        # ## NPY1R - CCKAR based
        # self.templates_to_skip = OrderedDict([('TM1',[['6BQG','6BQH'],[]]),('TM2',[[],[]]),('TM3',[[],[]]),('TM4',[[],[]]),
        #                                       ('TM5',[[],[]]),('TM6',[[],[]]),('TM7',[[],['5UNF','5UNG','5UNH']]),('H8',[[],[]])])
        self.debug = debug

    def find_ends(self, structure, protein_conformation):
        raw_res = Residue.objects.filter(protein_conformation=protein_conformation).exclude(
                                    protein_segment=None).order_by('protein_segment_id').distinct('protein_segment_id')
        raw_segs = [i.protein_segment for i in raw_res]
        ends = []
        for i in raw_segs:
            if i.slug[0] not in ['T','H']:
                continue
            end = SegmentEnds()
            end.start = list(Residue.objects.filter(protein_conformation=protein_conformation,
                                                    protein_segment__slug=i))[0].sequence_number
            end.end = list(Residue.objects.filter(protein_conformation=protein_conformation,
                                                  protein_segment__slug=i))[-1].sequence_number
            end.protein_segment = i
            ends.append([end, structure])
        return ends

    def fetch_struct_helix_ends_from_db(self, structure, H8_alt=None):
        ''' Returns structure's helix end generic numbers after updating them with annotated data.
        '''
        raw = self.find_ends(structure, structure.protein_conformation)
        anno_conf = ProteinConformation.objects.get(protein=structure.protein_conformation.protein.parent)
        annotated = self.find_ends(structure, anno_conf)

        if H8_alt!=None and H8_alt!=structure:
            H8_raw_conf = ProteinConformation.objects.get(protein=H8_alt.protein_conformation.protein.parent)
            if raw[-1][0].protein_segment.slug=='H8':
                try:
                    raw[-1] = [i for i in self.find_ends(H8_alt,H8_raw_conf) if i[0].protein_segment.slug=='H8'][0]
                except:
                    pass
            else:
                try:
                    raw.append([i for i in self.find_ends(H8_alt,H8_raw_conf) if i[0].protein_segment.slug=='H8'][0])
                except:
                    pass
            if annotated!=[]:
                if annotated[-1][0].protein_segment.slug=='H8':
                    try:
                        annotated[-1] = [i for i in self.find_ends(
                                            H8_alt,H8_alt.protein_conformation) if i[0].protein_segment.slug=='H8'][0]
                    except:
                        pass
                else:
                    try:
                        annotated.append([i for i in self.find_ends(
                                            H8_alt,H8_alt.protein_conformation) if i[0].protein_segment.slug=='H8'][0])
                    except:
                        pass
        ends = OrderedDict()
        for i in raw:
            if i[0].protein_segment.slug[0]=='T' or i[0].protein_segment.slug=='H8':
                if len(list(Residue.objects.filter(protein_conformation=i[1].protein_conformation,
                                                   protein_segment=i[0].protein_segment)))==0:
                    continue
                start_found = False
                break_point = 0
                while start_found==False:
                    if break_point==20:
                        i[0].start = None
                        break
                    try:
                        if Residue.objects.get(protein_conformation=i[1].protein_conformation,
                                               sequence_number=i[0].start).generic_number==None:
                            i[0].start+=1
                        else:
                            start_found = True
                    except:
                        break_point+=1
                        i[0].start+=1
                s = Residue.objects.get(protein_conformation=i[1].protein_conformation,sequence_number=i[0].start)
                end_found = False
                break_point = 0
                while end_found==False:
                    if break_point==20:
                        i[0].end = None
                        break
                    try:
                        if Residue.objects.get(protein_conformation=i[1].protein_conformation,
                                               sequence_number=i[0].end).generic_number==None:
                            i[0].end-=1
                        else:
                            end_found = True
                    except:
                        break_point+=1
                        i[0].end-=1
                e = Residue.objects.get(protein_conformation=i[1].protein_conformation,sequence_number=i[0].end)
                ends[s.protein_segment.slug] = [ggn(s.display_generic_number.label),ggn(e.display_generic_number.label)]
        for j in annotated:
            if j[0].protein_segment.slug[0]=='T' or j[0].protein_segment.slug=='H8':
                if j[0].start!=0:
                    found_start = False
                    break_point = 0
                    while found_start==False:
                        if break_point==20:
                            j[0].start = None
                            break
                        try:
                            if Residue.objects.get(protein_conformation=j[1].protein_conformation,
                                                   sequence_number=j[0].start).generic_number!=None:
                                found_start = True
                            else:
                                raise Exception()
                        except:
                            break_point+=1
                            j[0].start+=1
                    if j[0].start!=None:
                        sa = Residue.objects.get(protein_conformation=j[1].protein_conformation,sequence_number=j[0].start)
                        ends[j[0].protein_segment.slug][0] = ggn(sa.display_generic_number.label)
                if j[0].end!=0:
                    found_end = False
                    break_point = 0
                    while found_end==False:
                        if break_point==20:
                            j[0].end = None
                            break
                        try:
                            if Residue.objects.get(protein_conformation=j[1].protein_conformation,
                                                   sequence_number=j[0].end).generic_number!=None:
                                found_end = True
                            else:
                                raise Exception()
                        except:
                            break_point+=1
                            j[0].end-=1
                    if j[0].end!=None:
                        ea = Residue.objects.get(protein_conformation=j[1].protein_conformation,sequence_number=j[0].end)
                    try:
                        ends[j[0].protein_segment.slug][1] = ggn(ea.display_generic_number.label)
                    except:
                        pass
        return ends

    def fetch_struct_helix_ends_from_array(self, array):
        ''' Returns helix ends from structure array (GPCRDBParsingPDB.pdb_array_creator()).
        '''
        ends = OrderedDict()
        for seg_lab, seg in array.items():
            if seg_lab[0]=='T' or seg_lab=='H8':
                try:
                    ends[seg_lab] = [[i for i in seg.keys() if seg[i]!='x'][0].replace('.','x'),[i for i in seg.keys() if seg[i]!='x'][-1].replace('.','x')]
                except:
                    pass
        return ends

    def correct_helix_ends(self, main_structure, main_pdb_array, a, template_source, separate_H8=None):
        ''' Updates main template structure with annotated helix ends. If helix is too long, it removes residues, if it
            is too short, it superpositions residues from next closest template. Updates alignment with changes.

            @param separate_H8: if True, skip H8; if False, only H8; if None (default), both
        '''
        modifications = {'added':{'TM1':[[],[]],'TM2':[[],[]],'TM3':[[],[]],'TM4':[[],[]],'TM5':[[],[]],'TM6':[[],[]],
                                  'TM7':[[],[]], 'H8':[[],[]]},
                         'removed':{'TM1':[[],[]],'TM2':[[],[]],'TM3':[[],[]],'TM4':[[],[]],'TM5':[[],[]],'TM6':[[],[]],
                                    'TM7':[[],[]], 'H8':[[],[]]}}
        H8_alt = None
        if 'H8' in template_source:
            for i, j in template_source['H8'].items():
                if j[0]:
                    H8_alt = j[0]
                    break
        if separate_H8==True:
            H8_alt = None
        if self.debug:
            if not separate_H8 and not H8_alt:
                print('Error: Bad H8 alternative template')
        if separate_H8==False and H8_alt.pdb_code.index in structures_with_num_issues:
            ssno = StructureSeqNumOverwrite(H8_alt)
            ssno.seq_num_overwrite('pdb')
        raw_helix_ends = self.fetch_struct_helix_ends_from_array(main_pdb_array)
        anno_helix_ends = self.fetch_struct_helix_ends_from_db(main_structure, H8_alt)
        if self.debug:
            print('raw: ',raw_helix_ends)
            print('anno: ',anno_helix_ends)

        #### Exception for TM4 start of 5ZKP
        if main_structure.pdb_code.index in ['5ZKP']:
            for i in ['4x38','4x39','4x40','4x41','4x42','4x43','4x44','4x45','4x46','4x47','4x48']:
                a.template_dict['TM4'][i] = 'x'
        elif main_structure.pdb_code.index in ['5UZ7']:
            for i in list(range(35,50)):
                a.template_dict['TM6']['6x'+str(i)] = 'x'
        elif main_structure.pdb_code.index in ['6BQG', '6BQH'] and not self.revise_xtal:
            for i in list(range(32,40)):
                a.template_dict['TM1']['1x'+str(i)] = 'x'

        # Force active state with main template 5UNF, 5UNG or 5UNH to get new TM7 end
        skip_template = False
        if separate_H8 and main_structure.pdb_code.index in ['5UNF','5UNG','5UNH']:
            anno_helix_ends['TM7'][1] = '7x48'
            skip_template = True
        for lab,seg in a.template_dict.items():
            if separate_H8==True:
                if lab=='H8':
                    continue
            elif separate_H8==False:
                if lab!='H8':
                    continue
            for gn,res in seg.items():
                try:
                    if lab[0] in ['H']:
                        if res!='-':
                            r = Residue.objects.get(protein_conformation=H8_alt.protein_conformation,
                                                    display_generic_number__label=dgn(
                                                                            gn,H8_alt.protein_conformation))
                            if len(Rotamer.objects.filter(structure=H8_alt,residue=r))<1:
                                raise Exception()
                except:
                    a.template_dict[lab][gn] = 'x'
                    a.alignment_dict[lab][gn] = 'x'
        parser = GPCRDBParsingPDB()

        for raw_seg, anno_seg in zip(raw_helix_ends, anno_helix_ends):
            if separate_H8==True:
                if raw_seg=='H8':
                    continue
            elif separate_H8==False:
                if raw_seg!='H8':
                    continue
            if H8_alt!=None and H8_alt!=main_structure and raw_seg=='H8':
                template = H8_alt
            else:
                template = main_structure
            protein_conf = ProteinConformation.objects.get(protein=template.protein_conformation.protein.parent)
            try:
                s_dif = parser.gn_comparer(raw_helix_ends[raw_seg][0],anno_helix_ends[anno_seg][0],
                                           protein_conf)
            except:
                try:
                    s_dif = parser.gn_comparer(raw_helix_ends[raw_seg][0],anno_helix_ends[anno_seg][0],
                                               template.protein_conformation)
                    protein_conf = template.protein_conformation
                except:
                    for i in range(int(raw_helix_ends[raw_seg][0].split('x')[1]),
                                   int(anno_helix_ends[anno_seg][0].split('x')[1])):
                        a.template_dict[raw_seg]['8x{}'.format(str(i))]='x'
                        a.alignment_dict[raw_seg]['8x{}'.format(str(i))]='x'
                    s_dif=0
            if s_dif<0:
                s_gn = Residue.objects.get(protein_conformation=protein_conf,
                                           display_generic_number__label=dgn(raw_helix_ends[raw_seg][0],
                                                                             protein_conf))
                seq_nums = [i for i in range(s_gn.sequence_number,s_gn.sequence_number-s_dif)]
                gns = [ggn(j.display_generic_number.label) for j in list(Residue.objects.filter(
                            protein_conformation=protein_conf, sequence_number__in=seq_nums))]
                for gn in gns:
                    if gn in a.template_dict[raw_seg]:
                        a.template_dict[raw_seg][gn]='x'
                        a.alignment_dict[raw_seg][gn]='x'
                    else:
                        del main_pdb_array[raw_seg][gn.replace('x','.')]
                        modifications['removed'][raw_seg][0].append(gn)
            protein_conf = ProteinConformation.objects.get(protein=template.protein_conformation.protein.parent)
            try:
                e_dif = parser.gn_comparer(raw_helix_ends[raw_seg][1],anno_helix_ends[anno_seg][1],
                                           protein_conf)
            except:
                try:
                    e_dif = parser.gn_comparer(raw_helix_ends[raw_seg][1],anno_helix_ends[anno_seg][1],
                                               template.protein_conformation)
                    protein_conf = template.protein_conformation
                except:
                    for i in range(int(anno_helix_ends[anno_seg][1].split('x')[1])+1,
                                   int(raw_helix_ends[raw_seg][1].split('x')[1])+1):
                        a.template_dict[raw_seg]['8x{}'.format(str(i))]='x'
                        a.alignment_dict[raw_seg]['8x{}'.format(str(i))]='x'
                    e_dif = 0
            if e_dif>0:
                e_gn = Residue.objects.get(protein_conformation=protein_conf,
                                           display_generic_number__label=dgn(raw_helix_ends[raw_seg][1],
                                                                             protein_conf))
                seq_nums = [i for i in range(e_gn.sequence_number-e_dif+1,e_gn.sequence_number+1)]
                gns = [ggn(j.display_generic_number.label) for j in list(Residue.objects.filter(
                            protein_conformation=protein_conf, sequence_number__in=seq_nums))]
                for gn in gns:
                    a.template_dict[raw_seg][gn]='x'
                    a.alignment_dict[raw_seg][gn]='x'
                    try:
                        a.reference_dict[raw_seg][gn]
                    except:
                        a.reference_dict[raw_seg][gn]='x'
        self.helix_ends = raw_helix_ends

        protein_anomalies_gns = ProteinAnomaly.objects.all().values_list('generic_number__label', flat=True)

        for ref_seg, temp_seg, align_seg in zip(a.reference_dict, a.template_dict, a.alignment_dict):
            if separate_H8==True:
                if ref_seg=='H8':
                    continue
            elif separate_H8==False:
                if ref_seg!='H8':
                    continue
            if ref_seg=='H8' and H8_alt!=None:
                first_res = Residue.objects.get(protein_conformation=H8_alt.protein_conformation,
                                                display_generic_number__label=dgn(raw_helix_ends[ref_seg][0],
                                                                                  H8_alt.protein_conformation)).sequence_number
                for h in list(a.template_dict[temp_seg].keys())[::-1]:
                    if a.template_dict[temp_seg][h]!='x':
                        raw_helix_ends[ref_seg][1]=h
                        break
                last_res = Residue.objects.get(protein_conformation=H8_alt.protein_conformation,
                                               display_generic_number__label=dgn(raw_helix_ends[ref_seg][1],
                                                                                 H8_alt.protein_conformation)).sequence_number
                temp_seg_seq_len = len(list(Residue.objects.filter(protein_conformation=H8_alt.protein_conformation,
                                                                   sequence_number__in=range(first_res,last_res+1))))
                mid = temp_seg_seq_len/2

            elif ref_seg[0]=='T':
                first_res = Residue.objects.get(protein_conformation=main_structure.protein_conformation,
                                                display_generic_number__label=dgn(raw_helix_ends[ref_seg][0],
                                                                                  main_structure.protein_conformation)).sequence_number
                last_res = Residue.objects.get(protein_conformation=main_structure.protein_conformation,
                                               display_generic_number__label=dgn(raw_helix_ends[ref_seg][1],
                                                                                 main_structure.protein_conformation)).sequence_number
                temp_seg_seq_len = len(list(Residue.objects.filter(protein_conformation=main_structure.protein_conformation,
                                                                   sequence_number__in=range(first_res,last_res+1))))
                mid = temp_seg_seq_len/2
            if ref_seg[0] not in ['T','H']:
                continue
            if separate_H8==True:
                if ref_seg=='H8':
                    continue
            elif separate_H8==False:
                if ref_seg!='H8':
                    continue
            offset = 0
            increase_offset = True
            full_template_dict_seg = deepcopy(a.template_dict[temp_seg])

            delete_r = set()
            delete_t = set()
            delete_a = set()
            delete_ar = set()

            # Correct offset
            for ref_res in a.reference_dict[ref_seg]:
                if a.template_dict[temp_seg][ref_res]!='x' and ref_res.replace('x','.') not in main_pdb_array[temp_seg] and ref_res not in protein_anomalies_gns:
                    a.template_dict[temp_seg][ref_res] = 'x'

            for ref_res, temp_res, align_res in zip(a.reference_dict[ref_seg],a.template_dict[temp_seg],
                                                    a.alignment_dict[align_seg]):
                if a.template_dict[temp_seg][temp_res]=='x' and increase_offset==True:
                    offset+=1
                elif a.template_dict[temp_seg][temp_res]!='x':
                    increase_offset = False
                # Ignore in helix gaps due to anomaly differences
                if a.template_dict[temp_seg][temp_res]=='-':
                    continue
                if a.reference_dict[ref_seg][ref_res]=='x':
                    if list(full_template_dict_seg.keys()).index(ref_res)<mid+offset:
                        modifications['removed'][ref_seg][0].append(ref_res)
                    else:
                        modifications['removed'][ref_seg][1].append(ref_res)
                    delete_r.add((ref_seg,ref_res))
                    delete_t.add((temp_seg,temp_res))
                    delete_a.add((align_seg,align_res))
                    delete_ar.add((ref_seg, ref_res.replace('x','.')))
                elif a.template_dict[temp_seg][temp_res]=='x' or (temp_seg[0]=='T' and temp_res.replace('x','.') not in
                                                                                        list(main_pdb_array[temp_seg])):
                    if list(full_template_dict_seg.keys()).index(temp_res)<mid+offset:
                        modifications['added'][temp_seg][0].append(temp_res)
                    else:
                        modifications['added'][temp_seg][1].append(temp_res)

            for i,ii in delete_r:
                del a.reference_dict[i][ii]
            for i,ii in delete_t:
                del a.template_dict[i][ii]
            for i,ii in delete_a:
                del a.alignment_dict[i][ii]
            for i,ii in delete_ar:
                try:
                    del main_pdb_array[i][ii]
                except:
                    pass

            if ref_seg[0]=='T' or ref_seg=='H8':
                if len(modifications['added'][ref_seg][0])>0:
                    self.helix_ends[ref_seg][0] = modifications['added'][ref_seg][0][0]
                if len(modifications['added'][ref_seg][1])>0:
                    self.helix_ends[ref_seg][1] = modifications['added'][ref_seg][1][-1]
                if len(modifications['removed'][ref_seg][0])>0:
                    self.helix_ends[ref_seg][0] = parser.gn_indecer(modifications['removed'][ref_seg][0][-1], 'x', 1)
                if len(modifications['removed'][ref_seg][1])>0:
                    self.helix_ends[ref_seg][1] = parser.gn_indecer(modifications['removed'][ref_seg][1][0], 'x', -1)
                if len(modifications['added'][ref_seg][0])>0:
                    found_alt_start = False
                    for struct in list(self.similarity_table):
                        if struct!=main_structure:
                            try:
                                if skip_template and struct.pdb_code.index in self.templates_to_skip[ref_seg][0]:
                                    continue
                                alt_helix_ends = self.fetch_struct_helix_ends_from_db(struct)
                                protein_conf = ProteinConformation.objects.get(protein=struct.protein_conformation.protein.parent)
                                if parser.gn_comparer(alt_helix_ends[ref_seg][0],self.helix_ends[ref_seg][0],
                                                      protein_conf)<=0:
                                    all_keys = list(a.reference_dict[ref_seg].keys())[:len(modifications['added'][ref_seg][0])+4]
                                    ref_keys = [i for i in all_keys if i not in modifications['added'][ref_seg][0]]
                                    reference = parser.fetch_residues_from_array(main_pdb_array[ref_seg],ref_keys)
                                    template = parser.fetch_residues_from_pdb(struct,all_keys)
                                    superpose = sp.OneSidedSuperpose(reference,template,4,0)
                                    sup_residues = superpose.run()
                                    new_residues = OrderedDict()
                                    for gn, atoms in sup_residues.items():
                                        gn_ = gn.replace('.','x')
                                        if gn_ not in ref_keys:
                                            new_residues[gn] = atoms
                                            a.template_dict[temp_seg][gn_] = PDB.Polypeptide.three_to_one(
                                                                             atoms[0].get_parent().get_resname())
                                            if a.template_dict[temp_seg][gn_]==a.reference_dict[ref_seg][gn_]:
                                                a.alignment_dict[ref_seg][gn_] = a.reference_dict[ref_seg][gn_]
                                            else:
                                                a.alignment_dict[ref_seg][gn_] = '.'
                                    for gn, atoms in main_pdb_array[ref_seg].items():
                                        if gn not in new_residues:
                                            new_residues[gn] = atoms
                                    main_pdb_array[ref_seg] = new_residues
                                    self.template_source = update_template_source(self.template_source,modifications['added'][ref_seg][0],struct,ref_seg)
                                    found_alt_start = True
                                    break
                            except:
                                pass
                    if found_alt_start==False:
                        new_residues = OrderedDict()
                        for i in modifications['added'][ref_seg][0]:
                            new_residues[i.replace('x','.')] = 'x'
                            a.template_dict[ref_seg][i] = 'x'
                            a.alignment_dict[ref_seg][i] = 'x'
                        for i,j in main_pdb_array[ref_seg].items():
                            new_residues[i] = j
                        main_pdb_array[ref_seg] = new_residues
                if len(modifications['added'][ref_seg][1])>0:
                    found_alt_end = False
                    for struct in list(self.similarity_table):
                        if struct!=main_structure:
                            try:
                                if skip_template and struct.pdb_code.index in self.templates_to_skip[ref_seg][1]:
                                    continue
                                protein_conf = ProteinConformation.objects.get(protein=struct.protein_conformation.protein.parent)
                                alt_helix_ends = self.fetch_struct_helix_ends_from_db(struct)
                                if parser.gn_comparer(alt_helix_ends[ref_seg][1],self.helix_ends[ref_seg][1],
                                                      protein_conf)>=0:
                                    all_keys = list(a.reference_dict[ref_seg].keys())[-1*(len(modifications['added'][ref_seg][1])+4):]
                                    ref_keys = [i for i in all_keys if i not in modifications['added'][ref_seg][1]]
                                    reference = parser.fetch_residues_from_array(main_pdb_array[ref_seg],ref_keys)
                                    template = parser.fetch_residues_from_pdb(struct,all_keys)
                                    superpose = sp.OneSidedSuperpose(reference,template,4,1)
                                    sup_residues = superpose.run()
                                    new_residues = OrderedDict()
                                    for gn, atoms in sup_residues.items():
                                        if gn.replace('.','x') not in ref_keys:
                                            new_residues[gn]=atoms
                                    for gn, atoms in new_residues.items():
                                        gn_ = gn.replace('.','x')
                                        if gn_ in modifications['added'][ref_seg][1]:
                                            main_pdb_array[ref_seg][gn] = atoms
                                            a.template_dict[ref_seg][gn_] = PDB.Polypeptide.three_to_one(
                                                                            atoms[0].get_parent().get_resname())
                                            if a.template_dict[ref_seg][gn_]==a.reference_dict[ref_seg][gn_]:
                                                a.alignment_dict[ref_seg][gn_] = a.reference_dict[ref_seg][gn_]
                                            else:
                                                a.alignment_dict[ref_seg][gn_] = '.'
                                    self.template_source = update_template_source(self.template_source,modifications['added'][ref_seg][1],
                                                                struct,segment=ref_seg)
                                    found_alt_end = True
                                    break
                            except:
                                pass
                    if found_alt_end==False:
                        for i in modifications['added'][ref_seg][1]:
                            main_pdb_array[ref_seg][i.replace('x','.')] = 'x'
                            a.template_dict[ref_seg][i] = 'x'
                            a.alignment_dict[ref_seg][i] = 'x'
        self.helix_end_mods = modifications
        self.main_pdb_array = main_pdb_array
        self.alignment = a
        if separate_H8==False and H8_alt.pdb_code.index in structures_with_num_issues:
            ssno = StructureSeqNumOverwrite(H8_alt)
            ssno.seq_num_overwrite('wt')
        if self.debug:
            print(self.helix_end_mods)
        return main_pdb_array, a



class Loops(object):
    ''' Class to handle loops in GPCR structures.
    '''
    def __init__(self, reference_protein, loop_label, loop_template_structures, main_structure, helix_end_mods,
                 segment_order, revise_xtal):
        self.segment_order = segment_order
        if revise_xtal==True:
            ref_prot = reference_protein.parent
        else:
            ref_prot = reference_protein
        self.reference_protein = ref_prot
        self.prot_conf = ProteinConformation.objects.get(protein=ref_prot)
        self.loop_label = loop_label
        self.loop_template_structures = loop_template_structures
        self.main_structure = main_structure
        self.helix_end_mods = helix_end_mods
        self.loop_output_structure = None
        self.new_label = None
        self.aligned = False
        self.model_loop = False
        self.partialECL2_1 = False
        self.partialECL2_2 = False
        self.excluded_loops = {'ICL1':[],'ECL1':[],'ICL2':['5ZKP'],'ECL2':[],'ECL2_1':['5NM4'],'ECL2_mid':[],'ECL2_2':['6K41'],'ICL3':['3VW7'],'ECL3':['4DJH','6KJV','6KK1','6KK7','5VEW']}
        self.evade_chain_break = False

    def fetch_loop_residues(self, main_pdb_array, superpose_modded_loop=False):
        ''' Fetch list of Atom objects of the loop when there is an available template. Returns an OrderedDict().
        '''
        #################
        # d2 = Structure.objects.get(pdb_code__index='7L1V')
        # d3 = Structure.objects.get(pdb_code__index='3PBL')
        # if self.loop_label=='ECL2':
        #     self.loop_template_structures['ECL2_1'] = [d2]
        #     # self.loop_template_structures['ECL2_mid'] = [d2]
            # self.loop_template_structures['ECL2_2'] = [d2]
        # if self.loop_label in ['ICL1', 'ICL2', 'ECL1']:
        #     self.loop_template_structures = OrderedDict([(Structure.objects.get(pdb_code__index='7F9Y'), 100)])
        # if self.loop_label in ['ICL2', 'ICL1', 'ECL1']:
        #     self.loop_template_structures = OrderedDict([(d2, 100)])
        # if self.loop_label=='ECL3':
        #     self.loop_template_structures = OrderedDict([(d2, 100)])
        # if self.loop_label=='ICL2':
        #     self.loop_template_structures = OrderedDict([(self.main_structure, 100)])
        ## Kappa opioid
        # if self.loop_label in ['ICL1', 'ICL2', 'ICL3']:
        #     self.loop_template_structures = OrderedDict([(Structure.objects.get(pdb_code__index='6DDE'), 100)])
        #################

        if (self.loop_label=='ECL2' and (self.loop_template_structures==None or 'ECL2_mid' in
            self.loop_template_structures and self.loop_template_structures['ECL2_mid']==None)):
            return None
        if self.loop_template_structures!=None:
            ref_loop = list(Residue.objects.filter(protein_conformation=self.prot_conf,
                                                   protein_segment__slug=self.loop_label))
            parse = GPCRDBParsingPDB()
            seg_list = self.segment_order
            prev_seg = seg_list[seg_list.index(self.loop_label)-1]
            next_seg = seg_list[seg_list.index(self.loop_label)+1]
            if prev_seg=='C-term':
                orig_before_gns = []
            else:
                orig_before_gns = [i.replace('.','x') for i in list(main_pdb_array[prev_seg].keys())[-4:]]
            orig_after_gns = [j.replace('.','x') for j in list(main_pdb_array[next_seg].keys())[:4]]
            if len(orig_before_gns)==0:
                last_before_gn = None
            else:
                last_before_gn = orig_before_gns[-1]
            first_after_gn = orig_after_gns[0]
            if self.loop_label=='ECL2':
                try:
                    ref_res = Residue.objects.filter(protein_conformation__protein=self.reference_protein,
                                                     protein_segment__slug='ECL2')
                    r_first = list(ref_res)[0].sequence_number
                    r_last = list(ref_res)[-1].sequence_number
                    r_x50 = ref_res.get(display_generic_number__label='45.50x50').sequence_number
                except:
                    pass
            output = OrderedDict()
            if (self.loop_label=='ECL2' and 'ECL2_1' not in self.loop_template_structures) or self.loop_label!='ECL2' or superpose_modded_loop==True:
                for template in self.loop_template_structures:
                    if self.loop_label=='ICL2' and template!='aligned' and template.pdb_code.index=='2RH1' and self.reference_protein.entry_name=='adrb2_human':
                        continue
                    output = OrderedDict()
                    try:
                        if (template==self.main_structure or template=='aligned') and superpose_modded_loop==False and template.pdb_code.index not in self.excluded_loops[self.loop_label]:
                            if self.helix_end_mods!=None and (len(self.helix_end_mods['removed'][prev_seg][1])==0 and
                                                              len(self.helix_end_mods['removed'][next_seg][0])==0 and
                                                              len(self.helix_end_mods['added'][prev_seg][1])==0 and
                                                              len(self.helix_end_mods['added'][next_seg][0])==0):
                                if template=='aligned':
                                    self.aligned = True
                                else:
                                    self.aligned = False
                                try:
                                    l_res = self.compare_parent_loop_to_child(self.loop_label,template)
                                    if l_res==False:
                                        raise Exception()
                                    loop_res = [r.sequence_number for r in l_res[1]]
                                    at_least_one_gn = False
                                    x50_present, x50_temp_present = False, False
                                    for i in ref_loop:
                                        try:
                                            g = ggn(i.display_generic_number.label)
                                            at_least_one_gn = True
                                            if 'x50' in g:
                                                x50_present = True
                                                break
                                        except:
                                            pass
                                    for j in l_res[1]:
                                        try:
                                            g2 = ggn(j.display_generic_number.label)
                                            if 'x50' in g2:
                                                x50_temp_present = True
                                        except:
                                            pass
                                    if self.loop_template_structures[template]!=0:
                                        if x50_present==False and len(ref_loop)!=len(loop_res):
                                            continue
                                        elif x50_present==True and x50_temp_present==False and len(self.loop_template_structures)>1:
                                            continue
                                        partial = False
                                    else:
                                        partial = True
                                        if len(self.helix_end_mods['added'][prev_seg][1])!=0 or len(self.helix_end_mods['added'][next_seg][0])!=0:
                                            continue
                                    if at_least_one_gn==True:
                                        inter_array = parse.fetch_residues_from_pdb(self.main_structure,loop_res)
                                    else:
                                        inter_array = parse.fetch_residues_from_pdb(self.main_structure,loop_res,
                                                                                    just_nums=True)
                                    self.loop_output_structure = self.main_structure
                                    if partial==False:
                                        for id_, atoms in inter_array.items():
                                            output[str(id_)] = atoms
                                    else:
                                        p_c = ProteinConformation.objects.get(protein=self.main_structure.protein_conformation.protein.parent)
                                        p_loop_res = Residue.objects.filter(protein_conformation=p_c,
                                                                             protein_segment__slug=self.loop_label)
                                        for num in p_loop_res:
                                            try:
                                                output[str(num.sequence_number)] = inter_array[str(num.sequence_number)]
                                            except:
                                                output[str(num.sequence_number)] = '-'
                                    return output
                                except:
                                    self.aligned = False
                                    continue
                            else:
                                return self.fetch_loop_residues(main_pdb_array,superpose_modded_loop=True)
                        else:
                            if template.pdb_code.index in self.excluded_loops[self.loop_label]:
                                continue
                            if self.loop_label=='ICL4' and len(list(Residue.objects.filter(protein_conformation=self.prot_conf,protein_segment__slug='ICL4')))<3:
                                raise Exception()
                            if template=='aligned' or template==self.main_structure:
                                template = self.main_structure
                                self.aligned = True
                            if superpose_modded_loop==True:
                                self.model_loop = True
                                alt_last_before_gn = last_before_gn
                                b_num_found = False
                                break_count = 0
                                while b_num_found==False and break_count<30:
                                    try:
                                        b_num = Residue.objects.get(protein_conformation=template.protein_conformation,
                                                                    display_generic_number__label=dgn(alt_last_before_gn,
                                                                                                      template.protein_conformation)).sequence_number
                                        b_num_found = True
                                    except:
                                        alt_last_before_gn = parse.gn_indecer(alt_last_before_gn,'x',-1)
                                        break_count+=1
                                alt_first_after_gn = first_after_gn
                                a_num_found = False
                                break_count = 0
                                while a_num_found==False and break_count<30:
                                    try:
                                        a_num = Residue.objects.get(protein_conformation=template.protein_conformation,
                                                                    display_generic_number__label=dgn(alt_first_after_gn,
                                                                                                      template.protein_conformation)).sequence_number
                                        a_num_found = True
                                    except:
                                        alt_first_after_gn = parse.gn_indecer(alt_first_after_gn,'x',1)
                                        break_count+=1
                            else:
                                b_num = Residue.objects.get(protein_conformation=template.protein_conformation,
                                                            display_generic_number__label=dgn(last_before_gn,
                                                                    template.protein_conformation)).sequence_number
                                a_num = Residue.objects.get(protein_conformation=template.protein_conformation,
                                                            display_generic_number__label=dgn(first_after_gn,
                                                                    template.protein_conformation)).sequence_number
                            before4 = Residue.objects.filter(protein_conformation=template.protein_conformation,
                                                             sequence_number__in=[b_num,b_num-1,b_num-2,b_num-3])
                            after4 = Residue.objects.filter(protein_conformation=template.protein_conformation,
                                                             sequence_number__in=[a_num,a_num+1,a_num+2,a_num+3])
                            x50_present = False

                            for i in ref_loop:
                                try:
                                    if 'x50' in i.display_generic_number.label:
                                        x50_present = True
                                except:
                                    pass
                            if superpose_modded_loop==True and self.aligned==True:

                                loop_residues = Residue.objects.filter(protein_conformation=template.protein_conformation,
                                                                       protein_segment__slug=self.loop_label)
                                p_c_check = self.compare_parent_loop_to_child(self.loop_label,template)
                                if p_c_check==False:
                                    raise Exception()
                                else:
                                    del_from_loop = []
                                    for i in loop_residues:
                                        if i not in p_c_check[1]:
                                            del_from_loop.append(i)
                                    loop_residues = [i for i in loop_residues if i not in del_from_loop]
                                if (self.loop_label in ['ICL1','ECL1','ICL2'] and not x50_present and
                                    len(loop_residues)!=len(ref_loop)):
                                    raise Exception()
                            else:
                                loop_residues = Residue.objects.filter(protein_conformation=template.protein_conformation, protein_segment__slug=self.loop_label)
                                loop_residues_test = Residue.objects.filter(protein_conformation=template.protein_conformation,
                                                                            protein_segment__slug=self.loop_label)
                                p_c = ProteinConformation.objects.get(protein=template.protein_conformation.protein.parent)
                                loop_residues_test_parent = Residue.objects.filter(protein_conformation=p_c,
                                                                                   protein_segment__slug=self.loop_label)
                                gn_nums_loop = [i for i in loop_residues if i.generic_number!=None and i.protein_segment.slug[0] not in ['T','H']]
                                if self.loop_label in ['ICL1','ECL1','ICL2'] and x50_present==True:
                                    if len(gn_nums_loop)==0 or len(gn_nums_loop)<len([i for i in ref_loop if i.generic_number!=None and i.protein_segment.slug[0] not in ['T','H']]):
                                        continue
                                    else:
                                        loop_residues = loop_residues_test
                                if len(loop_residues_test)!=len(loop_residues_test_parent):
                                    continue
                            before_gns = [x.sequence_number for x in before4]
                            mid_nums = [x.sequence_number for x in loop_residues]
                            after_gns = [x.sequence_number for x in after4]
                            alt_residues = parse.fetch_residues_from_pdb(template, before_gns+mid_nums+after_gns)
                            orig_residues1 = parse.fetch_residues_from_array(main_pdb_array[prev_seg],orig_before_gns)
                            orig_residues2 = parse.fetch_residues_from_array(main_pdb_array[next_seg],orig_after_gns)
                            orig_residues = parse.add_two_ordereddict(orig_residues1,orig_residues2)
                            superpose = sp.LoopSuperpose(orig_residues, alt_residues)
                            new_residues = superpose.run()
                            key_list = list(new_residues.keys())[4:-4]
                            for key in key_list:
                                output[key] = new_residues[key]
                            self.loop_output_structure = template
                            return output
                    except:
                        self.aligned = False
                        continue
            else:
                output,ECL2_1,ECL2_mid,ECL2_2 = OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict()
                no_first_temp, no_second_temp = True,True
                main_temp_seq = Residue.objects.filter(protein_conformation=self.main_structure.protein_conformation,
                                                       protein_segment__slug=self.loop_label)
                cp = self.compare_parent_loop_to_child(self.loop_label, self.main_structure)
                try:
                    if cp[0]==True:
                        main_temp_seq = cp[1]
                except:
                    main_temp_seq = []
                for mid_template in self.loop_template_structures['ECL2_mid']:
                    ref_ECL2_res = Residue.objects.filter(protein_conformation=self.prot_conf, protein_segment__slug='ECL2')
                    ref_ECL2_45x50 = ref_ECL2_res.get(display_generic_number__label='45.50x50')
                    ref_ECL2_2_res = [i for i in ref_ECL2_res if i.sequence_number>ref_ECL2_45x50.sequence_number+2]
                    if type(main_temp_seq)==type([]):
                        if len(main_temp_seq)==0:
                            temp_ECL2_45x50 = []
                            temp_ECL2_2_res = []
                        else:
                            temp_ECL2_45x50 = [i for i in main_temp_seq if i.display_generic_number!=None and i.display_generic_number.label=='45.50x50']
                            if len(temp_ECL2_45x50)==0:
                                temp_ECL2_2_res = []
                            else:
                                temp_ECL2_45x50 = temp_ECL2_45x50[0]
                                temp_ECL2_2_res = [i for i in main_temp_seq if i.sequence_number>temp_ECL2_45x50.sequence_number+2]
                    else:
                        temp_ECL2_45x50 = main_temp_seq.get(display_generic_number__label='45.50x50')
                        temp_ECL2_2_res = [i for i in main_temp_seq if i.sequence_number>temp_ECL2_45x50.sequence_number+2]
                    if len(temp_ECL2_2_res)>len(ref_ECL2_2_res) and len(ref_ECL2_2_res)<2:
                        self.evade_chain_break = True
                    if mid_template==self.main_structure:
                        ECL2_mid = parse.fetch_residues_from_pdb(self.main_structure,['45x50','45x51','45x52'])
                        try:
                            x50 = [i.sequence_number for i in main_temp_seq if i.display_generic_number!=None and i.display_generic_number.label=='45.50x50'][0]
                        except:
                            x50 = main_temp_seq.get(display_generic_number__label='45.50x50').sequence_number
                        break
                    else:
                        try:
                            ECL2_mid = parse.fetch_residues_from_pdb(mid_template,[last_before_gn,first_after_gn,'3x25',
                                                                                   '45x50','45x51','45x52'])
                            ref_ECL2_mid1 = parse.fetch_residues_from_array(main_pdb_array['TM4'],[last_before_gn])
                            ref_ECL2_mid2 = parse.fetch_residues_from_array(main_pdb_array['TM5'],[first_after_gn])
                            ref_ECL2_mid3 = parse.fetch_residues_from_array(main_pdb_array['TM3'],['3x25'])
                            ref_ECL2_mid = parse.add_two_ordereddict(parse.add_two_ordereddict(ref_ECL2_mid1,
                                                                                               ref_ECL2_mid2),
                                                                                               ref_ECL2_mid3)
                            superpose = sp.ECL2MidSuperpose(ref_ECL2_mid,ECL2_mid)
                            new_mid_residues = superpose.run()
                            ECL2_mid = OrderedDict()
                            for i,j in new_mid_residues.items():
                                if i in ['45.50','45.51','45.52']:
                                    ECL2_mid[i] = j
                            break
                        except:
                            continue

                o1 = parse.fetch_residues_from_array(main_pdb_array[prev_seg],orig_before_gns)
                orig_residues1 = parse.add_two_ordereddict(o1,ECL2_mid)
                if self.loop_template_structures['ECL2_1']==None:
                    no_first_temp=True
                else:
                    for first_temp in self.loop_template_structures['ECL2_1']:
                        if first_temp.pdb_code.index in self.excluded_loops['ECL2_1']:
                            continue
                        if first_temp==self.main_structure:
                            try:
                                ECL2_1 = parse.fetch_residues_from_pdb(self.main_structure,
                                                                       list(range(list(main_temp_seq)[0].sequence_number,x50)))
                                no_first_temp=False
                                break
                            except:
                                try:
                                    partial_seq1 = Residue.objects.filter(protein_conformation=first_temp.protein_conformation,
                                                                          sequence_number__in=list(range(list(main_temp_seq)[0].sequence_number,x50)))
                                    partial_seq1_nums = [i.sequence_number for i in partial_seq1]
                                    ECL2_1 = parse.fetch_residues_from_pdb(first_temp, partial_seq1_nums)
                                    no_first_temp=False
                                    self.partialECL2_1 = True
                                    break
                                except:
                                    continue
                        else:
                            try:
                                b_num = Residue.objects.get(protein_conformation=first_temp.protein_conformation,
                                                            display_generic_number__label=dgn(last_before_gn,
                                                                                              first_temp.protein_conformation)).sequence_number
                                before4 = Residue.objects.filter(protein_conformation=first_temp.protein_conformation,
                                                                 sequence_number__in=[b_num,b_num-1,b_num-2,b_num-3])
                                alt_mid1 = Residue.objects.filter(protein_conformation=first_temp.protein_conformation,
                                                                  protein_segment__slug=self.loop_label,
                                                                  display_generic_number__label__in=['45.50x50','45.51x51','45.52x52'])
                                alt1_x50 = alt_mid1.get(display_generic_number__label='45.50x50').sequence_number
                                loop_res1 = Residue.objects.filter(protein_conformation=first_temp.protein_conformation,
                                                                   sequence_number__in=list(range(b_num+1, alt1_x50))).filter(protein_segment__slug=self.loop_label)
                                before_gns = [x.sequence_number for x in before4]
                                mid_gns1 = [x.sequence_number for x in loop_res1]
                                alt_residues1 = parse.fetch_residues_from_pdb(first_temp,before_gns+mid_gns1+['45x50','45x51','45x52'])
                                superpose = sp.LoopSuperpose(orig_residues1,alt_residues1,ECL2=True,part=1)
                                new_residues = superpose.run()
                                key_list = list(new_residues.keys())[4:-3]
                                for key in key_list:
                                    ECL2_1["1_"+key] = new_residues[key]
                                no_first_temp=False
                                break
                            except:
                                no_first_temp=True

                if no_first_temp==True:
                    for i in range(1,r_x50-r_first+1):
                        ECL2_1['1_'+str(i)]='x'
                    first_temp=None
                o2 = parse.fetch_residues_from_array(main_pdb_array[next_seg],orig_after_gns)
                orig_residues2 = parse.add_two_ordereddict(ECL2_mid,o2)

                if self.loop_template_structures['ECL2_2']==None:
                    no_second_temp=True
                else:
                    for second_temp in self.loop_template_structures['ECL2_2']:
                        if second_temp.pdb_code.index in self.excluded_loops['ECL2_2']:
                            continue
                        if second_temp==self.main_structure:
                            try:
                                ECL2_2 = parse.fetch_residues_from_pdb(self.main_structure,list(range(x50+3,list(main_temp_seq)[-1].sequence_number+1)))
                                no_second_temp=False
                                break
                            except:
                                try:
                                    partial_seq2 = Residue.objects.filter(protein_conformation=second_temp.protein_conformation,
                                                                          sequence_number__in=list(range(x50+3,list(main_temp_seq)[-1].sequence_number+1)))
                                    partial_seq2_nums = [i.sequence_number for i in partial_seq2]
                                    ECL2_2 = parse.fetch_residues_from_pdb(second_temp, partial_seq2_nums)
                                    no_second_temp=False
                                    self.partialECL2_2 = True
                                    break
                                except:
                                    continue
                        else:
                            try:
                                a_num = Residue.objects.get(protein_conformation=second_temp.protein_conformation,
                                                            display_generic_number__label=dgn(first_after_gn,
                                                                                              second_temp.protein_conformation)).sequence_number
                                after4 = Residue.objects.filter(protein_conformation=second_temp.protein_conformation,
                                                                sequence_number__in=[a_num,a_num+1,a_num+2,a_num+3])
                                alt_mid2 = Residue.objects.filter(protein_conformation=second_temp.protein_conformation,
                                                                  protein_segment__slug=self.loop_label,
                                                                  display_generic_number__label__in=['45.50x50','45.51x51','45.52x52'])
                                alt2_x50 = alt_mid2.get(display_generic_number__label='45.50x50').sequence_number
                                loop_res2 = Residue.objects.filter(protein_conformation=second_temp.protein_conformation,
                                                                   sequence_number__in=list(range(alt2_x50+3, a_num))).filter(protein_segment__slug=self.loop_label)
                                mid_gns2 = [x.sequence_number for x in loop_res2]
                                after_gns = [x.sequence_number for x in after4]
                                alt_residues2 = parse.fetch_residues_from_pdb(second_temp,['45x50','45x51','45x52']+mid_gns2+after_gns)
                                superpose = sp.LoopSuperpose(orig_residues2,alt_residues2,ECL2=True,part=2)
                                new_residues = superpose.run()
                                key_list = list(new_residues.keys())[3:-4]
                                for key in key_list:
                                    ECL2_2["2_"+key] = new_residues[key]
                                no_second_temp=False
                                break
                            except:
                                no_second_temp=True
                if no_second_temp==True:
                    for j in range(1,r_last-r_x50-1):
                        ECL2_2['2_'+str(j)]='x'
                    second_temp=None
                output['ECL2_1'] = ECL2_1
                output['ECL2_mid'] = ECL2_mid
                output['ECL2_2'] = ECL2_2
                self.loop_output_structure = [first_temp,mid_template,second_temp]
                return output
            if len(output.keys())==0:
                return None
        else:
            return None

    def insert_loop_to_arrays(self, loop_output_structure, main_pdb_array, loop_template, reference_dict,
                              template_dict, alignment_dict):
        ''' Updates the homology model with loop segments. Inserts previously fetched lists of loop Atom objects to
            the proper arrays, dictionaries.

            @param loop_output_structure: Structure object of loop template.
            @param main_pdb_array: nested OrderedDict(), output of GPCRDBParsingPDB().pdb_array_creator().
            @param loop_template: OrderedDict() of loop template with lists of Atom objects as values.
            @param reference_dict: reference dictionary of AlignedReferenceTemplate.
            @param template_dict: template dictionary of AlignedReferenceTemplate.
            @param alignment_dict: alignment dictionary of AlignedReferenceTemplate.
        '''
        shorter_ref, shorter_temp = False, False
        try:
            for r,t in zip(reference_dict[self.loop_label],template_dict[self.loop_label]):
                if reference_dict[self.loop_label][r] in ['-','x']:
                    shorter_ref = True
                    self.model_loop = True
                elif template_dict[self.loop_label][t] in ['-','x']:
                    shorter_temp = True
                    self.model_loop = True
        except:
            pass
        if loop_template!=None and loop_output_structure!=self.main_structure:
            loop_keys = list(loop_template.keys())[1:-1]
            continuous_loop = False
            self.main_pdb_array = self.discont_loop_insert_to_pdb(main_pdb_array, loop_template, loop_output_structure,
                                                                  temp_dict=template_dict)
        elif loop_template!=None and loop_output_structure==self.main_structure or self.aligned==True and (shorter_ref==True or shorter_temp==True):
            loop_keys = list(loop_template.keys())
            continuous_loop = True
            temporary_dict = OrderedDict()
            # correct for partial loops
            try:
                query_temp_dict = template_dict[self.loop_label]
            except:
                query_temp_dict = OrderedDict()
                parent_loop = Residue.objects.filter(protein_conformation__protein=loop_output_structure.protein_conformation.protein.parent, protein_segment__slug=self.loop_label)
                for i in parent_loop:
                    query_temp_dict[str(i.sequence_number)] = '-'
            try:
                if len(loop_keys)<len(query_temp_dict):
                    counter=0
                    for i in query_temp_dict:
                        if i.replace('x','.') in loop_keys:
                            temporary_dict[i.replace('x','.')] = loop_template[i.replace('x','.')]
                        else:
                            temporary_dict['gap{}'.format(str(counter))] = '-'
                        counter+=1
                    loop_template = temporary_dict
            except:
                pass
            self.main_pdb_array = self.cont_loop_insert_to_pdb(main_pdb_array, template_dict, loop_template)
        else:
            self.main_pdb_array = main_pdb_array

        if loop_template!=None:
            temp_ref_dict, temp_temp_dict, temp_aligned_dict = OrderedDict(),OrderedDict(),OrderedDict()
            if continuous_loop==True:
                if shorter_ref==True and shorter_temp==False:
                    ref_residues = list(reference_dict[self.loop_label].values())
                elif shorter_ref==True and shorter_temp==True:
                    ref_residues = list(reference_dict[self.loop_label].values())
                elif shorter_ref==False and shorter_temp==True:
                    ref_residues = list(reference_dict[self.loop_label].values())
                else:
                    ref_residues = [x.amino_acid for x in Residue.objects.filter(protein_conformation__protein=self.reference_protein,
                                                           protein_segment__slug=self.loop_label)]
            else:
                try:
                    ref_residues = list(reference_dict[self.loop_label].values())
                except:
                    ref_residues = list(Residue.objects.filter(protein_conformation__protein=self.reference_protein,
                                                               protein_segment__slug=self.loop_label))

            for ref_seg, temp_seg, aligned_seg in zip(reference_dict, template_dict, alignment_dict):
                if ref_seg[0]=='T' and self.segment_order.index(self.loop_label)-self.segment_order.index(ref_seg[:4])==1:
                    temp_ref_dict[ref_seg] = reference_dict[ref_seg]
                    temp_temp_dict[temp_seg] = template_dict[temp_seg]
                    temp_aligned_dict[aligned_seg] = alignment_dict[aligned_seg]
                    input_residues = list(loop_template.keys())
                    ref_loop_seg, temp_loop_seg, aligned_loop_seg = OrderedDict(),OrderedDict(),OrderedDict()
                    if continuous_loop==True:
                        l_res=0
                        for r_res, r_id in zip(ref_residues, input_residues):
                            l_res+=1
                            try:
                                loop_gn = ggn(Residue.objects.get(protein_conformation=self.main_structure.protein_conformation,
                                                                 display_generic_number__label=dgn(r_id.replace('.','x'),
                                                                                                   self.main_structure.protein_conformation)).display_generic_number.label)
                            except:
                                try:
                                    Residue.objects.get(protein_conformation=self.main_structure.protein_conformation,
                                                        sequence_number=r_id)
                                    loop_gn = self.loop_label+'|'+str(l_res)
                                except:
                                    loop_gn = self.loop_label+'?'+str(l_res)
                            ref_loop_seg[loop_gn] = r_res
                            try:
                                temp_loop_seg[loop_gn] = PDB.Polypeptide.three_to_one(loop_template[r_id][0].get_parent().get_resname())
                            except:
                                temp_loop_seg[loop_gn] = '-'
                            if ref_loop_seg[loop_gn]==temp_loop_seg[loop_gn]:
                                aligned_loop_seg[loop_gn] = ref_loop_seg[loop_gn]
                            elif ref_loop_seg[loop_gn]=='-' or temp_loop_seg[loop_gn]=='-':
                                aligned_loop_seg[loop_gn] = '-'
                            else:
                                aligned_loop_seg[loop_gn] = '.'
                        self.new_label = self.loop_label+'_cont'
                        temp_ref_dict[self.loop_label+'_cont'] = ref_loop_seg
                        temp_temp_dict[self.loop_label+'_cont'] = temp_loop_seg
                        temp_aligned_dict[self.loop_label+'_cont'] = aligned_loop_seg
                    else:
                        missing_indeces = []
                        try:
                            if len(list(template_dict[self.loop_label].keys()))<len(input_residues):
                                for i in input_residues:
                                    try:
                                        template_dict[self.loop_label][i.replace('.','x')]
                                    except:
                                        missing_indeces.append([i,input_residues.index(i)])
                                temp_ref_residues = []

                                for i in range(len(ref_residues)):
                                    for j in missing_indeces:
                                        if i==j[1]:
                                            temp_ref_residues.append(j[0])
                                    temp_ref_residues.append(ref_residues[i])
                                ref_residues = temp_ref_residues
                                if len(list(template_dict[self.loop_label].keys()))<len(input_residues):
                                    for j in missing_indeces:
                                        if j[1]>=len(ref_residues):
                                            ref_residues.append[j[0]]
                            elif len(list(template_dict[self.loop_label].keys()))>len(input_residues):
                                for i in list(template_dict[self.loop_label].keys()):
                                    try:
                                        input_residues.index(i.replace('x','.'))
                                    except:
                                        missing_indeces.append([i,list(template_dict[self.loop_label].keys()).index(i)])
                                temp_input_residues = []

                                if self.loop_label in ['ICL1', 'ECL1', 'ICL2']:
                                    for i in range(0,len(list(template_dict[self.loop_label].keys()))):
                                        gap_inserted = False
                                        for j in missing_indeces:
                                            if i==j[1]:
                                                temp_input_residues.append(j[0])
                                                gap_inserted = True
                                                break
                                        if gap_inserted==False:
                                            temp_input_residues.append(list(template_dict[self.loop_label])[i].replace('x','.'))
                                    input_residues = temp_input_residues
                                    if len(list(template_dict[self.loop_label].keys()))>len(input_residues):
                                        for j in missing_indeces:
                                            if j[1]>=len(input_residues):
                                                input_resi
                                else:
                                    for i in range(len(input_residues)):
                                        for j in missing_indeces:
                                            if i==j[1]:
                                                temp_input_residues.append(j[0])
                                        temp_input_residues.append(input_residues[i])
                                    input_residues = temp_input_residues
                                    if len(list(template_dict[self.loop_label].keys()))>len(input_residues):
                                        for j in missing_indeces:
                                            if j[1]>=len(input_residues):
                                                input_residues.append(j[0])
                        except:
                            pass
                        loop_ends = []
                        start_found, end_found = False, False
                        for i in range(len(ref_residues)):
                            if ref_residues[i]!='x' and start_found==False:
                                loop_ends.append(i)
                                start_found = True
                            if ref_residues[(i+1)*-1]!='x' and end_found==False:
                                loop_ends.append((i+1)*-1)
                                end_found = True
                            if start_found and end_found:
                                break
                        l_res=0
                        temp_pdb_dict = OrderedDict()
                        for r_res, r_id in zip(ref_residues, input_residues):
                            if l_res==loop_ends[0] or len(ref_residues)+loop_ends[1]==l_res:
                                try:
                                    ref_loop_seg[self.loop_label+'?'+str(l_res+1)] = ref_residues[l_res].amino_acid
                                except:
                                    ref_loop_seg[self.loop_label+'?'+str(l_res+1)] = ref_residues[l_res]
                                temp_loop_seg[self.loop_label+'?'+str(l_res+1)] = 'x'
                                aligned_loop_seg[self.loop_label+'?'+str(l_res+1)] = 'x'
                                temp_pdb_dict[self.loop_label+'?'+str(l_res+1)] = 'x'
                                l_res+=1
                                continue
                            l_res+=1
                            try:
                                try:
                                    loop_gn = ggn(Residue.objects.get(protein_conformation=loop_output_structure.protein_conformation,
                                                  display_generic_number__label=dgn(r_id.replace('.','x'),
                                                  loop_output_structure.protein_conformation)).display_generic_number.label)
                                    ggn(Residue.objects.get(protein_conformation=self.prot_conf,
                                                            display_generic_number__label=dgn(loop_gn,
                                                                                              self.prot_conf)).display_generic_number.label)
                                except:
                                    loop_gn = ggn(Residue.objects.get(protein_conformation=loop_output_structure.protein_conformation,
                                                  sequence_number=r_id).display_generic_number.label)
                                if len(loop_gn.split('x')[0])==1:
                                    raise Exception()
                            except:
                                try:
                                    if template_dict[self.loop_label][r_id.replace('.','x')] in ['x','-'] or r_res=='x':
                                        loop_gn = self.loop_label+'?'+str(l_res)
                                    else:
                                        raise Exception()
                                except:
                                    loop_gn = self.loop_label+'|'+str(l_res)
                            try:
                                ref_loop_seg[loop_gn] = r_res.amino_acid
                            except:
                                ref_loop_seg[loop_gn] = r_res
                            try:
                                temp_loop_seg[loop_gn] = PDB.Polypeptide.three_to_one(loop_template[r_id][0].get_parent().get_resname())
                            except:
                                temp_loop_seg[loop_gn] = '-'
                            if ref_loop_seg[loop_gn]==temp_loop_seg[loop_gn]:
                                aligned_loop_seg[loop_gn] = ref_loop_seg[loop_gn]
                            else:
                                aligned_loop_seg[loop_gn] = '.'
                            try:
                                temp_pdb_dict[loop_gn.replace('x','.')] = loop_template[r_id]
                            except:
                                temp_pdb_dict[loop_gn.replace('x','.')] = '-'
                        self.new_label = self.loop_label+'_dis'
                        temp_ref_dict[self.new_label] = ref_loop_seg
                        temp_temp_dict[self.new_label] = temp_loop_seg
                        temp_aligned_dict[self.new_label] = aligned_loop_seg
                        self.main_pdb_array[self.new_label] = temp_pdb_dict
                else:
                    temp_ref_dict[ref_seg] = reference_dict[ref_seg]
                    temp_temp_dict[temp_seg] = template_dict[temp_seg]
                    temp_aligned_dict[aligned_seg] = alignment_dict[aligned_seg]
            self.reference_dict = temp_ref_dict
            self.template_dict = temp_temp_dict
            self.alignment_dict = temp_aligned_dict
            try:
                del self.reference_dict[self.loop_label]
                del self.template_dict[self.loop_label]
                del self.alignment_dict[self.loop_label]
            except:
                pass
        else:
            self.reference_dict = reference_dict
            self.template_dict = template_dict
            self.alignment_dict = alignment_dict
            try:
                del self.reference_dict[self.loop_label]
                del self.template_dict[self.loop_label]
                del self.alignment_dict[self.loop_label]
            except:
                pass
        return self

    def insert_ECL2_to_arrays(self, loop_output_structure, main_pdb_array, loop_template, reference_dict,
                              template_dict, alignment_dict, partialECL2_1=False, partialECL2_2=False):
        temp_array = OrderedDict()
        parent = ProteinConformation.objects.get(protein=loop_output_structure[1].protein_conformation.protein.parent)
        seq = list(Residue.objects.filter(protein_conformation=parent, protein_segment__slug='ECL2'))
        x50 = [i for i in seq if i.generic_number!=None and i.generic_number.label=='45x50'][0]
        x50_i = seq.index(x50)

        # first part
        if loop_output_structure[0]!=None:
            if loop_output_structure[0]==self.main_structure:
                temp_array = self.cont_loop_insert_to_pdb(main_pdb_array, template_dict, loop_template['ECL2_1'],
                                                          ECL2='', x50_i=x50_i)
            else:
                ref_seq1 = list(Residue.objects.filter(protein_conformation__protein=self.reference_protein, protein_segment__slug='ECL2'))
                ref_x50 = [i for i in ref_seq1 if i.generic_number!=None and i.generic_number.label=='45x50'][0]
                ref_x50i = ref_seq1.index(ref_x50)
                if loop_output_structure[0]!=loop_output_structure[1]:
                    parent = ProteinConformation.objects.get(protein=loop_output_structure[0].protein_conformation.protein.parent)
                    seq = list(Residue.objects.filter(protein_conformation=parent, protein_segment__slug='ECL2'))
                t_dict1 = OrderedDict([('ECL2',OrderedDict())])
                for i in seq[:ref_x50i]:
                    # if i.sequence_number<x50.sequence_number:
                    t_dict1['ECL2']['1_'+str(i.sequence_number)] = '-'
                temp_array = self.discont_loop_insert_to_pdb(main_pdb_array, loop_template['ECL2_1'],
                                                             loop_output_structure, ECL2='', temp_dict=t_dict1)
        else:
            temp_array = self.gap_ECL2(main_pdb_array,loop_template['ECL2_1'],break_chain=True)
        # middle part
        for key, res in loop_template['ECL2_mid'].items():
            temp_array['ECL2'][key] = res
        # second part
        l_res = len(temp_array['ECL2'])
        if loop_output_structure[2]!=None:
            if loop_output_structure[2]==self.main_structure:
                if partialECL2_2==True:
                    for key in list(template_dict['ECL2'])[x50_i+3:]:
                        l_res+=1
                        if key in loop_template['ECL2_2']:
                            temp_array['ECL2'][self.loop_label+'|'+str(l_res)] = loop_template['ECL2_2'][key]
                        else:
                            temp_array['ECL2'][self.loop_label+'?'+str(l_res)] = '-'
                else:
                    for key, res in loop_template['ECL2_2'].items():
                        l_res+=1
                        if '.' in key:
                            temp_array['ECL2'][key] = res
                        else:
                            temp_array['ECL2'][self.loop_label+'|'+str(l_res)] = res
            else:
                loop_keys = list(loop_template['ECL2_2'].keys())[1:-1]
                temp_array['ECL2'][self.loop_label+'?'+str(l_res+1)] = 'x'
                l_res+=1
                if len(list(loop_template['ECL2_2'].keys()))>1:
                    for key in loop_keys:
                        l_res+=1
                        temp_array['ECL2'][self.loop_label+'|'+str(l_res)] = loop_template['ECL2_2'][key]
                    temp_array['ECL2'][self.loop_label+'?'+str(l_res+1)] = 'x'
        else:
            for key, res in loop_template['ECL2_2'].items():
                l_res+=1
                temp_array['ECL2'][self.loop_label+'?'+str(l_res)] = '-'

        self.main_pdb_array = temp_array
        temp_ref_dict, temp_temp_dict, temp_aligned_dict = OrderedDict(),OrderedDict(),OrderedDict()
        ref_residues = list(Residue.objects.filter(protein_conformation__protein=self.reference_protein,
                                                   protein_segment__slug='ECL2'))

        # correct for 1 res longer template
        ref_x50_i = ref_residues.index([i for i in ref_residues if i.generic_number!=None and ggn(i.display_generic_number.label)=='45x50'][0])
        if len(ref_residues)<len(self.main_pdb_array['ECL2']):
            if loop_output_structure[0]!=None:
                dif = len(list(self.main_pdb_array['ECL2'].keys())[:x50_i])-len(ref_residues[:ref_x50_i])
                ref_residues = ref_residues[:ref_x50_i] + list(dif*'-') + ref_residues[ref_x50_i:]
            if loop_output_structure[2]!=None:
                pdb_array_x50_i = list(self.main_pdb_array['ECL2'].keys()).index('45.50')
                dif = len(list(self.main_pdb_array['ECL2'].keys())[pdb_array_x50_i+3:])-len(ref_residues[ref_x50_i+3:])
                ref_residues = ref_residues[:ref_x50_i+3] + list(dif*'-') + ref_residues[ref_x50_i+3:]
        elif len(ref_residues)>len(self.main_pdb_array['ECL2']):
            if loop_output_structure[0]==self.main_structure:
                if len(seq)!=len(self.main_pdb_array['ECL2']):
                    x50_i = list(self.main_pdb_array['ECL2'].keys()).index('45.50')
                dif = len(ref_residues[:ref_x50_i])-len(list(self.main_pdb_array['ECL2'].keys())[:x50_i])
                temp_array = OrderedDict()
                for i in list(self.main_pdb_array['ECL2'])[:x50_i]:
                    temp_array[i] = self.main_pdb_array['ECL2'][i]
                if '|' in i:
                    s = i.split('|')
                elif '?':
                    s = i.split('?')
                num = int(s[1])
                for i in range(0,dif):
                    temp_array[s[0]+'?'+str(num+1)] = '-'
                    num+=1
                for i in list(self.main_pdb_array['ECL2'])[x50_i:x50_i+3]:
                    temp_array[i] = self.main_pdb_array['ECL2'][i]
                    num+=1
                if loop_output_structure[2]==self.main_structure:
                    dif = len(ref_residues[ref_x50_i+3:])-len(list(self.main_pdb_array['ECL2'].keys())[x50_i+3:])
                    for i in range(0,dif):
                        temp_array['ECL2'+'?'+str(num+1)] ='-'
                        num+=1
                for i in list(self.main_pdb_array['ECL2'])[x50_i+3:]:
                    temp_array['ECL2'+i[4]+str(num)] = self.main_pdb_array['ECL2'][i]
                    num+=1
                self.main_pdb_array['ECL2'] = temp_array

        for ref_seg, temp_seg, aligned_seg in zip(reference_dict, template_dict, alignment_dict):
            if ref_seg[0]=='T' and self.segment_order.index(self.loop_label)-self.segment_order.index(ref_seg[:4])==1:
                temp_ref_dict[ref_seg] = reference_dict[ref_seg]
                temp_temp_dict[temp_seg] = template_dict[temp_seg]
                temp_aligned_dict[aligned_seg] = alignment_dict[aligned_seg]
                temp_ref_dict['ECL2'],temp_temp_dict['ECL2'],temp_aligned_dict['ECL2'] = OrderedDict(),OrderedDict(),OrderedDict()

                for ref, key in zip(ref_residues, self.main_pdb_array['ECL2']):
                    try:
                        temp_ref_dict['ECL2'][key.replace('.','x')] = ref.amino_acid
                    except:
                        temp_ref_dict['ECL2'][key.replace('.','x')] = '-'
                    try:
                        temp_temp_dict['ECL2'][key.replace('.','x')] = PDB.Polypeptide.three_to_one(
                                                        self.main_pdb_array['ECL2'][key][0].get_parent().get_resname())
                    except:
                        temp_temp_dict['ECL2'][key.replace('.','x')] = self.main_pdb_array['ECL2'][key]
                    if temp_ref_dict['ECL2'][key.replace('.','x')]==temp_temp_dict['ECL2'][key.replace('.','x')]:
                        temp_aligned_dict['ECL2'][key.replace('.','x')] = temp_ref_dict['ECL2'][key.replace('.','x')]
                    elif temp_temp_dict['ECL2'][key.replace('.','x')]=='x':
                        temp_aligned_dict['ECL2'][key.replace('.','x')] = 'x'
                    elif temp_temp_dict['ECL2'][key.replace('.','x')]=='-':
                        temp_aligned_dict['ECL2'][key.replace('.','x')] = '-'
                    else:
                        temp_aligned_dict['ECL2'][key.replace('.','x')] = '.'
            else:
                if temp_seg=='ECL2':
                    continue
                temp_ref_dict[ref_seg] = reference_dict[ref_seg]
                temp_temp_dict[temp_seg] = template_dict[temp_seg]
                temp_aligned_dict[aligned_seg] = alignment_dict[aligned_seg]
        self.reference_dict = temp_ref_dict
        self.template_dict = temp_temp_dict
        self.alignment_dict = temp_aligned_dict
        return self

    def gap_ECL2(self, main_pdb_array, loop_template, break_chain=False):
        temp_array, temp_loop = OrderedDict(), OrderedDict()
        for seg_label, gns in main_pdb_array.items():
            if self.segment_order.index(self.loop_label)-self.segment_order.index(seg_label[:4])==1:
                temp_array[seg_label] = gns
                l_res = 0
                for key in loop_template:
                    l_res+=1
                    temp_loop[self.loop_label+'?'+str(l_res)] = '-'
                temp_array[self.loop_label] = temp_loop
            else:
                temp_array[seg_label] = gns
        return temp_array

    def compare_parent_loop_to_child(self, loop_label, structure, ignore_gns=False):
        l_res = list(Residue.objects.filter(protein_conformation=structure.protein_conformation,
                                            protein_segment__slug=loop_label))
        l_p_conf = ProteinConformation.objects.get(protein=structure.protein_conformation.protein.parent)
        parent_res = list(Residue.objects.filter(protein_conformation=l_p_conf,
                                                 protein_segment__slug=loop_label))
        parent_seq_nums = [i.sequence_number for i in parent_res]
        l_res_gn = [ggn(i.display_generic_number.label) for i in l_res if i.generic_number!=None]
        parent_res_gn = [ggn(i.display_generic_number.label) for i in parent_res if i.generic_number!=None]
        if len(l_res)!=len(parent_res) and l_res_gn==parent_res_gn:
            return True, [i for i in l_res if i.sequence_number in parent_seq_nums]
        elif l_res_gn!=parent_res_gn:
            return False
        else:
            return True, l_res

    def cont_loop_insert_to_pdb(self, main_pdb_array, template_dict, loop_template, ECL2=None, x50_i=None):
        temp_array, temp_loop = OrderedDict(), OrderedDict()
        for seg_label, gns in main_pdb_array.items():
            if self.segment_order.index(self.loop_label)-self.segment_order.index(seg_label[:4])==1:
                temp_array[seg_label] = gns
                l_res = 0
                if self.partialECL2_1==True:
                    for key in list(template_dict['ECL2'])[:x50_i]:
                        l_res+=1
                        if key in loop_template:
                            temp_loop[self.loop_label+'|'+str(l_res)] = loop_template[key]
                        else:
                            temp_loop[self.loop_label+'?'+str(l_res)] = '-'
                else:
                    for key in loop_template:
                        l_res+=1
                        if '.' in key:
                            temp_loop[key] = loop_template[key]
                        elif 'gap' in key:
                            temp_loop[self.loop_label+'?'+str(l_res)] = loop_template[key]
                        elif loop_template[key]=='-':
                            temp_loop[self.loop_label+'?'+str(l_res)] = loop_template[key]
                        else:
                            temp_loop[self.loop_label+'|'+str(l_res)] = loop_template[key]
                if ECL2!=None:
                    temp_array[self.loop_label] = temp_loop
                else:
                    temp_array[self.loop_label+'_cont'] = temp_loop
            else:
                temp_array[seg_label] = gns
        return temp_array

    def discont_loop_insert_to_pdb(self, main_pdb_array, loop_template, loop_output_structure, ECL2=None, temp_dict=None):
        temp_array, temp_loop = OrderedDict(), OrderedDict()
        loop_keys = list(loop_template.keys())[1:-1]
        for seg_label, gns in main_pdb_array.items():
            if self.segment_order.index(self.loop_label)-self.segment_order.index(seg_label[:4])==1:
                temp_array[seg_label] = gns
                l_res = 1
                temp_loop[self.loop_label+'?'+'1'] = 'x'
                if temp_dict!=None and self.loop_label in temp_dict:
                    iter_list = [i.replace('x','.') for i in temp_dict[self.loop_label]][1:-1]
                    if len(loop_keys)>len(iter_list):
                        iter_list = loop_keys
                else:
                    iter_list = loop_keys
                for key in iter_list:
                    l_res+=1
                    try:
                        try:
                            loop_gn = ggn(Residue.objects.get(protein_conformation=loop_output_structure.protein_conformation,
                                          display_generic_number__label=dgn(key.replace('.','x'),
                                          loop_output_structure.protein_conformation)).display_generic_number.label).replace('x','.')
                        except:
                            loop_gn = ggn(Residue.objects.get(protein_conformation=loop_output_structure.protein_conformation,
                                                             sequence_number=key).display_generic_number.label.replace('x','.'))
                        if len(loop_gn.split('.')[0])==1:
                            raise Exception()
                        if '.' in loop_gn:
                            Residue.objects.get(protein_conformation=self.prot_conf,
                                                display_generic_number__label=dgn(loop_gn.replace('.','x'),self.prot_conf))
                        temp_loop[loop_gn] = loop_template[key]
                    except:
                        try:
                            if loop_template[key]=='-':
                                raise Exception()
                            temp_loop[self.loop_label+'|'+str(l_res)] = loop_template[key]
                        except:
                            temp_loop[self.loop_label+'?'+str(l_res)] = '-'
                temp_loop[self.loop_label+'?'+str(l_res+1)] = 'x'
                if ECL2!=None:
                    temp_array[self.loop_label] = temp_loop
                else:
                    temp_array[self.loop_label+'_dis'] = temp_loop
            else:
                temp_array[seg_label] = gns
        return temp_array

    def insert_gaps_for_loops_to_arrays(self, main_pdb_array, reference_dict, template_dict, alignment_dict):
        ''' When there is no template for a loop region, this function inserts gaps for that region into the main
            template, fetches the reference residues and inserts these into the arrays. This allows for Modeller to
            freely model these loop regions.

            @param main_pdb_array: nested OrderedDict(), output of GPCRDBParsingPDB().pdb_array_creator().
            @param reference_dict: reference dictionary of AlignedReferenceTemplate.
            @param template_dict: template dictionary of AlignedReferenceTemplate.
            @param alignment_dict: alignment dictionary of AlignedReferenceTemplate.
        '''
        residues = Residue.objects.filter(protein_conformation__protein=self.reference_protein,
                                          protein_segment__slug=self.loop_label)
        temp_pdb_array = OrderedDict()
        for seg_id, seg in main_pdb_array.items():
            if self.segment_order.index(self.loop_label)-self.segment_order.index(seg_id[:4])==1:
                temp_loop = OrderedDict()
                count=0
                temp_pdb_array[seg_id] = seg
                for r in residues:
                    count+=1
                    temp_loop[self.loop_label+'?'+str(count)] = '-'
                temp_pdb_array[self.loop_label+'_free'] = temp_loop
                self.new_label = self.loop_label+'_free'
            else:
                temp_pdb_array[seg_id] = seg
        self.main_pdb_array = temp_pdb_array
        temp_ref_dict, temp_temp_dict, temp_aligned_dict = OrderedDict(), OrderedDict(), OrderedDict()
        for ref_seg, temp_seg, aligned_seg in zip(reference_dict, template_dict, alignment_dict):
            if ref_seg=='H8' and len(list(Residue.objects.filter(protein_conformation=self.prot_conf, protein_segment__slug='H8')))==0:
                continue
            if self.segment_order.index(self.loop_label)-self.segment_order.index(ref_seg[:4])==1:
                temp_ref_loop, temp_temp_loop, temp_aligned_loop = OrderedDict(), OrderedDict(), OrderedDict()
                temp_ref_dict[ref_seg] = reference_dict[ref_seg]
                temp_temp_dict[temp_seg] = template_dict[temp_seg]
                temp_aligned_dict[aligned_seg] = alignment_dict[aligned_seg]
                count=0
                for r in residues:
                    count+=1
                    temp_ref_loop[self.loop_label+'?'+str(count)] = r.amino_acid
                    temp_temp_loop[self.loop_label+'?'+str(count)] = '-'
                    temp_aligned_loop[self.loop_label+'?'+str(count)] = '.'
                temp_ref_dict[self.loop_label+'_free'] = temp_ref_loop
                temp_temp_dict[self.loop_label+'_free'] = temp_temp_loop
                temp_aligned_dict[self.loop_label+'_free'] = temp_aligned_loop
            else:
                temp_ref_dict[ref_seg] = reference_dict[ref_seg]
                temp_temp_dict[temp_seg] = template_dict[temp_seg]
                temp_aligned_dict[aligned_seg] = alignment_dict[aligned_seg]
        self.reference_dict = temp_ref_dict
        self.template_dict = temp_temp_dict
        self.alignment_dict = temp_aligned_dict


class Bulges(object):
    ''' Class to handle bulges in GPCR structures.
    '''
    def __init__(self, gn):
        self.gn = gn
        self.bulge_templates = []
        self.template = None

    def check_range(self, gn_list, protein_conformation, num):
        check = [dgn(i,protein_conformation) for i in gn_list]
        check_list = [i.sequence_number for i in list(Residue.objects.filter(protein_conformation=protein_conformation,
                                                                             display_generic_number__label__in=check))]
        ref_list = list(range(check_list[0],check_list[0]+num))
        if ref_list==check_list:
            return 1
        else:
            return 0

    def find_bulge_template(self, similarity_table, bulge_in_reference):
        ''' Searches for bulge template, returns residues of template (5 residues if the bulge is in the reference, 4
            residues if the bulge is in the template).

            @param gn: str, Generic number of bulge, e.g. 1x411 \n
            @param similarity_table: OrderedDict(), table of structures ordered by preference.
            Output of HomologyModeling().create_similarity_table(). \n
            @param bulge_in_reference: boolean, Set it to True if the bulge is in the reference, set it to False if the
            bulge is in the template.
        '''
        gn = self.gn
        parse = GPCRDBParsingPDB()
        for structure, value in similarity_table.items():
            anomaly_list = structure.protein_anomalies.all().prefetch_related()
            this_anomaly = ProteinAnomaly.objects.filter(generic_number__label=gn)
            if bulge_in_reference==True:
                try:
                    for anomaly in this_anomaly:
                        if anomaly in anomaly_list:
                            gn_list = [parse.gn_indecer(gn,'x',-2), parse.gn_indecer(gn,'x',-1), gn,
                                       parse.gn_indecer(gn,'x',+1), parse.gn_indecer(gn,'x',+2)]
                            if self.check_range(gn_list,structure.protein_conformation,5)==0:
                                raise Exception()
                            alt_bulge = parse.fetch_residues_from_pdb(structure, gn_list)
                            self.template = structure
                            return alt_bulge
                except:
                    pass
            elif bulge_in_reference==False:
                try:
                    suitable_temp = []
                    for anomaly in this_anomaly:
                        if anomaly not in anomaly_list:
                            pass
                        else:
                            suitable_temp.append('no')
                    if 'no' not in suitable_temp:
                        gn_list = [parse.gn_indecer(gn,'x',-2), parse.gn_indecer(gn,'x',-1),
                                   parse.gn_indecer(gn,'x',+1), parse.gn_indecer(gn,'x',+2)]
                        if self.check_range(gn_list,structure.protein_conformation,4)==0:
                            raise Exception()
                        alt_bulge = parse.fetch_residues_from_pdb(structure, gn_list)
                        self.template = structure
                        return alt_bulge
                except:
                    pass
        return None


class Constrictions(Bulges):
    ''' Class to handle constrictions in GPCRs.
    '''
    def __init__(self, gn):
        self.gn = gn
        self.constriction_templates = []
        self.template = None

    def find_constriction_template(self, similarity_table, constriction_in_reference):
        ''' Searches for constriction template, returns residues of template (4 residues if the constriction is in the
            reference, 5 residues if the constriction is in the template).

            @param gn: str, Generic number of constriction, e.g. 7x44 \n
            @param similarity_table: OrderedDict(), table of structures ordered by preference.
            Output of HomologyModeling().create_similarity_table(). \n
            @param constriction_in_reference: boolean, Set it to True if the constriction is in the reference, set it
            to False if the constriction is in the template.
        '''
        gn = self.gn
        parse = GPCRDBParsingPDB()
        for structure, value in similarity_table.items():
            this_anomaly = ProteinAnomaly.objects.filter(generic_number__label=gn)
            anomaly_list = structure.protein_anomalies.all().prefetch_related()
            if constriction_in_reference==True:
                try:
                    for anomaly in this_anomaly:
                        if anomaly in anomaly_list:
                            gn_list = [parse.gn_indecer(gn,'x',-2), parse.gn_indecer(gn,'x',-1),
                                       parse.gn_indecer(gn,'x',+1), parse.gn_indecer(gn,'x',+2)]
                            if self.check_range(gn_list,structure.protein_conformation,4)==0:
                                raise Exception()
                            alt_const = parse.fetch_residues_from_pdb(structure, gn_list)
                            self.template = structure
                            return alt_const
                except:
                    pass
            elif constriction_in_reference==False:
                try:
                    suitable_temp = []
                    for anomaly in this_anomaly:
                        if anomaly not in anomaly_list:
                            pass
                        else:
                            suitable_temp.append('no')
                    if 'no' not in suitable_temp:
                        gn_list = [parse.gn_indecer(gn,'x',-2), parse.gn_indecer(gn,'x',-1), gn,
                                   parse.gn_indecer(gn,'x',+1), parse.gn_indecer(gn,'x',+2)]
                        if self.check_range(gn_list,structure.protein_conformation,5)==0:
                            raise Exception()
                        alt_const = parse.fetch_residues_from_pdb(structure, gn_list)
                        self.template = structure
                        return alt_const
                except:
                    pass
        return None





class CreateStatistics(object):
    ''' Statistics dictionary for HomologyModeling.
    '''
    def __init__(self, reference):
        self.reference = reference
        self.info_dict = OrderedDict()

    def __repr__(self):
        return "<{} \n {} \n>".format(self.reference, self.info_dict)

    def items(self):
        ''' Returns the OrderedDict().items().
        '''
        return self.info_dict.items()

    def add_info(self, info_name, info):
        ''' Adds new information to the statistics dictionary.

            @param info_name: str, info name as dictionary key
            @param info: object, any object as value
        '''
        self.info_dict[info_name] = info
