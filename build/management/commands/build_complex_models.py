from html import parser

from django.conf import settings
from django.utils.text import slugify
from django.db import IntegrityError
from django.db.models import Q
from django.core.exceptions import FieldError

from build.management.commands.base_build import Command as BaseBuild
from protein.models import (Protein, ProteinConformation, ProteinState, ProteinSegment)
from residue.models import Residue
from common.models import WebLink, WebResource, Publication
from common.tools import test_model_updates
from common.definitions import G_PROTEIN_DISPLAY_NAME as g_prot_dict, ARRESTIN_DISPLAY_NAME as arr_dict
from structure.models import Structure, StructureType, PdbData, Rotamer, Fragment, StructureExtraProteins, StructureModelScores, StructureModelpLDDT
from construct.functions import *
from structure.management.commands.generate_complexes_to_model import get_stimulatory_peptide_like_ligand_AssayExperiment_obj, get_inhibitory_peptide_like_ligand_AssayExperiment_obj

from contactnetwork.models import *
from contactnetwork.cube import compute_interactions

from Bio.PDB import PDBParser, PPBuilder, PDBIO
from Bio import pairwise2

from structure.functions import ParseAFComplexModels
from ligand.models import Ligand, LigandPeptideStructure
from interaction.models import *
from interaction.views import regexaa, check_residue, extract_fragment_rotamer
from signprot.models import SignprotComplex
import structure.assign_generic_numbers_gpcr as as_gn

from structure.model_parsers.boltz_two import BoltzTwoComplexModelParserConfig, BoltzTwoComplexModelParser
from structure.model_parsers.alphafold_complex import AlphaFoldTwoComplexModelParser, AlphaFoldTwoComplexModelParserConfig
from structure.model_parsers.logging import ParserVerbosity

import django.apps
import logging
import os
import sys
import yaml
import time
import gc
from collections import OrderedDict
from datetime import datetime, date
import json
from io import StringIO
from Bio.PDB.Selection import *
import re

# import traceback

_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG
def dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))

def represent_ordereddict(dumper, data):
    value = []

    for item_key, item_value in data.items():
        node_key = dumper.represent_data(item_key)
        node_value = dumper.represent_data(item_value)

        value.append((node_key, node_value))

    return yaml.nodes.MappingNode(u'tag:yaml.org,2002:map', value)

yaml.add_representer(OrderedDict, represent_ordereddict)
yaml.add_constructor(_mapping_tag, dict_constructor)

class Command(BaseBuild):
    help = 'Reads source data and creates pdb structure records'

    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
            type=int,
            action='store',
            dest='proc',
            default=1,
            help='Number of processes to run')
        parser.add_argument('-s', '--structure',
            dest='structure',
            help='Structure to import (PDB ID)',
            nargs='+')
        parser.add_argument('-u', '--purge',
            action='store_true',
            dest='purge',
            default=False,
            help='Purge existing records')
        parser.add_argument('-m', '--model_set_name',
            dest='model_set_name',
            required=True,
            help='The name of the model set to process (i.e folder name in structure_data folder)')
        parser.add_argument('-r', '--parser',
            dest='parser',
            choices=['boltztwocomplex', 'alphafoldcomplex'],
            required=True,
            help='The parser to use for the model set (BoltzTwoComplex or AlphafoldComplex)')
        parser.add_argument('--cleaned_seq_csv',
            action='store',
            default=False,
            help='Load cleaned sequences from CSV (required for AlphafoldComplex parser)'),
        parser.add_argument('-y', '--parser_verbosity',
            choices=['silent', 'basic', 'everything'],
            default="basic",
            help='Set the verbosity level for the parser'),
        parser.add_argument('-e', '--error_handling',
            choices=["log", "raise", "log_then_raise", "log_with_trace", "log_with_trace_then_raise"],
            default="log_then_raise",
            help='Set the error handling strategy for the parser. "raise" will raise exceptions, "log" will log errors but suppress them, "log_then_raise" will log and then raise exceptions, "log_with_trace" will log errors with stack trace but suppress them, and "log_with_trace_then_raise" will log errors with stack trace and then raise exceptions.')

        self.tracker = {}
        self.all_models = django.apps.apps.get_models()[6:]
        test_model_updates(self.all_models, self.tracker, initialize=True)

        ### USE below to fix seg ends
        xtal_seg_end_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'mod_xtal_segends.yaml'])
        with open(xtal_seg_end_file, 'r') as f:
            self.xtal_seg_ends = yaml.load(f, Loader=yaml.Loader)

        xtal_anomalies_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'all_anomalies.yaml'])
        with open(xtal_anomalies_file, 'r') as f2:
            self.xtal_anomalies = yaml.load(f2, Loader=yaml.Loader)

        xtal_representatives_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'xtal_representatives.yaml'])
        with open(xtal_representatives_file, 'r') as f3:
            self.xtal_representatives = yaml.load(f3, Loader=yaml.Loader)

        non_xtal_seg_end_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'non_xtal_segends.yaml'])
        with open(non_xtal_seg_end_file, 'r') as f:
            self.non_xtal_seg_ends = yaml.load(f, Loader=yaml.Loader)

        self.s = ProteinSegment.objects.all()
        self.segments = {}
        for segment in self.s:
            self.segments[segment.slug] = segment

        self.parsed_pdb = None

        self.construct_errors, self.rotamer_errors, self.contactnetwork_errors, self.interaction_errors = [],[],[],[]

        with open(os.sep.join([settings.DATA_DIR, 'residue_data', 'unnatural_amino_acids.yaml']), 'r') as f_yaml:
            raw_uaa = yaml.safe_load(f_yaml)
            self.unnatural_amino_acids = {}
            for i, j in raw_uaa.items():
                self.unnatural_amino_acids[i] = j

    def handle(self, *args, **options):
        # Set verbosity level from integer enumeration based on command-line argument
        verbosity_level = ParserVerbosity.from_string_map.get(options['parser_verbosity'], ParserVerbosity.BASIC)

        # delete any existing structure data
        if options['purge']:
            try:
                self.purge_structures()
                self.tracker = {}
                test_model_updates(self.all_models, self.tracker, initialize=True)
            except Exception as msg:
                print(msg)
                self.logger.error(msg)

        #self.run_contactnetwork = not options['skip_cn']

        if options['parser'] == "alphafoldcomplex":
            config = AlphaFoldTwoComplexModelParserConfig(model_set_name=options['model_set_name'],
                                                        cleaned_seq_csv=options['cleaned_seq_csv'],
                                                        model_receptor_state="Active",
                                                        pdb_preferred_chain="A",
                                                        error_handling=options['error_handling'],
                                                        verbosity=verbosity_level)

            self.model_parser = AlphaFoldTwoComplexModelParser(config)
            self.model_parser.get_model_directories()
        elif options['parser'] == "boltztwocomplex":
            config = BoltzTwoComplexModelParserConfig(model_set_name=options['model_set_name'],
                                                        model_receptor_state="Active",
                                                        pdb_header_override={'deposition_date': '2026-03-01', 'release_date': '2026-03-01'},
                                                        default_model_version= {'default_version_number': '1', 'override': {'drd1_human-"zuclopenthixol"[5311507]': '2'}},
                                                        pdb_preferred_chain="A",
                                                        error_handling=options['error_handling'],
                                                        verbosity=verbosity_level)            
            self.model_parser = BoltzTwoComplexModelParser(config)
            self.model_parser.get_model_directories()

        if options['structure']:
            filtered_set = []
            for model_dir in self.model_parser.model_dirs:
                if model_dir in options['structure']:
                    filtered_set.append(model_dir)
            self.model_parser.model_dirs = filtered_set
            # self.parsed_structures.complexes = [i for i in self.parsed_structures.complexes if i in options['structure'] or i.lower() in options['structure']]

        # self.incremental_mode = options['incremental']        

        try:
            self.logger.info('CREATING STRUCTURES')
            self.prepare_input(options['proc'], self.model_parser.model_dirs)
            test_model_updates(self.all_models, self.tracker, check=True)
            self.logger.info('COMPLETED CREATING STRUCTURES')
        except Exception as msg:
            self.logger.error(msg)


    def purge_structures(self):
        
        alphafold_models_filter1 =  Q(structure__structure_type__slug__startswith='af-signprot') | \
                                   Q(structure__structure_type__slug__startswith='af-peptide') | \
                                   Q(structure__structure_type__slug__startswith='af-arrestin')
        
        alphafold_models_filter2 =  Q(structure_ligand_pair__structure__structure_type__slug__startswith='af-signprot') | \
                                   Q(structure_ligand_pair__structure__structure_type__slug__startswith='af-peptide') | \
                                   Q(structure_ligand_pair__structure__structure_type__slug__startswith='af-arrestin')               

        models = Structure.objects.filter(alphafold_models_filter1)

        for m in models:
            PdbData.objects.filter(pdb=m.pdb_data.pdb).delete()
            WebLink.objects.filter(index=m.pdb_code.index).delete()
        models.delete()
        
        rfi = ResidueFragmentInteraction.objects.filter(alphafold_models_filter2)
        rfi.delete()
        # ResidueFragmentInteractionType.objects.all().delete()
        sli = StructureLigandInteraction.objects.filter(alphafold_models_filter1)
        sli.delete()
        #Remove previous Rotamers/Residues to prepare repopulate
        f = Fragment.objects.filter(alphafold_models_filter1)
        f.delete()
        r = Rotamer.objects.filter(alphafold_models_filter1)
        r.delete()
        # PdbData.objects.all().delete()


    @staticmethod
    def get_peptide_ligand_effect_data():
        stimulatory_peptides = get_stimulatory_peptide_like_ligand_AssayExperiment_obj()
        inhibitory_peptides = get_inhibitory_peptide_like_ligand_AssayExperiment_obj()

        peptide_effects_dict = {}
        for effect,peptides in zip(['stimulatory', 'inhibitory'],[stimulatory_peptides, inhibitory_peptides]):
            peptide_effects_dict[effect] = {}
            for ep in peptides:
                if ep.ligand.sequence not in peptide_effects_dict[effect]:
                    peptide_effects_dict[effect][ep.ligand.sequence] = []
                try:
                    peptide_effects_dict[effect][ep.ligand.sequence].append(ep.ligand.gpcrdbid)
                except AttributeError:
                    peptide_effects_dict[effect][ep.ligand.sequence].append(ep.ligand.id)
        return peptide_effects_dict
    
    # @staticmethod
    # def parsecalculation(sd, data, molecule, ignore_ligand_preset=False):
    #     module_dir = '/tmp/interactions/'
    #     pdb_id = sd['pdb']
    #     pdb_name = sd['location'].split('/')[-1]
    #     complex_name = sd['location'].split('/')[-1].split('-r')[0]
    #     gpcrdb_id = molecule['gpcrdb id']
    #     pdb_location = module_dir + 'pdbs/' + complex_name + '/' + pdb_name
    #     web_resource = WebResource.objects.get(slug='pdb')
    #     web_link, _ = WebLink.objects.get_or_create(index=pdb_id)
    #     structure = Structure.objects.filter(pdb_code=web_link)
    #     if structure.exists():
    #         structure = Structure.objects.get(pdb_code=web_link)

    #         if structure.pdb_data is None:
    #             if os.path.isfile(pdb_location):
    #                 pdbdata, created = PdbData.objects.get_or_create(pdb=open(pdb_location, 'r').read())  # does this close the file?
    #             else:
    #                 print('quitting due to no pdb in filesystem')
    #                 quit()
    #             structure.pdb_data = pdbdata
    #             structure.save()

    #         protein = structure.protein_conformation
    #         lig_key = list(data.keys())[0]
    #         #/tmp/interactions/pdbs/oprd_mouse/oprd_mouse-1643-rank0.pdb
    #         prot_pep = sd['location'].split('/')[-1].split('-r')[0]
    #         # /tmp/interactions/results/ranked_0/interaction
    #         f = module_dir + "results/" + prot_pep + "/interaction/" + prot_pep + "_" + lig_key + ".pdb"
    #         print(f)

    #         if os.path.isfile(f):
    #             pdbdata, created = PdbData.objects.get_or_create(pdb=open(f, 'r').read())  # does this close the file?
    #         else:
    #             print('quitting due to no pdb for fragment in filesystem', f)
    #             quit()

    #         struct_lig_interactions = StructureLigandInteraction.objects.filter(pdb_reference=lig_key, ligand_id=gpcrdb_id, structure=structure, annotated=True) #, pdb_file=None
    #         if struct_lig_interactions.exists():  # if the annotated exists
    #             try:
    #                 struct_lig_interactions = struct_lig_interactions.get()
    #                 struct_lig_interactions.pdb_file = pdbdata
    #                 ligand = struct_lig_interactions.ligand
    #             except Exception as msg:
    #                 print('error with duplication structureligand',lig_key,msg)
    #                 quit() #not sure about this quit
    #         elif StructureLigandInteraction.objects.filter(pdb_reference=lig_key, structure=structure).exists():
    #             try:
    #                 struct_lig_interactions = StructureLigandInteraction.objects.filter(pdb_reference=lig_key, structure=structure).get()
    #                 struct_lig_interactions.pdb_file = pdbdata
    #             except StructureLigandInteraction.DoesNotExist: #already there
    #                 struct_lig_interactions = StructureLigandInteraction.objects.filter(pdb_reference=lig_key, structure=structure, pdb_file=pdbdata).get()
    #             ligand = struct_lig_interactions.ligand
    #         else:  # create ligand and pair
    #             print(pdb_id, "Skipping interactions with ", pdb_id)
    #             pass

    #         struct_lig_interactions.save()

    #         ResidueFragmentInteraction.objects.filter(structure_ligand_pair=struct_lig_interactions).delete()

    #         for interaction in data[lig_key]['interactions']:
    #             aa = interaction[0]
    #             if aa[-1] != structure.preferred_chain:
    #                 continue
    #             aa, pos, _ = regexaa(aa)
    #             residue = check_residue(protein, pos, aa)
    #             f = interaction[1]
    #             fragment, rotamer = extract_fragment_rotamer(f, residue, structure, ligand)
    #             if fragment is not None:
    #                 interaction_type, created = ResidueFragmentInteractionType.objects.get_or_create(
    #                                             slug=interaction[2],
    #                                             name=interaction[3],
    #                                             type=interaction[4], direction=interaction[5])
    #                 fragment_interaction, created = ResidueFragmentInteraction.objects.get_or_create(
    #                                                 structure_ligand_pair=struct_lig_interactions,
    #                                                 interaction_type=interaction_type,
    #                                                 fragment=fragment, rotamer=rotamer)
    #     else:
    #         print('Something went wrong and we passed')
    #         pass

    def main_func(self, positions, iteration, count, lock):
        self.model_parser.process_models(write = True, low_memory = True, offset_start = positions[0], offset_end = positions[1])
            
