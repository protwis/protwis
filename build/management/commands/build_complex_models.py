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


    def handle(self, *args, **options):
        tracker = {}
        all_models = django.apps.apps.get_models()[6:]
        test_model_updates(all_models, tracker, initialize=True)

        # Set verbosity level from integer enumeration based on command-line argument
        verbosity_level = ParserVerbosity.from_string_map.get(options['parser_verbosity'], ParserVerbosity.BASIC)

        if options['parser'] == "alphafoldcomplex":
            if not options['cleaned_seq_csv']:
                raise ValueError("The --cleaned_seq_csv argument is required for the 'alphafoldcomplex' parser.")
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

        try:
            self.logger.info('CREATING STRUCTURES')
            self.prepare_input(options['proc'], self.model_parser.model_dirs)
            test_model_updates(all_models, tracker, check=True)
            self.logger.info('COMPLETED CREATING STRUCTURES')
        except Exception as msg:
            self.logger.error(msg)

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

    def main_func(self, positions, iteration, count, lock):
        self.model_parser.process_models(write = True, low_memory = True, offset_start = positions[0], offset_end = positions[1])
            
