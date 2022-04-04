from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.utils.text import slugify
from django.db import IntegrityError

from protein.models import Protein, ProteinConformation
from residue.models import Residue
from structure.models import Structure
from structure.sequence_parser import SequenceParser
from construct.models import (Construct,Crystallization,CrystallizationLigandConc,ChemicalType,Chemical,ChemicalConc,ChemicalList,
CrystallizationMethods,CrystallizationTypes,ChemicalListName,ContributorInfo,ConstructMutation,ConstructInsertion,ConstructInsertionType,
ConstructDeletion,ConstructModification,CrystalInfo,ExpressionSystem,Solubilization,PurificationStep,Purification)

from Bio.PDB import parse_pdb_header

from collections import OrderedDict
import logging, os, json, datetime

class Command(BaseCommand):

    help = 'Parse construct structures and prepare auto_* .yaml files.'

    def add_arguments(self, parser):
        parser.add_argument('--filename', action='append', dest='filename',
            help='Structure to import (pdb format). Can be used multiple times')

    logger = logging.getLogger(__name__)

        # source file directory
    construct_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data','construct_data'])
    pdb_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'pdbs'])

    def handle(self, *args, **options):
        if options['filename']:
            filenames = options['filename']
        else:
            filenames = False

        try:
            self.analyze_construct(filenames)
        except Exception as msg:
            print(msg)
            self.logger.error(msg)


    def analyze_construct(self, filenames=None):
        self.logger.info("ANALYZING CONSTRUCT STRUCTURES")

        # read source files
        if not filenames:
            filenames = os.listdir(self.construct_data_dir)

        for filename in filenames:
            if filename[-3:]!='pdb' and filename[-3:]!='ent':
                continue
            root, ext = os.path.splitext(os.path.basename(filename))
            print(filename)
            print(root)
            filepath = os.sep.join([self.construct_data_dir, filename])
            self.logger.info("Working on a file: {}".format(filename))
            header = parse_pdb_header(filepath)
            parser = SequenceParser(filepath)


            json_data = OrderedDict()
            json_data["header"] = header
            json_data.update(parser.get_fusions())
            json_data.update(parser.get_mutations())
            json_data.update(parser.get_deletions())
            json.dump(json_data, open(os.sep.join([settings.DATA_DIR, "{}_auto.json".format(root)]), 'w'), indent=4, separators=(',', ': '))
