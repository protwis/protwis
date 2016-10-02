from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from structure.sequence_parser import SequenceParser

import Bio.PDB.Polypeptide as polypeptide
from Bio import SeqIO, pairwise2
from Bio.PDB import PDBParser, PPBuilder, parse_pdb_header


import logging, json, os

class Command(BaseCommand):

    help = "Function evaluating sequence compatibility with wildtype protein."

    logger = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument('pdb_file')

    def handle(self, *args, **options):
        root, ext = os.path.splitext(os.path.basename(options['pdb_file']))
        print("Working on file {}".format(options['pdb_file']))
        header = parse_pdb_header(options['pdb_file'])
        sp = SequenceParser(options['pdb_file'])
        print(sp.get_fusions())
        print(sp.get_mutations())
        print(sp.get_deletions())
        json_data = {}
        json_data["header"] = header
        json_data.update(sp.get_fusions())
        json_data.update(sp.get_mutations())
        json_data.update(sp.get_deletions())
        json.dump(json_data, open(os.sep.join([settings.DATA_DIR, "{}_auto.json".format(root)]), 'w'), indent=4, separators=(',', ': '))
        #json.dump(json_data, open("test.json", 'w'), indent=4, separators=(',', ': '))