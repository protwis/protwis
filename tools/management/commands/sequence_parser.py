from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from structure.sequence_parser import SequenceParser

import Bio.PDB.Polypeptide as polypeptide
from Bio import SeqIO, pairwise2
from Bio.PDB import PDBParser, PPBuilder, parse_pdb_header


import logging

class Command(BaseCommand):

    help = "Function evaluating sequence compatibility with wildtype protein."

    logger = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument('pdb_file')

    def handle(self, *args, **options):
        print("Working on file {}".format(options['pdb_file']))
        header = parse_pdb_header(options['pdb_file'])
        print(header['compound'])
        sp = SequenceParser(options['pdb_file'])
        c = list(sp.mapping.keys())[0]
        poly = sp.get_chain_peptides(c)
        for peptide in poly:
            print("Start: {} Stop: {} Len: {}".format(peptide[0].id[1], peptide[-1].id[1], len(peptide)))
            sp.map_to_wt_blast(c, peptide, None, int(peptide[0].id[1]))
        sp.map_seqres()
        sp.save_excel_report("test.xlsx")
        #sp.get_report()