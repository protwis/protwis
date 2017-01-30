from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from structure.functions import check_gn
from structure.assign_generic_numbers_gpcr import GenericNumbering
from structure.structural_superposition import FragmentSuperpose
#from structure.sequence_parser import SequenceParser

from Bio.PDB import *
from io import StringIO, BytesIO
import zipfile
import logging

class Command(BaseCommand):

    help = "Function evaluating sequence compatibility with wildtype protein."

    logger = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument('pdb_file')

    def handle(self, *args, **options):
        print("Working on file {}".format(options['pdb_file']))
        
        frag_sp = FragmentSuperpose(options['pdb_file'])
        superposed_fragments = frag_sp.superpose_fragments(use_similar=True)
        if superposed_fragments == []:
            print("No fragments were aligned.")
        else:
            print("{} fragments were aligned.".format(len(superposed_fragments)))
            io = PDBIO()
            zipf = zipfile.ZipFile(options['pdb_file'].replace('.pdb', '.zip'), 'a', zipfile.ZIP_DEFLATED)
            for fragment, pdb_data in superposed_fragments:
                io.set_structure(pdb_data)
                tmp = StringIO()
                io.save(tmp)
                zipf.writestr("all_fragments/{!s}".format(fragment.generate_filename()), tmp.getvalue())
            zipf.close()