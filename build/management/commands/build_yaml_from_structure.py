from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from protein.models import Protein
from residue.models import Residue
from common.models import WebLink, WebResource, Publication
from structure.models import Structure, StructureType, StructureStabilizingAgent
from ligand.models import Ligand, LigandRole

from optparse import make_option
from datetime import datetime
import logging, os


class Command(BaseCommand):
    
    def add_arguments(self, parser):
        parser.add_argument('filename', nargs='+', type=str)

    logger = logging.getLogger(__name__)
    structure_build_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data'])

    csv_fields = {
        'id' : 0, 
        'prot_name' : 1, 
        'class' : 2, 
        'pdb_code' : 3, 
        'endogenous_ligand' :4, 
        'resolution' : 5, 
        'xray_ligand' : 6, 
        'ligand_role' : 7, 
        'chain' : 8, 
        'pubmed_id' : 9, 
        'date' : 10, 
        'G_protein' : 11, 
        'stabilizing_agent' : 12, 
        'n-term' : 13, 
        'icl1' : 14, 
        'ecl1' : 15, 
        'icl2' : 16, 
        'ecl2.1' : 17, 
        'ecl2.2' : 18, 
        'icl3' : 19, 
        'ecl3m' : 20, 
        'c-term' : 21,
        }

    help = 'Creates structures from .csv formatted data strored in the directory {}'.format(structure_build_data_dir)

    def handle(self, *args, **options):
        # delete any existing protein data
        try:
            self.truncate_structure_tables()
        except Exception as msg:
            # print(msg)
            self.logger.error(msg)
        self.create_yamls(options['filename'])