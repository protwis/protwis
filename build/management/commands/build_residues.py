from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from protein.models import Protein
from residue.models import Residue

import logging


class Command(BaseCommand):
    help = 'Creates residues (without generic numbers) from protein records'

    logger = logging.getLogger(__name__)


    def handle(self, *args, **options):

        # delete any existing residue data
        try:
            self.truncate_residue_tables()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

        # create residue records for all proteins
        try:
            self.create_residues()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def truncate_residue_tables(self):
        cursor = connection.cursor()
        
        tables_to_truncate = [
            'residue_residue',
            'residue_residueset',
            'residue_residuenumber'
        ]

        for table in tables_to_truncate:
            cursor.execute("TRUNCATE TABLE " + table + " CASCADE")

    def create_residues(self):
        self.logger.info('CREATING RESIDUES')
        
        proteins = Protein.objects.all()

        for protein in proteins:
            for i, aa in enumerate(protein.sequence):
                r = Residue()
                r.protein = protein
                r.sequence_number = i
                r.amino_acid = aa

                try:
                    r.save()
                    self.logger.info('Created residue ' + aa + str(i) + ' for protein ' + protein.name)
                except Exception as msg:
                    print(msg)
                    self.logger.error('Failed creating residue ' + aa + str(i) + ' for protein ' + protein.name)

        self.logger.info('COMPLETED CREATING RESIDUES')