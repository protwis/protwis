from django.core.management.base import BaseCommand, CommandError
from django.db import connection
from django.db import IntegrityError

from ligand.models import Ligand, LigandFingerprint, LigandMol

import os
import django.apps
import logging

# CHECK POSTGRESQL LOG FILES FOR BUILDING ERRORS!!!
# e.g WARNING:  could not create molecule from SMILES '[N-]=[N]=N[C@H]1C(C(C)(C)C)C23C45C1OC(=O)[C@@]5(OC2OC(=O)[C@@H]3O)[C@@]1([C@H](C4O)OC(=O)C1C)O'

class Command(BaseCommand):
    help = 'Build cheminformatics molecule data and molecule fingerprints for database search. '

    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--verbose', default=False, action='store_true', help='Print progress in stdout.')

    logger = logging.getLogger(__name__)


    def handle(self, *args, **options):    
        with connection.cursor() as cursor:
            if options['verbose']: print('Truncating LigandMol...')
            LigandMol.custom_objects.truncate_table()
            self.logger.info('Truncated LigandMol.')
            if options['verbose']: print('Truncating LigandFingerprint...')
            LigandFingerprint.custom_objects.truncate_table()
            self.logger.info('Truncated LigandFingerprint.')
            ligand_column1_dbname = LigandMol._meta.get_field('ligand').get_attname_column()[1]
            mol_column_dbname = LigandMol._meta.get_field('molecule').get_attname_column()[1]
            ligand_pk_dbname = Ligand._meta.get_field(Ligand._meta.pk.name).get_attname_column()[1]
            ligand_smiles_column_dbname =  Ligand._meta.get_field('smiles').get_attname_column()[1]
            txt = 'Building cheminformatics molecule data...'
            self.logger.info(txt)
            if options['verbose']: print(txt)
            cursor.execute('INSERT INTO "%s" ("%s", "%s") SELECT * FROM (SELECT "%s", mol_from_smiles("%s"::cstring) AS mol FROM "%s") AS tmp WHERE mol IS NOT NULL;' % 
                (LigandMol.objects.model._meta.db_table,ligand_column1_dbname,mol_column_dbname,
                ligand_pk_dbname,ligand_smiles_column_dbname,Ligand.objects.model._meta.db_table))
            txt = 'COMPLETED building cheminformatics molecule data'
            self.logger.info(txt)
            if options['verbose']: print(txt)
            ligand_column2_dbname = LigandFingerprint._meta.get_field('ligand').get_attname_column()[1]
            mfp2_field_dbname = LigandFingerprint._meta.get_field('mfp2').get_attname_column()[1]
            txt = 'Building molecule fingerprints...'
            self.logger.info(txt)
            if options['verbose']: print(txt)
            cursor.execute('INSERT INTO "%s" ("%s","%s") SELECT * FROM (SELECT "%s",morganbv_fp("%s") AS mol FROM "%s") AS tmp;' % 
                (LigandFingerprint.objects.model._meta.db_table,ligand_column2_dbname,mfp2_field_dbname,
                ligand_column1_dbname,mol_column_dbname,LigandMol.objects.model._meta.db_table))
            txt = 'COMPLETED building molecule fingerprints'
            self.logger.info(txt)
            if options['verbose']: print(txt)



 