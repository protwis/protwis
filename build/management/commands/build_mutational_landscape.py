from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db import IntegrityError

from protein.models import (Protein, ProteinGProtein,ProteinGProteinPair, ProteinConformation, ProteinState, ProteinFamily, ProteinAlias,
        ProteinSequenceType, Species, Gene, ProteinSource, ProteinSegment)
from residue.models import (ResidueNumberingScheme, ResidueGenericNumber, Residue, ResidueGenericNumberEquivalent)
from mutational_landscape.models import NaturalMutations, CancerMutations, DiseaseMutations

import pandas as pd
import numpy as np
import math, os
import logging

class Command(BaseCommand):
    help = 'Build Mutational Landscape'

    # source file directory
    mutation_data_path = os.sep.join([settings.DATA_DIR, 'mutational_landscape'])

    logger = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument('--filename', action='append', dest='filename',
            help='Filename to import. Can be used multiple times')

    def handle(self, *args, **options):
        if options['filename']:
            filenames = options['filename']
        else:
            filenames = False

        try:
            self.purge_data()
            self.create_natural_mutations()
            # self.create_cancer_mutations()
            # self.create_disease_mutations()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def purge_data(self):
        try:
            NaturalMutations.objects.all().delete()
        except:
            self.logger.warning('Existing data cannot be deleted')

    def create_natural_mutations(self, filenames=False):
        self.logger.info('CREATING NATURAL MUTATIONS')

        # read source files
        if not filenames:
            filenames = [fn for fn in os.listdir(self.mutation_data_path) if fn.endswith('.csv')]

        for filename in filenames:
            filepath = os.sep.join([self.mutation_data_path, filename])

            snp_data =  pd.read_csv(filepath, low_memory=False)

            for index, entry in enumerate(snp_data.iterrows()):

                entry_name = snp_data[index:index+1]['EntryName'].values[0]
                sequence_number = snp_data[index:index+1]['Sequence Number'].values[0]
                residue = snp_data[index:index+1]['GPCRdb'].values[0]
                amino_acid = snp_data[index:index+1]['NMaa'].values[0]
                allele_frequency = float(snp_data[index:index+1]['Allele Frequency'].values[0])
                allele_count = int(snp_data[index:index+1]['Allele Count'].values[0])
                allele_number = int(snp_data[index:index+1]['Allele Number'].values[0])
                number_homozygotes = int(snp_data[index:index+1]['Number of Homozygotes'].values[0])

                try:
                    p = Protein.objects.get(entry_name=entry_name)
                except Protein.DoesNotExist:
                    self.logger.warning('Protein not found for entry_name {}'.format(entry_name))
                    continue

                try:
                    res=Residue.objects.get(protein_conformation__protein=p, sequence_number=sequence_number)
                except:
                    # self.logger.warning('No residue number (GAP - position) for', CGN, "in ", p.name, "")
                    continue

                if res:
                    try:
                        barcode, created = NaturalMutations.objects.get_or_create(protein=p, residue=res, amino_acid=amino_acid, allele_frequency=allele_frequency, allele_count=allele_count, allele_number=allele_number, number_homozygotes=number_homozygotes)
                        # if created:
                            # self.logger.info('Created SNP for ' + sequence_number + ' for protein ' + p.name)
                    except:
                        print(entry_name, sequence_number, allele_frequency, allele_count, allele_number, number_homozygotes)
                        # self.logger.error('Failed creating SNP for ' + sequence_number + ' for protein ' + p.name)

        self.logger.info('COMPLETED CREATING NATURAL MUTATIONS')

    def create_cancer_mutations(self, filenames=False):
        self.logger.info('CREATING CANCER MUTATIONS')

        # read source files
        if not filenames:
            filenames = [fn for fn in os.listdir(self.mutation_data_path) if fn.endswith('.csv')]

        for filename in filenames:
            filepath = os.sep.join([self.mutation_data_path, filename])

            cancer_data =  pd.read_csv(filepath, low_memory=False)

        self.logger.info('COMPLETED CREATING CANCER MUTATIONS')

    def create_disease_mutations(self, filenames=False):
        self.logger.info('CREATING DISEASE MUTATIONS')

        # read source files
        if not filenames:
            filenames = [fn for fn in os.listdir(self.mutation_data_path) if fn.endswith('.csv')]

        for filename in filenames:
            filepath = os.sep.join([self.mutation_data_path, filename])

            disease_data =  pd.read_csv(filepath, low_memory=False)

        self.logger.info('COMPLETED CREATING DISEASE MUTATIONS')


