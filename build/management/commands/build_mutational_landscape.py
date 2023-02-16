from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db import IntegrityError

from protein.models import (Protein, ProteinConformation, ProteinState, ProteinFamily, ProteinAlias,
        ProteinSequenceType, Species, Gene, ProteinSource, ProteinSegment)
from residue.models import (ResidueNumberingScheme, ResidueGenericNumber, Residue, ResidueGenericNumberEquivalent)
from mutational_landscape.models import NaturalMutations, PTMs
from common.tools import test_model_updates
import django.apps
import pandas as pd
import numpy as np
import math, os
import logging
import re
from decimal import *

getcontext().prec = 20

class Command(BaseCommand):
    help = 'Build Mutational Landscape'

    # source file directory
    mutation_data_path = os.sep.join([settings.DATA_DIR, 'mutational_landscape'])
    proteins_not_found = []

    logger = logging.getLogger(__name__)
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]

    def add_arguments(self, parser):
        parser.add_argument('--filename', action='append', dest='filename',
            help='Filename to import. Can be used multiple times')

    def handle(self, *args, **options):
        if options['filename']:
            filenames = options['filename']
        else:
            filenames = False

        self.proteins_not_found = []

        try:
            self.purge_data()
            test_model_updates(self.all_models, self.tracker, initialize=True)
            self.create_PTMs()
            self.create_natural_mutations()
            #here purge is not an argument, so we initialize the test after the purge
            test_model_updates(self.all_models, self.tracker, check=True)
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def purge_data(self):
        try:
            PTMs.objects.all().delete()
            NaturalMutations.objects.all().delete()
        except Exception as msg:
            print(msg)
            self.logger.warning('Existing data cannot be deleted')

    def create_natural_mutations(self, filenames=False):
        self.logger.info('CREATING NATURAL MUTATIONS')

        # read source files
        if not filenames:
            filenames = [fn for fn in os.listdir(self.mutation_data_path) if fn.endswith('exac.csv')]

        for filename in filenames:
            filepath = os.sep.join([self.mutation_data_path, filename])

            snp_data = pd.read_csv(filepath, low_memory=False)

            for index, entry in enumerate(snp_data.iterrows()):

                entry_name = snp_data[index:index + 1]['EntryName'].values[0]
                if entry_name in self.proteins_not_found:
                    continue
                sequence_number = snp_data[index:index + 1]['SequenceNumber'].values[0]
                allele_frequency = float(snp_data[index:index + 1]['af'].values[0]) # af/Allele Frequency
                allele_count = int(snp_data[index:index + 1]['ac'].values[0]) # ac/Allele Count
                allele_number = int(snp_data[index:index + 1]['an'].values[0]) # an/Allele Number
                number_homozygotes = int(snp_data[index:index + 1]['ac_hom'].values[0]) # ac_hm/ Number of Homozygotes
                type = snp_data[index:index + 1]['Type'].values[0]

                if type != 'missense':
                    prot_con = snp_data[index:index + 1]['Protein Consequence'].values[0]
                    splitterm = re.findall(r'\d+', prot_con)[0]
                    amino_acid = prot_con.split(splitterm)[1]
                    sift_score = None
                    polyphen_score = None
                else:
                    amino_acid = snp_data[index:index + 1]['NMaa'].values[0]
                    if 'SigProts' in filename:
                        sift_score = None
                        polyphen_score = None
                    else:
                        sift_score = float(snp_data[index:index + 1]['sift_score'].values[0])
                        polyphen_score = float(snp_data[index:index + 1]['polyphen_score'].values[0])


                try:
                    p = Protein.objects.get(entry_name=entry_name)
                except Protein.DoesNotExist:
                    self.proteins_not_found.append(entry_name)
                    self.logger.warning('Protein not found for entry_name {}'.format(entry_name))
                    continue

                try:
                    res=Residue.objects.get(protein_conformation__protein=p, sequence_number=sequence_number)
                except:
                    # self.logger.warning('No residue number (GAP - position) for', sequence_number, "in ", p.name, "")
                    continue
                if res:
                    # try:
                    snp, created = NaturalMutations.objects.get_or_create(protein=p, residue=res, amino_acid=amino_acid, allele_frequency=allele_frequency, allele_count=allele_count, allele_number=allele_number, number_homozygotes=number_homozygotes,
                    sift_score=sift_score, type=type, polyphen_score=polyphen_score) #
                        # if created:
                            # self.logger.info('Created SNP for ' + str(sequence_number) + ' for protein ' + str(p.name))
                    # except:
                        # print(entry_name, sequence_number, allele_frequency, allele_count, allele_number, number_homozygotes, type)
                        # self.logger.error('Failed creating SNP for ' + sequence_number + ' for protein ' + p.name)

        self.logger.info('COMPLETED CREATING NATURAL MUTATIONS')

    def create_PTMs(self, filenames=False):
        self.logger.info('CREATING PTM SITES')

        # read source files
        if not filenames:
            filenames = [fn for fn in os.listdir(self.mutation_data_path) if fn.endswith('ptms.csv')]

        for filename in filenames:
            filepath = os.sep.join([self.mutation_data_path, filename])

            ptm_data =  pd.read_csv(filepath, low_memory=False)

            for index, entry in enumerate(ptm_data.iterrows()):

                entry_name = ptm_data[index:index+1]['EntryName'].values[0]
                if entry_name in self.proteins_not_found:
                    continue
                sequence_number = ptm_data[index:index+1]['SequenceNumber'].values[0]
                modification = ptm_data[index:index+1]['Type'].values[0]
                # source = ptm_data[index:index+1]['Source'].values[0]

                try:
                    p = Protein.objects.get(entry_name=entry_name)
                except Protein.DoesNotExist:
                    self.proteins_not_found.append(entry_name)
                    self.logger.warning('Protein not found for entry_name {}'.format(entry_name))
                    continue

                try:
                    res=Residue.objects.get(protein_conformation__protein=p, sequence_number=sequence_number)
                except:
                    continue
                if res:
                    # g = PTMsType.objects.get_or_create(modification=modification)
                    snp, created = PTMs.objects.get_or_create(protein=p, residue=res, modification=modification) #

        self.logger.info('COMPLETED CREATING PTM SITES')
