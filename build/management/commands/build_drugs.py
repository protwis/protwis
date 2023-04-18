from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from common.models import WebResource, WebLink
from protein.models import Protein
from drugs.models import Drugs
from common.tools import test_model_updates
from optparse import make_option
import logging
import csv
import os
import pandas as pd
import django.apps

class Command(BaseCommand):
    help = 'Build Drug Data'

    def add_arguments(self, parser):
        parser.add_argument('--filename', action='append', dest='filename',
            help='Filename to import. Can be used multiple times')

    logger = logging.getLogger(__name__)

        # source file directory
    drugdata_data_dir = os.sep.join([settings.DATA_DIR, 'drug_data'])
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]

    def handle(self, *args, **options):
        if options['filename']:
            filenames = options['filename']
        else:
            filenames = False

        try:
            self.purge_drugs()
            test_model_updates(self.all_models, self.tracker, initialize=True)
            self.create_drug_data(filenames)
            test_model_updates(self.all_models, self.tracker, check=True)
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def purge_drugs(self):
        try:
            Drugs.objects.all().delete()
        except Drugs.DoesNotExist:
            self.logger.warning('Drugs mod not found: nothing to delete.')

    def create_drug_data(self, filenames=False):
        self.logger.info('CREATING DRUGDATA')

        # read source files
        if not filenames:
            filenames = [fn for fn in os.listdir(self.drugdata_data_dir) if fn.endswith('drug_data.csv')]

        for filename in filenames:

            filepath = os.sep.join([self.drugdata_data_dir, filename])

            data = pd.read_csv(filepath, low_memory=False, encoding = "ISO-8859-1")

            for index, row in enumerate(data.iterrows()):

                drugname = data[index:index+1]['Drug Name'].values[0]
                trialname = data[index:index+1]['Trial name'].values[0]
                drugalias_raw = data[index:index+1]['DrugAliases'].values[0]
                drugalias = ['' if str(drugalias_raw) == 'nan' else ', '+str(drugalias_raw)]
                # trialadd = ['' if str(trialname) == drugname else ' ('+str(trialname)+')']
                drugname = drugname + drugalias[0]

                entry_name = data[index:index+1]['EntryName'].values[0]

                phase = data[index:index+1]['Phase'].values[0]
                PhaseDate = data[index:index+1]['PhaseDate'].values[0]
                ClinicalStatus = data[index:index+1]['ClinicalStatus'].values[0]
                moa = data[index:index+1]['ModeOfAction'].values[0]
                targetlevel = data[index:index+1]['TargetCategory'].values[0]

                drugtype = data[index:index+1]['Drug Class'].values[0]
                indication = data[index:index+1]['Indication(s)'].values[0]
                novelty = data[index:index+1]['Target_novelty'].values[0]
                approval = data[index:index+1]['Approval'].values[0]
                status = data[index:index+1]['Status'].values[0]

                references = data[index:index+1]['PMID'].values[0]

                drug, created = Drugs.objects.get_or_create(name=drugname, synonym=', '.join(drugalias), drugtype=drugtype, indication=indication, novelty=novelty, approval=approval, phase=phase, phasedate=PhaseDate, clinicalstatus=ClinicalStatus, moa=moa, status=status, targetlevel=targetlevel,references=references)

                # fetch protein
                try:
                    p = Protein.objects.get(entry_name=entry_name)
                    drug.target.add(p)
                    drug.save()
                except Protein.DoesNotExist:

                    self.logger.error('Protein not found for entry_name {}'.format(entry_name))
                    continue


                # target_list = drug.target.all()

        self.logger.info('COMPLETED CREATING DRUGDATA')
