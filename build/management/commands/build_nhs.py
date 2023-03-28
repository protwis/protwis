from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db import IntegrityError

from drugs.models import Drugs
from mutational_landscape.models import NHSPrescribings

import pandas as pd
import numpy as np
import math, os
import logging
import re
from decimal import *
getcontext().prec = 20

class Command(BaseCommand):
    help = 'Build NHS data'

    # source file directory
    drug_data_path =  os.sep.join([settings.DATA_DIR, 'drug_data'])

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
            self.create_NHS()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def purge_data(self):
        try:
            NHSPrescribings.objects.all().delete()
        except Exception as msg:
            print(msg)
            self.logger.warning('Existing data cannot be deleted')

    def create_NHS(self, filenames=False):
        print('Creating NHS')
        self.logger.info('CREATING NHS PRESCRIBINGS')

        #read source files
        if not filenames:
            filenames = [fn for fn in os.listdir(self.drug_data_path) if fn.endswith('nhs.csv')]

        for filename in filenames:
            filepath = os.sep.join([self.drug_data_path, filename])
            nhs_data = pd.read_csv(filepath, low_memory=False)

            for _, entry in nhs_data.iterrows():

                date = entry['date']
                quantity = entry['quantity']
                items = entry['items']
                actual_cost = entry['actual_cost']
                drugCode = entry['drugCode']
                op_name = entry['drugName']
                bnf_section_raw = entry['section']
                bnf_section_name = bnf_section_raw.split(': ')[1]
                bnf_section_id = bnf_section_raw.split(': ')[0]
                drugNameQuery = entry['drugNameQuery']

                try:
                    drugname = Drugs.objects.filter(name=drugNameQuery)[0]

                except:
                    self.logger.warning('Drug not found for chemical {}'.format(drugNameQuery))
                    continue

                nhs, created = NHSPrescribings.objects.get_or_create(date=date, quantity=quantity, items=items, actual_cost=actual_cost, drugCode=drugCode, op_name=op_name, bnf_section=bnf_section_name, drugname=drugname)

        self.logger.info('COMPLETED CREATING NHS PRESCRIBINGS')
