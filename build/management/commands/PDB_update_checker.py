from django.core.management.base import BaseCommand, CommandError
from django.core.exceptions import MultipleObjectsReturned
from django.conf import settings
from django.db import connection
from django.utils.text import slugify
from django.db import IntegrityError
from structure.models import Structure
import time
import pandas as pd
import numpy as np
import requests
from bs4 import BeautifulSoup
from requests.exceptions import ConnectionError
import logging
import os
from urllib.request import urlopen



class Command(BaseCommand):
    help = 'Check latest PDBs versions'
    old_pdbs = []
    logger = logging.getLogger(__name__)
    pdb_codes = list(Structure.objects.all().values_list('pdb_code__index', flat=True).distinct())
    to_skip = ['6ORV', '6Z4Q', '6Z4S', '6Z4V', '6Z66', '6YVR', '6Z8N', '6ZA8', '6ZIN', '7B6W', '7PP1', '7F1T', '7XBX', '7F1R', '7F1Q', '7TMW'] #don't download these
    URL = 'http://rcsb.org/versions/{}'
    pdb_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'pdbs'])

    def handle(self, *args, **options):
        # try:
            self.update_pdbs()
        # except Exception as msg:
        #     print(msg)
        #     self.logger.error(msg)

    @staticmethod
    def get_soup(URL, id):
        r = requests.get(URL.format(id))
        soup = BeautifulSoup(r.text, "html.parser")
        return soup

    @staticmethod
    def date_conversion(date):
        months = {'JAN': '01', 'APR': '04', 'JUL': '07', 'OCT': '10',
                  'FEB': '02', 'MAY': '05', 'AUG': '08', 'NOV': '11',
                  'MAR': '03', 'JUN': '06', 'SEP': '09', 'DEC': '12'}
        day = date.split('-')[0]
        mon = date.split('-')[1]
        year = date.split('-')[2]
        output = '-'.join(['20'+year, months[mon], day])
        return output

    def scrapping_data(self, counter, to_download, pdb_codes):
        for code in pdb_codes[counter:]:
            try:
                counter += 1
                revs = []
                if code.startswith('AFM'):
                    continue

                code = code.lower()
                soup = self.get_soup(self.URL, code)
                table = soup.findAll('table', { 'class' : 'table-bordered' })
                release_table = table[0].findAll('tr')

                if len(release_table) == 2:
                    continue

                latest_release = release_table[-1].findAll('td')
                latest_date = latest_release[1].text
                db_data = Structure.objects.filter(pdb_code__index=code.upper())[0]
                pdb_data = db_data.pdb_data.pdb.split('\n')

                for line in pdb_data:
                    if line.startswith('REVDAT'):
                        revs.append(line.split()[2])

                if len(revs) == 0:
                    converted = db_data.publication_date.strftime('%Y-%m-%d')
                else:
                    converted = self.date_conversion(revs[0])

                if latest_date != converted:
                    to_download.append(code.upper())

            except ConnectionError:
                print('Connection Error detected, now sleeping and restarting...')
                time.sleep(20)
                scrapping_data(counter, to_download, pdb_codes)
        return to_download

    def update_pdbs(self):
        keyword = 'REVDAT'
        counter = 0
        start = time.time()
        old_pdbs = self.scrapping_data(counter, self.old_pdbs, self.pdb_codes)
        end = round(time.time() - start, 1)
        print('Calculation time: {}'.format(end))
        print('Total PDBs to be updated: {}'.format(len(old_pdbs)))
        for pdb in old_pdbs:
            if pdb in self.to_skip:
                self.logger.info('There is a new release of PDB {}'.format(pdb))
                print('There is a new release of PDB {}'.format(pdb))
                continue
            pdb_path = os.sep.join([self.pdb_data_dir, pdb + '.pdb'])
            self.logger.info('Fetching PDB file {}'.format(pdb))
            print('Fetching PDB file {}'.format(pdb))
            url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdb
            pdbdata_raw = urlopen(url).read().decode('utf-8')
            with open(pdb_path, 'w') as f:
                f.write(pdbdata_raw)
