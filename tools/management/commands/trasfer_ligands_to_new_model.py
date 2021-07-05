from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from protein.models import Protein, ProteinGProteinPair
from ligand.models import *
from common.models import WebLink, Publication
import pandas as pd

import logging
import os
from decimal import Decimal
import gzip, json
import datetime

class Command(BaseCommand):
    mylog = logging.getLogger(__name__)
    mylog.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    file_handler = logging.FileHandler('Selectivity.log')
    file_handler.setLevel(logging.ERROR)
    file_handler.setFormatter(formatter)
    mylog.addHandler(file_handler)

    help = 'Calcualtes ligand/receptor selectivity'
    # source file directory
    # structure_data_dir = os.sep.join([settings.EXCEL_DATA, 'ligand_data', 'bias'])

    publication_cache = {}
    ligand_cache = int()
    f_receptor_count = dict()
    b_receptor_count = dict()

    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
                            type=int,
                            action='store',
                            dest='proc',
                            default=1,
                            help='Number of processes to run')
        parser.add_argument('-f', '--filename',
                            action='append',
                            dest='filename',
                            help='Filename to import. Can be used multiple times')
        parser.add_argument('-u', '--purge',
                            action='store_true',
                            dest='purge',
                            default=False,
                            help='Purge existing bias records')
        parser.add_argument('--test_run', action='store_true', help='Skip this during a test run',
                            default=False)

    def handle(self, *args, **options):
        if options['test_run']:
            print('Skipping in test run')
            return
        if options['purge']:
            try:
                print('Started purging bias data')
                self.purge_bias_data()
                print('Ended purging bias data')
            except Exception as msg:
                print(msg)
        self.calculate_selectivity()
        try:
            print('transfer all ligands')

        except Exception as msg:
            print('--error--', msg, '\n')

    def purge_bias_data(self):
        delete_bias_experiment = Ligand_v2.objects.all()
        delete_bias_experiment.delete()



    def calculate_selectivity(self):
        # get unique ligands
        assays = self.get_ligands()
        print('#Step 1 - Done')
        # iterate throgu assayexperiments using ligand ids
        ligands_with_data = self.get_data(assays)
        print('#Step 2 - Done')
        #process ligand assay queryset
        processed_ligand_assays = self.process_assays(ligands_with_data)
        print('#Step 3 - Done')
        self.prepare_to_save(processed_ligand_assays)
        print('#Step 4 - Done')

    def get_ligands(self):
        #Getting ligands from the model
        try:
            content = Ligand.objects.all().order_by(
                'id')
        except Ligand.DoesNotExist:
            content = None
        return content

    def get_data(self, ligands):
        ligand_list=list()
        #Getting data from the model for a ligand\n##limiting only by EC50 | IC50 (values)'
        for ligand in ligands:
            ligand_data = dict()
            # try:
            ligand_data['name'] = ligand.name
            ligand_data['old_ligand_id'] = ligand.id
            ligand_data['pdbe'] = ligand.pdbe
            ligand_data['smiles'] = ligand.properities.smiles
            ligand_data['inchikey'] = ligand.properities.inchikey
            ligand_data['sequence'] = ligand.properities.sequence
            ligand_data['mw'] = ligand.properities.mw
            ligand_data['rotatable_bonds'] = ligand.properities.rotatable_bonds
            ligand_data['hacc'] = ligand.properities.hacc
            ligand_data['hdon'] = ligand.properities.hdon
            ligand_data['logp'] = ligand.properities.logp
            ligand_data['type'] = ligand.properities.ligand_type
            ligand_data['web_links'] = ligand.properities.web_links.all()
            ligand_data['vendors'] = ligand.properities.vendors.all()
            # except:
            #     ligand_data = None
            ligand_list.append(ligand_data)
        return ligand_list

    def process_assays(self, assays):
        context = dict()
        print('*****len before merge******', len(assays))
        for ligand in assays:
            name = None
            if ligand['smiles']:
                name = str(ligand['smiles'])
            if name is None and ligand['inchikey']:
                name = str(ligand['inchikey'])
            if name is None and ligand['sequence']:
                name = str(ligand['sequence'])
            if name is None:
                self.ligand_cache+=1
            if name is not None:
                if name in context:
                    context[name].append(ligand)
                else:
                    context[name] = list()
                    context[name].append(ligand)
        print('*****len after merge******', len(context))
        print('*****len self.ligand_cache******', self.ligand_cache)
        return context

    def prepare_to_save(self, ligands):
        for i in ligands.items():
            for ligand_data in i[1]:
                try:
                    type = ligand_data['type'].name
                except:
                    type = 5
                ligand_v2 = Ligand_v2(
                    pubchem_id = ligand_data['name'],
                    default_name = ligand_data['name'],
                    smiles = ligand_data['smiles'],
                    inchikey = ligand_data['inchikey'],
                    pdb = ligand_data['pdbe'],
                    old_ligand_id = ligand_data['old_ligand_id'],
                    sequence = ligand_data['sequence'],
                    ligand_type = type
                    )
                ligand_v2.save()
                ligand_v2_physchem = Ligand_v2_Physchem(
                    ligand = ligand_v2,
                    mw = ligand_data['mw'],
                    rotatable_bonds = ligand_data['rotatable_bonds'],
                    hacc = ligand_data['hacc'],
                    hdon = ligand_data['hdon'],
                    logp = ligand_data['logp'],
                )
                ligand_v2_physchem.save()
                for link in ligand_data['web_links']:
                    synonyms = Synonyms(
                        ligand = ligand_v2,
                        name = ligand_data['name'],
                        resource = link.web_resource.name,
                        link = link,
                        specialty = 'ligand',
                    )
                    synonyms.save()
                for link in ligand_data['vendors']:

                    synonyms = Synonyms(
                        ligand = ligand_v2,
                        name = link.sid,
                        resource = link.vendor.name,
                        link = link.url,
                        specialty = 'vendor',
                    )
                    synonyms.save()

    def sort_assay(self, assays):
        return sorted(assays, key=lambda i: i['pchembl_value'], reverse=True)
