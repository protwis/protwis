import math
import logging
import time
from build.management.commands.base_build import Command as BaseBuild
from django.db import models
from protein.models import Protein
from ligand.models import *

MISSING_PROTEINS = {}
SKIPPED = 0

class Command(BaseBuild):
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
    ligand_cache = {}
    data_all = []

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

        parser.add_argument('--test_run', action='store_true', help='Skip this during a test run',
                            default=False)

    def handle(self, *args, **options):
        if options['test_run']:
            print('Skipping in test run')
            return
        try:
            print('Calculating ligand/receptor selectivities')
            print(options['filename'])
            self.calculate_selectivity()
            self.logger.info('COMPLETED CREATING BIAS')

        except Exception as msg:
            print('--error--', msg, '\n')
            self.logger.info("The error appeared in def handle")

# pylint: disable=E0602
    def calculate_selectivity(self):
        start = time.time()
        # get unique ligands
        ligands = self.get_ligands()
        # iterate throgu assayexperiments using ligand ids
        for ligand in ligands:
            # for every ligand compare assays and find most potent receptor
            assay_raw = self.get_data(ligand.ligand)
            # pass assays to compare potencies for different receptor__species
            # list can be empty, has only 1 element or multiple
            try:
                assay_f = assay_raw.filter(Q(standard_type='Potency') |
                                                Q(standard_type='EC50'))
                assay_b = assay_raw.filter(standard_type='IC50')
            except assay_raw.DoesNotExist:
                print('no data for ligand', ligand)
            self.process_data(assay_f, "F")
            try:
                self.process_data(assay_f, "F")
            except:
                print("f selectivity error")
                continue
            try:
                self.process_data(assay_b, "B")
            except:
                print("b selectivity error")
                continue

        end = time.time()
        print('---temp_increment time---', end - start)

    def process_data(self, assay, type_d):
        assay_list = list()
        assay_list = self.process_assays(assay)
        try:
            sorted_assay_list = self.sort_assay(assay_list)
        except:
            print('sorting error')
        self.analyze_assay(sorted_assay_list, type_d)

# pylint: disable=R0201
    def get_ligands(self):
        #Getting ligands from the model
        try:
            content = AssayExperiment.objects.all().order_by(
                'ligand').distinct('ligand').only('ligand')
        except AssayExperiment.DoesNotExist:
            content = None
        return content

# pylint: disable=R0201
    def get_data(self, ligand_name):
        #Getting data from the model for a ligand\n##limiting only by EC50 | IC50 (values)'
        try:
            content = AssayExperiment.objects.filter(ligand=ligand_name
                                                     ).filter(Q(assay_type='F') | Q(assay_type='B')
                                                              ).filter(Q(standard_type='IC50') |
                                                                       Q(standard_type='EC50') |
                                                                       Q(standard_type='Potency')).prefetch_related(
                                                                           'ligand', 'protein'
                                                                           ).only('protein', 'ligand', 'standard_type', 'standard_value', 'assay_type'
                                                                           ).order_by('ligand')
        except AssayExperiment.DoesNotExist:
            content = None
        return content

# pylint: disable=R0201
    def process_assays(self, assays):
        processed_data = list()
        for i in assays:
            try:
                assay_data = dict()
                assay_data["protein"] = i.protein
                assay_data["ligand"] = i.ligand
                assay_data["assay_type"] = i.assay_type
                assay_data["standard_type"] = i.standard_type
                assay_data["standard_value"] = float(i.standard_value)
                assay_data['reference_protein'] = Protein()
                processed_data.append(assay_data)
            except:
                print('process data', i)
                continue

        return processed_data

# pylint: disable=R0201
    def sort_assay(self, assays):
        return sorted(assays, key=lambda i: i['standard_value'], reverse=True)

# pylint: disable=R0201
    def analyze_assay(self, assays, type_d):
        # select most potent if more than 10 folds

        try:
            most_potent = min(assays, key=lambda x: x['standard_value'])
            # most_potent = min(item['standard_value'] for item in assays if item != None)
        except:
            most_potent = {}
            return most_potent
        # print(most_potent)
        if most_potent and most_potent != None:
            for i in assays:
                try:
                    if (most_potent['standard_value']*10) < (i['standard_value']):
                        if most_potent['protein'] != i['protein']:
                            most_potent['reference_protein'] = i['protein']
                            mp_value = self.fetch_measurements(most_potent['standard_value'], most_potent['standard_type'], 'nm')
                            i_value = self.fetch_measurements(i['standard_value'], i['standard_type'], 'nm')
                            most_potent['value'] = round(i_value - mp_value, 3)
                            # print('\n stabdard value and protein of I',most_potent,'\nidata:', i )
                            self.save_data(most_potent, type_d)
                except:
                    pass

        return most_potent
        # for assay in assays:
        #     if most_potent['pchembl_value']

    def save_data(self, final_assay, type_d):
        #saving assay ---', final_assay
        if self.check_dublicate(final_assay) == False:
            save_assay = LigandReceptorStatistics(
                ligand=final_assay['ligand'],
                protein=final_assay['protein'],
                type=type_d,
                value=final_assay['value'],
                reference_protein=final_assay['reference_protein']
            )
            save_assay.save()

# pylint: disable=R0201
    def check_dublicate(self, final_assay):
        try:
            experiment = LigandReceptorStatistics.objects.filter(
                protein=final_assay['protein'],
                ligand=final_assay['ligand'],
                reference_protein=final_assay['reference_protein'],
                value=final_assay['value'],
                type=final_assay['assay_type']
                )
            experiment = experiment.get()
            return True
        except LigandReceptorStatistics.DoesNotExist:
            return False

    def fetch_measurements(self, potency, p_type, unit):
        if p_type.lower() == 'pec50':
            potency = 10**(potency*(-1))
            # pp = (-1)*log(potency)
            p_type = 'EC50'
        elif p_type.lower() == 'logec50':
            potency = 10**(potency)
            p_type = 'EC50'
        elif p_type.lower() == 'pic50':
            potency = 10**(potency*(-1))
            p_type = 'IC50'
        elif p_type.lower() == 'logic50':
            potency = 10**(potency)
            p_type = 'IC50'

        if p_type.lower() == 'ec50':
            if unit.lower() == 'nm':
                potency = potency* 10**(-9)
            elif unit.lower() == 'µm':
                potency = potency* 10**(-9)
            elif unit.lower() == 'pm':
                potency = potency* 10**(-12)
            elif unit.lower() == 'mm':
                potency = potency* 10**(-6)
            else:
                pass
        if p_type.lower() == 'ic50':
            if unit.lower() == 'nm':
                potency = potency* 10**(-9)
            elif unit.lower() == 'µm':
                potency = potency* 10**(-9)
            elif unit.lower() == 'pm':
                potency = potency* 10**(-12)
            elif unit.lower() == 'mm':
                potency = potency* 10**(-6)
            else:
                pass
        if potency:
            potency = math.log10(potency)
        return potency
