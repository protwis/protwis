from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from residue.models import Residue
from protein.models import Protein, ProteinGProteinPair
from ligand.models import *
from common.models import WebLink, WebResource, Publication
from django.db.models import Q, Count
import logging
from datetime import datetime
import time

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
                self.purge_data()
                print('Ended purging bias data')
            except Exception as msg:
                print(msg)
        self.calculate_selectivity()
        try:
            print('Calculating ligand/receptor selectivities')
            print(options['filename'])
            # self.calculate_selectivity()
            self.logger.info('COMPLETED CREATING BIAS')

        except Exception as msg:
            print('--error--', msg, '\n')
            self.logger.info("The error appeared in def handle")

    def purge_data(self):
        delete_bias_experiment = LigandReceptorStatistics.objects.all()
        delete_bias_experiment.delete()
        self.logger.info('Data is purged')


    def calculate_selectivity(self):
        # TODO: Discuss with David/Albert : how to chose selectivity/ compare among what?
        start = time.time()
        # get unique ligands
        assays = self.get_ligands()
        print(type(assays), len(assays))
        # iterate throgu assayexperiments using ligand ids
        ligands_with_data = self.get_data(assays)
        #process ligand assay queryset
        processed_ligand_assays = self.process_assays(ligands_with_data)
        self.save_data()

        end = time.time()
        print('---temp_increment time---', end - start)

    def get_ligands(self):
        #Getting ligands from the model
        try:
            content = AssayExperiment.objects.all().order_by(
                'ligand').distinct('ligand').only('ligand')
        except AssayExperiment.DoesNotExist:
            content = None
        return content

    def get_data(self, assays):
        ligand_list=list()
        #Getting data from the model for a ligand\n##limiting only by EC50 | IC50 (values)'
        for assay in assays:
            ligand_data = dict()
            try:
                content = AssayExperiment.objects.filter(ligand=assay.ligand
                                                         ).filter(Q(assay_type='F') | Q(assay_type='B')
                                                                  ).filter(Q(standard_type='IC50') |
                                                                           Q(standard_type='EC50') |
                                                                           Q(standard_type='Potency')).prefetch_related(
                    'ligand', 'protein'
                ).only('protein', 'ligand', 'standard_type', 'standard_value', 'assay_type', 'pchembl_value'
                       ).order_by('ligand')

                ligand_data['ligand'] = assay.ligand
                ligand_data['assays'] = content
            except AssayExperiment.DoesNotExist:
                ligand_data = None
            ligand_list.append(ligand_data)
        return ligand_list

    def process_assays(self, assays):
        for assay in assays:
            assay = self.process_every_queryset(assay)
            self.analyze_f_b(assay)

    def analyze_f_b(self, assay):
        assay_f=None
        assay_b=None
        if len(assay['assay_f'])>1:
            assay_f=assay['assay_f']
        if len(assay['assay_b'])>1:
            assay_b=assay['assay_b']
        if assay_b is not None:
            try:
                if float(assay_b[0]['pchembl_value'])-1 > float(assay_b[1]['pchembl_value']):
                    if assay_b[0]['protein'] is not assay_b[1]['protein']:
                        if assay_b[0]['protein'] in self.b_receptor_count:
                            self.b_receptor_count[assay_b[0]['protein']] = self.b_receptor_count[assay_b[0]['protein']]+1
                        else:
                            self.b_receptor_count[assay_b[0]['protein']] = 1
                        # import pdb; pdb.set_trace()
            except:
                pass
        if assay_f is not None:
            try:
                if float(assay_f[0]['pchembl_value'])-1 > float(assay_f[1]['pchembl_value']):
                    if assay_f[0]['protein'] is not assay_f[1]['protein']:
                        if assay_f[0]['protein'] in self.f_receptor_count:
                            self.f_receptor_count[assay_f[0]['protein']] = self.f_receptor_count[assay_f[0]['protein']]+1
                        else:
                            self.f_receptor_count[assay_f[0]['protein']] = 1
                        # import pdb; pdb.set_trace()
            except:
                pass

    def process_every_queryset(self, assays):
        processed_data = dict()
        ligand = assays['ligand']
        assay_b = assays['assays'].filter(standard_type='IC50')
        assay_f = assays['assays'].filter(Q(standard_type='Potency') |
                             Q(standard_type='EC50'))

        assay_b = self.sort_assay(self.process_querysets(assay_b))
        assay_f = self.sort_assay(self.process_querysets(assay_f))
        processed_data['ligand'] = ligand
        processed_data['assay_b'] = assay_b
        processed_data['assay_f'] = assay_f
        return processed_data

    def sort_assay(self, assays):
        return sorted(assays, key=lambda i: i['pchembl_value'], reverse=True)

    def process_querysets(self,querysets):
        processed_data = list()
        for i in querysets:
            try:
                assay_data = dict()
                assay_data["protein"] = i.protein
                assay_data["ligand"] = i.ligand
                assay_data["pchembl_value"] = i.pchembl_value
                assay_data["assay_type"] = i.assay_type
                assay_data["standard_type"] = i.standard_type
                assay_data["standard_value"] = i.standard_value
                assay_data['reference_protein'] = Protein()
                processed_data.append(assay_data)
            except:
                print('process data fail', i)
                continue

        return processed_data

    def sort_assay(self, assays):
        return sorted(assays, key=lambda i: i['pchembl_value'], reverse=True)

    def fetch_receptor_trunsducers(self, receptor):
        primary = set()
        temp = str()
        temp1 = str()
        secondary = set()
        try:
            gprotein = ProteinGProteinPair.objects.filter(protein=receptor)
            for x in gprotein:
                if x.transduction and x.transduction == 'primary':
                    primary.add(x.g_protein.name)
                elif x.transduction and x.transduction == 'secondary':
                    secondary.add(x.g_protein.name)
            for i in primary:
                temp += str(i.replace(' family', '')) + str(', ')

            for i in secondary:
                temp1 += str(i.replace('family', '')) + str(', ')
            return temp, temp1
        except:
            self.logger.info('receptor not found error')
            return None, None


    def save_data(self):
        #saving assay ---', final_assay
        for protein, counter in self.f_receptor_count.items():
            primary, secondary = self.fetch_receptor_trunsducers(
                protein)
            save_assay = LigandReceptorStatistics(
                protein=protein,
                type='f',
                value=counter,
                primary=primary,
                secondary=secondary)
            save_assay.save()

        for protein, counter in self.b_receptor_count.items():
            primary, secondary = self.fetch_receptor_trunsducers(
                protein)
            save_assay = LigandReceptorStatistics(
                protein=protein,
                type='b',
                value=counter,
                primary=primary,
                secondary=secondary)
            save_assay.save()

        print(self.f_receptor_count)
        print(self.b_receptor_count)
