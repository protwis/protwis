from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db import IntegrityError
from django.utils.text import slugify
from django.http import HttpResponse, JsonResponse

from build.management.commands.base_build import Command as BaseBuild
from common.tools import fetch_from_cache, save_to_cache, fetch_from_web_api
from residue.models import Residue
from protein.models import Protein, ProteinGProteinPair
from ligand.models import *
from mutation.models import Mutation
from ligand.functions import get_or_make_ligand
from common.models import WebLink, WebResource, Publication

import logging
import os
from datetime import datetime
import xlrd
import operator
import traceback
import time
import math
import json


MISSING_PROTEINS = {}
SKIPPED = 0


class Command(BaseBuild):
    mylog = logging.getLogger(__name__)
    mylog.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    file_handler = logging.FileHandler('biasDataTest.log')
    file_handler.setLevel(logging.ERROR)
    file_handler.setFormatter(formatter)
    mylog.addHandler(file_handler)

    help = 'Reads bias data and imports it'
    # source file directory
    # structure_data_dir = os.sep.join([settings.EXCEL_DATA, 'ligand_data', 'bias'])
    structure_data_dir = '/protwis/sites/protwis/excel/'
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
        # delete any existing structure data
        if options['purge']:
            try:
                self.purge_bias_data()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
        # import the structure data
        self.bias_list()
        try:
            print('CREATING BIAS DATA')
            print(options['filename'])
            # self.bias_list()
            self.logger.info('COMPLETED CREATING BIAS')

        except Exception as msg:
            print('--error--', msg, '\n')
            self.logger.info("The error appeared in def handle")

    def fetch_experiment(self, publication, ligand, receptor, residue, mutation, source):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        try:
            experiment = AnalyzedExperiment.objects.filter(
                publication=publication, ligand=ligand, receptor=receptor, source=source)
            experiment = experiment.get()
            return True
        except Exception as msg:
            experiment = None
            self.mylog.exception(
                "Protein AnalyzedExperiment error | module: AnalyzedExperiment.")
            return False

    def process_data(self, content):
        '''
        Merge BiasedExperiment with its children
        and pass it back to loop through dublicates
        '''
        rd = []
        counter = 0
        for instance in enumerate(content):
            temp_obj = []
            fin_obj = {}
            fin_obj['main'] = (instance[1])
            rd.append(fin_obj)

        print('---counter of assays at process data---', counter)
        return rd

    def fetch_protein(self, protein_from_excel):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        test = None
        if Protein.objects.get(id=protein_from_excel):
            protein = Protein.objects.get(id=protein_from_excel)

        return protein

    def fetch_ligand(self, ligand_id):
        """
        fetch ligands with Ligand model
        requires: ligand id
        """
        l = None
        if str(ligand_id) in self.ligand_cache:
            if ligand_id in self.ligand_cache[str(ligand_id)]:
                l = self.ligand_cache[str(ligand_id)][ligand_id]
        else:
            self.ligand_cache[str(ligand_id)] = {}

        if not l:
            l = Ligand.objects.get(id = ligand_id)
        return l

    #revise
    def fetch_ligand_cache(self, ligand_cache):
        """
        fetch ligands with Ligand model
        requires: ligand id
        """
        ligand_list = list();
        for ligand in ligand_cache:
            ligand_list.append(ligand)

        # l = None
        # if str(ligand_id) in self.ligand_cache:
        #     if ligand_id in self.ligand_cache[str(ligand_id)]:
        #         l = self.ligand_cache[str(ligand_id)][ligand_id]
        # else:
        #     self.ligand_cache[str(ligand_id)] = {}
        #
        # if not l:
        #     l = Ligand.objects.get(id = ligand_id)
        return ligand_list

    def fetch_vendor(self, ligand):

        links = LigandVendorLink.objects.filter(lp=ligand.properities.id)
        vendor_count = 0
        for x in links:

            vendor_count = vendor_count + 1

        return vendor_count

    def define_protein(self, description):
        # if description exists
        signalling_protein = None
        if description:
            #to lowercase

            description = description.lower()

        # define keywords
            G12 = 0
            Gio = 0
            Gq = 0
            GS = 0
            Barr = 0
            # separate description
            for word in description.split():
                if word.lower == 'gi' or word.lower == 'gio':
                    Gio = Gio+100
                if word.lower == 'g12' or word == 'G12':
                    G12 = G12+100
                if word.lower == 'gq' or word =='Gq':
                    Gq = Gq+100
                if word.lower == 'gs' or word == 'GS':
                    GS = GS+100
                if word.lower == 'arrestin' or word == 'barrestin':
                    Barr = Barr+100

                if word == 'camp':
                    for key in description.split():
                        if key == 'inhibition' or key == 'decrease' or key == 'forskolin-stimulated':
                            Gio = Gio+10
                            continue
                        elif (key == 'accumulation' or key == 'formation' or
                              key == 'release' or key == 'production' or
                              key == 'generation' or key == 'creb'):
                            GS = GS+10

                elif (word == 'ca2+' or word == 'calcium or ip-1' or
                      word == 'inositol' or word == 'phosphate'):
                    for key in description.split():
                        if (key == 'accumulation' or key == 'formation' or
                              key == 'release' or key == 'flux' or
                              key == 'production' or key == 'mobilization' or key == 'generation'):
                            Gq = Gq+10
                elif (word == 'arrestin'):
                    Barr = Barr+100
            signal_prot = list()
            signal_prot.append(G12, Gio ,Gq ,GS ,Barr)
            signalling_protein = str(max(signal_prot))

        return signalling_protein

    def fetch_receptor_trunsducers(self, receptor):
        primary = ""
        secondary = ""
        gprotein = ProteinGProteinPair.objects.filter(protein=receptor)
        for x in gprotein:
            if x.transduction and x.transduction == 'primary':
                primary += str(x.g_protein.name)
            elif x.transduction and x.transduction == 'secondary':
                secondary += str(x.g_protein.name)
        return primary, secondary

    def change(self, rd, ligand_list):
        '''
        Merge bias experminet data with assay data
        Define family in accordance with subfamily
        '''
        test = list()
        send = dict()
        increment = 0
        counter = 0

        for j in rd:
            temp_dict = dict()
            temp = dict()
            doubles = []

            temp['publication'] = None
            temp['ligand'] = self.fetch_ligand( j['main'].ligand_id)
            if temp['ligand'].id in ligand_list:
                print('\n\n---yes---')
                continue
            temp['vendors'] = self.fetch_vendor(temp['ligand'])
            temp['receptor'] = self.fetch_protein(j['main'].protein_id)
            temp['residue'] = None
            temp['mutation'] = None
            if temp['receptor']:
                temp['species'] = temp['receptor'].species.common_name
            else:
                temp['species'] = None
            temp_dict['signalling_protein'] = self.define_protein(j['main'].assay_description)
            if temp_dict['signalling_protein'] == None:
                continue
            temp_dict['family'] = temp_dict['signalling_protein']
            temp_dict['assay'] = j['main'].assay_id
            temp_dict['assay_description'] = j['main'].assay_description
            temp_dict['quantitive_measure_type'] = j['main'].standard_type
            temp_dict['quantitive_activity'] = float(j['main'].standard_value)
            temp_dict['quantitive_unit'] = j['main'].standard_units
            temp_dict['assay_type'] = j['main'].assay_type
            temp_dict['potency'] = None

            if not isinstance(temp_dict['quantitive_activity'], (int, float)):
                temp_dict['quantitive_activity'] = None
            else:
                if temp_dict['quantitive_unit'].lower() == 'nm':
                    temp_dict['quantitive_activity'] = temp_dict['quantitive_activity']* 10**(-9)
                elif temp_dict['quantitive_unit'].lower() == 'µm':
                    temp_dict['quantitive_activity'] = temp_dict['quantitive_activity']* 10**(-9)
                elif temp_dict['quantitive_unit'].lower() == 'pm':
                    temp_dict['quantitive_activity'] = temp_dict['quantitive_activity']* 10**(-12)
                elif temp_dict['quantitive_unit'].lower() == 'mm':
                    temp_dict['quantitive_activity'] = temp_dict['quantitive_activity']* 10**(-6)
                else:
                    pass
            doubles.append(temp_dict)
            temp['assay'] = doubles
            send[increment] = temp
            increment = increment + 1
            counter += 1

        print('---counter of assays at change---', counter)
        return send

    def process_group(self, send):
        '''
        Recieve data from "process_data"
        search for objects with same publication ligand receptor
        Create new object of biased_experiment and all children
        '''
        doubles = list()
        result = list()
        context = dict()
        counter = 0
        recep_residue = dict()

        for j in send.items():

            name = str(str(j[1]['ligand']) + '/' + str(j[1]['receptor']))
            temp_obj = list()

            if(name in context):
                temp_obj = context[name]['assay']

            for i in j[1]['assay']:
                if i not in temp_obj:
                    temp_obj.append(i)

            context[name] = j[1]
            context[name]['assay'] = temp_obj
        keys = [k for k, v in context.items() if len(v['assay']) <= 1]
        for x in keys:
            del context[x]

        self.process_calculation(context)
        for i in context.items():
            temp_obj = dict()
            name = str(i[1]['ligand'])
            if(name in recep_residue):
                test1=0
                test2=0
                # compare order_no 0 for receptor and remain lowest ec50
                for qa in recep_residue[name][1]['biasdata']:
                    if(qa['order_no'] == 0):
                        test1 = qa['quantitive_activity']
                for qa in i[1]['biasdata']:
                    if(qa['order_no'] == 0):
                        test2 = qa['quantitive_activity']

                if(test1>test2):
                    recep_residue[name] = i[1]
            temp_obj = i
            recep_residue[name] = temp_obj
        return recep_residue

    def process_calculation(self, context):
        countq = 0
        counter = 0
        counter1 = 0
        for i in context.items():
            test = dict()
            temp_obj = list()
            # checking for dublicates
            for j in i[1]['assay']:
                if j not in temp_obj:
                    temp_obj.append(j)
                else:
                    print('passing dublicate___-')
            i[1]['assay'] = temp_obj
            test = sorted(i[1]['assay'], key=lambda k: k['quantitive_activity']
                          if k['quantitive_activity'] else 999999,  reverse=False)

            for x in enumerate(test):
                x[1]['order_no'] = x[0]
            i[1]['biasdata'] = test
            i[1].pop('assay')
            self.calc_potency(i[1]['biasdata'])

    def calc_potency(self, biasdata):
        count = 0
        most_potent = dict()
        for i in biasdata:
            count+=1
            if i['order_no'] == 0:
                most_potent = i
                i['potency'] = None
        for i in biasdata:
            try:
                if i['order_no'] > 0:
                    if i['quantitive_measure_type'].lower() == 'ec50':
                        if i['quantitive_activity'] is not None and most_potent['quantitive_activity'] is not None:
                            i['potency'] = round(
                                i['quantitive_activity'] / most_potent['quantitive_activity'], 1)
                    elif i['quantitive_measure_type'] == 'IC50':
                        # TODO: round
                        i['potency'] = i['quantitive_activity'] - most_potent['quantitive_activity']
                    else:
                        i['potency'] = None
            except:
                i['potency'] = None

        # print('---counter of calc_potency---', count)

    def limit_family(self,send):
        result_list = list()
        gs_list = list()
        gs = dict()
        gio_list = list()
        gio = dict()
        barr_list = list()
        barr = dict()
        gq_list = list()
        gq = dict()
        g12_list = list()
        g12 = dict()
        for i in send:
            ## TODO: add standard deviation
            if i['signalling_protein'] == 'GS':
                if i['signalling_protein'] not in gs_list:
                    gs = i
                    print('\n---gs1---', gs)
                    gs_list.append(i['signalling_protein'])
                else:
                    if gs['potency'] is not None and i['potency'] < gs['potency']:
                        gs = i
                        print('\n---gs2---', gs)
                    else:
                        pass
            if i['signalling_protein'] == 'Gio':
                if i['signalling_protein'] not in gio_list:
                    gio = i
                    gio_list.append(i['signalling_protein'])
                else:
                    if i['potency'] < gi['potency']:
                        gio = i
                    else:
                        pass
            if i['signalling_protein'] == 'Barr':
                if i['signalling_protein'] not in barr_list:
                    barr = i
                    barr_list.append(i['signalling_protein'])
                else:
                    if i['potency'] < barr['potency']:
                        barr = i
                    else:
                        pass
            if i['signalling_protein'] == 'Gq':
                if i['signalling_protein'] not in gq_list:
                    gq = i
                    gq_list.append(i['signalling_protein'])
                else:
                    if gq['potency'] is not None and i['potency'] < gq['potency']:
                        gq = i
                    else:
                        pass
            if i['signalling_protein'] == 'G12':
                if i['signalling_protein'] not in g12_list:
                    g12 = i
                    g12_list.append(i['signalling_protein'])
                else:
                    if g12['potency'] is not None and i['potency'] < g12['potency']:
                        g12 = i
                    else:
                        pass
        if len(gs) > 0:
            result_list.append(gs)
        if len(gio) > 0:
            result_list.append(gio)
        if len(barr) > 0:
            result_list.append(barr)
        if len(gq) > 0:
            result_list.append(gq)
        if len(g12) > 0:
            result_list.append(g12)
        print('\n---result_list---', result_list)
        return result_list

    def fetch_chembl(self,ligand):
        temp = ligand
        chembl_id = None
        links = temp.properities.web_links.all()
        # print('\n----link id---', links)
        for x in links:
            if x.web_resource.slug=='chembl_ligand':
                chembl_id = [x for x in links if x.web_resource.slug=='chembl_ligand'][0].index
        print('\n----chembl id---', chembl_id)
        return chembl_id

    def bias_list(self):
        print('i am in')
        context = {}
        ligand_cache = AnalyzedExperiment.objects.values_list('ligand', flat=True).filter(source='different_publication')
        ligand_list = self.fetch_ligand_cache(ligand_cache)
        content = AssayExperiment.objects.filter(assay_type = 'F').exclude(standard_type = 'Ki').exclude(standard_type = 'AC50').exclude(standard_type = 'Kd').prefetch_related(
            'assay', 'ligand', 'protein'
            ).order_by('protein').order_by('ligand')

        # merge children
        pre_data = self.process_data(content)
        # # transform to combined dictionary
        combined = self.change(pre_data, ligand_list)
        context.update({'data': self.process_group(combined)})
        for i in context['data'].items():
            try:
                # i[1].pop('reference')
                # i[1].pop('biasdata')
                    # TODO: move by one tab when uncomment try catch
                source = 'chembl3'
                assays = dict()

                if self.fetch_experiment(i[1][1]['publication'], i[1][1]['ligand'], i[1][1]['receptor'], i[1][1]['residue'], i[1][1]['mutation'], source) == False:
                    primary, secondary = self.fetch_receptor_trunsducers(i[1][1]['receptor'])
                    # i[1][1]['biasdata'] = self.limit_family(i[1][1]['biasdata'])
                    chembl = None
                    chembl = self.fetch_chembl(i[1][1]['ligand'])
                    experiment_entry = AnalyzedExperiment(publication=i[1][1]['publication'],
                                                          ligand=i[1][1]['ligand'],
                                                          receptor=i[1][1]['receptor'],
                                                          mutation=i[1][1]['mutation'],
                                                          residue=i[1][1]['residue'],
                                                          chembl=chembl,
                                                          vendor_quantity = i[1][1]['vendors'],
                                                          primary=primary,
                                                          secondary=secondary,
                                                          source = source
                                                          )
                    experiment_entry.save()

                    for ex in i[1][1]['biasdata']:
                        # print('--saving---', '\n')

                        experiment_assay = AnalyzedAssay(experiment=experiment_entry,
                                                         order_no=ex['order_no'],
                                                         signalling_protein= ex['signalling_protein'],
                                                         quantitive_measure_type=ex['quantitive_measure_type'],
                                                         quantitive_activity=ex['quantitive_activity'],
                                                         quantitive_unit=ex['quantitive_unit'],
                                                         potency=ex['potency'],
                                                         assay_description = ex['assay_description']
                                                         )
                        experiment_assay.save()
                        # print('saved')
                else:
                    pass
                    print("already defined")
            except:
                #print('\n---saving error---', msg)
                continue
