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

    def fetch_vendor(self, ligand):

        links = LigandVendorLink.objects.filter(lp=ligand.properities.id)
        vendor_count = 0
        for x in links:

            vendor_count = vendor_count + 1

        return vendor_count



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

    def limit_family(self,send):

        assay = list()
        families = list()
        G12 = dict()
        Gio = dict()
        Gq = dict()
        GS = dict()
        Barr = dict()
        for i in send:
            if i['signalling_protein'] == 'Barr':
                if bool(Barr) == False:
                    Barr = i
                else:
                    if i['quantitive_activity'] < Barr['quantitive_activity']:
                        Barr = i

            if i['signalling_protein'] == 'G12':
                if bool(G12) == False:
                    G12 = i
                else:
                    if i['quantitive_activity'] < G12['quantitive_activity']:
                        G12 = i

            if i['signalling_protein'] == 'Gio':
                if bool(Gio) == False:
                    Gio = i
                else:
                    if i['quantitive_activity'] < Gio['quantitive_activity']:
                        Gio = i

            if i['signalling_protein'] == 'Gq':
                if bool(Gq) == False:
                    Gq = i
                else:
                    if i['quantitive_activity'] < Gq['quantitive_activity']:
                        Gq = i

            if i['signalling_protein'] == 'GS':
                if bool(GS) == False:
                    GS = i
                else:
                    if i['quantitive_activity'] < GS['quantitive_activity']:
                        GS = i

        if Barr:
            families.append(Barr)
        if G12:
            families.append(G12)
        if Gio:
            families.append(Gio)
        if Gq:
            families.append(Gq)
        if GS:
            families.append(GS)
        return families


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

        content = AssayExperiment.objects.filter(standard_units = 'nM').exclude(standard_value__isnull=True).prefetch_related(
            'publication', 'ligand', 'receptor'
            ).order_by('protein').order_by('ligand')
        print(len(content))
        # print('content', content)
        # # merge children
        pre_data = self.process_data(content)
        # # # transform to combined dictionary
        combined = self.change(pre_data)

        combined = self.process_group(combined)
        print('len1 - combined',len(combined))
        self.family_reduction(combined)
        keys = [k for k, v in combined.items() if len(v['biasdata']) <= 1]
        for x in keys:
            del combined[x]
        print('len2 - combined',len(combined))
        context.update({'data': combined})
        print('len3 - context',len(context['data']))
        #
        #
        for i in context['data'].items():

            # try:
            # i[1].pop('reference')
            # i[1].pop('biasdata')
                # TODO: move by one tab when uncomment try catch
            source = 'chembl3'
            assays = dict()

            # if self.fetch_experiment(i[1]['publication'], i[1]['ligand'], i[1]['receptor'], i[1]['residue'], i[1]['mutation'], source) == False:
            #


            primary, secondary = self.fetch_receptor_trunsducers(i[1]['receptor'])
            experiment_entry = AnalyzedExperiment(publication=i[1]['publication'],
                                                  ligand=i[1]['ligand'],
                                                  receptor=i[1]['receptor'],
                                                  chembl=i[1]['chembl'],
                                                  vendor_quantity = 1,
                                                  primary=primary,
                                                  secondary=secondary,
                                                  source = source
                                                  )
            experiment_entry.save()

            for ex in i[1]['biasdata']:
                # print('--saving---', '\n')

                experiment_assay = AnalyzedAssay(experiment=experiment_entry,
                                                 order_no=ex['order_no'],
                                                 assay_type= ex['assay_type'],
                                                 signalling_protein= ex['signalling_protein'],
                                                 quantitive_measure_type=ex['quantitive_measure_type'],
                                                 quantitive_activity=ex['quantitive_activity'],
                                                 quantitive_unit=ex['quantitive_unit'],
                                                 potency=ex['potency'],
                                                 assay_description = ex['assay_description']
                                                 )
                experiment_assay.save()
                    # print('saved')
            # else:
            #     pass
            #     print("already defined")
            # # except:
            #     print('\n---saving error---')
            #     continue

    def family_reduction(self, combined):
        for i in combined.items():
            i[1]['biasdata'] = self.limit_family(i[1]['biasdata'])


    def process_data(self, content):
        '''
        Merge BiasedExperiment with its children
        and pass it back to loop through dublicates
        '''
        rd = []
        counter = 0

        for instance in content:
            temp_obj = []
            fin_obj = {}
            fin_obj['main'] = instance
            rd.append(fin_obj)
        return rd

    def change(self, rd):
        '''
        Merge bias experminet data with assay data
        Define family in accordance with subfamily
        '''
        send = dict()
        increment = 0

        for j in rd:

            temp_dict = dict()
            temp = dict()
            doubles = []
            temp['publication'] = j['main'].publication
            temp['ligand'] =  j['main'].ligand
            temp['vendors'] = None
            temp['chembl'] = j['main'].chembl
            temp['receptor'] = j['main'].receptor
            if temp['receptor']:
                temp['species'] = temp['receptor'].species.common_name
            else:
                temp['species'] = None
            temp_dict['signalling_protein'] = self.define_protein(j['main'].assay_description)
            if temp_dict['signalling_protein'] == None:
                continue
            temp_dict['family'] = temp_dict['signalling_protein']
            temp_dict['assay_description'] = j['main'].assay_description
            temp_dict['quantitive_measure_type'] = j['main'].standard_type
            temp_dict['quantitive_unit'] = j['main'].standard_units
            if j['main'].publication:
                temp_dict['assay_type'] = j['main'].publication.web_link.index
            temp_dict['potency'] = None
            if j['main'].standard_value:
                temp_dict['quantitive_activity'] = float(j['main'].standard_value)
                if temp_dict['quantitive_unit'].lower() == 'nm':
                    temp_dict['quantitive_activity'] = temp_dict['quantitive_activity']* 10**(-9)
                elif temp_dict['quantitive_unit'].lower() == 'µm':
                    temp_dict['quantitive_activity'] = temp_dict['quantitive_activity']* 10**(-9)
                elif temp_dict['quantitive_unit'].lower() == 'pm':
                    temp_dict['quantitive_activity'] = temp_dict['quantitive_activity']* 10**(-12)
                elif temp_dict['quantitive_unit'].lower() == 'mm':
                    temp_dict['quantitive_activity'] = temp_dict['quantitive_activity']* 10**(-6)
            else:
                temp_dict['quantitive_activity'] = None

            doubles.append(temp_dict)
            temp['assay'] = doubles
            send[increment] = temp
            increment = increment + 1
        return send

    def define_protein(self, description):
        # if description exists
        signalling_protein = None
        if description:
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
                elif word.lower == 'g12' or word == 'G12':
                    G12 = G12+100
                elif word.lower == 'gq' or word =='Gq':
                    Gq = Gq+100
                elif word.lower == 'gs' or word == 'GS':
                    GS = GS+100
                elif word.lower == 'arrestin' or word == 'barrestin':
                    Barr = Barr+100

                elif word == 'camp':
                    for key in description.split():
                        if key == 'inhibition' or key == 'decrease' or key == 'forskolin-stimulated':
                            Gio = Gio+10
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

            if G12 > max(Gio,Gq,GS,Barr):
                signalling_protein = 'G12'
            elif Gio > max(G12,Gq,GS,Barr):
                signalling_protein = 'Gio'
            elif Gq > max(G12,Gio,GS,Barr):
                signalling_protein = 'Gq'
            elif GS > max(G12,Gq,Gio,Barr):
                signalling_protein = 'GS'
            elif Barr > max(G12,Gq,GS,Gio):
                signalling_protein = 'Barr'
        return signalling_protein

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

        self.process_calculation(context)
        # for i in context.items():
        #     print('---send[increment]---', i[1])
        print('len',len(context))

        return context

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
            # TODO: Change 9999999 to normal method that skips None value
            # self.convert_activity(i[1]['assay'])
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
