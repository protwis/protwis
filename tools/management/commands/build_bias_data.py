from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db import IntegrityError
from django.utils.text import slugify
from django.http import HttpResponse, JsonResponse
from decimal import Decimal
from build.management.commands.base_build import Command as BaseBuild
from common.tools import fetch_from_cache, save_to_cache, fetch_from_web_api
from residue.models import Residue
from protein.models import Protein, ProteinGProteinPair
from ligand.models import *
from mutation.models import Mutation
from ligand.functions import get_or_make_ligand
from common.models import WebLink, WebResource, Publication
from django.db import connection
import logging
import os
import queue
import logging
import os
from datetime import datetime
import xlrd
import operator
import traceback
import time
import math
import json
import threading
import concurrent.futures
import pytz


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
    structure_data_dir = 'excel/'
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
                print('Started purging bias data')
                self.purge_bias_data()
                print('Ended purging bias data')
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
        # import the structure data
        try:
            print('CREATING BIAS DATA')
            print(options['filename'])
            self.build_bias_data()
            self.build_bias_data_subtypes()
            self.logger.info('COMPLETED CREATING BIAS')
        except Exception as msg:
            print('--error--', msg, '\n')
            self.logger.info("The error appeared in def handle")

    def purge_bias_data(self):
        
        delete_bias_experiment = AnalyzedExperiment.objects.all()
        delete_bias_experiment.delete()

    def fetch_publication_authors(self, publication, experiment_assay):
        counter = 0

        author_list = list()
        if publication.authors != None:

            for authors in publication.authors.split(','):
                author_list.append(authors.strip())

            author_list.reverse()
            for i in author_list:
                if counter < 3:
                    assay_author = ExperimentAssayAuthors(experiment=experiment_assay,
                                                          author=i)
                    assay_author.save()
                    counter = counter + 1
            # assay_author = ExperimentAssayAuthors(experiment = experiment_assay,

    def fetch_measurements(self, potency, p_type, unit):
        if p_type.lower() == 'pec50':
            potency = 10**(potency * (-1))
            # pp = (-1)*log(potency)
            p_type = 'EC50'
        elif p_type.lower() == 'logec50':
            potency = 10**(potency)
            p_type = 'EC50'
        elif p_type.lower() == 'pic50':
            potency = 10**(potency * (-1))
            p_type = 'IC50'
        elif p_type.lower() == 'logic50':
            potency = 10**(potency)
            p_type = 'IC50'

        if p_type.lower() == 'ec50':
            if unit.lower() == 'nm':
                potency = potency * 10**(-9)
            elif unit.lower() == 'µm':
                potency = potency * 10**(-9)
            elif unit.lower() == 'pm':
                potency = potency * 10**(-12)
            elif unit.lower() == 'mm':
                potency = potency * 10**(-6)
            else:
                pass
        if p_type.lower() == 'ic50':
            if unit.lower() == 'nm':
                potency = potency * 10**(-9)
            elif unit.lower() == 'µm':
                potency = potency * 10**(-9)
            elif unit.lower() == 'pm':
                potency = potency * 10**(-12)
            elif unit.lower() == 'mm':
                potency = potency * 10**(-6)
            else:
                pass
        # if potency:
        #     potency = (-1)*math.log10(potency)
        #     p_type = 'pec50'
        #     # potency = "{:.2E}".format(Decimal(potency))
        return potency, p_type

    def define_g_family(self, protein, assay_type):

        family = None
        if (protein == 'b-arrestin' or
            protein == 'b-arrestin-1 (non-visual arrestin-2)' or
                protein == 'b-arrestin-2 (non-visual arrestin-3)'):
            family = 'B-arr'

        elif (protein == 'gi/o-family' or
                protein == 'gai1' or
                protein == 'gai2' or
                protein == 'gai3' or
                protein == 'gao' or
                protein == 'gaoA' or
                protein == 'gai' or
                protein == 'gai1' or
                protein == 'gai2' or
                protein == 'gai3' or
                protein == 'gai1/2' or
                protein == 'gao' or
                protein == 'gaoA' or
                protein == 'gaoB' or
                protein == 'gao1' or
                protein == 'gat1' or
                protein == 'gat2' or
                protein == 'gat3' or
                protein == 'gaz' or
                protein == 'gaoB'):
            family = 'Gi/o'

        elif (protein == 'gq-family' or
                protein == 'ga12' or
                protein == ' gaq' or
                protein == 'gaq/11' or
                protein == 'gaq/14' or
                protein == 'gaq/15' or
                protein == 'gaq/16'):
            family = 'Gq/11'

        elif (protein == 'g12/13-family' or
                protein == 'ga12' or
                protein == 'ga13'):
            family = 'G12/13'

        elif (protein == 'gs-family' or
              protein == 'gas' or
              protein == 'gaolf'):
            family = 'Gs'

        elif (protein == '' or
              protein == None):
            if assay_type == 'pERK1/2 activation' or assay_type == "pERK1-2":
                family = 'pERK1-2'
        else:
            family == protein

        return family

    def fetch_endogenous(self, protein):
        try:
            with connection.cursor() as cursor:
                cursor.execute(
                    "SELECT * FROM protein_endogenous_ligands WHERE protein_id =%s", [protein.pk])
                row = cursor.fetchone()
                end_ligand = Ligand.objects.filter(id=row[2])
                test = end_ligand.get()

            return test
        except:
            return None

    def fetch_vendor(self, ligand, experiment_entry):
        temp = ligand
        links = LigandVendorLink.objects.filter(lp=ligand.properities.id)
        # vendor_count = 0
        for x in links:
            if x.vendor.name not in ['ZINC', 'ChEMBL', 'BindingDB', 'SureChEMBL', 'eMolecules', 'MolPort', 'PubChem']:
                ligand_vendor = BiasedExperimentVendors(experiment=experiment_entry,
                                                        vendor=x)
                ligand_vendor.save()
            # vendor_count = vendor_count + 1

        # return vendor_count

    def fetch_chembl(self, ligand):
        temp = ligand
        chembl_id = None
        links = temp.properities.web_links.all()

        for x in links:
            if x.web_resource.slug == 'chembl_ligand':
                chembl_id = [
                    x for x in links if x.web_resource.slug == 'chembl_ligand'][0].index
        return chembl_id

    def fetch_protein(self, protein_from_excel, source):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        test = None
        if Protein.objects.filter(entry_name=protein_from_excel):
            protein = Protein.objects.filter(entry_name=protein_from_excel)
            test = protein.get()
        elif Protein.objects.filter(web_links__index=protein_from_excel, web_links__web_resource__slug='uniprot'):
            protein1 = Protein.objects.filter(
                web_links__index=protein_from_excel, web_links__web_resource__slug='uniprot')
            test = protein1[0]
        # if test == None:
        #     print('---protein error---',protein_from_excel,source )
        return test

    def fetch_ligand(self, ligand_id, ligand_type, ligand_name, source_file):
        """
        fetch ligands with Ligand model
        requires: ligand id, ligand id type, ligand name
        requires: source_file name
        """
        l = None

        try:
            if ligand_id in self.ligand_cache:
                l = self.ligand_cache[ligand_id]
            else:
                l = get_or_make_ligand(ligand_id, ligand_type, ligand_name)
                self.ligand_cache[ligand_id] = l
        except Exception as msg:
            l = None
            # print('ligand_id---',l,'\n end')
        return l

    def fetch_publication(self, publication_doi):
        """
        fetch publication with Publication model
        requires: publication doi or pmid
        """
        try:
            float(publication_doi)
            publication_doi = str(int(publication_doi))
        except ValueError:
            pass

        if publication_doi.isdigit():  # assume pubmed
            pub_type = 'pubmed'
        else:  # assume doi
            pub_type = 'doi'
        if publication_doi not in self.publication_cache:
            try:
                wl = WebLink.objects.get(
                    index=publication_doi, web_resource__slug=pub_type)
            except WebLink.DoesNotExist:

                try:
                    wl = WebLink.objects.create(index=publication_doi,
                                                web_resource=WebResource.objects.get(slug=pub_type))
                except IntegrityError:

                    wl = WebLink.objects.get(
                        index=publication_doi, web_resource__slug=pub_type)


            try:
                pub = Publication.objects.get(web_link=wl)
            except Publication.DoesNotExist:
                pub = Publication()
                try:
                    pub.web_link = wl
                    pub.save()
                except IntegrityError:
                    pub = Publication.objects.get(web_link=wl)

                if pub_type == 'doi':
                    pub.update_from_doi(doi=publication_doi)
                elif pub_type == 'pubmed':
                    pub.update_from_pubmed_data(index=publication_doi)
                try:
                    pub.save()
                except:
                    self.mylog.debug(
                        "publication fetching error | module: fetch_publication. Row # is : " + str(publication_doi) + ' ' + pub_type)
                    # if something off with publication, skip.
            self.publication_cache[publication_doi] = pub
        else:
            pub = self.publication_cache[publication_doi]

        return pub

    def fetch_experiment(self, publication, ligand, receptor, source):
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
                "Experiment AnalyzedExperiment error | module: AnalyzedExperiment.")
            return False



    # Bias data block #
    def build_bias_data(self):
        print('Build bias data gproteins')
        context = dict()
        context_subtypes = dict()
        # get data_all
        content = self.get_from_model()
        print('stage # 1 : Getting data finished')
        # process data_all
        content_with_children = self.process_data(content)
        print('stage # 2: Processing children in queryset finished')
        # change data
        changed_data = self.queryset_to_dict(content_with_children)
        print('stage # 3: Converting queryset into dict finished')
        # process references
        data_with_references = self.process_references(changed_data)
        print('stage # 4: Separating reference assays is finished')
        # link reference assays
        send = self.process_reference_assays(data_with_references)
        send_subtype = send
        print('stage # 5: affiliating test assays with references finished')
        # merging experiments for ligand/receptor/publication,
        # calculation bias factors / labs / articles
        results_temp = self.process_group(send)
        print('stage # 6: Merging assays with same ligand/receptor/publication is finished')
        # prepare dataset for saving to model
        context.update({'data': results_temp})
        print('stage # 7: combining data into common dict is finished')
        # save dataset to model
        self.save_data_to_model(context, 'different_family')
        print('stage # 8: saving data to model is finished')

    def build_bias_data_subtypes(self):
        print('Build bias data for subtypes')
        context = dict()
        context_subtypes = dict()
        # get data_all
        content = self.get_from_model()
        print('stage # 1 : Getting data finished')
        # process data_all
        content_with_children = self.process_data(content)
        print('stage # 2: Processing children in queryset finished')
        # change data
        changed_data = self.queryset_to_dict(content_with_children)
        print('stage # 3: Converting queryset into dict finished')
        # process references
        data_with_references = self.process_references(changed_data)
        print('stage # 4: Separating reference assays is finished')
        # link reference assays
        send = self.process_reference_assays(data_with_references)
        send_subtype = send
        print('stage # 5: affiliating test assays with references finished')
        # merging experiments for ligand/receptor/publication,
        # calculation bias factors / labs / articles
        results_temp = self.process_group_same_family(send)
        print('stage # 6: Merging assays with same ligand/receptor/publication is finished')
        # prepare dataset for saving to model
        context.update({'data': results_temp})
        print('stage # 7: combining data into common dict is finished')
        # save dataset to model
        self.save_data_to_model(context, 'same_family')
        print('stage # 8: saving data to model is finished')


    def save_data_to_model(self, context, source):
        for i in context['data'].items():

            if self.fetch_experiment(i[1]['publication'], i[1]['ligand'], i[1]['receptor'], source) == False:
                primary, secondary = self.fetch_receptor_trunsducers(
                    i[1]['receptor'])
                # primary = None
                # secondary = None

                experiment_entry = AnalyzedExperiment(publication=i[1]['publication'],
                                                      ligand=i[1]['ligand'],
                                                      receptor=i[1]['receptor'],
                                                      chembl=i[1]['chembl'],
                                                      source=source,
                                                      endogenous_ligand=i[1]['endogenous_ligand'],
                                                      vendor_quantity=i[1]['vendor_counter'],
                                                      reference_ligand=i[1]['reference_ligand'],
                                                      primary=primary,
                                                      secondary=secondary,
                                                      article_quantity=i[1]['article_quantity'],
                                                      labs_quantity=i[1]['labs'],
                                                      )
                experiment_entry.save()
                for ex in i[1]['biasdata']:
                    if ex['quantitive_activity'] is not None:
                        ex['quantitive_activity'] = '%.2E' % Decimal(
                            ex['quantitive_activity'])
                    # print('--saving---', '\n')
                    emax_ligand = ex['emax_reference_ligand']
                    experiment_assay = AnalyzedAssay(experiment=experiment_entry,
                                                     family=ex['family'],
                                                     order_no=ex['order_no'],
                                                     signalling_protein=ex['signalling_protein'],
                                                     cell_line=ex['cell_line'],
                                                     assay_type=ex['assay_type'],
                                                     assay_measure=ex['assay_measure_method'],
                                                     assay_time_resolved=ex['assay_time_resolved'],
                                                     ligand_function=ex['ligand_function'],
                                                     quantitive_measure_type=ex['quantitive_measure_type'],
                                                     quantitive_activity=ex['quantitive_activity'],
                                                     quantitive_activity_initial=ex['quantitive_activity_initial'],
                                                     quantitive_unit=ex['quantitive_unit'],
                                                     qualitative_activity=ex['qualitative_activity'],
                                                     quantitive_efficacy=ex['quantitive_efficacy'],
                                                     efficacy_measure_type=ex['efficacy_measure_type'],
                                                     efficacy_unit=ex['efficacy_unit'],
                                                     potency=ex['potency'],
                                                     t_coefficient=ex['t_coefficient'],
                                                     t_value=ex['t_coefficient_initial'],
                                                     t_factor=ex['t_factor'],
                                                     log_bias_factor=ex['log_bias_factor'],
                                                     emax_ligand_reference=emax_ligand
                                                     )
                    experiment_assay.save()

            else:
                pass

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
                temp += str(i) + str(', ')

            for i in secondary:
                temp1 += str(i) + str(', ')
            return temp, temp1
        except:
            print('receptor not found', receptor)
            return None,None

    def fetch_experiment(self, publication, ligand, receptor, source):
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

    def get_from_model(self):
        try:
            content = BiasedExperiment.objects.all().prefetch_related(
                'experiment_data', 'ligand', 'receptor', 'publication', 'publication__web_link', 'experiment_data__emax_ligand_reference',
            ).order_by('publication', 'receptor', 'ligand')
        except EmptyResultSet:
            print('no data returned')
            content = None
        return content

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

            vendor_counter = 0
            for i in instance[1].experiment_data_vendors.all():
                vendor_counter = vendor_counter + 1

            for entry in instance[1].experiment_data.all():
                author_list = list()
                for author in entry.experiment_data_authors.all():
                    author_list.append(author.author)

                temp_obj.append(entry)
                counter += 1
            fin_obj['authors'] = author_list
            fin_obj['children'] = temp_obj
            fin_obj['vendor_counter'] = vendor_counter
            rd.append(fin_obj)
        return rd

    def queryset_to_dict(self, rd):
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
            temp['publication'] = j['main'].publication
            if j['main'].receptor:
                temp['species'] = j['main'].receptor.species.common_name
            else:
                temp['species'] = None
            temp['ligand'] = j['main'].ligand
            temp['endogenous_ligand'] = j['main'].endogenous_ligand
            temp['chembl'] = j['main'].chembl
            temp['receptor'] = j['main'].receptor
            temp['residue'] = j['main'].residue
            temp['mutation'] = j['main'].mutation
            temp['assay'] = dict()
            temp['vendor_counter'] = j['vendor_counter']
            temp['reference'] = list()
            temp['authors'] = j['authors']
            temp['ref_ligand_experiment'] = dict()
            temp['article_quantity'] = 0
            temp['labs_quantity'] = 0
            temp['reference_ligand'] = None
            if not j['children']:
                continue
            else:
                temp_dict['potency'] = ''
                temp_dict['t_factor'] = ''
                temp_dict['log_bias_factor'] = ''
                temp_dict['order_no'] = 0
                temp_dict['reference_ligand'] = None

                temp_dict['signalling_protein'] = j['children'][0].signalling_protein
                temp_dict['cell_line'] = j['children'][0].cell_line
                temp_dict['family'] = j['children'][0].family
                temp_dict['assay_type'] = j['children'][0].assay_type
                temp_dict['assay_measure_method'] = j['children'][0].assay_measure
                temp_dict['assay_time_resolved'] = j['children'][0].assay_time_resolved
                if j['children'][0].quantitive_activity:
                    temp_dict['quantitive_activity'] = j['children'][0].quantitive_activity
                    temp_dict['quantitive_activity_initial'] = j['children'][0].quantitive_activity
                else:
                    temp_dict['quantitive_activity'] = None
                    temp_dict['quantitive_activity_initial'] = None
                temp_dict['qualitative_activity'] = j['children'][0].qualitative_activity
                temp_dict['quantitive_unit'] = j['children'][0].quantitive_unit
                temp_dict['quantitive_efficacy'] = j['children'][0].quantitive_efficacy
                temp_dict['efficacy_unit'] = j['children'][0].efficacy_unit
                temp_dict['quantitive_measure_type'] = j['children'][0].quantitive_measure_type
                temp_dict['efficacy_measure_type'] = j['children'][0].efficacy_measure_type
                temp_dict['t_coefficient'] = j['children'][0].bias_value
                temp_dict['t_coefficient_initial'] = j['children'][0].bias_value_initial

                temp_dict['bias_reference'] = j['children'][0].bias_reference
                temp_dict['emax_reference_ligand'] = j['children'][0].emax_ligand_reference
                temp_dict['ligand_function'] = j['children'][0].ligand_function
                temp_dict['ligand'] = j['main'].ligand

                if temp_dict['quantitive_activity_initial']:
                    temp_dict['quantitive_activity_initial'] = (-1) * math.log10(
                        temp_dict['quantitive_activity_initial'])
                    temp_dict['quantitive_activity_initial'] = "{:.2F}".format(
                        Decimal(temp_dict['quantitive_activity_initial']))

                temp['ref_ligand_experiment'] = j['children'][0].emax_ligand_reference

                doubles.append(temp_dict)
                temp['assay'] = doubles
                send[increment] = temp
                increment = increment + 1
                counter += 1

        return send

    def process_references(self, dictionary):
        context = dict()
        send = dict()
        counter = 0
        for j in dictionary.items():
            name = str(j[1]['publication'])
            temp_obj = list()

            if name in context:
                counter += 1
                temp_obj = context[name]

            temp_obj.append(j[1])
            context[name] = temp_obj
            context[name][0]['publication'] == name

        return context

    def process_reference_assays(self, context):
        '''
        Recieve data from "data_with_references"
        search for objects with same publication ligand receptor
        Create new object of biased_experiment and all children
        '''
        send = dict()
        increment = 0

        for j in context.items():

            ref = list()
            for data in j[1]:
                if data['assay'][0]['bias_reference'].lower() != "" and data['assay'][0]['bias_reference'] == 'Reference':
                    if data in ref:
                        pass
                    else:
                        ref.append(data)

            for data in j[1]:
                for i in ref:

                    if (data['receptor'] == i['receptor'] and
                        data['species'] == i['species'] and
                        data['assay'][0]['assay_type'] == i['assay'][0]['assay_type'] and
                        data['assay'][0]['signalling_protein'] == i['assay'][0]['signalling_protein'] and
                        data['assay'][0]['cell_line'] == i['assay'][0]['cell_line'] and
                        data['assay'][0]['assay_measure_method'] == i['assay'][0]['assay_measure_method'] and
                            data['assay'][0]['bias_reference'] != 'Reference'):
                        data['assay'][0]['reference_quantitive_activity'] = i['assay'][0]['quantitive_activity']
                        data['assay'][0]['reference_quantitive_efficacy'] = i['assay'][0]['quantitive_efficacy']
                        data['assay'][0]['reference_t_coefficient_initial'] = i['assay'][0]['t_coefficient_initial']
                        data['assay'][0]['reference_ligand'] = i['assay'][0]['ligand']
                        data['reference_ligand'] = i['assay'][0]['ligand']
                        data['assay'][0]['reference_measure_type'] = i['assay'][0]['quantitive_measure_type']

            for data in j[1]:
                if data['assay'][0]['bias_reference'] != 'Reference':
                    send[increment] = data
                    increment += 1

        return send

    def process_group(self, send):
        '''
        Recieve data from "process_data"
        search for objects with same publication ligand receptor
        Create new object of biased_experiment and all children
        '''
        context = dict()

        for j in send.items():

            name = str(j[1]['publication']) + \
                '/' + str(j[1]['ligand']) + '/' + str(j[1]['receptor'])
            temp_obj = list()

            if(name in context):
                temp_obj = context[name]['assay']
                reference_obj = context[name]['reference']

            for i in j[1]['assay']:
                if i not in temp_obj:
                    temp_obj.append(i)

            context[name] = j[1]
            context[name]['assay'] = temp_obj
            # select most potent assay for every subtype if >1 assays exists
            context[name]['assay'] = self.limit_family(context[name]['assay'])

        self.process_calculation(context)
        self.count_publications(context)
        return context

    def process_group_same_family(self, send):
        '''
        Recieve data from "process_data"
        search for objects with same publication ligand receptor
        Create new object of biased_experiment and all children
        '''
        context = dict()

        for j in send.items():

            name = str(j[1]['publication']) + \
                '/' + str(j[1]['ligand']) + '/' + str(j[1]['receptor'])
            temp_obj = list()

            if(name in context):
                temp_obj = context[name]['assay']
                reference_obj = context[name]['reference']

            for i in j[1]['assay']:
                if i not in temp_obj:
                    temp_obj.append(i)

            context[name] = j[1]
            context[name]['assay'] = temp_obj
            # select most potent assay for every subtype if >1 assays exists

        self.process_calculation(context)
        self.count_publications(context)
        return context

    def limit_family(self, send):
        assay = list()
        families = list()
        G12 = dict()
        Gio = dict()
        Gq = dict()
        GS = dict()
        Barr = dict()
        for i in send:
            try:
                if i['family'] == 'B-arr':
                    if bool(Barr) == False:
                        Barr = i
                    else:
                        if i['quantitive_activity'] < Barr['quantitive_activity']:
                            Barr = i
            except:
                continue
            try:
                if i['family'] == 'G12/13':
                    if bool(G12) == False:
                        G12 = i
                    else:
                        if i['quantitive_activity'] < G12['quantitive_activity']:
                            G12 = i
            except:
                continue

            try:
                if i['family'] == 'Gi/o':
                    if bool(Gio) == False:
                        Gio = i
                    else:
                        if i['quantitive_activity'] < Gio['quantitive_activity']:
                            Gio = i
            except:
                continue
            try:
                if i['family'] == 'Gq/11':
                    if bool(Gq) == False:
                        Gq = i
                    else:
                        if i['quantitive_activity'] < Gq['quantitive_activity']:
                            Gq = i
            except:
                continue
            try:
                if i['family'] == 'Gs':
                    if bool(GS) == False:
                        GS = i
                    else:
                        if i['quantitive_activity'] < GS['quantitive_activity']:
                            GS = i
            except:
                continue

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

    def process_calculation(self, context):
        for i in context.items():
            test = dict()
            lgb_refine = dict()
            temp_obj = list()
            for j in i[1]['assay']:
                if j not in temp_obj:
                    temp_obj.append(j)
                else:
                    pass
            i[1]['assay'] = temp_obj
            test = sorted(i[1]['assay'], key=lambda k: k['quantitive_activity']
                          if k['quantitive_activity'] else 999999,  reverse=False)

            for x in enumerate(test):
                x[1]['order_no'] = x[0]

            i[1]['biasdata'] = test
            i[1].pop('assay')
            # calculate log bias
            self.calc_bias_factor(i[1]['biasdata'])
            # recalculates lbf if it is negative
            i[1]['biasdata'] = self.validate_lbf(i)

            self.calc_potency_and_transduction(i[1]['biasdata'])

    def validate_lbf(self, i):
        for x in i[1]['biasdata']:
            if x['log_bias_factor'] and isinstance(x['log_bias_factor'], int) and x['log_bias_factor'] < 0:
                for j in i[1]['biasdata']:
                    if j['order_no'] == x['order_no'] + 1:
                        j['order_no'], x['order_no'] = x['order_no'], j['order_no']

                self.calc_bias_factor(i[1]['biasdata'])
                validate_lbf(i[1]['biasdata'])
        return i[1]['biasdata']

    def calc_bias_factor(self, biasdata):
        most_reference = dict()
        most_potent = dict()
        for i in biasdata:
            if i['order_no'] == 0:
                most_potent = i
                i['log_bias_factor'] = None

        for i in biasdata:
            temp_reference = dict()
            try:
                if i['order_no'] != 0:
                    if (i['quantitive_measure_type'].lower() == 'ec50' and i['reference_measure_type'].lower() == 'ec50' and
                            most_potent['quantitive_measure_type'].lower() == 'ec50' and most_potent['reference_measure_type'].lower() == 'ec50'):
                        a = 0
                        b = 0
                        c = 0
                        d = 0
                        a = math.log10(
                            most_potent['quantitive_efficacy'] / most_potent['quantitive_activity'])
                        b = math.log10(
                            most_potent['reference_quantitive_efficacy'] / most_potent['reference_quantitive_activity'])
                        c = math.log10(
                            i['quantitive_efficacy'] / i['quantitive_activity'])
                        d = math.log10(
                            i['reference_quantitive_efficacy'] / i['reference_quantitive_activity'])
                        temp_calculation = self.caclulate_bias_factor_variables(
                            a, b, c, d)
                        i['log_bias_factor'] = round(temp_calculation, 1)
                    elif (i['quantitive_measure_type'].lower() == 'ic50' and temp_reference['quantitive_measure_type'].lower() == 'ic50'):
                        i['log_bias_factor'] = 'Only agonist in main pathway'

                    elif (i['qualitative_activity'] == 'No activity'):
                        i['log_bias_factor'] = "Full Bias"

                    elif (i['qualitative_activity'] == 'Low activity'):
                        i['log_bias_factor'] = "High Bias"

                    elif (i['qualitative_activity'] == 'High activity'):
                        i['log_bias_factor'] = "Low Bias"

                else:
                    i['log_bias_factor'] = None
            except:
                i['log_bias_factor'] = None

    def caclulate_bias_factor_variables(self, a, b, c, d):
        """calculations for log bias factor inputs"""
        lgb = 0
        try:
            lgb = (a - b) - (c - d)
        except:
            lbg = 0
        return lgb

    def calc_potency_and_transduction(self, biasdata):
        count = 0
        most_potent = dict()
        for i in biasdata:
            count += 1
            if i['order_no'] == 0:
                most_potent = i
        # T_factor -- bias factor
        for i in biasdata:
            if i['order_no'] > 0:
                if i['quantitive_measure_type'].lower() == 'ec50' or i['quantitive_measure_type'].lower() == 'ic50':

                    if i['quantitive_activity'] is not None and i['quantitive_activity'] != 0 and most_potent['quantitive_activity'] is not None:
                        i['potency'] = round(
                            i['quantitive_activity'] / most_potent['quantitive_activity'], 1)
                    elif i['quantitive_measure_type'].lower() == 'pec50' or i['quantitive_measure_type'].lower() == 'pic50':
                        # TODO: round
                        i['potency'] = round(
                            most_potent['quantitive_activity'] - i['quantitive_activity'], 1)
                    else:
                        i['potency'] = None

                if i['t_coefficient'] is not None and most_potent['t_coefficient'] is not None:
                    # TODO: validate if difference is non negative
                    i['t_factor'] = round(
                        most_potent['t_coefficient'] - i['t_coefficient'], 1)
                else:
                    i['t_factor'] = None

        # print('---counter of calc_potency---', count)

    def count_publications(self, context):
        temp = dict()
        for i in context.items():
            labs = list()
            i[1]['labs'] = 0
            labs.append(i[1]['publication'])
            lab_counter = 1
            for j in context.items():
                if j[1]['publication'] not in labs:
                    if set(i[1]['authors']) & set(j[1]['authors']):
                        lab_counter += 1
                        labs.append(j[1]['publication'])
                        i[1]['labs'] = lab_counter

            temp_obj = 1
            name = str(i[1]['ref_ligand_experiment']) + \
                '/' + str(i[1]['ligand']) + '/' + str(i[1]['receptor'])
            if(name in temp):
                for assays in i[1]['biasdata']:
                    if assays['order_no'] > 0:
                        if assays['log_bias_factor'] != None and assays['log_bias_factor'] != '' or assays['t_factor'] != None and assays['t_factor'] != '':
                            temp_obj = temp[name] + 1

            temp[name] = temp_obj

        for i in context.items():
            temp_obj = 0
            name = str(i[1]['ref_ligand_experiment']) + \
                '/' + str(i[1]['ligand']) + '/' + str(i[1]['receptor'])
            counter = 0
            if(name in temp):
                i[1]['article_quantity'] = temp[name]
