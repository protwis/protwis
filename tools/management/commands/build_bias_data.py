from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db import IntegrityError
from django.utils.text import slugify
from django.http import HttpResponse, JsonResponse

from build.management.commands.base_build import Command as BaseBuild
from common.tools import fetch_from_cache, save_to_cache, fetch_from_web_api
from residue.models import Residue
from protein.models import Protein
from ligand.models import Ligand, LigandProperities, LigandRole, LigandType, BiasedExperiment, ExperimentAssay, AnalyzedExperiment, AnalyzedAssay
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
            #self.prepare_all_data(options['filename'])


            self.logger.info('COMPLETED CREATING BIAS')

        except Exception as msg:
            print('--error--', msg, '\n')
            self.logger.info("The error appeared in def handle")

    def loaddatafromexcel(self, excelpath):
        """
        Reads excel file (require specific excel sheet)
        """
        num_rows = 0
        try:
            workbook = xlrd.open_workbook(excelpath)
            worksheets = workbook.sheet_names()
            temp = []

            for worksheet_name in worksheets:
                worksheet = workbook.sheet_by_name(worksheet_name)
                num_rows = worksheet.nrows - 1
                num_cells = worksheet.ncols - 1
                curr_row = 0  # skip first, otherwise -1
                while curr_row < num_rows:
                    curr_row += 1
                    row = worksheet.row(curr_row)
                    curr_cell = -1
                    temprow = []
                    while curr_cell < num_cells:
                        curr_cell += 1
                        cell_value = worksheet.cell_value(curr_row, curr_cell)
                        cell_type = worksheet.cell_type(curr_row, curr_cell)

                        # fix wrong spaced cells
                        if cell_value == " ":
                            cell_value = ""
                        temprow.append(cell_value)
                    temp.append(temprow)
                    # if curr_row>10: break
            return temp
        except:
            self.logger.info(
                "The error appeared during reading the excel", num_rows)

    def analyse_rows(self, rows, source_file):
        """
        Reads excel rows one by one
        Fetch data to models
        Saves to DB
        """
        skipped = 0
        # Analyse the rows from excel and assign the right headers
        temp = []
        print("start")
        for i, r in enumerate(rows, 1):
            # code to skip rows in excel for faster testing
            # if i < 0:
            #     continue
            # if i > 1:
            #     break
            print(i)
            d = {}
            if r[4] != '':  # checks if the ligand field exists
                try:
                    res = ''
                    mut = ''
                    d['submitting_group'] = r[0]
                    d['reference'] = r[1]

                    d['ligand_name'] = r[4]
                    d['ligand_type'] = r[5]
                    d['ligand_id'] = r[6]
                    d['bias_reference'] = r[7]
                    d['bias_ligand_name'] = r[8]
                    d['bias_ligand_type'] = r[9]
                    d['bias_ligand_id'] = r[10]

                    d['receptor'] = r[11].lower().strip()
                    d['receptor_mutant_aa_no'] = r[12]
                    d['receptor_wt'] = r[13].strip()
                    d['receptor_mut_aa'] = r[14]
                    d['cell_line'] = r[15]
                    d['protein'] = r[16].lower().strip()
                    d['protein_assay'] = r[17].strip()
                    d['protein_assay_method'] = r[18]
                    d['protein_time_resolved'] = r[19]

                    d['protein_ligand_function'] = r[20]
                    d['protein_mtype'] = r[21]
                    d['protein_activity_equation'] = r[22]
                    d['protein_activity_quantity'] = r[23]
                    d['protein_activity_quantity_unit'] = r[24]
                    d['protein_activity_quality'] = r[25]
                    d['protein_efficacy_measure'] = r[26]
                    d['protein_efficacy_equation'] = r[27]
                    d['protein_efficacy_quantity'] = r[28]
                    d['protein_efficacy_quantity_unit'] = r[29]
                    d['pathway_bias_initial'] = r[30]
                    d['pathway_bias'] = r[31]

                    d['source_file'] = source_file + str(i)

                    if not isinstance(d['ligand_id'], str):
                        d['ligand_id'] = int(d['ligand_id'])
                    if not isinstance(d['bias_ligand_id'], str):
                        d['bias_ligand_id'] = int(d['bias_ligand_id'])
                    # coverts string to object
                    if d['protein_activity_quantity'] == "":
                        d['protein_activity_quantity'] = None
                    if d['protein_efficacy_quantity'] == "":
                        d['protein_efficacy_quantity'] = None

                    if not isinstance(d['pathway_bias'], float):
                        d['pathway_bias'] = None
                    if not isinstance(d['pathway_bias_initial'], float):
                        d['pathway_bias_initial'] = None
                    # fetch publicaition
                    pub = self.fetch_publication(d['reference'])
                    # fetch main ligand
                    l = self.fetch_ligand(
                        d['ligand_id'], d['ligand_type'], d['ligand_name'], d['source_file'])
                    # fetch reference_ligand
                    reference_ligand = self.fetch_ligand(
                        d['bias_ligand_id'], d['bias_ligand_type'], d['bias_ligand_name'], d['source_file'])
                    # fetch protein
                    protein = self.fetch_protein(
                        d['receptor'], d['source_file'])
                    # fetch protein Residue
                    if d['receptor_wt'] is "":
                        res = None
                    else:
                        res = self.fetch_residue(
                            protein, d['receptor_mutant_aa_no'], d['receptor_wt'])
                    # fetch protein mutation
                    if d['receptor_mut_aa'] is "":
                        mut = None
                    else:
                        mut = self.fetch_mutation(
                            protein, res, d['receptor_mut_aa'], d['source_file'])

                    experiment_entry = BiasedExperiment(submission_author=d['submitting_group'],
                                                        publication=pub,
                                                        ligand=l,
                                                        receptor=protein,
                                                        mutation=mut,
                                                        residue=res
                                                        )
                    experiment_entry.save()
                    print('---bias---', d['bias_reference'], '\n')
                    experiment_assay = ExperimentAssay(biased_experiment=experiment_entry,
                                                       signalling_protein=d['protein'],
                                                       cell_line=d['cell_line'],
                                                       assay_type=d['protein_assay'],
                                                       assay_measure=d['protein_assay_method'],
                                                       assay_time_resolved=d['protein_time_resolved'],
                                                       ligand_function=d['protein_ligand_function'],
                                                       quantitive_measure_type=d['protein_mtype'],
                                                       quantitive_activity=d['protein_activity_quantity'],
                                                       quantitive_activity_sign=d['protein_activity_equation'],
                                                       quantitive_unit=d['protein_activity_quantity_unit'],
                                                       qualitative_activity=d['protein_activity_quality'],
                                                       quantitive_efficacy=d['protein_efficacy_quantity'],
                                                       efficacy_measure_type=d['protein_efficacy_measure'],
                                                       efficacy_sign=d['protein_efficacy_equation'],
                                                       efficacy_unit=d['protein_efficacy_quantity_unit'],
                                                       bias_reference=d['bias_reference'],
                                                       bias_value=d['pathway_bias'],
                                                       bias_value_initial=d['pathway_bias_initial'],
                                                       emax_ligand_reference=reference_ligand
                                                       )

                    experiment_assay.save()

                except Exception as msg:
                    print(msg)
                    self.mylog.exception(
                        "Experiment save error. Row: " + d['source_file'], msg)
                    print(d['source_file'], msg)
                    continue

            temp.append(d)
        return temp

    def fetch_mutation(self, protein, res, amino_acid, source):
        """
        fetch mutant with Mutation model
        required: receptorID
        required: residue
        required: mutant amino_acid
        """
        try:
            print("protein is {}, residue = {} and aa is {}".format(
                {protein}, {res}, {amino_acid}))
            print(protein.pk, res.pk)
            mut = Mutation.objects.filter(protein=protein,
                                          residue=res,
                                          amino_acid=amino_acid)
            mut = mut.get()

        except Exception as msg:
            # print("error message from mutation: " + source, msg)
            new_mut = Mutation(protein=protein,
                               residue=res,
                               amino_acid=amino_acid)
            new_mut.save()
            print("Mutation created")
            mut = self.fetch_mutation(protein, res, amino_acid, source)
            self.mylog.exception(
                "Mutation fetching error | module: fetch_mutation. Row # is : ", source)
        return mut

    def fetch_residue(self, protein, sequence_number, amino_acid):
        """"
        fetch residue with Residue model
        required: receptorID
        required: sequence_number, amino_acid
        """
        try:
            res = Residue.objects.filter(protein_conformation__protein=protein,
                                         sequence_number=int(sequence_number),
                                         amino_acid=amino_acid)
            res = res.get()
        except Exception as msg:
            print("error message from residue: ", msg)
            self.mylog.exception(
                "residue fetching error | module: fetch_residue. " + source)
            res = None
        return res

    def fetch_protein(self, protein_from_excel, source):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        try:
            protein = Protein.objects.filter(entry_name=protein_from_excel)
            protein = protein.get()
        except Exception as msg:
            protein = None
            self.mylog.exception(
                "Protein fetching error | module: fetch_protein. Row # is : ", source)
        return protein

    def fetch_ligand(self, ligand_id, ligand_type, ligand_name, source_file):
        """
        fetch ligands with Ligand model
        requires: ligand id, ligand id type, ligand name
        requires: source_file name
        """
        l = None
        if str(ligand_id) in self.ligand_cache:
            if ligand_id in self.ligand_cache[str(ligand_id)]:
                l = self.ligand_cache[str(ligand_id)][ligand_id]
        else:
            self.ligand_cache[str(ligand_id)] = {}

        if not l:
            try:
                l = get_or_make_ligand(
                    ligand_id, ligand_type, str(ligand_name))
            except Exception as msg:
                l = None
                self.ligand_cache[str(ligand_name), ligand_id] = l
                self.mylog.exception("Protein fetching error | module: fetch_ligand. Row # is : ",
                                     ligand_name, ligand_type, ligand_id, source_file)
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

    def fetch_experiment(self, publication, ligand, receptor):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        try:
            experiment = AnalyzedExperiment.objects.filter(
                publication=publication, ligand=ligand, receptor=receptor)
            experiment = experiment.get()
            return True
        except Exception as msg:
            experiment = None
            self.mylog.exception(
                "Protein AnalyzedExperiment error | module: AnalyzedExperiment.")
            return False

    def prepare_all_data(self, filenames):
        if not filenames:
            filenames = os.listdir(self.structure_data_dir)
        for source_file in filenames:
            # print("source_file " + str(source_file))
            source_file_path = os.sep.join(
                [self.structure_data_dir, source_file]).replace('//', '/')
            # print("source_file_path " + str(source_file_path))
            if os.path.isfile(source_file_path) and source_file[0] != '.':
                self.logger.info('Reading file {}'.format(source_file_path))
                print('Reading file {}'.format(source_file_path))
                # read the yaml file
                rows = []
                if source_file[-4:] == 'xlsx' or source_file[-3:] == 'xls':
                    if "~$" in source_file:
                        # ignore open excel files
                        continue
                    rows = self.loaddatafromexcel(source_file_path)
                    rows = self.analyse_rows(rows, source_file)
                else:
                    self.mylog.debug('unknown format'.source_file)
                    continue

                self.data_all += rows
        print(len(self.data_all), " total data points")
        print("Finished")

    def process_data(self, content):
        '''
        Merge BiasedExperiment with its children
        and pass it back to loop through dublicates
        '''
        rd = []
        for instance in enumerate(content):
            temp_obj = []
            fin_obj = {}
            fin_obj['main'] = (instance[1])
            for entry in instance[1].experiment_data.all():
                temp_obj.append(entry)
            fin_obj['children'] = temp_obj
            rd.append(fin_obj)

        return rd

    def change(self, rd):
        '''
        Merge bias experminet data with assay data
        Define family in accordance with subfamily
        '''
        test = list()
        send = dict()
        increment = 0

        for j in rd:

            temp_dict = dict()
            temp = dict()
            doubles = []
            temp['publication'] = j['main'].publication
            temp['ligand'] = j['main'].ligand
            temp['receptor'] = j['main'].receptor
            if not j['children']:
                continue
            else:
                if (j['children'][0].signalling_protein == 'β-arrestin' or
                    j['children'][0].signalling_protein == 'β-arrestin-1 (non-visual arrestin-2)' or
                        j['children'][0].signalling_protein == 'β-arrestin-2 (non-visual arrestin-3)'):
                    temp_dict['pathway'] = 'β-arrestin'

                elif (j['children'][0].signalling_protein == 'gi/o-family' or
                      j['children'][0].signalling_protein == 'gαi1' or
                      j['children'][0].signalling_protein == 'gαi2' or
                      j['children'][0].signalling_protein == 'gαi3' or
                      j['children'][0].signalling_protein == 'gαo' or
                      j['children'][0].signalling_protein == 'gαoA' or
                      j['children'][0].signalling_protein == 'gαoB'):
                    temp_dict['pathway'] = 'Gi/o-family'

                elif (j['children'][0].signalling_protein == 'gq-family' or
                        j['children'][0].signalling_protein == 'gαq' or
                        j['children'][0].signalling_protein == 'gαq11' or
                        j['children'][0].signalling_protein == 'gαq14' or
                        j['children'][0].signalling_protein == 'gαq14' or
                        j['children'][0].signalling_protein == 'gαq16' or
                        j['children'][0].signalling_protein == 'gαq14 (gαq16)'):
                    temp_dict['pathway'] = 'Gq-family'

                elif (j['children'][0].signalling_protein == 'g12/13-family' or
                      j['children'][0].signalling_protein == 'gα12' or
                      j['children'][0].signalling_protein == 'gα13'):
                    temp_dict['pathway'] = 'G12/13-family'

                elif (j['children'][0].signalling_protein == 'gs-family' or
                      j['children'][0].signalling_protein == 'gαs' or
                      j['children'][0].signalling_protein == 'gαolf'):
                    temp_dict['pathway'] = 'Gs-family'
                else:
                    temp_dict['pathway'] = 'No data'

                temp_dict['signalling_protein'] = j['children'][0].signalling_protein
                temp_dict['cell_line'] = j['children'][0].cell_line
                temp_dict['assay_type'] = j['children'][0].assay_type
                temp_dict['assay_measure_method'] = j['children'][0].assay_measure
                temp_dict['assay_time_resolved'] = j['children'][0].assay_time_resolved
                temp_dict['quantitive_activity'] = j['children'][0].quantitive_activity
                temp_dict['qualitative_activity'] = j['children'][0].qualitative_activity
                temp_dict['quantitive_unit'] = j['children'][0].quantitive_unit
                temp_dict['quantitive_efficacy'] = j['children'][0].quantitive_efficacy
                temp_dict['efficacy_unit'] = j['children'][0].efficacy_unit
                temp_dict['quantitive_measure_type'] = j['children'][0].quantitive_measure_type
                temp_dict['efficacy_measure_type'] = j['children'][0].efficacy_measure_type
                temp_dict['t_coefficient'] = j['children'][0].bias_value
                temp_dict['t_coefficient_initial'] = j['children'][0].bias_value_initial
                temp_dict['potency'] = 'Not yet'
                temp_dict['bias_reference'] = j['children'][0].bias_reference
                temp_dict['bias_value_initial'] = j['children'][0].bias_value_initial
                temp_dict['emax_reference_ligand'] = j['children'][0].emax_ligand_reference
                temp_dict['ligand_function'] = j['children'][0].ligand_function
                temp_dict['ligand'] = j['main'].ligand
                temp_dict['order_no'] = 0
                doubles.append(temp_dict)
                temp['assay'] = doubles
                send[increment] = temp
                increment = increment + 1

        return send

    def process_references(self, dictionary):
        '''
        1) Separate references from experiment data_all
        2) Assign references to experiments
        3) Saves refernces to DB (ReferenceLigand model)
        '''
        references = list()
        for item in dictionary.items():
            if item[1]['assay'][0]['bias_reference'].lower() != '':
                references.append(item[1])
                #del item

        for item in dictionary.items():
            for i in references:
                if (i['publication'] == item[1]['publication'] and
                        i['receptor'] == item[1]['receptor']):
                    if (item[1]['assay'][0]['signalling_protein'] == i['assay'][0]['signalling_protein'] and
                        item[1]['assay'][0]['assay_type'] == i['assay'][0]['assay_type'] and
                        item[1]['assay'][0]['cell_line'] == i['assay'][0]['cell_line'] and
                            item[1]['assay'][0]['assay_measure_method'] == i['assay'][0]['assay_measure_method']):
                        item[1]['reference'] = i['assay']
                    else:
                        continue

    def process_dublicates(self, dictionary):
        '''
        Recieve data from "process_data"
        search for objects with same publication ligand receptor
        Create new object of biased_experiment and all children
        '''
        doubles = []
        result = []
        context = {}
        send = {}
        for j in dictionary.items():
            name = str(j[1]['publication']) + \
                '/' + str(j[1]['ligand']) + '/' + str(j[1]['receptor'])
            temp_obj = list()
            reference_obj = list()
            if(name in context):
                temp_obj = context[name]['sub']
                reference_obj = context[name]['reference']

                for i in j[1]['assay']:
                    temp_obj = temp_obj + j[1]['assay']
                if 'reference' not  in j[1]:
                    continue
                else:
                    for i in j[1]['reference']:
                        reference_obj = reference_obj + j[1]['reference']
            else:
                for entry in j[1]['assay']:
                    temp_obj.append(entry)
                if 'reference' not  in j[1]:
                    continue
                else:
                    for entry in j[1]['reference']:
                        reference_obj.append(entry)

            context[name] = j[1]
            context[name].pop('assay')
            context[name]['sub'] = temp_obj
            context[name]['reference'] = reference_obj
            send[name] = self.parse_children(context[name])
            send[name]['reference'] = reference_obj

        self.group_family(send)
        for i in send.items():
            test = dict()
            i[1].pop('assay')
            # family sorting
            try:
                test = sorted(i[1]['family'].items(), key=lambda x: x[1]['quantitive_activity']
                              if x[1]['quantitive_activity'] else 0,  reverse=False)
            except Exception as msg:
                print('error---', msg, i[1],'\n')
                continue
            for x in enumerate(test):
                x[1][1]['order_no'] = x[0]

            i[1]['biasdata'] = test

            self.verify_calc(i[1]['biasdata'], i[1]['reference'])
            self.calc_potency(i[1]['biasdata'])
            self.calc_bias_factor(i[1]['biasdata'], i[1]['reference'])

            for j in i[1]['family'].items():
                if(j[1]['quantitive_activity'] is not None and
                   j[1]['quantitive_efficacy'] is not None):
                    j[1]['quantitive_activity'] = round(
                    j[1]['quantitive_activity'], 1)
                    j[1]['quantitive_efficacy'] = int(
                        j[1]['quantitive_efficacy'])
        # self.assay_five(send)
        return send

    def verify_calc(self, biasdata, reference):
        '''
        check for  transduction t_coefficient
        if exists, then calcualte Relative transduction t_coefficient
        '''
        for item in enumerate(biasdata):

            if (item[1][1]['t_coefficient_initial'] is not None and
                    item[1][1]['t_coefficient'] is None):
                for j in reference:
                    if (item[1][1]['assay_type'] == j['assay_type'] and
                        item[1][1]['cell_line'] == j['cell_line'] and
                            item[1][1]['signalling_protein'] == j['signalling_protein']):
                        item[1][1]['t_coefficient'] = round(math.log10(
                            item[1][1]['t_coefficient_initial']) - math.log10(j['t_coefficient_initial']),1)
                    if j['emax_reference_ligand'] is None:
                        j['emax_reference_ligand'] = j['ligand']

    def assay_five(self, send):
        '''
        Used to fill empty template (5 in total, 5 families) spaces
        to hide columns im template
        '''
        for i in send.items():
            temp_dict = dict()
            bias_dict = 0

            if len(i[1]['family']) < 6:
                length = 5 - len(i[1]['family'])
                temp_dict['pathway'] = ''
                temp_dict['bias'] = ''
                temp_dict['cell_line'] = ''
                temp_dict['assay_type'] = ''
                temp_dict['log_bias_factor'] = ''
                temp_dict['t_factor'] = ''
                temp_dict['ligand_function'] = ''

                for j in range(length, 6):
                    test = ('No_data', temp_dict)
                    i[1]['family'][str(length)] = test
                    length += 1

    def parse_children(self, context):

        temp = dict()
        test = list()
        assay_list = list()

        for k in context.items():

            if k[0] == 'ligand':
                temp['ligand'] = k[1]
            elif k[0] == 'receptor':
                temp['receptor'] = k[1]
            elif k[0] == 'publication':
                temp['publication'] = k[1]

            elif k[0] == 'sub':
                counter = 0
                for i in k[1]:

                    temp_dict = dict()
                    # add logic of  families

                    if (i['signalling_protein'] == 'β-arrestin' or
                        i['signalling_protein'] == 'β-arrestin-1 (non-visual arrestin-2)' or
                            i['signalling_protein'] == 'β-arrestin-2 (non-visual arrestin-3)'):
                        temp_dict['pathway'] = 'β-arrestin'

                    elif (i['signalling_protein'] == 'gi/o-family' or
                          i['signalling_protein'] == 'gαi1' or
                          i['signalling_protein'] == 'gαi2' or
                          i['signalling_protein'] == 'gαi3' or
                          i['signalling_protein'] == 'gαo' or
                          i['signalling_protein'] == 'gαoA' or
                          i['signalling_protein'] == 'gαoB'):
                        temp_dict['pathway'] = 'Gi/o-family'

                    elif (i['signalling_protein'] == 'gq-family' or
                            i['signalling_protein'] == 'gαq' or
                            i['signalling_protein'] == 'gαq11' or
                            i['signalling_protein'] == 'gαq14' or
                            i['signalling_protein'] == 'gαq14' or
                            i['signalling_protein'] == 'gαq16' or
                            i['signalling_protein'] == 'gαq14 (gαq16)'):
                        temp_dict['pathway'] = 'Gq-family'

                    elif (i['signalling_protein'] == 'g12/13-family' or
                          i['signalling_protein'] == 'gα12' or
                          i['signalling_protein'] == 'gα13'):
                        temp_dict['pathway'] = 'G12/13-family'

                    elif (i['signalling_protein'] == 'gs-family' or
                          i['signalling_protein'] == 'gαs' or
                          i['signalling_protein'] == 'gαolf'):
                        temp_dict['pathway'] = 'Gs-family'
                    else:
                        temp_dict['pathway'] = 'No data'

                    temp_dict['signalling_protein'] = i['signalling_protein']
                    temp_dict['cell_line'] = i['cell_line']
                    temp_dict['assay_type'] = i['assay_type']
                    temp_dict['assay_measure_method'] = i['assay_measure_method']
                    temp_dict['assay_time_resolved'] = i['assay_time_resolved']
                    temp_dict['quantitive_activity'] = i['quantitive_activity']
                    temp_dict['qualitative_activity'] = i['qualitative_activity']
                    temp_dict['quantitive_unit'] = i['quantitive_unit']
                    temp_dict['quantitive_efficacy'] = i['quantitive_efficacy']
                    temp_dict['efficacy_unit'] = i['efficacy_unit']
                    temp_dict['quantitive_measure_type'] = i['quantitive_measure_type']
                    temp_dict['efficacy_measure_type'] = i['efficacy_measure_type']
                    temp_dict['t_coefficient'] = i['t_coefficient']
                    temp_dict['t_coefficient_initial'] = i['t_coefficient_initial']
                    temp_dict['potency'] = ''
                    temp_dict['t_factor'] = ''
                    temp_dict['log_bias_factor'] = ''
                    temp_dict['bias_reference'] = i['bias_reference']
                    temp_dict['bias_value_initial'] = i['bias_value_initial']
                    temp_dict['emax_reference_ligand'] = i['emax_reference_ligand']
                    temp_dict['ligand_function'] = i['ligand_function']

                    temp_dict['order_no'] = 0
                    assay_list.append(temp_dict)
                temp['assay'] = assay_list

        context = temp
        return context

    def group_family(self, send):
        '''
        reformat data by families
        to correctly display it
        '''
        for item in send.items():

            context = dict()
            for i in item[1]['assay']:

                if i['pathway'] == 'β-arrestin':
                    context['br'] = i
                    #item[1]['family']= context
                elif i['pathway'] == 'Gi/o-family':
                    context['gi'] = i
                    #item[1]['gi']= i
                elif i['pathway'] == 'Gq-family':
                    context['gq'] = i
                    #item[1]['gq']= i
                elif i['pathway'] == 'G12/13-family':
                    context['g12'] = i
                    #item[1]['g12']= i
                elif i['pathway'] == 'Gs-family':
                    context['gs'] = i
                    #item[1]['gs']= i
                else:
                    continue
                item[1]['family'] = context

    def calc_potency(self, biasdata):
        most_potent = dict()
        for i in enumerate(biasdata):
            if i[0] == 0:
                most_potent = i[1][1]
            else:
                if i[1][1]['quantitive_activity'] is not None and most_potent['quantitive_activity'] is not None:
                    i[1][1]['potency'] = round(i[1][1]['quantitive_activity'] / most_potent['quantitive_activity'], 1)

                if i[1][1]['t_coefficient'] is not None and most_potent['t_coefficient'] is not None:
                    i[1][1]['t_factor'] = round(math.pow(10, (most_potent['t_coefficient'] - i[1][1]['t_coefficient'])), 1)
                else:
                    i[1][1]['t_factor'] = None

    def calc_bias_factor(self, biasdata, reference):
        most_reference = dict()
        most_potent = dict()
        for i in enumerate(biasdata):
            temp_reference = dict()
            for j in reference:
                if i[0] == 0:
                    most_potent = i[1][1]
                    if (i[1][1]['assay_type'] == j['assay_type'] and
                        i[1][1]['cell_line'] == j['cell_line'] and
                            i[1][1]['signalling_protein'] == j['signalling_protein']):
                        most_reference = j
                else:
                    if (i[1][1]['quantitive_activity'] is not None and most_potent['quantitive_activity'] is not None and
                        i[1][1]['quantitive_efficacy'] is not None and most_potent['quantitive_efficacy'] is not None and
                            j['quantitive_efficacy'] is not None and j['quantitive_efficacy'] is not None):
                        if (i[1][1]['assay_type'] == j['assay_type'] and
                            i[1][1]['cell_line'] == j['cell_line'] and
                                i[1][1]['signalling_protein'] == j['signalling_protein'] and
                                most_reference['quantitive_activity'] is not None):
                            temp_reference = j
                            #print('---experiment---', type(most_potent['quantitive_efficacy']), type( most_potent['quantitive_activity']),type(most_reference['quantitive_efficacy']),type(most_reference['quantitive_activity']),type(temp_reference['quantitive_efficacy']),type(temp_reference['quantitive_activity']))
                            i[1][1]['log_bias_factor'] = round((most_potent['quantitive_efficacy'] / most_potent['quantitive_activity']) - (most_reference['quantitive_efficacy'] / most_reference['quantitive_activity']) - (
                                (i[1][1]['quantitive_efficacy'] / i[1][1]['quantitive_activity']) - (temp_reference['quantitive_efficacy'] / temp_reference['quantitive_activity'])), 1)


                        else:
                            i[1][1]['log_bias_factor'] = ''

    def bias_list(self):
        print('i am in')
        context = {}
        content = BiasedExperiment.objects.all().prefetch_related(
            'experiment_data', 'ligand', 'receptor', 'publication', 'publication__web_link', 'experiment_data__emax_ligand_reference')
        # merge children
        pre_data = self.process_data(content)
        # transform to combined dictionary
        combined = self.change(pre_data)
        # combine references
        self.process_references(combined)
        context.update({'data': self.process_dublicates(combined)})
        for i in context['data'].items():
            try:
                i[1].pop('reference')
                i[1].pop('biasdata')
            except Exception as msg:
                print('\n---saving error---',msg)
                continue
            if self.fetch_experiment(i[1]['publication'], i[1]['ligand'], i[1]['receptor']) == False:
                experiment_entry = AnalyzedExperiment(publication=i[1]['publication'],
                                                      ligand=i[1]['ligand'],
                                                      receptor=i[1]['receptor']
                                                      )
                experiment_entry.save()
                for ex in i[1]['family'].items():
                    emax_ligand = ex[1]['emax_reference_ligand']
                    experiment_assay = AnalyzedAssay(experiment=experiment_entry,
                                                     family=ex[0],
                                                     order_no=ex[1]['order_no'],
                                                     signalling_protein=ex[1]['signalling_protein'],
                                                     cell_line=ex[1]['cell_line'],
                                                     assay_type=ex[1]['assay_type'],
                                                     assay_measure=ex[1]['assay_measure_method'],
                                                     assay_time_resolved=ex[1]['assay_time_resolved'],
                                                     ligand_function=ex[1]['ligand_function'],
                                                     quantitive_measure_type=ex[1]['quantitive_measure_type'],
                                                     quantitive_activity=ex[1]['quantitive_activity'],
                                                     # quantitive_activity_sign=ex[1]['quantitive_activity_sign'],
                                                     quantitive_unit=ex[1]['quantitive_unit'],
                                                     qualitative_activity=ex[1]['qualitative_activity'],
                                                     quantitive_efficacy=ex[1]['quantitive_efficacy'],
                                                     efficacy_measure_type=ex[1]['efficacy_measure_type'],
                                                     # efficacy_sign=ex[1]['efficacy_sign'],
                                                     efficacy_unit=ex[1]['efficacy_unit'],
                                                     potency=ex[1]['potency'],
                                                     t_coefficient=ex[1]['t_coefficient'],
                                                     t_value=ex[1]['t_coefficient'],
                                                     t_factor=ex[1]['t_factor'],
                                                     log_bias_factor=ex[1]['log_bias_factor'],
                                                     emax_ligand_reference=emax_ligand
                                                     )
                    experiment_assay.save()
                    print('saved')
            else:
                print("already defined")
