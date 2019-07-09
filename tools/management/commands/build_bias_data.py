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
            # self.prepare_all_data(options['filename'])
            # self.bias_list()
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
            if i < 0:
                continue
            if i > 1000:
                break
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

    def fetch_experiment(self, publication, ligand, receptor, residue, mutation):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        try:
            experiment = AnalyzedExperiment.objects.filter(
                publication=publication, ligand=ligand, receptor=receptor, residue=residue, mutation=mutation)
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
            temp['residue'] = j['main'].residue
            temp['mutation'] = j['main'].mutation

            if not j['children']:
                continue
            else:
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
            if (item[1]['assay'][0]['bias_reference'].lower() != '' and
                    item[1]['assay'][0]['bias_reference'] not in references):
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
                '/' + str(j[1]['ligand']) + '/' + str(j[1]['receptor']
                                                      ) + str(j[1]['residue']) + str(j[1]['mutation'])
            temp_obj = list()
            reference_obj = list()
            if(name in context):

                temp_obj = context[name]['sub']
                reference_obj = context[name]['reference']

                for i in j[1]['assay']:
                    temp_obj = temp_obj + j[1]['assay']
                    temp_obj = [i for n, i in enumerate(
                        temp_obj) if i not in temp_obj[n + 1:]]

                if 'reference' not in j[1]:
                    continue
                else:
                    for i in j[1]['reference']:
                        reference_obj = reference_obj + j[1]['reference']
                        reference_obj = [i for n, i in enumerate(
                            reference_obj) if i not in reference_obj[n + 1:]]

            else:
                for entry in j[1]['assay']:
                    temp_obj.append(entry)
                if 'reference' not in j[1]:
                    continue
                else:
                    for entry in j[1]['reference']:
                        reference_obj.append(entry)
                        reference_obj = [i for n, i in enumerate(
                            reference_obj) if i not in reference_obj[n + 1:]]

            context[name] = j[1]

            context[name].pop('assay')
            context[name]['sub'] = temp_obj
            context[name]['reference'] = reference_obj
            send[name] = self.parse_children(context[name])
            send[name]['reference'] = reference_obj

        for i in send.items():
            test = dict()
            # try:
            i[1]['assay'][:] = [d for d in i[1]['assay']
                                if d.get('bias_reference').lower() != 'yes']
            # TODO: Change 9999999 to normal method that skips None value
            test = sorted(i[1]['assay'], key=lambda x: x['quantitive_activity']
                          if x['quantitive_activity'] else 999999,  reverse=False)

            # except Exception as msg:
            #     print('error---', msg, i[1],'\n')
            #     continue
            for x in enumerate(test):
                x[1]['order_no'] = x[0]

            i[1]['biasdata'] = test
            i[1].pop('assay')

            # TODO: CHECK VARIABLES FOR CALCULATIONS
            # TODO: ASSIGN QUALITATIVE VALUE FOR IC50 OR PIC50
            # TODO: Change COLOR DEPENDING IC50
            self.calc_potency(i[1]['biasdata'])
            self.verify_calc(i[1]['biasdata'], i[1]['reference'])
            self.calc_bias_factor(i[1]['biasdata'], i[1]['reference'])

            for j in i[1]['biasdata']:
                if(j['quantitive_activity'] is not None and
                   j['quantitive_efficacy'] is not None):
                    j['quantitive_activity'] = round(
                        j['quantitive_activity'], 1)
                    j['quantitive_efficacy'] = int(
                        j['quantitive_efficacy'])
        # self.assay_five(send)
        return send

    def verify_calc(self, biasdata, reference):
        '''
        check for  transduction t_coefficient_initial both in assay and reference
        if exists, then calcualte Relative transduction t_coefficient
        '''
        for item in enumerate(biasdata):

            if not isinstance(item[1]['t_coefficient_initial'], (int, float)):
                item[1]['t_coefficient_initial'] = None
            if not isinstance(item[1]['quantitive_activity'], (int, float)):
                item[1]['quantitive_activity'] = None
        # TODO: verify calc with Kasper or Sahar
            else:
                if item[1]['quantitive_measure_type'].lower() == 'ec50':
                    continue
                elif item[1]['quantitive_measure_type'].lower() == 'pec50':
                    item[1]['quantitive_activity'] = 10**(
                        item[1]['quantitive_activity'] * (-1))
                    item[1]['quantitive_measure_type'] = 'EC50'
                elif item[1]['quantitive_measure_type'].lower() == 'logec50':
                    item[1]['quantitive_activity'] = 10**(
                        item[1]['quantitive_activity'])
                    item[1]['quantitive_measure_type'] = 'EC50'
                elif item[1]['quantitive_measure_type'].lower() == 'pic50':
                    item[1]['quantitive_activity'] = 10**(
                        item[1]['quantitive_activity'] * (-1))
                    item[1]['quantitive_measure_type'] = 'IC50'
                else:
                    item[1]['quantitive_activity'] = None

                    # TODO: continue ec50 types
            if not isinstance(item[1]['quantitive_efficacy'], (int, float)):
                item[1]['quantitive_efficacy'] = None
            if not isinstance(item[1]['t_coefficient'], (int, float)):
                item[1]['t_coefficient'] = None

            for j in reference:
                if not isinstance(j['t_coefficient_initial'], (int, float)):
                    j['t_coefficient_initial'] = None
                if not isinstance(j['quantitive_activity'], (int, float)):
                    j['quantitive_activity'] = None
                if not isinstance(j['quantitive_efficacy'], (int, float)):
                    j['quantitive_efficacy'] = None
                if not isinstance(j['t_coefficient'], (int, float)):
                    j['t_coefficient'] = None
            for j in reference:
                if (j['assay_type'] == item[1]['assay_type'] and
                        j['cell_line'] == item[1]['cell_line'] and
                        j['signalling_protein'] == item[1]['signalling_protein'] and
                        item[1]['t_coefficient'] == None and
                        isinstance(item[1]['t_coefficient_initial'], (int, float)) and
                        isinstance(j['t_coefficient_initial'], (int, float))):
                    item[1]['t_coefficient'] = round(math.log10(
                        item[1]['t_coefficient_initial']) - math.log10(j['t_coefficient_initial']), 1)
                else:
                    continue

    def calc_potency(self, biasdata):
        most_potent = dict()
        for i in biasdata:

            if i['order_no'] == 0:
                most_potent = i
            else:
                # TODO: change color, differentiate according pc50 ic50 logpc50 logic50
                if i['quantitive_measure_type'].lower() == 'ec50':
                    if i['quantitive_activity'] is not None and most_potent['quantitive_activity'] is not None:
                        i['potency'] = round(
                            i['quantitive_activity'] / most_potent['quantitive_activity'], 1)
                    elif i['quantitive_measure_type'].lower() == 'ic50':
                        i['potency'] = round(
                            i['quantitive_activity'] - most_potent['quantitive_activity'], 1)
                    else:
                        i['potency'] = None

                if i['t_coefficient'] is not None and most_potent['t_coefficient'] is not None:
                    # TODO: validate if difference is non negative
                    i['t_factor'] = round(
                        math.pow(10, (most_potent['t_coefficient'] - i['t_coefficient'])), 1)
                else:
                    i['t_factor'] = None

    def calc_bias_factor(self, biasdata, reference):
        most_reference = dict()
        most_potent = dict()
        for i in biasdata:

            temp_reference = dict()
            if i['order_no'] == 0:
                most_potent = i
                for j in reference:
                    if (i['assay_type'] == j['assay_type'] and
                        i['cell_line'] == j['cell_line'] and
                            i['signalling_protein'] == j['signalling_protein']):
                        most_reference = j
            else:
                for j in reference:
                    if (i['assay_type'] == j['assay_type'] and
                        i['cell_line'] == j['cell_line'] and
                            i['signalling_protein'] == j['signalling_protein']):
                        temp_reference = j
                if (i['quantitive_activity'] is not None and most_potent['quantitive_activity'] is not None and
                    i['quantitive_efficacy'] is not None and most_potent['quantitive_efficacy'] is not None and
                    j['quantitive_activity'] is not None and j['quantitive_efficacy'] is not None and
                        temp_reference['quantitive_activity'] is not None and temp_reference['quantitive_efficacy'] is not None):
                    if (i['quantitive_measure_type'].lower() == 'ec50' and j['quantitive_measure_type'].lower() == 'ec50' and
                            most_potent['quantitive_measure_type'].lower() == 'ec50' and temp_reference['quantitive_measure_type'].lower() == 'ec50'):
                        temp_calculation = math.log10(most_potent['quantitive_efficacy'] / most_potent['quantitive_activity']) - math.log10(most_reference['quantitive_efficacy'] / most_reference['quantitive_activity']) - (
                            math.log10(i['quantitive_efficacy'] / i['quantitive_activity']) - math.log10(temp_reference['quantitive_efficacy'] / temp_reference['quantitive_activity']))
                        i['log_bias_factor'] = round(temp_calculation, 1)
                    elif (i['quantitive_measure_type'].lower() == 'ic50' and j['quantitive_measure_type'].lower() == 'ic50' and
                          most_potent['quantitive_measure_type'].lower() == 'ic50' and temp_reference['quantitive_measure_type'].lower() == 'ic50'):
                        i['log_bias_factor'] = 'Only agonism in main pathway'
                    else:
                        i['log_bias_factor'] = None
                else:
                    i['log_bias_factor'] = None

    def assay_five(self, send):
        '''
        Used to fill empty template (5 in total, 5 families) spaces
        to hide columns im template
        '''
        for i in send.items():
            temp_dict = dict()
            bias_dict = 0

            if len(i[1]['biasdata']) < 5:
                length = 5 - len(i[1]['biasdata'])
                temp_dict['pathway'] = ''
                temp_dict['bias'] = ''
                temp_dict['cell_line'] = ''
                temp_dict['assay_type'] = ''
                temp_dict['log_bias_factor'] = ''
                temp_dict['t_factor'] = ''
                temp_dict['ligand_function'] = ''

                for j in range(length, 6):
                    i[1]['biasdata'].append(temp_dict)
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
            elif k[0] == 'residue':
                temp['residue'] = k[1]
            elif k[0] == 'mutation':
                temp['mutation'] = k[1]

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

    def bias_list(self):
        print('i am in')
        context = {}
        content = BiasedExperiment.objects.all().prefetch_related(
            'experiment_data', 'ligand', 'receptor', 'publication', 'publication__web_link', 'experiment_data__emax_ligand_reference', 'mutation', 'mutation__residue')
        # merge children
        pre_data = self.process_data(content)
        # transform to combined dictionary
        combined = self.change(pre_data)
        # combine references
        self.process_references(combined)
        context.update({'data': self.process_dublicates(combined)})
        for i in context['data'].items():
            # try:
                # i[1].pop('reference')
                # i[1].pop('biasdata')
                # TODO: move by one tab when uncomment try catch
            if self.fetch_experiment(i[1]['publication'], i[1]['ligand'], i[1]['receptor'], i[1]['residue'], i[1]['mutation']) == False:
                experiment_entry = AnalyzedExperiment(publication=i[1]['publication'],
                                                      ligand=i[1]['ligand'],
                                                      receptor=i[1]['receptor'],
                                                      mutation=i[1]['mutation'],
                                                      residue=i[1]['residue']
                                                      )
                experiment_entry.save()
                for ex in i[1]['biasdata']:
                    print('--saving---', ex, '\n')
                    emax_ligand = ex['emax_reference_ligand']
                    experiment_assay = AnalyzedAssay(experiment=experiment_entry,
                                                     family=ex['pathway'],
                                                     order_no=ex['order_no'],
                                                     signalling_protein=ex['signalling_protein'],
                                                     cell_line=ex['cell_line'],
                                                     assay_type=ex['assay_type'],
                                                     assay_measure=ex['assay_measure_method'],
                                                     assay_time_resolved=ex['assay_time_resolved'],
                                                     ligand_function=ex['ligand_function'],
                                                     quantitive_measure_type=ex['quantitive_measure_type'],
                                                     quantitive_activity=ex['quantitive_activity'],
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
                    print('saved')
            else:
                print("already defined")
            # except Exception as msg:
            #     print('\n---saving error---', msg)
            #     continue
