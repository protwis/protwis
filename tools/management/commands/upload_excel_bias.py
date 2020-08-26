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
from protein.models import Protein
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
            self.prepare_all_data(options['filename'])
            self.logger.info('COMPLETED CREATING BIAS')
        except Exception as msg:
            print('--error--', msg, '\n')
            self.logger.info("The error appeared in def handle")

    def purge_bias_data(self):
        delete_bias_excel = BiasedExperiment.objects.all()
        delete_bias_excel.delete()
        delete_bias_experiment = AnalyzedExperiment.objects.all()
        delete_bias_experiment.delete()

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
        for i, r in enumerate(rows, 1):
            # code to skip rows in excel for faster testing
            # if i < 15:
            #     continue
            # if i > 1200:
            #     break
            if i % 100 == 0:
                print(i)
            d = dict()

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

            d['cell_line'] = r[12]
            d['protein'] = r[13].strip().replace('α','a').replace('β','B').replace('g','G').lower()
            d['protein_assay'] = r[14].strip()
            d['protein_assay_method'] = r[15]
            d['protein_time_resolved'] = r[16]

            d['protein_ligand_function'] = r[17]
            d['protein_mtype'] = r[18]
            d['protein_activity_equation'] = r[19]
            d['protein_activity_quantity'] = r[20]
            d['protein_activity_quantity_unit'] = r[21]
            d['protein_activity_quality'] = r[22]
            d['protein_efficacy_measure'] = r[23]
            d['protein_efficacy_equation'] = r[24]
            d['protein_efficacy_quantity'] = r[25]
            d['protein_efficacy_quantity_unit'] = r[26]
            d['pathway_bias_initial'] = r[27]
            d['pathway_bias'] = r[28]

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
            elif d['protein_efficacy_quantity'] !=None:
                d['protein_efficacy_quantity'] = round(d['protein_efficacy_quantity'],0)
            if not isinstance(d['pathway_bias'], float):
                d['pathway_bias'] = None
            if not isinstance(d['pathway_bias_initial'], float):
                d['pathway_bias_initial'] = None

            #fetch bias measurements
            if not isinstance(d['protein_activity_quantity'], (int, float)):
                d['protein_activity_quantity'] = None
            else:
                d['protein_activity_quantity'],d['protein_mtype'] = self.fetch_measurements(d['protein_activity_quantity'],
																	     d['protein_mtype'],
																	     d['protein_activity_quantity_unit'])

            if not isinstance(d['pathway_bias_initial'], (int, float)):
                d['pathway_bias_initial'] = None
            if not isinstance(d['pathway_bias'], (int, float)):
                bias_value=d['pathway_bias'] = None


            with concurrent.futures.ThreadPoolExecutor() as executor:
                family_future = executor.submit(self.define_g_family,d['protein'], d['protein_assay'])
                family = family_future.result()
                pub_future = executor.submit(self.fetch_publication,d['reference'])
                pub = pub_future.result()


            # fetch main ligand
            l = self.fetch_ligand(
                d['ligand_id'], d['ligand_type'], d['ligand_name'], d['source_file'])
            if not l:
                continue
            #fetch endogenous ligand
            protein = self.fetch_protein(d['receptor'], d['source_file'])

            # fetch reference_ligand
            reference_ligand = self.fetch_ligand(
                d['bias_ligand_id'], d['bias_ligand_type'], d['bias_ligand_name'], d['source_file'])


            #fetch ChEMBL
            chembl = None
            chembl = self.fetch_chembl(l)

            # fetch protein
            protein = self.fetch_protein(d['receptor'], d['source_file'])
            if protein == None:
                continue
            end_ligand  = self.fetch_endogenous(protein)


## TODO:  check if it was already uploaded
            experiment_entry = BiasedExperiment(submission_author=d['submitting_group'],
                                                publication=pub,
                                                ligand=l,
                                                receptor=protein,
                                                chembl = chembl,
                                                endogenous_ligand = end_ligand

                                                )
            experiment_entry.save()
            self.fetch_vendor(l,experiment_entry)

            experiment_assay = ExperimentAssay(biased_experiment=experiment_entry,
                                               signalling_protein=d['protein'],
                                               family = family,
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
            #fetch authors
            self.fetch_publication_authors(pub,experiment_assay)


            temp.append(d)
        return temp


    def fetch_publication_authors(self,publication, experiment_assay):
        counter = 0

        author_list = list()
        if publication.authors != None:

            for authors in publication.authors.split(','):
                author_list.append(authors.strip())

            author_list.reverse()
            for i in author_list:
                if counter < 3:
                    assay_author = ExperimentAssayAuthors(experiment = experiment_assay,
                    author=i)
                    assay_author.save()
                    counter=counter+1
            # assay_author = ExperimentAssayAuthors(experiment = experiment_assay,



    def fetch_measurements(self, potency, p_type, unit):
        if p_type.lower()  == 'pec50':
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

        if p_type.lower()  == 'ec50':
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
        if p_type.lower()  == 'ic50':
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
        # if potency:
        #     potency = (-1)*math.log10(potency)
        #     p_type = 'pec50'
        #     # potency = "{:.2E}".format(Decimal(potency))
        return potency,p_type

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
                protein==' gaq' or
                protein=='gaq/11' or
                protein=='gaq/14' or
                protein=='gaq/15' or
                protein=='gaq/16'):
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
            if assay_type == 'pERK1/2 activation' or assay_type =="pERK1-2":
                family = 'pERK1-2'
        else:
            family == protein

        return family

    def fetch_endogenous(self, protein):
        try:
            with connection.cursor() as cursor:
                cursor.execute("SELECT * FROM protein_endogenous_ligands WHERE protein_id =%s", [protein.pk])
                row = cursor.fetchone()
                end_ligand = Ligand.objects.filter(id=row[2])
                test = end_ligand.get()

            return test
        except:
            return None

    def fetch_vendor(self, ligand,experiment_entry):
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

    def fetch_chembl(self,ligand):
        temp = ligand
        chembl_id = None
        links = temp.properities.web_links.all()

        for x in links:
            if x.web_resource.slug=='chembl_ligand':
                chembl_id = [x for x in links if x.web_resource.slug=='chembl_ligand'][0].index
        return chembl_id

    def fetch_protein(self,protein_from_excel, source):
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

    def prepare_all_data(self, filenames):
        if not filenames:
            filenames = os.listdir(self.structure_data_dir)
        for source_file in filenames:

            source_file_path = os.sep.join(
                [self.structure_data_dir, source_file]).replace('//', '/')

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
