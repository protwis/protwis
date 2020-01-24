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
from ligand.models import *
from mutation.models import Mutation
from ligand.functions import get_or_make_ligand
from common.models import WebLink, WebResource, Publication
from django.db import connection
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
    structure_data_dir = 'excel/pathways/'
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
        #self.bias_list()
        try:
            print('CREATING BIAS PATHWAYS DATA')
            print(options['filename'])
            self.prepare_all_data(options['filename'])
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
            # if i < 58:
            #     continue
            # if i > 58:
            #     break
            # print(r)
            d = {}
            if r[4] != '':  # checks if the receptor field exists
                # try:
                res = ''
                mut = ''
                d['submitting_group'] = r[0]
                d['relevance'] = r[15]
                # doi
                d['reference'] = r[1]
                # protein
                d['receptor'] = r[4].lower()
                d['signalling_protein'] = r[5].lower().strip()
                # ligand
                d['ligand_name'] = r[6]
                d['ligand_type'] = r[7]
                d['ligand_id'] = r[8]
                # pathway
                d['pathway_outcome'] = r[14]
                d['pathway_summary'] = r[13]
                d['pathway_detail'] = r[12]
                #experiment
                d['experiment_disctinction'] = r[9]
                d['experiment_system'] = r[10]
                d['experiment_method'] = r[11]

                d['source_file'] = source_file + str(i)

                if not isinstance(d['ligand_id'], str):
                    d['ligand_id'] = int(d['ligand_id'])

                #define G family
                family = self.define_g_family(d['signalling_protein'])

                # fetch publicaition
                pub = self.fetch_publication(d['reference'])

                # fetch main ligand
                l = self.fetch_ligand(
                    d['ligand_id'], d['ligand_type'], d['ligand_name'], d['source_file'])

                #fetch ChEMBL
                chembl = None
                chembl = self.fetch_chembl(l)

                # fetch protein
                protein = self.fetch_protein(d['receptor'], d['source_file'])
                if protein == None:
                    print('уккщк')
                    continue

## TODO:  check if it was already uploaded
                experiment_entry = BiasedPathways(submission_author=d['submitting_group'],
                                                    publication=pub,
                                                    ligand=l,
                                                    receptor=protein,
                                                    chembl = chembl,
                                                    relevance = d['relevance'],
                                                    signalling_protein = d['signalling_protein']

                                                    )
                experiment_entry.save()

                experiment_assay = BiasedPathwaysAssay(biased_pathway=experiment_entry,
                                                  pathway_outcome_high = d['pathway_outcome'],
                                                  pathway_outcome_summary = d['pathway_summary'],
                                                  pathway_outcome_detail  = d['pathway_detail'],
                                                  experiment_pathway_distinction = d['experiment_disctinction'],
                                                  experiment_system = d['experiment_system'],
                                                  experiment_outcome_method= d['experiment_method']
                                                   )

                experiment_assay.save()

                # except Exception as msg:
                #     print(d['source_file'], msg)
                #     continue

            temp.append(d)
        return temp

    def define_g_family(self, protein):
        if (protein == 'β-arrestin' or
            protein == 'β-arrestin-1 (non-visual arrestin-2)' or
                protein == 'β-arrestin-2 (non-visual arrestin-3)'):
            family = 'B-arr'

        elif (protein == 'gi/o-family' or
              protein == 'gαi1' or
              protein == 'gαi2' or
              protein == 'gαi3' or
              protein == 'gαo' or
              protein == 'gαoA' or
              protein == 'gαoB'):
            family = 'Gi/o'

        elif (protein == 'gq-family' or
                protein == 'gαq' or
                protein == 'gαq11' or
                protein == 'gαq14' or
                protein == 'gαq14' or
                protein == 'gαq16' or
                protein == 'gαq14 (gαq16)'):
            family = 'Gq/11'

        elif (protein == 'g12/13-family' or
              protein == 'gα12' or
              protein == 'gα13'):
            family = 'G12/13'

        elif (protein == 'gs-family' or
              protein == 'gαs' or
              protein == 'gαolf'):
            family = 'Gs'
        else:
            family = 'No data'

        return family

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

        return test

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
