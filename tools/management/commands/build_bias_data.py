from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db import IntegrityError
from django.utils.text import slugify


from build.management.commands.base_build import Command as BaseBuild
from common.tools import fetch_from_cache, save_to_cache, fetch_from_web_api
from residue.models import Residue
from protein.models import Protein
from ligand.models import Ligand, LigandProperities, LigandRole, LigandType, BiasedExperiment, ExperimentAssay
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

        try:

            print('CREATING BIAS DATA')
            print(options['filename'])
            self.prepare_all_data(options['filename'])
            self.logger.info('COMPLETED CREATING BIAS')

        except Exception as msg:
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
            #code to skip rows in excel for faster testing
            # if i < 0:
            #     continue
            # if i > 100:
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
                                                       bias_reference = d['bias_reference'],
                                                       bias_value=d['pathway_bias'],
                                                       bias_value_initial=d['pathway_bias_initial'],
                                                       bias_ligand_reference=reference_ligand
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
            mut=Mutation.objects.filter(protein=protein,
                                          residue=res,
                                          amino_acid=amino_acid)
            mut=mut.get()

        except Exception as msg:
            # print("error message from mutation: " + source, msg)
            new_mut=Mutation(protein=protein,
                               residue=res,
                               amino_acid=amino_acid)
            new_mut.save()
            print("Mutation created")
            mut=self.fetch_mutation(protein, res, amino_acid, source)
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
            res=Residue.objects.filter(protein_conformation__protein=protein,
                                         sequence_number=int(sequence_number),
                                         amino_acid=amino_acid)
            res=res.get()
        except Exception as msg:
            print("error message from residue: ", msg)
            self.mylog.exception(
                "residue fetching error | module: fetch_residue. " + source)
            res=None
        return res

    def fetch_protein(self, protein_from_excel, source):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        try:
            protein=Protein.objects.filter(entry_name=protein_from_excel)
            protein=protein.get()
        except Exception as msg:
            protein=None
            self.mylog.exception(
                "Protein fetching error | module: fetch_protein. Row # is : ", source)
        return protein

    def fetch_ligand(self, ligand_id, ligand_type, ligand_name, source_file):
        """
        fetch ligands with Ligand model
        requires: ligand id, ligand id type, ligand name
        requires: source_file name
        """
        l=None
        if str(ligand_id) in self.ligand_cache:
            if ligand_id in self.ligand_cache[str(ligand_id)]:
                l=self.ligand_cache[str(ligand_id)][ligand_id]
        else:
            self.ligand_cache[str(ligand_id)]={}

        if not l:
            try:
                l=get_or_make_ligand(
                    ligand_id, ligand_type, str(ligand_name))
            except Exception as msg:
                l=None
                self.ligand_cache[str(ligand_name), ligand_id]=l
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
            publication_doi=str(int(publication_doi))
        except ValueError:
            pass

        if publication_doi.isdigit():  # assume pubmed
            pub_type='pubmed'
        else:  # assume doi
            pub_type='doi'
        if publication_doi not in self.publication_cache:
            try:
                wl=WebLink.objects.get(
                    index=publication_doi, web_resource__slug=pub_type)
            except WebLink.DoesNotExist:
                try:
                    wl=WebLink.objects.create(index=publication_doi,
                                                web_resource=WebResource.objects.get(slug=pub_type))
                except IntegrityError:
                    wl=WebLink.objects.get(
                        index=publication_doi, web_resource__slug=pub_type)

            try:
                pub=Publication.objects.get(web_link=wl)
            except Publication.DoesNotExist:
                pub=Publication()
                try:
                    pub.web_link=wl
                    pub.save()
                except IntegrityError:
                    pub=Publication.objects.get(web_link=wl)

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
            self.publication_cache[publication_doi]=pub
        else:
            pub=self.publication_cache[publication_doi]

        return pub

    def prepare_all_data(self, filenames):
        if not filenames:
            filenames=os.listdir(self.structure_data_dir)
        for source_file in filenames:
            # print("source_file " + str(source_file))
            source_file_path=os.sep.join(
                [self.structure_data_dir, source_file]).replace('//', '/')
            # print("source_file_path " + str(source_file_path))
            if os.path.isfile(source_file_path) and source_file[0] != '.':
                self.logger.info('Reading file {}'.format(source_file_path))
                print('Reading file {}'.format(source_file_path))
                # read the yaml file
                rows=[]
                if source_file[-4:] == 'xlsx' or source_file[-3:] == 'xls':
                    if "~$" in source_file:
                        # ignore open excel files
                        continue
                    rows=self.loaddatafromexcel(source_file_path)
                    rows=self.analyse_rows(rows, source_file)
                else:
                    self.mylog.debug('unknown format'.source_file)
                    continue

                self.data_all += rows
        print(len(self.data_all), " total data points")
        print("Finished")
        self.mylog.info('finished')
