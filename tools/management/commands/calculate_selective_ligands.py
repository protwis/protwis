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
from django.db.models import Q, Count

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
        # self.bias_list()
        try:
            print('CREATING BIAS PATHWAYS DATA')
            print(options['filename'])
            self.calculate_selectivity()
            self.logger.info('COMPLETED CREATING BIAS')

        except Exception as msg:
            print('--error--', msg, '\n')
            self.logger.info("The error appeared in def handle")

    def calculate_selectivity(self):
        # TODO: Discuss with David/Albert : how to chose selectivity/ compare among what?
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
                print('no data for ligand',ligand)
            try:
                self.process_data(assay_f,"F")
            except:
                print("f selectivity errir")
                continue
            try:
                self.process_data(assay_b,"B")
            except:
                print("b selectivity errir")
                continue

        end = time.time()
        print('---temp_increment time---', end - start)

    def process_data(self, assay,type):
        assay_list = list()
        assay_list = self.process_assays(assay)
        try:
            sorted_assay_list = self.sort_assay(assay_list)
        except:
            print('sorted error')

        # compare assays by standard value, leave only ones with 1p fold selectivity
        final_assay = self.analyze_assay(sorted_assay_list)
        # if final_assay, then save it to db
        if final_assay:
            self.save_data(final_assay, type)

    def get_ligands(self):
        #Getting ligands from the model
        try:
            content = AssayExperiment.objects.all().order_by(
                'ligand').distinct('ligand').only('ligand')
        except AssayExperiment.DoesNotExist:
            content = None
        return content

    def get_data(self, ligand_name):
        #Getting data from the model for a ligand\n##limiting only by EC50 | IC50 (values)'
        try:

            content = AssayExperiment.objects.filter(ligand=ligand_name
                                                     ).filter(Q(assay_type='F') | Q(assay_type='B')
                                                              ).filter(Q(standard_type='IC50') |
                                                                       Q(standard_type='EC50') |
                                                                       Q(standard_type='Potency')).prefetch_related(
                'ligand', 'protein'
            ).only('protein', 'ligand', 'standard_type', 'pchembl_value', 'assay_type'
                   ).order_by('ligand')
        except AssayExperiment.DoesNotExist:
            content = None
        return content

    def process_assays(self, assays):
        processed_data = list()
        for i in assays:
            try:
                assay_data = dict()
                assay_data["protein"] = i.protein
                assay_data["ligand"] = i.ligand
                assay_data["assay_type"] = i.assay_type
                assay_data["standard_type"] = i.standard_type
                assay_data["pchembl_value"] = float(i.standard_value)
                processed_data.append(assay_data)
            except:
                print('process data', i)
                break

        return processed_data

    def sort_assay(self, assays):
        return sorted(assays, key=lambda i: i['pchembl_value'], reverse=True)

    def analyze_assay(self, assays):
        # select most potent if more than 10 folds
        most_potent = next(iter(assays or []), None)
        if most_potent and most_potent != None:
            for i in assays:
                try:
                    if most_potent['pchembl_value'] > (i['pchembl_value'] + 1):
                        return most_potent
                        break
                except:
                    print('analyze rows',i)

        # for assay in assays:
        #     if most_potent['pchembl_value']

    def save_data(self, final_assay,type):
        #saving assay ---', final_assay
        save_assay = LigandReceptorStatistics(
            ligand=final_assay['ligand'],
            protein=final_assay['protein'],
            type=type,
            potency=final_assay['pchembl_value']
        )
        save_assay.save()
