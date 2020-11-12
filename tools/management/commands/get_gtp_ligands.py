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
from multiprocessing.pool import ThreadPool
from chembl_webresource_client.new_client import new_client
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
import requests
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

    help = 'Updates ChEMBL data and imports it'
    publication_cache = {}
    ligand_cache = {}
    data_all = []
    my_queue = queue.Queue()

    def storeInQueue(f):
        def wrapper(*args):
            my_queue.put(f(*args))
        return wrapper

    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
                            type=int,
                            action='store',
                            dest='proc',
                            default=3,
                            help='Number of processes to run')
        parser.add_argument('-f', '--filename',
                            action='append',
                            dest='filename',
                            help='Filename to import. Can be used multiple times')

        parser.add_argument('--test_run', action='store_true', help='Skip this during a test run',
                            default=False)

    def handle(self, *args, **options):
        if options['test_run']:
            print('Skipping in test run')
            return
        try:
            print('Updatind ligand data from GuideToPharma')
            self.analyse_rows()
            self.logger.info('COMPLETED updating GuideToPharma Data')
        except Exception as msg:
            print('--error--', msg, '\n')
            self.logger.info("The error appeared in def handle")

    def purge_bias_data(self):
        pass

    def get_endogenous(self, targets):
        for target in targets:
            protein = self.fetch_protein(target)
            response = requests.get(
                "https://www.guidetopharmacology.org/services/targets/" + str(target) + "/naturalLigands")
            data = response.json()
            for i in data:
                try:
                    ligand = self.fetch_ligand(
                        data[0]['ligandId'], data[0]['type'])
                    if ligand and protein:
                        protein.endogenous_ligands.add(ligand)
                        lig, created = Ligand.objects.update_or_create(
                            id=ligand.id,
                            defaults={'endogenous': True},
                        )
                    else:
                        pass
                except:
                    continue

    def save_ligand_copy(self, ligand):
        """
        fetch ligands with Ligand model
        requires: ligand id, ligand id type, ligand name
        requires: source_file name
        """
        lig = Ligand()
        l = lig.load_by_gtop_id(
            ligand['name'], ligand['entry'], ligand['type'])
        if l != None:
            self.ligand_cache[ligand['entry']] = l

    def get_ligands(self):
        ligands = list()
        response = requests.get(
            "https://www.guidetopharmacology.org/services/ligands")
        print('total ligands:', len(response.json()))
        for entry in response.json():
            try:
                temp = dict()
                temp['entry'] = entry['ligandId']
                temp['name'] = entry['name']
                temp['type'] = entry['type']
                self.save_ligand_copy(temp)
            except:
                pass

    def get_gpcrs(self):
        target_list = list()
        response = requests.get(
            "https://www.guidetopharmacology.org/services/targets/families")
        for entry in response.json():
            try:
                if entry['parentFamilyIds'] != None and entry['parentFamilyIds'][0] == 694:
                    target_list.extend(entry['targetIds'])
            except:
                pass
                print('error entry', entry)
        return target_list

    def get_ligand_assays(self, targets):
        assay_list = list()
        response = requests.get(
            "https://www.guidetopharmacology.org/services/interactions")
        for entry in response.json():
            try:
                if entry['targetId'] != None and entry['targetId'] in targets:
                    assay_list.append(entry)
            except:
                pass
        print('\ninteractions', assay_list[:2])
        # print(len(assay_list))
        return assay_list

    def process_ligand_assays(self, assays):
        for i in assays:
            temp_dict = dict()
            temp_dict['protein'] = self.fetch_protein(i['targetId'])
            temp_dict['ligand'] = self.fetch_ligand(i['ligandId'], i['type'])
            print(i['refIds'], 'ligand g: ', temp_dict['ligand'],
                  ' ligand fra', i['ligandId'],)
            if i['refIds']:
                temp_dict['doi'] = self.fetch_publication(i['refIds'])
            else:
                temp_dict['doi'] = None
            if temp_dict['protein'] == None:
                continue
            if temp_dict['ligand'] == None:
                continue
            temp_dict['standard_type'] = i['affinityParameter']
            temp_dict['standard_value'] = i['affinity']
            temp_dict['assay_description'] = i['ligandContext']
            if temp_dict['assay_description'] == None:
                temp_dict['assay_description'] = "No data available"
            self.upload_to_db(temp_dict)

    def analyse_rows(self):
        """
        Fetch data to models
        Saves to DB
        """
        print('---Starting---\n')
        start = time.time()
        print("\n#1 Get GPCR ids from target families")
        target_list = self.get_gpcrs()
        print("\n#2 Get Ligands")
        self.get_ligands()
        print("\n#3 Get Ligand assays", self.ligand_cache)
        assays = self.get_ligand_assays(target_list)
        print("\n#4 Process Ligand assays")
        self.process_ligand_assays(assays)
        print("\n#5 Get Endogeneous ligands")
        self.get_endogenous(target_list)
        print('\n\n---Finished---')
        # range should be set to number of total objects/20 (API batch size)
        # 555200 is the last id saved before session was aborted

    def upload_to_db(self, i):
        # saves data
        chembl_data = AssayExperiment(ligand=i["ligand"],
                                      publication=i["doi"],
                                      protein=i["protein"],
                                      standard_type=i["standard_type"],
                                      standard_value=i["standard_value"],
                                      assay_description=i["assay_description"],
                                      )
        chembl_data.save()

    def fetch_protein(self, target):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        try:
            test = None
            test = Protein.objects.filter(
                web_links__index=target, web_links__web_resource__slug='gtop').first()
            return test
        except:
            return None

    def fetch_ligand(self, ligand_id, type):
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
                l = Ligand.objects.filter(
                    properities__web_links__index=ligand_id).first()
                if l:
                    cid = l.properities.web_links.filter(
                        web_resource__slug='gtoplig').first()
                    if cid:
                        return l
                    else:
                        l = None
                else:
                    lig = Ligand()
                    l = lig.load_by_gtop_id('', ligand_id, ligand['type'])
        except Exception as msg:
            l = None
            # print('ligand_id---',l,'\n end')
        return l

    def fetch_publication(self, refs):
        """
        fetch publication with Publication model
        requires: publication doi or pmid
        """
        publication_doi = None
        reference_response = requests.get(
            "https://www.guidetopharmacology.org/services/refs/" + str(refs[0]))

        if reference_response.status_code == 200:
            ligand_data = reference_response.json()
            get_pubmed = ligand_data['pmid']
            get_doi = ligand_data['doi']
            if get_pubmed != None:
                pub_type = 'pubmed'
                publication_doi = get_pubmed
            elif get_doi != None:
                pub_type = 'doi'
                publication_doi = get_doi
            else:
                return None
        else:
            return None
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
