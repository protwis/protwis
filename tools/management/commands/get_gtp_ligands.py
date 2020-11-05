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
        # self.analyse_rows()
        try:
            print('Updatind ChEMBL data')
            self.analyse_rows()
            self.logger.info('COMPLETED updating ChEMBL Data')
        except Exception as msg:
            print('--error--', msg, '\n')
            self.logger.info("The error appeared in def handle")

    def purge_bias_data(self):
        pass

    def get_ligands(self):
        ligands = list()
        response = requests.get("https://www.guidetopharmacology.org/services/ligands/")
        print('\nTotal ligands', len(response.json()))
        for entry in response.json():
            ligands.append(entry)
        print(ligands[:10])
        return ligands

    def get_gpcrs(self):
        target_list = list()
        response = requests.get("https://www.guidetopharmacology.org/services/targets/families")
        for entry in response.json():
            try:
                if entry['parentFamilyIds'] != None and entry['parentFamilyIds'][0] == 694 :
                    target_list.extend(entry['targetIds'])
            except:
                pass
                # print('error entry', entry)
        return target_list

    def get_ligand_assays(self, targets):
        assay_list = list()
        print("\n Targets", targets)
        response = requests.get("https://www.guidetopharmacology.org/services/interactions")
        print('\n\ntotal', len(response.json()))
        for entry in response.json():
            try:
                if entry['targetId'] != None and entry['targetId'] in targets:
                    assay_list.append(entry)
            except:
                pass
        # print(assay_list[:10])
        # print(len(assay_list))
        return assay_list


    def get_dois(self, dci, q):
        # gets references for assays from ChEMBL (DOI)
        pubs = new_client.document.filter(document_chembl_id=dci).only('doi')
        if len(pubs) > 0:
            doi = pubs[0]['doi']
            q.put(doi)
        else:
            q.put(None)

    def get_cell_line(self, assay_id, q):
        # gets cell line info for assays from ChEMBL
        new_test = new_client.assay.filter(
            assay_chembl_id=assay_id).only('assay_cell_type')
        if len(new_test) > 0:
            cell_line = new_test[0]['assay_cell_type']
            q.put(cell_line)
        else:
            q.put('no data')

    def valdiate_data(self, i):
        # validates ChEMBL assays in accordance with addtional Tsonko rules
        result = False
        if(i['standard_units'] == 'nM' or i['standard_units'] == 'um' or i['standard_units'] == 'M'
           or i['standard_units'] == 'pmol' or i['standard_units'] == 'mM' or i['standard_units'] == 'fmol'
           or i['standard_units'] == 'pM' or i['standard_units'] == 'nmol' or i['standard_units'] == 'fM'):
            if (i['assay_type'] != 'U' or i['assay_type'] != 'A'):
                if(i['activity_comment'] != 'inconclusive' or i['activity_comment'] != 'Inconclusive'):
                    result = True
        return result

    def process_chembl(self, chembl_assays, temp_increment):
        # Loop through API results (20 objects per batch)
        chembl_data = dict()
        main_dict = dict()
        increment = 0
        for i in chembl_assays:
            temp_increment = temp_increment + 1
            if self.valdiate_data(i) == False:
                continue
            temp_dict = dict()
            temp_dict['protein'] = self.fetch_protein(i['target_chembl_id'])
            temp_dict['doi'] = None
            if temp_dict['protein'] == None:
                continue

            temp_dict['smiles'] = i['canonical_smiles']
            temp_dict['ligand'] = self.fetch_ligand(
                i['molecule_chembl_id'], i['canonical_smiles'])
            if temp_dict['ligand'] == None:
                continue

            if(self.check_dublicates(ligand=temp_dict["ligand"], protein=temp_dict["protein"], assay_description=i["assay_description"],
                                     document=i['document_chembl_id'],
                                     standard_value=i["standard_value"],
                                     activity=i["activity_id"],
                                     pchembl_value=i["pchembl_value"]) == True):

                continue
                # q = queue.Queue()
                # x=threading.Thread(target=self.get_cell_line, args=(i['assay_chembl_id'], q)).start()
            cell_line = None
            #
            # pub_q = queue.Queue
            # y=threading.Thread(target=self.get_dois, args=(i['document_chembl_id'], q)).start()
            # pub = q.get()
            # if pub is not None:
            #     temp_dict['doi'] = self.fetch_publication(pub)
            temp_dict['activity_id'] = i['activity_id']
            temp_dict['standard_type'] = i['standard_type']
            temp_dict['standard_value'] = i['standard_value']
            temp_dict['standard_units'] = i['standard_units']
            temp_dict['standard_relation'] = i['standard_relation']
            temp_dict['assay_description'] = i['assay_description']
            temp_dict['assay_type'] = i['assay_type']
            temp_dict['cell_line'] = cell_line
            temp_dict['pchembl_value'] = i['pchembl_value']
            temp_dict['document_chembl_id'] = i['document_chembl_id']
            temp_dict['chembl_id'] = i['molecule_chembl_id']
            temp_dict['assay_id'] = i['assay_chembl_id']

            self.upload_to_db(temp_dict)

    def analyse_rows(self):
        """
        Fetch data to models
        Saves to DB
        """
        print('---Starting---')

          #       "targetId" : 1,
          # "ligandAsTargetId" : 0,
          # "targetSpecies" : "Human",
          # "primaryTarget" : false,
          # "ligandId" : 10454,
          # "ligandContext" : null,
          # "type" : "Agonist",
          # "action" : "Agonist",
          # "affinity" : "5.6",
          # "affinityParameter" : "pIC50",
          # "dataPointInteractionIds" : [ 86671 ],
          # "refIds" : [ 38175 ]

        start = time.time()
        print("#1 Get GPCR ids from target families")
        target_list = self.get_gpcrs()
        print("#2 Get Ligands")
        ligand_list = self.get_ligands()
        print("#3 Get Ligand assays ")
        assays = self.get_ligand_assays(target_list)

        print('---process_chembl---')
        # range should be set to number of total objects/20 (API batch size)
        # 555200 is the last id saved before session was aborted


    def check_dublicates(self, ligand, protein, assay_description, acitivity, standard_value, document, pchembl_value):
        # Checks if assay experiment is already saved
        try:
            experiment = AssayExperiment.objects.filter(
                ligand=ligand, protein=protein, assay_description=assay_description,
                standard_value=standard_value, pchembl_value=pchembl_value, acitivity=acitivity, document_chembl_id=document)
            experiment = experiment.get()
            if experiment:
                return True
            else:
                return False
        except Exception as msg:
            experiment = None
            self.mylog.exception(
                "Experiment AnalyzedExperiment error | module: ASSAYEXPERIMENT.")
            return True

    def upload_to_db(self, i):
        # saves data

        chembl_data = AssayExperiment(ligand=i["ligand"],
                                      publication=i["doi"],
                                      protein=i["protein"],
                                      chembl=i["chembl_id"],
                                      smiles=i["smiles"],
                                      cell_line=i['cell_line'],
                                      activity=i["activity_id"],
                                      standard_type=i["standard_type"],
                                      standard_value=i["standard_value"],
                                      standard_units=i["standard_units"],
                                      standard_relation=i["standard_relation"],
                                      assay_description=i["assay_description"],
                                      assay_type=i["assay_type"],
                                      pchembl_value=i["pchembl_value"],
                                      document_chembl_id=i["document_chembl_id"],
                                      )
        chembl_data.save()
        # print('--saved---')

    def fetch_protein(self, target):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        test = None
        test = Protein.objects.filter(
            web_links__index=target, web_links__web_resource__slug='chembl').first()
        return test

    def fetch_ligand(self, ligand_id, smiles):
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
                        web_resource__slug='pubchem').first()
                    if cid:
                        cid = cid.index
                    else:
                        l = None
                else:
                    l = get_or_make_ligand(smiles, 'SMILES', ligand_id,)
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
