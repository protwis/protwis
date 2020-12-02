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
        # self.analyse_rows()
        try:
            print('Updatind ChEMBL data')
            self.analyse_rows()
            self.logger.info('COMPLETED updating ChEMBL Data')
        except Exception as msg:
            print('--error--', msg, '\n')
            self.logger.info("The error appeared in def handle")


    def purge_bias_data(self):
        delete_bias_excel = BiasedExperiment.objects.all()
        delete_bias_excel.delete()
        delete_bias_experiment = AnalyzedExperiment.objects.all()
        delete_bias_experiment.delete()


    def get_gpcrs(self):
        print('---get_gpcrs from ChEMBL---')
        g_family_ids = [ 1,147,165,166,202,281,407,435,446,460,468,479,480,483,484,486,487,491,499,500,501,502,503,504,506,507,508,509,
        510,515,516,517,518,528,533,534,535,540,541,542,544,547,548,549,550,551,554,555,
        556,558,559,561,562,563,565,566,567,568,569,570,571,573,574,603,604,605,606,607,
        608,609,610,611,612,613,614,615,616,617,618,619,620,621,830,1020,1021,1022,1038,
        1082,1083,1088,1089,1251,1253,1256,1259,1265,1266,1267,1268,1269,1270,1271,1272,
        1273,1274,1275,1277,1278]
        target_list = set()
        for item in g_family_ids:
            gproteins = new_client.target_component.filter(protein_classifications__protein_classification_id=item).only('targets')
            for protein in gproteins:
                for target in protein['targets']:
                    if target['target_chembl_id'] not in target_list:
                        target_list.add(target['target_chembl_id'])
                    else:
                        pass
        print('GPCRS ready')
        return target_list


    def get_chembl_assay(self,targets, prev_id, current_id):
        # gets assays from ChEMBL (by batch sie of 20). Filters using Tsonko rules, GPCRs.
        # as arguments takes GPCR ids(targets), batch start size (prev_id) and end size(current_id)
        new_test = new_client.activity.filter(pchembl_value__isnull=False).filter(data_validity_comment__isnull=True
                                        ).filter(standard_value__isnull = False
                                        ).filter(standard_units__isnull = False
                                        ).filter(target_chembl_id__in = targets
                                        ).only(['molecule_chembl_id', 'target_chembl_id' ,'standard_type',
                                        'standard_value','standard_units','standard_relation','activity_comment',
                                        'assay_description','assay_type',
                                        'document_chembl_id','pchembl_value',
                                        'activity_id','canonical_smiles','assay_chembl_id'])[prev_id:current_id]
        return new_test

    def get_dois(self, dci, q):
        # gets references for assays from ChEMBL (DOI)
        pubs = new_client.document.filter(document_chembl_id = dci).only('doi')
        if len(pubs) > 0:
            doi = pubs[0]['doi']
            q.put(doi)
        else:
            q.put(None)

    def get_cell_line(self, assay_id, q):
        # gets cell line info for assays from ChEMBL
        new_test = new_client.assay.filter(assay_chembl_id = assay_id).only('assay_cell_type')
        if len(new_test) > 0:
            cell_line = new_test[0]['assay_cell_type']
            q.put(cell_line)
        else:
            q.put('no data')

    def valdiate_data(self, i):
        #validates ChEMBL assays in accordance with addtional Tsonko rules
        result = False
        if(i['standard_units'] == 'nM' or i['standard_units'] == 'um' or i['standard_units'] == 'M'
        or i['standard_units'] == 'pmol' or i['standard_units'] == 'mM' or i['standard_units'] == 'fmol'
        or i['standard_units'] == 'pM' or i['standard_units'] == 'nmol' or i['standard_units'] == 'fM'):
            if ( i['assay_type'] != 'U' or i['assay_type'] != 'A'):
                if( i['activity_comment'] != 'inconclusive'  or i['activity_comment'] != 'Inconclusive'):
                    result = True
        return result

    def process_chembl(self,chembl_assays, temp_increment):
        #Loop through API results (20 objects per batch)
        chembl_data = dict()
        main_dict = dict()
        increment = 0
        for i in chembl_assays:
            temp_increment = temp_increment+1
            if self.valdiate_data(i) == False:
                continue
            temp_dict = dict()
            temp_dict['protein'] = self.fetch_protein( i['target_chembl_id'])
            temp_dict['doi']=None
            if temp_dict['protein'] == None:
                continue

            temp_dict['smiles'] = i['canonical_smiles']
            temp_dict['ligand'] = self.fetch_ligand(i['molecule_chembl_id'],i['canonical_smiles'])
            if temp_dict['ligand'] == None:
                continue

            if( self.check_dublicates(temp_dict["ligand"], temp_dict["protein"], i["assay_description"],
                           i["molecule_chembl_id"],
                           i["standard_value"],i["standard_units"],
                           i["pchembl_value"]) == True):

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
            chembl_data[increment] = temp_dict
            increment=increment+1
            self.upload_to_db(chembl_data)

    def analyse_rows(self):
        """
        Fetch data to models
        Saves to DB
        """
        print('---Starting---')
        current_id = 0
        prev_id = 0
        start = time.time()
        target_list = self.get_gpcrs()
        target_list_list = list(target_list)
        start = time.time()
        chembl_assays = None
        print('---process_chembl---')
        #range should be set to number of total objects/20 (API batch size)
        #555200 is the last id saved before session was aborted
        for i in range(30578):
            current_id =  591900 + ((i+1) * 20)
            prev_id =  591900 + (i *20)
            chembl_assays = self.get_chembl_assay(target_list_list, prev_id, current_id)
            chembl_assays=list(chembl_assays)
            self.process_chembl(chembl_assays,current_id)
            # control the flow
            if(current_id%100==0):
                end = time.time()
                print('---temp_increment time---',current_id,  end - start)


    def check_dublicates(self, ligand, protein, assay_description, chembl,standard_value,standard_units, pchembl_value ):
        # Checks if assay experiment is already saved
        try:
            experiment = AssayExperiment.objects.filter(
                ligand=ligand, protein=protein, assay_description=assay_description,
                chembl=chembl,standard_value=standard_value,standard_units=standard_units,pchembl_value=pchembl_value )
            experiment = experiment.get()
            if experiment:
                return True
            else:
                return False
        except Exception as msg:
            experiment = None
            self.mylog.exception(
                "Experiment AnalyzedExperiment error | module: AnalyzedExperiment.")
            return False

    def upload_to_db(self, chembl):
        # saves data
        for i in chembl.items():
            chembl_data = AssayExperiment(ligand = i[1]["ligand"],
                                        publication = i[1]["doi"],
                                        protein = i[1]["protein"],
                                        chembl = i[1]["chembl_id"],
                                        smiles = i[1]["smiles"],
                                        cell_line = i[1]['cell_line'],
                                        activity = i[1]["activity_id"],
                                        standard_type  = i[1]["standard_type"],
                                        standard_value = i[1]["standard_value"],
                                        standard_units = i[1]["standard_units"],
                                        standard_relation = i[1]["standard_relation"],
                                        assay_description = i[1]["assay_description"],
                                        assay_type = i[1]["assay_type"],
                                        pchembl_value = i[1]["pchembl_value"],
                                        document_chembl_id = i[1]["document_chembl_id"],
                                                  )
            chembl_data.save()
            # print('--saved---')

    def fetch_measurements(self, potency, p_type, unit):
        # it was used for bias prediction build. Temporarily unused
        if p_type.lower()  == 'pec50':
            potency = 10**(potency*(-1))
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
        elif p_type.lower()  == 'ec50':
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
        if potency:
            potency = "{:.2E}".format(Decimal(potency))
        return potency,p_type

    def fetch_protein(self,target):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        test = None
        test = Protein.objects.filter(web_links__index = target, web_links__web_resource__slug = 'chembl').first()
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
                l = Ligand.objects.filter(properities__web_links__index=ligand_id).first()
                if l:
                    cid = l.properities.web_links.filter(web_resource__slug = 'pubchem').first()
                    if cid:
                        cid = cid.index
                    else:
                        l = None
                else:
                    l = get_or_make_ligand(smiles, 'SMILES', ligand_id,  )
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
