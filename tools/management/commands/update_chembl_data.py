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

    help = 'Reads bias data and imports it'
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
        print('---get_gpcrs---')
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


    def get_chembl_assay(self):
        print('---get_chembl_assay---')

                        # new_test = new_client.activity.filter(standard_units = 'nM' ).filter(standard_units = 'um' ).filter(standard_units = 'M'
                        # ).filter(standard_units = 'pmol' ).filter(standard_units = 'mM' ).filter(standard_units = 'fmol'
                        # ).filter(standard_units = 'pM' ).filter(standard_units = 'nmol' ).filter(standard_units = 'fM'
                        # ).filter(pchembl_value__isnull=False
                        # ).filter(data_validity_comment__isnull=True
                        # ).filter(standard_value__isnull = False
                        # ).only(['molecule_chembl_id', 'target_chembl_id' ,'standard_type',
                        #                                 'standard_value','standard_units','standard_relation','activity_comment',
                        #                                 'assay_description','assay_type',
                        #                                 'document_chembl_id','pchembl_value',
                        #                                 'activity_id','canonical_smiles','assay_chembl_id'])[:10]

        new_test = new_client.activity.filter(pchembl_value__isnull=False).filter(standard_value__isnull = False).filter(data_validity_comment__isnull=True).only(['molecule_chembl_id', 'target_chembl_id','standard_type',
                                        'standard_value','standard_units','standard_relation',
                                        'assay_description','assay_type',
                                        'document_chembl_id','pchembl_value',
                                        'activity_id','canonical_smiles','assay_chembl_id'])[:2]



        # print('--cheml len---', len(new_test))
        print('--cheml len---', new_test)
        return new_test


    def get_dois(self, dci, q):
        pubs = new_client.document.filter(document_chembl_id = dci).only('doi')
        if len(pubs) > 0:
            doi = pubs[0]['doi']
            q.put(doi)
        else:
            q.put(None)


    def get_cell_line(self, assay_id, q):
        new_test = new_client.assay.filter(assay_chembl_id = assay_id).only('assay_cell_type')
        if len(new_test) > 0:
            cell_line = new_test[0]['assay_cell_type']
            q.put(cell_line)
        else:
            q.put('no data')


    def process_chembl(self,chembl_assays, targets, start_time):
        chembl_data = dict()
        main_dict = dict()
        increment = 0
        temp_increment = 0
        target_list = targets
        for i in chembl_assays:
            temp_increment = temp_increment+1
            if(i['target_chembl_id'] in target_list):
                temp_dict = dict()
                temp_dict['receptor'] = self.fetch_protein( i['target_chembl_id'])
                temp_dict['doi']=None
                res = []


                if temp_dict['receptor'] and temp_dict['receptor'] != None:
                    temp_dict['smiles'] = i['canonical_smiles']
                    temp_dict['ligand'] = self.fetch_ligand(i['molecule_chembl_id'],i['canonical_smiles'])
                    if temp_dict['ligand']:
                        q = queue.Queue()
                        x=threading.Thread(target=self.get_cell_line, args=(i['assay_chembl_id'], q)).start()
                        cell_line = q.get()

                        pub_q = queue.Queue
                        y=threading.Thread(target=self.get_dois, args=(i['document_chembl_id'], q)).start()
                        pub = q.get()
                        if pub is not None:
                            temp_dict['doi'] = self.fetch_publication(pub)
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
                        if x is not None:
                            x.join()
                        if y is not None:
                            y.join()
                        increment=increment+1
                        if increment%100==0:
                            self.upload_to_db(chembl_data)
                            # x = threading.Thread(target=self.upload_to_db, args=(chembl_data,))
                            # x.start()
                            # x.join()

                            print('--status--', temp_increment)
                            increment = 0
                            chembl_data = dict()

                    else:
                        pass
                else:
                    pass
            else:
                pass
        return None

# 2-35 = 12-18
    def analyse_rows(self):
        """
        Fetch data to models
        Saves to DB
        """
        start = time.time()
        print('---Starting---')
        target_list = None
        chembl_assays = None
        with concurrent.futures.ThreadPoolExecutor() as executor:
            target_future = executor.submit(self.get_gpcrs,)
            target_list = target_future.result()
            chembl_assays_future = executor.submit(self.get_chembl_assay,)
            chembl_assays = chembl_assays_future.result()
        ## TODO: add timescale for operation
        print('---process_chembl---')

        chembl_data = self.process_chembl(chembl_assays,target_list,start)

        end = time.time()

        print('\n\n---Total time---', end - start)


    def check_dublicates(self, ligand, pub, receptor, assay_description, chembl,standard_value,standard_units, pchembl_value ):
        try:
            experiment = ChemblAssays.objects.filter(
                publication=pub, ligand=ligand, receptor=receptor, assay_description=assay_description,
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
        for i in chembl.items():
            if( self.check_dublicates(i[1]["ligand"], i[1]["doi"], i[1]["receptor"], i[1]["assay_description"],
                                   i[1]["chembl_id"],i[1]["standard_value"],i[1]["standard_units"], i[1]["pchembl_value"]) == False):
                chembl_data = ChemblAssays(ligand = i[1]["ligand"],
                                            publication = i[1]["doi"],
                                            receptor = i[1]["receptor"],
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
            else:

                pass


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
