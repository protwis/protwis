from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Pool
from django.db import IntegrityError
from build.management.commands.base_build import Command as BaseBuild
from protein.models import Protein
from ligand.models import *
from common.models import WebLink, WebResource, Publication
import logging
import time
import requests
from multiprocessing.dummy import Pool as ThreadPool

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
        if options['purge']:
            try:
                self.purge_bias_data()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
        self.analyse_rows()
        try:
            print('Updatind ligand data from GuideToPharma')

            print('COMPLETED updating GuideToPharma Data')
        except Exception as msg:
            print('--error--', msg, '\n')
            self.logger.info("The error appeared in def handle")

# pylint: disable=R0201
    def purge_bias_data(self):
        print("# Purging data")

        # delete_bias_experiment = Ligand.objects.all()
        # delete_bias_experiment.delete()
        # delete_ligand_props = LigandProperities.objects.all()
        # delete_ligand_props.delete()
        print("# Data purged")

    def analyse_rows(self):
        """
        Fetch data to models
        Saves to DB
        """
        print('---Starting---\n')
        print("\n#1 Get GPCR ids from target families")
        target_list = self.get_gpcrs()
        print("\n#2 Get Ligands")
        self.get_ligands()
        print("\n#3 Get Ligand assays")
        assays = self.get_ligand_assays(target_list)
        print("\n#4 Process Ligand assays", len(assays), ' assays')
        self.process_ligand_assays(assays)
        print("\n#5 Get Endogeneous ligands")
        # self.get_endogenous(target_list)
        print('\n\n---Finished---')

    def get_endogenous(self, targets):
        for target in targets:
            protein = self.fetch_protein(target)
            response = requests.get(
                "https://www.guidetopharmacology.org/services/targets/" + str(target) + "/naturalLigands")
            if response.status_code == 200:
                data = response.json()
                for i in data:
                    try:
                        ligand_name = str()
                        try:
                            ligand_name = self.ligand_cache[i['targetId']]
                        except:
                            ligand_name = ""
                        ligand = self.fetch_ligand(
                            data[0]['ligandId'], data[0]['type'], ligand_name)
                        if ligand and protein:
                            protein.endogenous_ligands.add(ligand)
                            Ligand.objects.update_or_create(
                                id=ligand.id,
                                defaults={'endogenous': True},
                            )
                        else:
                            pass
                    except:
                        pass

    def upload_to_db(self, i):
        # saves data
        chembl_data = AssayExperiment(source=i["source"],
                                      ligand=i["ligand"],
                                      publication=i["doi"],
                                      protein=i["protein"],
                                      published_type=i["standard_type"],
                                      published_value=i["standard_value"],
                                      assay_description=i["assay_description"],
                                      )
        chembl_data.save()

# pylint: disable=R0201
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

# pylint: disable=R0201
    def fetch_ligand(self, ligand_id, ligand_type, name):
        """
        fetch ligands with Ligand model
        requires: ligand id, ligand id type, ligand name
        requires: source_file name
        """
        l = None

        try:
            l = Ligand.objects.filter(
                properities__web_links__index=ligand_id).first()
            if l:
                cid = l.properities.web_links.filter(
                    web_resource__slug='gtoplig').first()
                if cid:
                    return l
                else:
                    lig = Ligand()
                    l = lig.load_by_gtop_id(name, ligand_id, ligand_type)
            else:
                lig = Ligand()
                l = lig.load_by_gtop_id(name, ligand_id, ligand_type)
        except Exception:
            l = None
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

    def save_ligand_copy(self, ligand):
        """
        fetch ligands with Ligand model
        requires: ligand id, ligand id type, ligand name
        requires: source_file name
        """
        self.ligand_cache.update({ligand['entry']: ligand['name']})

    def get_ligands(self):
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

# pylint: disable=R0201
    def get_gpcrs(self):
        target_list = list()
        url = "https://www.guidetopharmacology.org/services/targets/families"
        response = ''
        while response == '':
            try:
                response = requests.get(url)
                break
            except:
                print("Connection refused by the server..")
                print("Let me sleep for 1 second")
                print("ZZzzzz...")
                time.sleep(1)
                print("Was a nice sleep, now let me continue...")
                continue

        for entry in response.json():
            try:
                if entry['parentFamilyIds'] != None:
                    if entry['parentFamilyIds'][0] == 694 or entry['parentFamilyIds'][0] == 115:
                        target_list.extend(entry['targetIds'])
            except:
                pass
        return target_list

# pylint: disable=R0201
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
        return assay_list

    def create_source(self, ligand_id):
        source = AssayExperimentSource()
        web_resource = WebResource.objects.get(slug='gtoplig')
        try:
            wl, created = WebLink.objects.get_or_create(index=ligand_id,
                                                        web_resource=web_resource)
        except IntegrityError:
            wl = Weblink.objects.get(
                index=ligand_id, web_resource=web_resource)

        source.database = 'GuideToPharmacology'
        source.database_id = ligand_id
        try:
            source.save()
            source.web_links.add(wl)
        except IntegrityError:
            return AssayExperimentSource.objects.get(web_links=wl,
                                                    database='GuideToPharmacology',
                                                    database_id=ligand_id)
        return source

    def process_ligand_assays(self, assays):
        start = time.time()
        data = assays
        print('type assay', type(data))

        pool = ThreadPool(4)
        pool.map(self.process_save, data)
        pool.close()
        pool.join()
        print('5 process/thread total time: ', time.time() - start, '\n\n')

    def do_something(self, data):
        print(data)

    def process_save(self, i):
        temp_dict = dict()
        temp_dict['protein'] = self.fetch_protein(i['targetId'])
        ligand_name = str()
        try:
            ligand_name = self.ligand_cache[i['targetId']]
        except:
            ligand_name = i['ligandId']
        temp_dict['ligand'] = self.fetch_ligand(
            i['ligandId'], i['type'], ligand_name)
        if i['refIds']:
            try:
                temp_dict['doi'] = self.fetch_publication(i['refIds'])
            except:
                temp_dict['doi'] = None
        else:
            temp_dict['doi'] = None
        if temp_dict['protein'] is not None:
            if temp_dict['ligand'] is not None:
                temp_dict['source'] = self.create_source(i['ligandId'])
                temp_dict['standard_type'] = i['affinityParameter']
                temp_dict['standard_value'] = i['affinity']
                temp_dict['assay_description'] = i['ligandContext']
                if temp_dict['assay_description'] == None:
                    temp_dict['assay_description'] = "No data available"
                self.upload_to_db(temp_dict)
