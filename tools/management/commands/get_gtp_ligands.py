from django.db import IntegrityError
from build.management.commands.base_build import Command as BaseBuild
from protein.models import Protein
from ligand.models import GTP_endogenous_ligand, Ligand
from common.models import WebLink, WebResource, Publication
import logging
import time
import requests

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

    help = 'Updates GuideToPharma data and imports it'
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


# pylint: disable=R0201
    def purge_bias_data(self):
        print("# Purging data")
        delete_bias_experiment = GTP_endogenous_ligand.objects.all()
        delete_bias_experiment.delete()

    def analyse_rows(self):
        """
        Fetch data to models
        Saves to DB
        """
        print('---Starting---\n')
        print("\n#1 Get GPCR ids from target families")
        target_list = self.get_gpcrs()
        print("\n#2 Get Endogeneous ligands")
        self.get_endogenous(target_list)
        print("\n#3 Get endogenous data from GPCRDb")
        endogenous_ligands_from_db = self.endogenous_ligands_from_db()
        print("\n#4 Convert_query_to_dict" )
        endogenous_list = self.convert_query_to_dict(endogenous_ligands_from_db)
        print("\n#5 Combine same assays")
        combined_assays = self.process_assays(endogenous_list)
        print("\n#6 Prepare to calculate averages and save")
        prepared_data = self.prepare_calculate_averages_and_save(combined_assays)
        print("\n#7 Calculate averages")
        averaged_list = self.calculate_averages(prepared_data)
        print("\n#8 Save")
        self.save_data(averaged_list)
        print('\n\n---Finished---')

    def get_endogenous(self, targets):
        for target in targets:
            protein = self.fetch_protein(target)
            response = ''
            while response == '':
                try:
                    response = requests.get(
                        "https://www.guidetopharmacology.org/services/targets/" + str(target) + "/naturalLigands")
                    if response.status_code == 200:
                        data = response.json()
                        for i in data:
                            # try:
                            # import pdb; pdb.set_trace()
                            ligand_name = str()
                            try:
                                ligand_name = self.ligand_cache[i['targetId']]
                            except:
                                ligand_name = ""

                            ligand = self.fetch_ligand(
                                i['ligandId'], i['type'], ligand_name)
                            if ligand and protein:
                                # if target == 1:

                                self.get_ligand_interactions(target=target,ligand=ligand, ligand_id_gtp=data[0]['ligandId'], ligand_type=data[0]['type'],receptor=protein)
                            else:
                                pass
                            # except:
                            #     pass
                except:
                    print("Connection refused by the server..")
                    time.sleep(1)
                    response == ''


    def get_ligand_interactions(self, target, ligand, ligand_id_gtp, ligand_type, receptor):
        response = requests.get(
            "https://www.guidetopharmacology.org/services/ligands/" + str(ligand_id_gtp) + "/interactions")
        if response.status_code == 200:
            # import pdb; pdb.set_trace()
            data = response.json()
            for interaction in data:
                if interaction['targetId'] == target:
                    # import pdb; pdb.set_trace()
                    emin=emax=eavg=kmin=kmax=kavg = None
                    if interaction['endogenous'] == True:
                        # import pdb; pdb.set_trace()
                        if interaction['affinityParameter'] == 'pKi':
                            try:
                                temp = interaction['affinity'].strip().split("-")
                                kmin = float(temp[0])
                                kmax = float(temp[1])
                                kavg = (kmin + kmax) / 2
                            except:
                                kavg = float(interaction['affinity'])

                        if interaction['affinityParameter'] == 'pEC50':
                            try:
                                temp = interaction['affinity'].strip().split("-")
                                emin = float(temp[0])
                                emax = float(temp[1])
                                eavg = (emin + emax) / 2
                            except:
                                eavg = float(interaction['affinity'])
                        if kavg or eavg:

                            try:
                                ligand_type = ligand.properities.ligand_type.name
                            except:
                                ligand_type = None
                            try:
                                endo_ligand_type = interaction['type']
                            except:
                                endo_ligand_type = None

                            link = "https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId="+str(ligand_id_gtp)
                            if self.fetch_experiment(ligand=ligand, receptor=receptor, pavg=kavg, eavg=eavg)==False:
                                gtp_data = GTP_endogenous_ligand(
                                        ligand = ligand,
                                        ligand_type = ligand_type,
                                        endogenous_princip = endo_ligand_type,
                                        receptor = receptor,
                                        pec50_avg = eavg,
                                        pec50_min = emin,
                                        pec50_max = emax,
                                        pKi_avg = kavg,
                                        pKi_min = kmin,
                                        pKi_max = kmax,
                                        gpt_link = link,
                                )
                                gtp_data.save()
                            try:
                                for reference in interaction['refs']:
                                    publication = self.fetch_publication(reference['pmid'])
                                    gtp_data.publication.add(publication)
                            except:
                                publication= None
                        else:
                            pass

    def fetch_experiment(self, ligand, receptor, pavg, eavg):
        '''
        fetch receptor with Protein model
        requires: protein id, source
        '''
        try:
            experiment = GTP_endogenous_ligand.objects.filter(
                 ligand=ligand, receptor=receptor, pKi_avg=pavg,pec50_avg=eavg)
            experiment = experiment.get()
            return True
        except Exception:
            self.logger.info('fetch_experiment error')
            experiment = None
            return False

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
                self.mylog.debug(
                    "ligand fetching error | module: fetch_publication. Row # is :")


    def get_gpcrs(self):
        target_list = list()
        url = "https://www.guidetopharmacology.org/services/targets/families"
        response = ''
        while response == '':
            try:
                response = requests.get(url)

            except:
                print("Connection refused by the server..")
                print("Let me sleep for 1 second")
                print("ZZzzzz...")
                time.sleep(1)
                print("Was a nice sleep, now let me continue...")
                response == ''

        for entry in response.json():
            try:
                if entry['parentFamilyIds'] != None:
                    if entry['parentFamilyIds'][0] == 694 or entry['parentFamilyIds'][0] == 115:
                        target_list.extend(entry['targetIds'])
            except:
                print("Was a nice sleep, now let me continue...")
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
        return assay_list

    def endogenous_ligands_from_db(self):
        return GTP_endogenous_ligand.objects.all()

    def convert_query_to_dict(self, endogenous_queryset):
        endogenous_list = list()
        for assay in endogenous_queryset:
            single_dict = dict()
            # import pdb; pdb.set_trace()
            single_dict['id'] = assay
            single_dict['ligand'] = assay.ligand
            single_dict['ligand_type'] = assay.ligand_type
            single_dict['endogenous_princip'] = assay.endogenous_princip
            single_dict['receptor'] = assay.receptor
            single_dict['pec50_avg'] = assay.pec50_avg
            single_dict['pec50_min'] = assay.pec50_min
            single_dict['pec50_max'] = assay.pec50_max
            single_dict['pKi_avg'] = assay.pKi_avg
            single_dict['pKi_min'] = assay.pKi_min
            single_dict['pKi_max'] = assay.pKi_max
            single_dict['gpt_link'] = assay.gpt_link
            single_dict['references'] = list()
            for publication in assay.publication.all():
                single_dict['references'].append(publication)
            endogenous_list.append(single_dict)
        return endogenous_list

    def process_assays(self, endogenous_list):
        context = dict()
        for assay in endogenous_list:
            name = str(assay['ligand'].id) + \
                '/' + str(assay['receptor'].id)
            if name not in context:
                context[name] = dict()
                context[name]['ligand'] = assay['ligand']
                context[name]['ligand_type'] = assay['ligand_type']
                context[name]['endogenous_princip'] = assay['endogenous_princip']
                context[name]['receptor'] = assay['receptor']
                context[name]['pec50_avg'] = list()
                context[name]['pec50_min'] = list()
                context[name]['pec50_max'] = list()
                context[name]['pKi_avg'] = list()
                context[name]['pKi_min'] = list()
                context[name]['pKi_max'] = list()
                context[name]['gpt_link'] = assay['gpt_link']
                context[name]['references'] = list()

            if assay['pec50_avg'] is not None and isinstance(assay['pec50_avg'], float):
                context[name]['pec50_avg'].append(assay['pec50_avg'])
            if assay['pec50_min'] is not None and isinstance(assay['pec50_min'], float):
                context[name]['pec50_min'].append(assay['pec50_min'])
            if assay['pec50_max'] is not None and isinstance(assay['pec50_max'], float):
                context[name]['pec50_max'].append(assay['pec50_max'])
            if assay['pKi_avg'] is not None and isinstance(assay['pKi_avg'], float):
                context[name]['pKi_avg'].append(assay['pKi_avg'])
            if assay['pKi_min'] is not None and isinstance(assay['pKi_min'], float):
                context[name]['pKi_min'].append(assay['pKi_min'])
            if assay['pKi_max'] is not None and isinstance(assay['pKi_max'], float):
                context[name]['pKi_max'].append(assay['pKi_max'])
            context[name]['references'].extend(assay['references'])
        return context

    def prepare_calculate_averages_and_save(self, context):
        prepared_data = list()
        for i in context.items():
            prepared_data.append(i[1])
        return prepared_data

    def calculate_averages(self, prepared_data):
        for assay in prepared_data:

            try:
                assay['pKi_avg'] = round(sum(assay['pKi_avg'])/len(assay['pKi_avg']),2)
            except:
                if len(assay['pKi_avg']) < 1:
                    assay['pKi_avg'] = None
                # import pdb; pdb.set_trace()
            try:
                assay['pKi_min'] = round(sum(assay['pKi_min'])/len(assay['pKi_min']),2)
            except:
                if len(assay['pKi_min']) < 1:
                    assay['pKi_min'] = None
                # import pdb; pdb.set_trace()
            try:
                assay['pKi_max'] = round(sum(assay['pKi_max'])/len(assay['pKi_max']),2)
            except:
                if len(assay['pKi_max']) < 1:
                    assay['pKi_max'] = None
                # import pdb; pdb.set_trace()
            try:
                assay['pec50_avg'] = round(sum(assay['pec50_avg'])/len(assay['pec50_avg']),2)
            except:
                if len(assay['pec50_avg']) < 1:
                    assay['pec50_avg'] = None
                # import pdb; pdb.set_trace()
            try:
                assay['pec50_min'] = round(sum(assay['pec50_min'])/len(assay['pec50_min']),2)
            except:
                if len(assay['pec50_min']) < 1:
                    assay['pec50_min'] = None
                # import pdb; pdb.set_trace()
            try:
                assay['pec50_max'] = round(sum(assay['pec50_max'])/len(assay['pec50_max']),2)
            except:
                if len(assay['pec50_max']) < 1:
                    assay['pec50_max'] = None
                # import pdb; pdb.set_trace()
        return prepared_data

    def save_data(self, save_data):
        for assay in save_data:
            gtp_data = GTP_endogenous_ligand(
                    ligand = assay['ligand'],
                    ligand_type = assay['ligand_type'],
                    endogenous_princip = assay['endogenous_princip'],
                    receptor = assay['receptor'],
                    pec50_avg = assay['pec50_avg'],
                    pec50_min = assay['pec50_min'],
                    pec50_max = assay['pec50_max'],
                    pKi_avg = assay['pKi_avg'],
                    pKi_min = assay['pKi_min'],
                    pKi_max = assay['pKi_max'],
                    gpt_link = "GPCRDb",
            )
            gtp_data.save()
