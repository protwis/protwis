from django.core.management.base import BaseCommand, CommandError
from django.db import IntegrityError
from django.conf import settings
from build.management.commands.base_build import Command as BaseBuild
from protein.models import Protein, ProteinGProteinPair
from ligand.models import *
from ligand.functions import get_or_make_ligand
from common.models import WebLink, WebResource, Publication
from multiprocessing.pool import ThreadPool
import queue
import logging
import xlrd
from datetime import datetime
import operator
import json
import requests
import concurrent.futures
import os
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
    structure_data_dir = os.sep.join([settings.DATA_DIR, 'ligand_data', 'pubchem'])

    help = 'Updates PubChem data and imports it'
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
            print('Updatind ligand data from PubChem')
            # self.analyse_rows()
            print('COMPLETED updating PubChem Data')
        except Exception as msg:
            print('--error--', msg, '\n')
            self.logger.info("The error appeared in def handle")

    def analyse_rows(self):
        """
        Fetch data to models
        Saves to DB
        """
        aids = list()
        start = time.time()
        print('---Starting---\n')
        # print("\n#1 Get pubchem gene ids from target aceesions")
        # TODO: uncomment if gene data required
        # target_list = self.read_csv_genes()
        # final_target_list = self.filter_genes_by_gpcrdb_accessions(target_list,gpcrdb_prots)
        # TODO: finish uncomment Part

        # TODO: test part - delete
        # self.prosess_assay('asd')
        self.prosess_assay(259977)

        # TODO: end of test_run

        print("\n#1 Get Proteins from GPCRDb")
        gpcrdb_prots = self.get_proteins_from_gpcrdb()
        print("\n#2 Get assays for gpcrs")
        # test = gpcrdb_prots[:10]
        # aids = self.process_accessions(test)
        # self.process_aids(aids)
        print("\n#3 Get assays")
        # assays = self.get_ligand_assays(target_list)
        # print("\n#4 Process Ligand assays", len(assays), ' assays')
        # self.process_ligand_assays(assays)
        print("\n#5 Get Endogeneous ligands")
        # self.get_endogenous(target_list)
        end = time.time()

        print('\n\n---Finished---',end - start)
    #processing excel part
    def read_csv_genes(self):
        filenames = os.listdir(self.structure_data_dir)
        print('***Working destination: ',filenames)
        for source_file in filenames:
            source_file_path = os.sep.join(
                [self.structure_data_dir, source_file]).replace('//', '/')
            if os.path.isfile(source_file_path) and source_file[0] != '.':
                self.logger.info('Reading file {}'.format(source_file_path))
                # read the yaml file
                rows = []
                if source_file[-4:] == 'xlsx' or source_file[-3:] == 'xls':
                    if "~$" in source_file:
                        # ignore open excel files
                        continue
                    rows = self.loaddatafromexcel(source_file_path)
                    rows = self.analyse_excel_rows(rows, source_file)
                else:
                    self.mylog.debug('unknown format'.source_file)
                    continue

                self.data_all += rows
        print("***Total data points", len(self.data_all))
        print("***Finished reading excel")
        return self.data_all

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

    def analyse_excel_rows(self, rows, source_file):
        """
        Reads excel rows one by one
        Fetch data to models
        Saves to DB
        """
        print('***Accessing gen_id for every accession code')
        skipped = 0
        # Analyse the rows from excel and assign the right headers
        temp = list()
        for i, r in enumerate(rows, 1):
            if( isinstance(r[2], str) and len(r[2].split('|')) > 1):
                for i in r[2].split('|'):
                    temp.append({i:int(r[0])})
            else:
                temp.append({r[2]:int(r[0])})
        return temp
    #end of processing excel part

    #retrieve proteins from GPCRDb and assign genes by accession
    def get_proteins_from_gpcrdb(self):
        proteins = list()
        test = None
        test = Protein.objects.values_list('accession', flat=True).distinct()
        print('***Total unique accession codes',len(test))
        return test

    def filter_genes_by_gpcrdb_accessions(self, pubchem_data, gpcrdb_data):
        print('***Comparing lists')
        process_set = list()
        increment_1 = 0
        test = gpcrdb_data[:10]

        if 'A0A096N764' in test:
            print(test)
            print(type(test))
        print('***pubchem and gpcrdb', len(pubchem_data),len(gpcrdb_data))
        for entry in pubchem_data:
            increment_1 = increment_1+1
            a = list(entry.keys())[0]
            if a in gpcrdb_data:
                process_set.append(entry)
        print(len(process_set),increment_1)
        return process_set
    #end of retrieve proteins from GPCRDb and assign genes by accession

    #get aids from pubchem
    def process_accessions(self, gpcrdb_prots):
        aids = list()
        for gpcr in gpcrdb_prots:
            try:
                local_aid = self.get_aids_from_pubchem(gpcr)
                aids.extend(local_aid)
            except:
                pass
        return aids

    #get sid and assay data_all
    def process_aids(self, aids):
        assays = list()
        for aid in aids:
            try:
                local_aid = self.prosess_assay(aid)
                assays.extend(local_aid)
            except:
                pass
        return assays

    def prosess_assay(self,aid):
        assay = dict()
        assay_response = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/"+str(aid)+"/json")
        if assay_response.status_code == 200:
            try:
                assay_info = dict()
                assay_subs = list()
                assay_data = assay_response.json()
                assay_info['assay_id'] = assay_data['PC_AssaySubmit']['assay']['descr']['aid']['id']
                assay_info['aid_source_db'] = assay_data['PC_AssaySubmit']['assay']['descr']['aid_source']['db']['name']
                assay_info['aid_source_db_id'] = assay_data['PC_AssaySubmit']['assay']['descr']['aid_source']['db']['source_id']['str']
                assay_info['assay_description'] = assay_data['PC_AssaySubmit']['assay']['descr']['name']
                temp_dois = assay_data['PC_AssaySubmit']['assay']['descr']['xref']
                for i in temp_dois:
                    if any('pmid' in d for d in i['xref']) == True:
                        assay_info['DOI/pubmed']= i['xref']['pmid']
                assay_quantitive = assay_data['PC_AssaySubmit']['data']
                for substance in assay_quantitive:
                    temp_assay= dict()
                    temp_assay['sid'] = substance['sid']
                    for i in substance['data']:
                        if i['tid'] == 2:
                            temp_assay['Standard_type'] = i['value']['sval']
                        if i['tid'] == 3:
                            temp_assay['Standard_relation'] = i['value']['sval']
                        if i['tid'] == 4:
                            temp_assay['Standard_value'] = i['value']['fval']
                        if i['tid'] == 5:
                            temp_assay['Standard-unit'] = i['value']['sval']

                    assay_subs.append(temp_assay)
                assay_info['assays'] = assay_subs
                cid= self.get_cid(substance['sid'])
                assay_info['ligand'] = self.get_ligand_or_create(cid)
                print(assay_info)
                return assay_info
            except:
                return None

        return assay

    def get_cid(self,sid):

        substance_response = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/"+str(sid)+"/json")
        if substance_response.status_code == 200:
            # TODO: try except
            substance_data = substance_response.json()
            cid = substance_data['PC_Substances'][0]['compound'][1]['id']['id']['cid']
            return cid

    def get_ligand_or_create(self,cid):
        ligand_name = str()
        properties = dict()
        ligand_name_response = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"+str(cid)+"/synonyms/json")
        if ligand_name_response.status_code == 200:
            ligand_name = ligand_name_response.json()
            ligand_name = ligand_name['InformationList']['Information'][0]['Synonym'][0]
        compound_response = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"+str(cid)+"/property/CanonicalSMILES,InChIKey,MolecularWeight,HBondDonorCount,HBondAcceptorCount,XLogP,RotatableBondCount/json")
        if compound_response.status_code == 200:
            # TODO: try except
            compound_data = compound_response.json()
            pubchem = compound_data
            if pubchem['PropertyTable']['Properties'][0]:
                if 'HBondAcceptorCount' in pubchem['PropertyTable']['Properties'][0] :
                    properties['hacc'] =  pubchem['PropertyTable']['Properties'][0]['HBondAcceptorCount']
                if 'HBondDonorCount' in pubchem['PropertyTable']['Properties'][0] :
                    properties['hdon'] =  pubchem['PropertyTable']['Properties'][0]['HBondDonorCount']
                if 'XLogP' in pubchem['PropertyTable']['Properties'][0] :
                    properties['logp'] =  pubchem['PropertyTable']['Properties'][0]['XLogP']
                if 'RotatableBondCount' in pubchem['PropertyTable']['Properties'][0] :
                    properties['rotatable_bonds'] =  pubchem['PropertyTable']['Properties'][0]['RotatableBondCount']
                if 'MolecularWeight' in pubchem['PropertyTable']['Properties'][0] :
                    properties['mw'] = pubchem['PropertyTable']['Properties'][0]['MolecularWeight']
            try:
                properties['smiles'] =  pubchem['PropertyTable']['Properties'][0]['CanonicalSMILES']
                properties['inchikey'] =  pubchem['PropertyTable']['Properties'][0]['InChIKey']
            except:
                return None
            lp = self.create_ligand_properties(cid,properties)
            ligand = self.create_ligand(lp, ligand_name)

            return ligand

    def create_ligand(self, lp, ligand_name):
        try:
            existing_ligand = Ligand.objects.get(name=ligand_name, canonical=True)
            return existing_ligand
        except Ligand.DoesNotExist:
            ligand = Ligand()
            ligand.properities = lp
            ligand.name = ligand_name
            ligand.canonical = True
            ligand.ambigious_alias = False
            ligand.pdbe = None
            try:
                ligand.save()
            except IntegrityError:
                return Ligand.objects.get(name=ligand_name, canonical=True)
            return ligand

    def create_ligand_properties(self, cid, structure):
        web_resource = WebResource.objects.get(slug='pubchem')
        try:
            wl, created = WebLink.objects.get_or_create(index=cid, web_resource=web_resource)
        except IntegrityError:
            wl = Weblink.objects.get(index=cid, web_resource=web_resource)
        lp = LigandProperities()
        try:
            lt = LigandType.objects.filter(name = ligand_type)[0]
            lp.ligand_type = lt
        except :
            lt =  LigandType.objects.filter(name = 'small molecule')[0]
            lp.ligand_type = lt
        lp.smiles = structure['smiles']
        lp.inchikey = structure['inchikey']
        # lp.sequence= structure['sequence']
        lp.mw = structure['mw']
        lp.rotatable_bonds = structure['rotatable_bonds']
        lp.hacc = structure['hacc']
        lp.hdon = structure['hdon']
        lp.logp = structure['logp']
        try:
            lp.save()
            lp.web_links.add(wl)
        except IntegrityError:
            lp = LigandProperities.objects.get(inchikey=structure['inchikey'])
        return lp

    def get_aids_from_pubchem(self, accession):
        structure_response = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/accession/"+accession+"/aids/json")
        if structure_response.status_code == 200:
            try:
                ligand_data = structure_response.json()
                temp=ligand_data['IdentifierList']
                return temp['AID']
            except:
                return None

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
                            data[0]['ligandId'], data[0]['type'],ligand_name)
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

    def upload_to_db(self, i):
        # saves data
        print('data saved')
        chembl_data = AssayExperiment(ligand=i["ligand"],
                                      publication=i["doi"],
                                      protein=i["protein"],
                                      published_type=i["standard_type"],
                                      published_value=i["standard_value"],
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

    def fetch_ligand(self, ligand_id, type, name):
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
                    l = lig.load_by_gtop_id(name, ligand_id, type)
            else:
                lig = Ligand()
                l = lig.load_by_gtop_id(name, ligand_id, type)
        except Exception:
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

    def save_ligand_copy(self, ligand):
        """
        fetch ligands with Ligand model
        requires: ligand id, ligand id type, ligand name
        requires: source_file name
        """
        self.ligand_cache.update({ligand['entry']:ligand['name'] })

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
                if entry['parentFamilyIds'] != None:
                    if entry['parentFamilyIds'][0] == 694 or entry['parentFamilyIds'][0] == 115:
                        target_list.extend(entry['targetIds'])
            except:
                pass
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

    def process_ligand_assays(self, assays):

        for i in assays:
            temp_dict = dict()
            temp_dict['protein'] = self.fetch_protein(i['targetId'])
            ligand_name = str()
            try:
                ligand_name = self.ligand_cache[i['targetId']]
            except:
                ligand_name = ""
            temp_dict['ligand'] = self.fetch_ligand(i['ligandId'], i['type'], ligand_name)
            if i['refIds']:
                try:
                    temp_dict['doi'] = self.fetch_publication(i['refIds'])
                except:
                    temp_dict['doi'] = None
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
