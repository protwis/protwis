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
    data_point = 0

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
        # self.prosess_assay(259977,'H0XT53')

        # TODO: end of test_run

        print("\n#1 Get Proteins from GPCRDb")
        gpcrdb_prots = self.get_proteins_from_gpcrdb()
        print("\n#2 Get assays for gpcrs")
        test = gpcrdb_prots[:10]
        aids = self.process_accessions(gpcrdb_prots)
        print('*** len aids', len(aids))
        print("\n#3 Process assays")
        self.process_aids(aids)
        # self.get_endogenous(target_list)
        end = time.time()

        print('\n\n---Finished---',end - start)


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
            # try:
            if gpcr:
                local_aid = self.get_aids_from_pubchem_orig(gpcr)
                if local_aid and local_aid != None:
                    self.process_aids(local_aid)
                aids.append(local_aid)
            # except:
            #     pass
        return aids

    #get sid and assay data_all
    def process_aids(self, aids):
        assays = list()
        # for aid in aids:
        # try:
        print('\n***New accession', aids)
        accession = list(aids.keys())[0]
        localaid = aids[accession]
        for t_aid in localaid:
            self.prosess_assay(t_aid, accession)
        # except:
        #     pass

    def prosess_assay(self,aid, accession):
        print("\n***Getting compound data for Pubchem assays:", aid )
        assay_response = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/"+str(aid)+"/json")
        if assay_response.status_code == 200:
            assay_data = assay_response.json()
            assay_quantitive = assay_data['PC_AssaySubmit']['data']
            for substance in assay_quantitive:
                self.process_assay_and_save(substance, aid, accession)


            # assay_info['assay_id'] = assay_data['PC_AssaySubmit']['assay']['descr']['aid']['id']
            # assay_info['aid_source_db'] = assay_data['PC_AssaySubmit']['assay']['descr']['aid_source']['db']['name']
            # assay_info['aid_source_db_id'] = assay_data['PC_AssaySubmit']['assay']['descr']['aid_source']['db']['source_id']['str']
            # assay_info['assay_name'] = assay_data['PC_AssaySubmit']['assay']['descr']['name']
            # temp_dois = assay_data['PC_AssaySubmit']['assay']['descr']['xref']
            # for i in temp_dois:
            #     if any('pmid' in d for d in i['xref']) == True:
            #         assay_info['DOI/pubmed']= i['xref']['pmid']

            #     assay_subs.append(temp_assay)
            # assay_info['assays'] = assay_subs


    def process_assay_and_save(self, substance,aid, accession):

        temp_assay= dict()
        temp_assay['sid'] = substance['sid']
        for i in substance['data']:
            if i['tid'] == 2:
                try:
                    temp_assay['Standard_type'] = i['value']['sval']
                except:
                    print('error apperead', i)
                    temp_assay['Standard_type'] = None
            if i['tid'] == 3:
                try:
                    temp_assay['Standard_relation'] = i['value']['sval']
                except:
                    temp_assay['Standard_relation'] = '='
            if i['tid'] == 4:
                try:
                    temp_assay['Standard_value'] = i['value']['fval']
                except:
                    print('error apperead', i)
                    temp_assay['Standard_value'] = None
            if i['tid'] == 5:
                try:
                    temp_assay['Standard_unit'] = i['value']['sval']
                except:
                    temp_assay['Standard_unit'] = 'nM'

        cid, temp_assay['assay_description'], temp_assay['reference'] = self.get_compound_data(substance['sid'],aid)
        temp_assay['ligand'] = self.get_ligand_or_create(cid)
        temp_assay['compound'] = cid
        temp_assay['protein'] = self.fetch_protein(accession)
        temp_assay['publication'] = self.fetch_publication(temp_assay['reference'])
        self.upload_to_db(temp_assay)

    def get_compound_data(self,sid,aid):
        substance_response = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/"+str(sid)+"/assaysummary/json")
        assay_description = str()
        pubmed = str()
        compound_id = str()
        if substance_response.status_code == 200:
            # TODO: try except
            substance_data = substance_response.json()
            cid = substance_data['Table']['Row']
            for i in cid:
                temp_aid = i['Cell'][0]
                temp_sid = i['Cell'][2]
                if str(aid) == str(temp_aid) and str(sid) ==str(temp_sid):
                    assay_description = i['Cell'][9]
                    pubmed = i['Cell'][11]
                    compound_id = i['Cell'][3]
                    return compound_id, assay_description, pubmed
        return compound_id, assay_description, pubmed

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
        try:
            lp.smiles = structure['smiles']
            lp.inchikey = structure['inchikey']
            lp.mw = structure['mw']
            lp.rotatable_bonds = structure['rotatable_bonds']
            lp.hacc = structure['hacc']
            lp.hdon = structure['hdon']
            lp.logp = structure['logp']
        except:
            lp.logp = 0.0
        try:
            lp.save()
            lp.web_links.add(wl)
        except IntegrityError:
            lp = LigandProperities.objects.get(inchikey=structure['inchikey'])
        return lp

    def get_aids_from_pubchem(self, accession):
        result = {accession:'259977'}
        return result

    def get_aids_from_pubchem_orig(self, accession):
        structure_response = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/accession/"+accession+"/aids/json")
        if structure_response.status_code == 200:
            try:
                ligand_data = structure_response.json()
                temp=ligand_data['IdentifierList']
                result = {accession:temp['AID']}
                return result
            except:
                return None

    #end of ligand Part
    def upload_to_db(self, i):
        # saves data
        try:
            chembl_data = AssayExperiment(ligand=i["ligand"],
                                          publication=i["publication"],
                                          protein=i["protein"],
                                          standard_value=i["Standard_value"],
                                          standard_relation=i["Standard_relation"],
                                          standard_type=i["Standard_type"],
                                          standard_units=i["Standard_unit"],
                                          assay_description=i["assay_description"],
                                          document_chembl_id=i['compound'],
                                          )
            chembl_data.save()
            self.data_point += 1
            if self.data_point%10==0:
                print('\n\ndata saved',self.data_point)
        except:
            print('data unsaved', i)


    def fetch_protein(self, target):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        try:
            test = None
            test = Protein.objects.filter(accession=target).first()
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
