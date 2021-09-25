from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db import IntegrityError
from django.utils.text import slugify
from django.http import HttpResponse, JsonResponse
from decimal import Decimal
from build.management.commands.base_build import Command as BaseBuild
from common.tools import fetch_from_cache, save_to_cache, fetch_from_web_api
from protein.models import Protein
from ligand.models import BiasedExperiment, BiasedExperimentVendors,AnalyzedExperiment, BiasedExperimentAssay, ExperimentAssayAuthors, Ligand, LigandProperities, LigandType, LigandVendorLink

from ligand.functions import get_or_make_ligand
from common.models import WebLink, WebResource, Publication
import logging
import math
import pandas as pd
import os
import traceback
import time
import requests
import timeit

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
    structure_data_dir = os.sep.join([settings.DATA_DIR, 'ligand_data', 'bias_data'])
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
                print('Started purging bias data')
                Command.purge_bias_data()
                print('Ended purging bias data')
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
        # import the structure data
        Command.prepare_all_data()

    @staticmethod
    def purge_bias_data():
        delete_bias_excel = BiasedExperiment.objects.all()
        delete_bias_excel.delete()
        delete_bias_experiment = AnalyzedExperiment.objects.all()
        delete_bias_experiment.delete()


    @staticmethod
    def prepare_all_data():
        start = timeit.default_timer()
        print('**** Stage #1 : read Excel  ****')
        df = Command.read_excel_pandas()
        print('**** Stage #2 : parse Excel  ****')
        df_from_excel = Command.convert_df_to_dict(df)
        # df_from_excel = Command.initialise_data_dict(df)
        print('**** Stage #3 : save Excel  ****')
        Command.main_process(df_from_excel)
        print('**** Stage #4 : finish Excel  ****')
        stop = timeit.default_timer()
        print('Total Time:', stop - start)

    @staticmethod
    def convert_df_to_dict(df):
        #cast everything to str
        df = df.astype(str)
        #cast NaN to none
        for column in df:
            df[column] = df[column].replace({'nan':None})
        # df = df.replace({'nan':None})
        #convert pandas df into list of dicts
        return_list_of_dicts = df.to_dict('records')
        return return_list_of_dicts

    @staticmethod
    def read_excel_pandas():
        source_file_path = os.sep.join(
            [Command.structure_data_dir, 'Biased_ligand_single_pathway_data.xlsx']).replace('//', '/')
        xls = pd.ExcelFile(source_file_path)
        df = pd.read_excel(xls, 'Data')
        return df

    @staticmethod
    def main_process(df_from_excel):
        # row_counter = 0
        for d in df_from_excel:
            # row_counter = row_counter + 1
            # if(row_counter < 434):
            # #     continue
            # if(row_counter > 2):
            #     break
            #     import pdb; pdb.set_trace()
            try:
                 d['Alt 1)\nQuantitative activity'] = float( d['Alt 1)\nQuantitative activity'])
            except:
                 d['Alt 1)\nQuantitative activity'] =  d['Alt 1)\nQuantitative activity']
            try:
                d['Alt 1)\nQuantitative efficacy'] = float(d['Alt 1)\nQuantitative efficacy'])
            except:
                d['Alt 1)\nQuantitative efficacy'] = d['Alt 1)\nQuantitative efficacy']
            try:
                d['Transduction Coefficient [log(τ/KA)]'] = float(d['Transduction Coefficient [log(τ/KA)]'])
            except:
                d['Transduction Coefficient [log(τ/KA)]'] = None
            try:
                d['Relative Transduction Coefficient [Δlog(τ/KA)]'] = float(d['Relative Transduction Coefficient [Δlog(τ/KA)]'])
            except:
                d['Relative Transduction Coefficient [Δlog(τ/KA)]'] = None
            try:
                d['Alt 1)\nQuantitative activity'] = float(d['Alt 1)\nQuantitative activity'])
            except:
                d['Alt 1)\nQuantitative activity'] = None
            try:
                d['Alt 1)\nQuantitative efficacy'] = float(d['Alt 1)\nQuantitative efficacy'])
            except:
                d['Alt 1)\nQuantitative efficacy'] = None

            try:
                if d['Alt 1)\nQuantitative activity'] < 5 and d['Measure type'] == 'pEC50' and d['Alt 1)\nQuantitative efficacy']>0.0:
                        d['Alt 1)\nQuantitative activity'] = 4.9
            except:
                d['Alt 1)\nQuantitative activity'] = d['Alt 1)\nQuantitative activity']
            try:
                if d['Alt 2)\nQualitative activity'].lower() == 'low activity':
                    if d['Alt 1)\nQuantitative efficacy'] == None or d['Alt 1)\nQuantitative efficacy']==0.0:
                        d['Alt 1)\nQuantitative activity'] = 4.9
                        d['Measure type'] = 'pEC50'
                        d['Alt 1)\nQuantitative efficacy'] = 20
                        d['Alt 2)\nQualitative activity'] = None
                    else:
                        d['Alt 1)\nQuantitative activity'] = 4.9
                        d['Measure type'] = 'pEC50'
                        d['Alt 2)\nQualitative activity'] = None
            except:
                d['Alt 2)\nQualitative activity'] = d['Alt 2)\nQualitative activity']
            try:
                d['Unit'] = str(d['Unit'])
            except:
                d['Unit'] = d['Unit']
            d['Alt 1)\nQuantitative activity'], d['Measure type'] = Command.fetch_measurements(potency=d['Alt 1)\nQuantitative activity'],
                                                                        p_type= d['Measure type'], unit = d['Unit'])
            protein = Command.fetch_protein(d['Receptor\nUniProt entry name or code'].lower())
            # family = self.define_g_family(d['Primary effector subtype'].lower(), d['assay_type'], protein )
            pub = Command.fetch_publication(d['Reference\nDOI or PMID'])
            l = Command.fetch_ligand(
                d['ID'], d['ID type'], d['Ligand tested for bias or func. Sel.\nName'])

            # fetch reference_ligand
            if (d['Emax reference ligand\nName'] is not None):
                reference_ligand = Command.fetch_ligand(
                    d['ID.1'], d['ID type.1'], d['Emax reference ligand\nName'])
            else:
                reference_ligand = None
            # fetch protein
            if protein == None:
                return None
            end_ligand  = Command.fetch_endogenous(protein)
            signalling_protein = d['Primary effector subtype']
            try:
                signalling_protein = signalling_protein.strip().replace('α','a').replace('β','B').replace('g','G').lower()
            except:
                signalling_protein = None

            experiment_entry = BiasedExperiment(submission_author=d['Data submitting group leader'],
                                                publication=pub,
                                                ligand=l,
                                                receptor=protein,
                                                auxiliary_protein = d['Auxiliary protein\nUniProt entry name or code'],
                                                endogenous_ligand = end_ligand,
                                                ligand_source_id = d['ID'],
                                                ligand_source_type = d['ID type'],
                                                receptor_isoform = d['UniProt identifier code (isoform)'],
                                                receptor_gtpo = d['GtoP receptor name']
                                                )
            # try:
            experiment_entry.save()
            Command.fetch_vendor(l,experiment_entry)
            experiment_assay = BiasedExperimentAssay(biased_experiment=experiment_entry,
                                                   signalling_protein=signalling_protein,
                                                   family = d['Primary effector family'],
                                                   cell_line=d['Cell line'],
                                                   assay_type=d['Assay type'],
                                                   molecule_1=d['Measured molecule 1'],
                                                   molecule_2=d['Measured molecule 2'],
                                                   pathway_level = d['Pathway level'],
                                                   measured_biological_process=d['Measured process'],
                                                   signal_detection_tecnique=d['Signal detection technique'],
                                                   assay_time_resolved=d['Time resolved'],
                                                   ligand_function=d['Signaling protein\nligand activity\nLigand modality'],
                                                   quantitive_measure_type=d['Measure type'],
                                                   quantitive_activity=d['Alt 1)\nQuantitative activity'],
                                                   quantitive_activity_sign=d['>\n<\n=\n~'],
                                                   quantitive_unit=d['Unit'],
                                                   qualitative_activity=d['Alt 2)\nQualitative activity'],
                                                   quantitive_efficacy=d['Alt 1)\nQuantitative efficacy'],
                                                   efficacy_measure_type=d['Measure type.1'],
                                                   efficacy_sign=d['>\n<\n=\n~.1'],
                                                   efficacy_unit=d['Unit.1'],
                                                   bias_reference=d['Endogenous and/or reference ligand'],
                                                   transduction_coef=d['Transduction Coefficient [log(τ/KA)]'],
                                                   relative_transduction_coef=d['Relative Transduction Coefficient [Δlog(τ/KA)]'],
                                                   emax_ligand_reference=reference_ligand,
                                                   )
            experiment_assay.save()
                #fetch authors
            Command.fetch_publication_authors(pub,experiment_assay)
            # return d

    def analyse_rows(self, rows):
        """
        Reads excel rows one by one
        """

        # Analyse the rows from excel and assign the right headers
        temp = []
        start = time.time()
        print('1 process/thread start')

        # pool = Pool(4)
        # pool.map(self.main_process, rows)

        for i, r in enumerate(rows, 1):
            if i%100==0:
                print(i)
            d = self.main_process(r)
            temp.append(d)

        print('1 process/thread total time: ', time.time() - start, '\n\n')
        return temp

    @staticmethod
    def fetch_publication_authors(publication, experiment_assay):
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

    @staticmethod
    def fetch_measurements(potency, p_type, unit):
        if potency is not None:
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

        if potency is not None:
            if unit:
                if p_type.lower()  == 'ec50':
                    if unit.lower() == 'nm':
                        potency = potency* 10**(-9)
                    elif unit.lower() == 'µm':
                        potency = potency* 10**(-6)
                    elif unit.lower() == 'pm':
                        potency = potency* 10**(-12)
                    elif unit.lower() == 'mm':
                        potency = potency* 10**(-3)
                if p_type.lower()  == 'ic50':
                    if unit.lower() == 'nm':
                        potency = potency* 10**(-9)
                    elif unit.lower() == 'µm':
                        potency = potency* 10**(-6)
                    elif unit.lower() == 'pm':
                        potency = potency* 10**(-12)
                    elif unit.lower() == 'mm':
                        potency = potency* 10**(-3)
                return potency,p_type
        else:
            return None, None

    @staticmethod
    def fetch_endogenous(protein):
        try:
            with connection.cursor() as cursor:
                cursor.execute("SELECT * FROM protein_endogenous_ligands WHERE protein_id =%s", [protein.pk])
                row = cursor.fetchone()
                end_ligand = Ligand.objects.filter(id=row[2])
                test = end_ligand.get()
            return test
        except:
            return None

    @staticmethod
    def fetch_vendor(ligand,experiment_entry):
        temp = ligand
        links = LigandVendorLink.objects.filter(lp=ligand.properities.id)
        # vendor_count = 0
        for x in links:
            if x.vendor.name not in ['ZINC', 'ChEMBL', 'BindingDB', 'SureChEMBL', 'eMolecules', 'MolPort', 'PubChem']:
                ligand_vendor = BiasedExperimentVendors(experiment=experiment_entry,
                                                        vendor=x)
                ligand_vendor.save()


    @staticmethod
    def fetch_protein(protein_from_excel):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        try:
            test = None
            if Protein.objects.filter(entry_name=protein_from_excel):
                protein = Protein.objects.filter(entry_name=protein_from_excel)
                test = protein.get()
            elif Protein.objects.filter(web_links__index=protein_from_excel, web_links__web_resource__slug='uniprot'):
                protein1 = Protein.objects.filter(
                    web_links__index=protein_from_excel, web_links__web_resource__slug='uniprot')
                test = protein1[0]
            return test
        except:
            return None

    @staticmethod
    def fetch_ligand(ligand_id, ligand_type, ligand_name):
        """
        fetch ligands with Ligand model
        requires: ligand id, ligand id type, ligand name
        requires: source_file name
        """
        l = None
        try:
            if ligand_id in Command.ligand_cache:
                l = Command.ligand_cache[ligand_id]
            else:
                if l == None:
                    l = get_or_make_ligand(ligand_id, ligand_type, ligand_name)
                    Command.ligand_cache[ligand_id] = l
            if l == None:
                l = Command.create_empty_ligand(ligand_name)
        except:
            web_resource = WebResource.objects.get(slug='pubchem')
            try:
                l = Ligand.objects.get(properities__web_links__web_resource=web_resource,
                properities__web_links__index=ligand_id)
            except:
                l = Command.create_empty_ligand(ligand_name)
                # print('null ligand', l)
        return l

    @staticmethod
    def fetch_publication(publication_doi):
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
        try:
            if publication_doi not in Command.publication_cache:
                pub = False
                if pub_type == 'doi':
                    pub = Publication.get_or_create_from_doi(publication_doi)
                elif pub_type == 'pubmed':
                    pub = Publication.get_or_create_from_pubmed(publication_doi)
                Command.publication_cache[publication_doi] = pub
            else:
                pub = Command.publication_cache[publication_doi]
        except:
            pub = Publication.objects.filter(web_link__index = publication_doi).first()
        return pub

    def fetch_experiment(self, publication, ligand, receptor, source):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        try:
            experiment = AnalyzedExperiment.objects.filter(
                publication=publication, ligand=ligand, receptor=receptor, source=source)
            experiment = experiment.get()
            return True
        except Exception as msg:
            experiment = None
            self.mylog.exception(
                "Experiment AnalyzedExperiment error | module: AnalyzedExperiment.")
            return False

    @staticmethod
    def create_empty_ligand(ligand_name):
        # gtoplig webresource
        lp = Command.build_ligand_properties()
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

    @staticmethod
    def build_ligand_properties():
        lp = LigandProperities()
        lt =  LigandType.objects.get(name = 'small molecule')
        lp.ligand_type = lt
        lp.smiles = None
        lp.inchikey = None
        lp.sequence= None
        lp.mw = None
        lp.rotatable_bonds = None
        lp.hacc = None
        lp.hdon = None
        lp.logp = None
        lp.save()
        return lp

    @staticmethod
    def get_ligand_name(cid):
        ligand_name = None
        ligand_name_response = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"+str(cid)+"/synonyms/json")
        if ligand_name_response.status_code == 200:
            try:
                ligand_name = ligand_name_response.json()
                ligand_name = ligand_name['InformationList']['Information'][0]['Synonym'][0]
            except:
                ligand_name = None
        return ligand_name

    @staticmethod
    def get_ligand_properties(cid):
        properties = dict()
        compound_response = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"+str(cid)+"/property/CanonicalSMILES,InChIKey,MolecularWeight,HBondDonorCount,HBondAcceptorCount,XLogP,RotatableBondCount/json")
        if compound_response.status_code == 200:
            # TODO: try except
            compound_data = compound_response.json()
            pubchem = compound_data
            if pubchem['PropertyTable']['Properties'][0]:
                try:
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
                    properties['smiles'] =  pubchem['PropertyTable']['Properties'][0]['CanonicalSMILES']
                    properties['inchikey'] =  pubchem['PropertyTable']['Properties'][0]['InChIKey']
                except:
                    properties = dict()
        return properties

    @staticmethod
    def get_ligand_or_create(cid):
        ligand_name = str()
        ligand_name = Command.get_ligand_name(cid)
        properties = Command.get_ligand_properties(cid)
        lp = Command.create_ligand_properties(cid,properties)
        ligand = Command.create_ligand(lp, ligand_name)
        return ligand

    @staticmethod
    def create_ligand_properties(cid, structure):
        web_resource = WebResource.objects.get(slug='pubchem')
        try:
            wl = WebLink.objects.get_or_create(index=cid, web_resource=web_resource)
        except IntegrityError:
            wl = WebLink.objects.get(index=cid, web_resource=web_resource)
        lp = LigandProperities()
        try:
            lt = LigandType.objects.filter(name = 'ID type')[0]
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

    @staticmethod
    def create_ligand(lp, ligand_name):
        try:
            existing_ligand = Ligand.objects.get(name=ligand_name, canonical=True)
            return existing_ligand
        except Ligand.DoesNotExist:
            try:
                ligand = Ligand()
                ligand.properities = lp
                ligand.name = ligand_name
                ligand.canonical = True
                ligand.ambigious_alias = False
                ligand.pdbe = None
                ligand.save()
            except:
                ligand = None
            return ligand
