from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db import IntegrityError
from django.utils.text import slugify
from django.http import HttpResponse, JsonResponse
from decimal import Decimal
from build.management.commands.base_build import Command as BaseBuild
from common.tools import fetch_from_cache, save_to_cache, fetch_from_web_api
from protein.models import Protein, ProteinCouplings
from ligand.models import  Ligand, LigandProperities, LigandType, Endogenous_GTP, BiasedData

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

# The pEC50 is defined as the negative logarithm of the EC50

MISSING_PROTEINS = {}
SKIPPED = 0
#Copying code structure from Alibek upload_excel_bias data script
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
    endogenous_data_dir = os.sep.join([settings.DATA_DIR, 'ligand_data', 'endogenous_data'])
    cell_data_dir = os.sep.join([settings.DATA_DIR, 'ligand_data', 'cell_line'])
    publication_cache = {}
    cell_cache = {}
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
        delete_bias_excel = BiasedData.objects.all() #New Model Biased Data
        delete_bias_excel.delete()

    @staticmethod
    def prepare_all_data():
        start = timeit.default_timer()
        print('**** Stage #1 : reading Excel  ****')
        bias_data = Command.read_excel_pandas(Command.structure_data_dir, 'Biased_ligand_single_pathway_data.xlsx')
        endo_data = Command.read_excel_pandas(Command.endogenous_data_dir, 'GtoP_Endogenous_Full_Data.xlsx')
        cell_data = Command.read_excel_pandas(Command.cell_data_dir, 'Cell_lines.xlsx')
        print('**** Stage #2 : parsing Excel  ****')
        df_from_excel = Command.convert_df_to_dict(bias_data)
        print('**** Stage #3 : processing data  ****')
        Command.main_process(df_from_excel)
        print('**** Stage #4 : uploading data  ****')
        stop = timeit.default_timer() #TO DEFINE
        print('Total Time:', stop - start)

    @staticmethod
    def read_excel_pandas(datadir, filename):
        source_file_path = os.sep.join([datadir, filename]).replace('//', '/')
        xls = pd.ExcelFile(source_file_path)
        df = pd.read_excel(xls, 'Data')
        return df

    @staticmethod
    def convert_df_to_dict(df):
        #Remove excel rows where no Emax value and EC50 has “=” sign and unit is “nM” (38 rows)
        df.drop(df[(df['Alt 1)\nQuantitative efficacy'].isna()) & (df['>\n<\n=\n~'] == '=') & (df['Unit'] == 'nM')].index, inplace=True)
        #cast everything to str
        df = df.astype(str)
        #cast NaN to none
        for column in df:
            df[column] = df[column].replace({'nan':None})
        #convert pandas df into list of dicts
        return_list_of_dicts = df.to_dict('records')
        return return_list_of_dicts

    @staticmethod
    def main_process(df_from_excel):
        for d in df_from_excel:
            #checking data values: float, string and low_activity checks
            d = Command.data_checks(d)
            #converting pEC50 to EC50. Normalizing to Molar
            d['Alt 1)\nQuantitative activity'], d['Measure type'], d['>\n<\n=\n~'] = Command.fetch_measurements(potency=d['Alt 1)\nQuantitative activity'], p_type= d['Measure type'], unit = d['Unit'], sign = d['>\n<\n=\n~'])
            #fetch protein name
            protein = Command.fetch_protein(d['Receptor\nUniProt entry name or code'].lower())
            #re-labeling "G protein" and "Gq/11 or Gi/o" based on data from the GProtein Couplings db
            if d['Primary effector family'] == 'G protein':
                d['Primary effector family'] = Command.fetch_g_protein(protein.id)
            #fetching publication info
            pub = Command.fetch_publication(d['Reference\nDOI or PMID'])
            #fetching the tissue and specie info from the excel sheet
            specie, tissue = Command.fetch_cell_line_info(d['Cell line'])
            #fetching ligand information
            l = Command.fetch_ligand(d['ID'], d['ID type'], d['Ligand tested for bias or func. Sel.\nName'])
            #assessing protein
            if protein == None:
                return None
            #Fetching from the endogenous excel datasheet (will need to be updated with GtoP data)
            #also need to be added some parameter to check onto (not all our data has PubChem CIDs)
            #translate SMILES into CID?
            # if d['ID type'] == 'PubChem CID':
            endogenous_status  = Command.fetch_endogenous(protein.id, l.id)

            signalling_protein = d['Primary effector subtype']
            try:
                signalling_protein = signalling_protein.strip().replace('α','a').replace('β','B').replace('g','G').lower()
            except:
                signalling_protein = None

            #assessing EC50 value (No value in case of IC50):
            EC50 = None
            if d['Measure type'] != 'IC50':
                EC50 = d['Alt 1)\nQuantitative activity']

            experiment_data= BiasedData(
                                        ligand = l,
                                        publication = pub,
                                        experiment = d['Fig./table No. '], #to be added by Kasper (TBA)
                                        endogenous_status = endogenous_status, #need to fetch from endogenous ligand browser now fetching from the excel datasheet
                                        receptor = protein,
                                        receptor_isoform = d['UniProt identifier code (isoform)'],
                                        active_receptor_complex = d['GtoP receptor name'],
                                        cell_line = d['Cell line'],
                                        tissue = tissue,
                                        specie = specie,
                                        primary_effector_family = d['Primary effector family'],
                                        primary_effector_subtype = signalling_protein,
                                        molecule_1 = d['Measured molecule 1'],
                                        molecule_2 = d['Measured molecule 2'],
                                        pathway_level = d['Pathway level'],
                                        assay_type = d['Assay type'],
                                        EC50 = EC50,
                                        EC50_sign = d['>\n<\n=\n~'],
                                        qualitative_activity=d['Alt 2)\nQualitative activity'],
                                        Emax = d['Alt 1)\nQuantitative efficacy'],
                                        Emax_sign = d['>\n<\n=\n~.1'],
                                        Tau_KA=d['Transduction Coefficient [log(τ/KA)]'],
                                        delta_Tau_KA=d['Relative Transduction Coefficient [Δlog(τ/KA)]'],
                                        )

            experiment_data.save()

    @staticmethod
    def data_checks(data):
        #floats check
        floaters = ['Alt 1)\nQuantitative activity', 'Alt 1)\nQuantitative efficacy', 'Transduction Coefficient [log(τ/KA)]', 'Relative Transduction Coefficient [Δlog(τ/KA)]']
        for key in floaters:
            try:
                data[key] = float(data[key])
            except (TypeError, ValueError):
                data[key] = None
        try: ## should be 166 rows
            if data['Alt 1)\nQuantitative activity'] < 5 and data['Measure type'] == 'pEC50' and data['Alt 1)\nQuantitative efficacy'] > 0.0:
                    data['Alt 1)\nQuantitative activity'] = 4.9
        except TypeError:
            pass
        try: ## relabeling qualitative activity when EC50 but no Emax/log(Tau/KA)
            if data['Measure type'] == 'EC50' or data['Measure type'] == 'pEC50':
                if (data['Alt 1)\nQuantitative efficacy'] ==  None) and (data['Transduction Coefficient [log(τ/KA)]'] ==  None) and (data['Alt 1)\nQuantitative activity'] is not None):
                    data['Alt 2)\nQualitative activity'] = 'No Activity'
        except TypeError:
            pass
        try: ## relabeling qualitative activity when IC50 or pIC50
            if data['Measure type'] == 'pIC50' or data['Measure type'] == 'IC50':
                    data['Alt 2)\nQualitative activity'] = 'No Activity'
        except TypeError:
            pass
        #low activity check
        try:
            if data['Alt 2)\nQualitative activity'].lower() == 'low activity': ###131 total rows
                if data['Alt 1)\nQuantitative efficacy'] == None or data['Alt 1)\nQuantitative efficacy'] == 0.0: ##9 rows
                    data['Alt 1)\nQuantitative activity'] = 4.9
                    data['Measure type'] = 'pEC50'
                    data['Alt 1)\nQuantitative efficacy'] = 20
                    data['Alt 2)\nQualitative activity'] = None
                else: ##122 rows, changing a lot
                    data['Alt 1)\nQuantitative activity'] = 4.9
                    data['Measure type'] = 'pEC50'
                    data['Alt 2)\nQualitative activity'] = None
        except AttributeError:
            pass
        #string check
        try:
            data['Unit'] = str(data['Unit'])
        except ValueError:
            pass
        return data

    @staticmethod
    def fetch_measurements(potency, p_type, unit, sign):
        if potency is not None:
            if p_type.lower()  == 'pec50':
                potency = 10**(potency*(-1))
                p_type = 'EC50'
                if sign == '<':
                    sign = '>'
                elif sign == '>':
                    sign = '<'
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
                return potency, p_type, sign
        else:
            return None, None, None

    @staticmethod
    def fetch_cell_line_info(cell_line):
        if cell_line in Command.cell_cache:
            specie = Command.cell_cache[cell_line]["Species"]
            tissue = Command.cell_cache[cell_line]["Tissue/organ"]
        else:
            specie = cell_line
            tissue = cell_line
        return specie, tissue

    @staticmethod
    def fetch_endogenous(protein, ligand):
        try:
            data = Endogenous_GTP.objects.filter(receptor=protein, ligand=ligand).values_list("endogenous_status")
            if data.count() > 0:
                return data.first()[0]
            else:
                return None
        except:
            return None


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
                protein1 = Protein.objects.filter(web_links__index=protein_from_excel, web_links__web_resource__slug='uniprot')
                test = protein1[0]
            return test
        except:
            return None

    @staticmethod
    def fetch_g_protein(id):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        try:
            #first check Bouvier data
            if ProteinCouplings.objects.filter(protein_id=id, source='Bouvier'):
                data = ProteinCouplings.objects.filter(protein_id=id, source='Bouvier').values_list("g_protein_id__name", "logmaxec50",)
                gprot = data.order_by('-logmaxec50')[0][0] #whichever value is needed
            elif ProteinCouplings.objects.filter(protein_id=id, source='GuideToPharma', transduction='primary'):
                data = ProteinCouplings.objects.filter(protein_id=id, source='GuideToPharma', transduction='primary').values_list('g_protein_id__name')
                gprot = data[0][0] #whichever value is needed
            return gprot
        except:
            gprot = 'Gq/11'
            return gprot

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
