from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db import IntegrityError
from django.utils.text import slugify
from django.http import HttpResponse, JsonResponse
from decimal import Decimal
from build.management.commands.base_build import Command as BaseBuild
from build.management.commands.build_ligand_functions import *

from common.tools import fetch_from_cache, save_to_cache, fetch_from_web_api, test_model_updates
from common.models import WebLink, WebResource, Publication
from protein.models import Protein, ProteinCouplings
from ligand.models import  Ligand, LigandType, Endogenous_GTP, BiasedData, BalancedLigands
from ligand.functions import OnTheFly

import logging
import math
import pandas as pd
import numpy as np
import os
import traceback
import time
import requests
import timeit
import django.apps


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
    cell_data_dir = os.sep.join([settings.DATA_DIR, 'ligand_data', 'cell_line'])
    publication_cache = {}
    cell_cache = {}
    ligand_cache = {}
    data_all = []
    #Setting the variables for the test tracking of the model upadates
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

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
                print('Started purging bias and balanced data')
                Command.purge_bias_data()
                self.tracker = {}
                test_model_updates(self.all_models, self.tracker, initialize=True)
                print('Ended purging bias and balanced data')
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
        # import the structure data
        self.prepare_all_data()

    @staticmethod
    def purge_bias_data():
        delete_bias_excel = BiasedData.objects.all() #New Model Biased Data
        delete_balanced_ligands = BalancedLigands.objects.all() #New Model Balanced Ligands
        delete_bias_excel.delete()
        delete_balanced_ligands.delete()

    def prepare_all_data(self):
        start = timeit.default_timer()
        print('**** Stage #1: reading Excel  ****')
        bias_data = Command.read_excel_pandas(Command.structure_data_dir, 'Biased_ligand_single_pathway_data.xlsx')
        cell_data = Command.read_excel_pandas(Command.cell_data_dir, 'Cell_lines.xlsx')
        print('**** Stage #2: parsing Excel  ****')
        df_from_excel = Command.convert_df_to_dict(bias_data)
        print('**** Stage #3: processing & uploading data  ****')
        Command.main_process(df_from_excel, cell_data)
        print('**** Checkin models update ****')
        test_model_updates(self.all_models, self.tracker, check=True)
        print('**** Stage #4: updating physiology biased ligands')
        Command.update_biased_columns()
        print('**** Stage #5: updating physiology biased (subtypes) ligands')
        Command.update_biased_columns(subtype=True)
        print('**** Stage #6: updating pathway preference')
        Command.update_biased_columns(pathway=True)
        print('**** Stage #7: creating Balanced Ligand database')
        Command.model_assemble()
        print('**** Checkin models update ****')
        test_model_updates(self.all_models, self.tracker, check=True)
        print('**** Stage #8: updating biased ligands using balanced reference')
        Command.update_biased_columns(balanced=True)
        print('**** Stage #9: updating biased ligands using balanced reference (subtype)')
        Command.update_biased_columns(subtype=True, balanced=True)
        stop = timeit.default_timer()
        print('Total Time:', stop - start)

    @staticmethod
    def read_excel_pandas(datadir, filename):
        source_file_path = os.sep.join([datadir, filename]).replace('//', '/')
        xls = pd.ExcelFile(source_file_path)
        df = pd.read_excel(xls, 'Data', dtype=str)
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
    def main_process(df_from_excel, cell):
        prot_dict = {}
        gprot_dict = {}
        lig_dict = {}
        for d in df_from_excel:
            #checking data values: float, string and low_activity checks
            d = Command.data_checks(d)

            #converting pEC50 to EC50. Normalizing to Molar
            d['Alt 1)\nQuantitative activity'], d['Measure type'], d['>\n<\n=\n~'] = Command.fetch_measurements(potency=d['Alt 1)\nQuantitative activity'], p_type= d['Measure type'], unit = d['Unit'], sign = d['>\n<\n=\n~'])

            #fetch protein name - check for empty lines
            if d['Receptor\nUniProt entry name or code'] == None:
                continue
            else:
                prot_code = d['Receptor\nUniProt entry name or code'].lower()
                if prot_code not in prot_dict:
                    prot_dict[prot_code] = Command.fetch_protein(prot_code)
                protein = prot_dict[prot_code]
                if protein == None:
                    continue # Skip if protein not found

            #re-labeling "G protein" and "Gq/11 or Gi/o" based on data from the GProtein Couplings db
            if d['Primary effector family'] == 'G protein':
                if protein.id not in gprot_dict:
                    gprot_dict[protein.id] = Command.fetch_g_protein(protein.id)
                d['Primary effector family'] = gprot_dict[protein.id]

            #fetching publication info
            pub = Command.fetch_publication(d['Reference\nDOI or PMID'])

            #fetching the tissue and specie info from the excel sheet
            species, tissue = Command.fetch_cell_line_info(d['Cell line'], cell)

            #fetching ligand information
            types = {"PubChem CID":"pubchem", "SMILES": "smiles", "IUPHAR/BPS Guide to pharmacology": "gtoplig"}
            if d['ID'] != None:
                if d['ID type'] != None:
                    key = d['ID'] + "|" + d['ID type']
                else:
                    key = d['ID'] + "|None"
                if key in lig_dict:
                    l = lig_dict[key]
                else:
                    ids = {}
                    if d['ID type'] in types:
                        if isinstance(d['ID'], list):
                            ids[types[d['ID type']]] = d['ID'][0]
                        else:
                            ids[types[d['ID type']]] = d['ID']
                    elif d['ID type'] == "PubChem SID":
                        # Try to resolve SID to CID
                        cid = resolve_pubchem_SID(d['ID'])
                        if cid != None:
                            ids["pubchem"] = cid

                    l = get_or_create_ligand(d['Ligand tested for bias or func. Sel.\nName'], ids, "small-molecule", False, True)
                    lig_dict[key] = l

            # What about the other ligand => use as reference?
            # if d['ID.1'] != None:
            #     ids = {}
            #     if d['ID type.1'] in types:
            #         ids[types[d['ID type.1']]] = d['ID.1']
            #     ligand = get_or_create_ligand(d['Emax reference ligand\nName'], ids)

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
                                        species = species,
                                        primary_effector_family = d['Primary effector family'],
                                        primary_effector_subtype = signalling_protein,
                                        molecule_1 = d['Measured molecule 1'],
                                        molecule_2 = d['Measured molecule 2'],
                                        measured_process = d['Measured process'],
                                        pathway_level = d['Pathway level'],
                                        assay_type = d['Assay type'],
                                        EC50 = EC50,
                                        EC50_sign = d['>\n<\n=\n~'],
                                        qualitative_activity=d['Alt 2)\nQualitative activity'],
                                        Emax = d['Alt 1)\nQuantitative efficacy'],
                                        Emax_sign = d['>\n<\n=\n~.1'],
                                        Tau_KA=d['Transduction Coefficient [log(τ/KA)]'],
                                        delta_Tau_KA=d['Relative Transduction Coefficient [Δlog(τ/KA)]'],
                                        time_resolved=d['Time resolved'],
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
                    data['Alt 2)\nQualitative activity'] = 'No activity'
        except TypeError:
            pass
        try: ## relabeling qualitative activity when IC50 or pIC50
            if data['Measure type'] == 'pIC50' or data['Measure type'] == 'IC50':
                    data['Alt 2)\nQualitative activity'] = 'No activity'
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
    def fetch_cell_line_info(cell_line, cell_data):
        if cell_line in Command.cell_cache:
            species = Command.cell_cache[cell_line]["Species"]
            tissue = Command.cell_cache[cell_line]["Tissue/organ"]
        else:
            if cell_line in list(cell_data['Cell_line_name'].unique()):
                species = cell_data.loc[cell_data['Cell_line_name'] == cell_line, 'Species'].item()
                tissue = cell_data.loc[cell_data['Cell_line_name'] == cell_line, 'Tissue/organ'].item()
                Command.cell_cache[cell_line] ={"Species": species, "Tissue/organ": tissue}
            else:
                species = cell_line
                tissue = cell_line
        return species, tissue

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
    def update_biased_columns(subtype=False, pathway=False, balanced=False):

        receptors = list(BiasedData.objects.all().values_list("receptor_id", flat=True).distinct())

        key = 'physiology_biased'

        if subtype:
            if balanced:
                key = 'pathway_subtype_biased'
            else:
                key =  'subtype_biased'
        elif pathway:
            key = 'pathway_preferred'
        elif balanced:
            key = 'pathway_biased'


        for protein in receptors:
            data = OnTheFly(protein, subtype=subtype, pathway=pathway, balanced=balanced)
            for publication in data.keys():
                for row in data[publication]:
                    if pathway:
                        if 'P1' in data[publication][row].keys():
                            BiasedData.objects.filter(ligand_id=data[publication][row]['ligand_id'],
                                                      publication_id=publication,
                                                      receptor_id=protein).update(**{key: data[publication][row]['P1']})
                        elif 'Pathway Rank' in data[publication][row].keys():
                            BiasedData.objects.filter(ligand_id=data[publication][row]['ligand_id'],
                                                      publication_id=publication,
                                                      receptor_id=protein).update(**{key: data[publication][row]['primary_effector_family']})
                        else:
                            continue

                    else:
                        try:
                            if (data[publication][row]['Bias factor'] is None) or (data[publication][row]['Bias factor'] >= 5):
                                BiasedData.objects.filter(ligand_id=data[publication][row]['ligand_id'],
                                                          publication_id=publication,
                                                          receptor_id=protein).update(**{key: data[publication][row]['P1']})
                        except (KeyError, TypeError):
                            continue

###############################################################################
############## IMPLEMENTING THE BALANCEDLIGANDS BUILD IN HERE #################
###############################################################################

    @staticmethod
    def model_assemble():
        """
        Fetch data to models
        Saves to DB
        """
        print('**** Stage #7.1 Parsing data from BiasedData model')
        receptors = list(BiasedData.objects.all().values_list("receptor_id", flat=True).distinct())
        queried_db = Command.parse_data(receptors)
        print('**** Stage #7.22 Generating the reference balanced ligands dataframe')
        balanced_db = Command.creating_balanced_database(queried_db)
        print('**** Stage #7.3 Pushing the generated dataframe to the model')
        Command.create_model(balanced_db)

    @staticmethod
    def define_ligand_pathways(master):
        ligands = {}
        for key in master:
            if master[key]['ligand_id'] not in ligands.keys():
                ligands[master[key]['ligand_id']] = []
            ligands[master[key]['ligand_id']].append(key)
        return ligands

    @staticmethod
    def assess_pathway_preferences(comparisons, tested, subtype=False):
        pathway_preference = {}
        #calculate values for ranking (or replace with qualitative activity when missing)
        for assay in comparisons.keys():
            for test in comparisons[assay]:
                if tested[test]['ligand_id'] not in pathway_preference.keys():
                    pathway_preference[tested[test]['ligand_id']] = {}
                if subtype:
                    path_label = str(tested[test]['primary_effector_family']) + ' - ' + str(tested[test]['primary_effector_subtype'])
                else:
                    path_label = tested[test]['primary_effector_family']
                Tau_KA = tested[test]['Tau_KA']
                try:
                    Emax_EC50 = round(math.log((tested[test]['Emax']/tested[test]['EC50']),10), 3)
                except (TypeError, ValueError):
                    Emax_EC50 = None
                pathway_preference[tested[test]['ligand_id']][path_label] = [Tau_KA, Emax_EC50]
        #ranking accordingly to Log(Tau/KA) or Log(Emax/EC50) (depending on how many missing values are)
        for val in pathway_preference.keys():
            temp = []
            none_tau = len([pathway_preference[val][key][0] for key in pathway_preference[val] if pathway_preference[val][key][0] is None])
            none_emax = len([pathway_preference[val][key][1] for key in pathway_preference[val] if pathway_preference[val][key][1] is None])
            if none_tau <= none_emax:
                pathway_preference[val] = list(dict(sorted(pathway_preference[val].items(), key=lambda item: -1000 if item[1][0] == None else item[1][0], reverse=True)).keys())
            else:
                pathway_preference[val] = list(dict(sorted(pathway_preference[val].items(), key=lambda item: -1000 if item[1][1] == None else item[1][1], reverse=True)).keys())
        #provide ranked keys
            for path in pathway_preference[val]:
                if subtype:
                    if path.split(' - ')[1] == 'None':
                        temp.append([key for key in tested.keys() if (tested[key]['ligand_id'] == val and tested[key]['primary_effector_family'] == path.split(' - ')[0] and tested[key]['primary_effector_subtype'] ==  None)][0])
                    else:
                        temp.append([key for key in tested.keys() if (tested[key]['ligand_id'] == val and tested[key]['primary_effector_family'] == path.split(' - ')[0] and tested[key]['primary_effector_subtype'] == path.split(' - ')[1])][0])
                else:
                    temp.append([key for key in tested.keys() if (tested[key]['ligand_id'] == val and tested[key]['primary_effector_family'] == path)][0])
            pathway_preference[val] = temp
        return pathway_preference

    @staticmethod
    def calculate_delta(comparisons, tested, subtype=False):
        ranking = Command.assess_pathway_preferences(comparisons, tested, subtype)
        #STEPS
        calculated_output = {}
        for drug in ranking.keys():
            #perform analysis only when we have multiple pathways to actually compare
            if len(ranking[drug]) > 1:
            #Set reference pathway (first in list)
                path_count = 0
                while path_count < len(ranking[drug]) - 1:
                    reference = ranking[drug][path_count]
                    for test in ranking[drug][path_count+1:]:
                        path_count +=1
                        #Match pathway levels + skip matching if Arrestin involved
                        if (tested[test]['pathway_level'] == tested[reference]['pathway_level']) or ('Arrestin' in [tested[reference]['primary_effector_family'], tested[test]['primary_effector_family']]):
                            if subtype:
                                #make a check for excluding missing subtype info,
                                #because that way we are not investigating subtype properly
                                if (tested[reference]['primary_effector_subtype'] == None) or (tested[test]['primary_effector_subtype'] == None):
                                    continue
                                Pathway_Pair = tested[reference]['primary_effector_subtype'] + ' - ' + tested[test]['primary_effector_subtype']
                            else:
                                Pathway_Pair = tested[reference]['primary_effector_family'] + ' - ' + tested[test]['primary_effector_family']

                            ID = str(drug) + '_' + Pathway_Pair
                            calculated_output[ID] = {}
                            calculated_output[ID]['receptor_id'] = tested[test]['receptor_id']
                            calculated_output[ID]['ligand_id'] = tested[test]['ligand_id']
                            calculated_output[ID]['doi'] = tested[test]['doi']
                            calculated_output[ID]['Comparison'] = Pathway_Pair
                            try:
                                delta_logtauka = round(tested[reference]['Tau_KA'] - tested[test]['Tau_KA'], 3)
                            except TypeError:
                                delta_logtauka = None
                            try:
                                tested[reference]['Log(Emax/EC50)'] = round(math.log((tested[reference]['Emax']/tested[reference]['EC50']),10), 3)
                            except (TypeError, ValueError):
                                tested[reference]['Log(Emax/EC50)'] = None
                            try:
                                tested[test]['Log(Emax/EC50)'] = round(math.log((tested[test]['Emax']/tested[test]['EC50']),10), 3)
                            except (TypeError, ValueError):
                                tested[test]['Log(Emax/EC50)'] = None
                            try:
                                delta_logemaxec50 = round(tested[reference]['Log(Emax/EC50)'] - tested[test]['Log(Emax/EC50)'], 3)
                            except TypeError:
                                delta_logemaxec50 = None

                            calculated_output[ID]['Delta Log(Tau/KA)'] = delta_logtauka
                            calculated_output[ID]['Delta Log(Emax/EC50)'] = delta_logemaxec50
                            if subtype:
                                calculated_output[ID]['subtype'] = 'YES'
                            else:
                                calculated_output[ID]['subtype'] = 'NO'

        return calculated_output

    @staticmethod
    def AdjustedFly(receptor_id, subtype=False):
        #fetching data given the receptor id
        test_data = BiasedData.objects.filter(receptor=receptor_id, Emax__gte=90)
        pub_ids = list(BiasedData.objects.filter(receptor=receptor_id).values_list("publication", flat=True).distinct())

        # Performance: first collect all publication and ligand data
        pub_objs = Publication.objects.filter(id__in=pub_ids).values_list("id", "web_link_id__index", "year", "journal_id__name", "authors")
        pub_objs_dict = {pub_obj[0]:pub_obj[1:] for pub_obj in pub_objs}

        publications = {}
        for entry in test_data:
            if entry.publication_id not in publications.keys():
                publications[entry.publication_id] = {}
            if entry.id not in publications[entry.publication_id].keys():
                if entry.publication_id in pub_objs_dict:
                    pub_data = pub_objs_dict[entry.publication_id]
                else:
                    pub_data = Publication.objects.filter(id=entry.publication_id).values_list("web_link_id__index")

                publications[entry.publication_id][entry.id] = {'doi': pub_data[0],
                                                                'receptor_id': entry.receptor_id,
                                                                'ligand_id': entry.ligand_id,
                                                                'primary_effector_family': entry.primary_effector_family,
                                                                'primary_effector_subtype': entry.primary_effector_subtype,
                                                                'pathway_level': entry.pathway_level,
                                                                'EC50': entry.EC50,
                                                                'Emax': entry.Emax,
                                                                'Tau_KA': entry.Tau_KA,}

        for pub in list(publications.keys()):
            ligands = Command.define_ligand_pathways(publications[pub])
            publications[pub] = Command.calculate_delta(ligands, publications[pub], subtype)

        return publications

    @staticmethod
    def parse_data(receptors):
        pathway_dump = pd.DataFrame()
        for protein in receptors:
            data = Command.AdjustedFly(protein)
            data_subtype = Command.AdjustedFly(protein, True)
            for publication in data.keys():
                for row in data[publication]:
                    pathway_dump = pathway_dump.append(data[publication][row], ignore_index=True)
            for publication in data_subtype.keys():
                for row in data_subtype[publication]:
                    pathway_dump = pathway_dump.append(data_subtype[publication][row], ignore_index=True)
        return pathway_dump

    @staticmethod
    def creating_balanced_database(parsed_data):
        parsed_data['preferred_metric'] = parsed_data['Delta Log(Tau/KA)'].combine_first(parsed_data['Delta Log(Emax/EC50)'])
        filtered = parsed_data.loc[~parsed_data['preferred_metric'].isnull()]

        test = filtered[filtered['preferred_metric'].between(-0.2, 0.2, inclusive="neither")]

        balanced_db = test[['receptor_id', 'ligand_id', 'doi', 'Comparison', 'Delta Log(Tau/KA)', 'Delta Log(Emax/EC50)', 'preferred_metric', 'subtype']].copy()
        balanced_db[['Pathway 1', 'Pathway 2']] = balanced_db['Comparison'].str.split(' - ', 1, expand=True)
        #this step is needed to undiscriminate order of pathways in comparisons
        #while assessing the balance reference ligand for that specific pathway pair
        balanced_db['sorter'] = ''
        for index, row in balanced_db.iterrows():
            row['sorter'] = ''.join(list(set([row['Pathway 1'], row['Pathway 2']])))

        balanced_db['Ranking'] = np.nan
        to_be_parsed = balanced_db.groupby(['receptor_id', 'sorter', 'doi']).size().reset_index().rename(columns={0:'count'})
        the_tuples = list(zip(to_be_parsed.receptor_id, to_be_parsed.sorter, to_be_parsed.doi))

        for pair in the_tuples:
            piece = balanced_db.loc[(balanced_db['receptor_id'] == pair[0]) & (balanced_db['sorter'] == pair[1]) & (balanced_db['doi'] == pair[2])]
            if len(piece) != 1:
                ranked_dict = dict(zip(piece.ligand_id, abs(piece.preferred_metric)))
                ranked_dict = {ligand: diff_value for ligand, diff_value in sorted(ranked_dict.items(), key=lambda item: item[1])}
                Rank = 1
                for key in ranked_dict.keys():
                    balanced_db.loc[(balanced_db['receptor_id'] == pair[0]) & (balanced_db['sorter'] == pair[1]) & (balanced_db['doi'] == pair[2]) & (balanced_db['ligand_id'] == key), 'Ranking'] =  Rank
                    Rank += 1
            else:
                balanced_db.loc[(balanced_db['receptor_id'] == pair[0]) & (balanced_db['sorter'] == pair[1]) & (balanced_db['doi'] == pair[2]), 'Ranking'] =  1

        return balanced_db

    @staticmethod
    def create_model(balanced_db):
        reference_balanced_db = balanced_db.loc[balanced_db['Ranking'] == 1].drop_duplicates()
        for index, row in reference_balanced_db.iterrows():
            p = Command.fetch_protein_from_id(row['receptor_id'])
            l = Command.fetch_ligand(row['ligand_id'])
            pub = Command.fetch_publication(row['doi'])
            subtype = True if row['subtype'] == 'YES' else False
            balanced_data = BalancedLigands(
                        ligand = l,   #link to ligand model
                        receptor = p, #link to protein model
                        first_pathway = row['Pathway 1'],
                        second_pathway = row['Pathway 2'],
                        delta_logEmaxEC50 = row['Delta Log(Emax/EC50)'],
                        delta_logTauKA = row['Delta Log(Tau/KA)'],
                        publication = pub,
                        subtype_balanced = subtype,
                        )
            balanced_data.save()


    @staticmethod
    def fetch_protein_from_id(protein_id):
        """
        fetch receptor with Protein model
        requires: protein id.
        Should NEVER return None given the
        logic of this build script
        """
        try:
            test = None
            if Protein.objects.filter(id=protein_id):
                protein = Protein.objects.filter(id=protein_id)
                test = protein.get()
            return test
        except:
            return None

    @staticmethod
    def fetch_ligand(ligand_id):
        """
        fetch ligands with Ligand model
        requires: ligand id.
        Should NEVER return None given the
        logic of this build script
        """
        try:
            test = None
            if Ligand.objects.filter(id=ligand_id):
                drug = Ligand.objects.filter(id=ligand_id)
                test = drug.get()
            return test
        except:
            return None
