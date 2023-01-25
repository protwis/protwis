from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import IntegrityError
from django.utils.text import slugify

from build.management.commands.base_build import Command as BaseBuild
from build.management.commands.build_ligand_functions import *

from mutation.models import *
from common.tools import fetch_from_web_api, test_model_updates
from residue.models import Residue
from protein.models import Protein
from ligand.models import Ligand, LigandRole, LigandType
from common.models import WebLink, WebResource, Publication

import json
import yaml
import logging
import os
import random
import re
from datetime import datetime
import math
import xlrd
import operator
import traceback
import time
import django.apps

class Command(BaseBuild):
    help = 'Reads source data and creates pdb structure records'

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
            help='Purge existing mutations records')
        parser.add_argument('--test_run', action='store_true', help='Skip this during a test run', default=False)

    logger = logging.getLogger(__name__)

    # source file directory
    structure_data_dir = os.sep.join([settings.DATA_DIR, 'mutant_data'])
    #Setting the variables for the test tracking of the model upadates
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    publication_cache = {}
    ligand_cache = {}
    ref_ligand_cache = {}
    data_all = []
    mol_types = {"PubChem CID": "pubchem",
                "ChEMBL Compound ID": "chembl_ligand",
                "IUPHAR/BPS Guide to pharmacology": "gtoplig",
                "SMILES" : "smiles",
                "FASTA sequence (peptide)": "sequence",
                "UniProt Entry Code (peptide)" : "uniprot"}

    def handle(self, *args, **options):
        if options['test_run']:
            print('Skipping in test run')
            return

        # delete any existing mutant data
        if options['purge']:
            try:
                self.purge_mutants()
                self.tracker = {}
                test_model_updates(self.all_models, self.tracker, initialize=True)
            except Exception as msg:
                print(msg)
                self.logger.error(msg)

        # import the mutant data
        try:
            self.logger.info('CREATING MUTANT DATA')
            print('Preparing data')
            self.prepare_all_data(options['filename'])
            # random.shuffle(self.data_all)
            print('Preparing input')
            self.prepare_input(options['proc'], self.data_all)
            print('Performing check')
            test_model_updates(self.all_models, self.tracker, check=True)
            self.logger.info('COMPLETED CREATING MUTANTS')

        except Exception as msg:
            print(msg)
            traceback.print_exc()
            self.logger.error(msg)

    def purge_mutants(self):
        Mutation.objects.all().delete()
        MutationRaw.objects.all().delete()
        MutationExperiment.objects.all().delete()

    def loaddatafromexcel(self,excelpath):
        workbook = xlrd.open_workbook(excelpath)
        worksheets = workbook.sheet_names()
        temp = []
        old = 0
        for worksheet_name in worksheets:
            worksheet = workbook.sheet_by_name(worksheet_name)
            try:
                if worksheet.cell_value(0, 0) == "REFERENCE \nDOI (or PMID)": #old format FIXME
                    pass
                elif worksheet.cell_value(0, 0) == "REFERENCE \nDOI or PMID": #new format
                    pass
                    old = 1
                elif worksheet.cell_value(0, 0) == "REFERENCE \nDOI or PMID (original)": #newest format
                    pass
                elif worksheet.cell_value(0, 1) == "REFERENCE \nDOI or PMID (original)": #newest format
                    pass
                else: #skip non-matching xls files
                    continue
            except:
                continue
            num_rows = worksheet.nrows - 1
            num_cells = worksheet.ncols - 1
            curr_row = 0 #skip first, otherwise -1
            while curr_row < num_rows:
                curr_row += 1
                row = worksheet.row(curr_row)
                curr_cell = -1
                temprow = []
                if not old and worksheet.cell_value(curr_row, 1) == '': #if empty reference
                    continue
                elif old and worksheet.cell_value(curr_row, 0) == '': #if empty reference
                    continue
                while curr_cell < num_cells:
                    curr_cell += 1
                    cell_type = worksheet.cell_type(curr_row, curr_cell)
                    cell_value = worksheet.cell_value(curr_row, curr_cell)

                    # fix wrong spaced cells
                    if cell_value==" ":
                        cell_value = ""

                    temprow.append(cell_value)
                temp.append(temprow)
                #if curr_row>10: break
        return [temp, old]

    def analyse_rows(self,rows,source_file, old):
        # Analyse the rows from excel and assign the right headers
        temp = []
        for i,r in enumerate(rows,1):
            d = {}
            if not old and r[9]!='':
                # if multi mutant group skip it
                self.logger.info('Skipped row due to being a multi group ' + source_file + "_" + str(i))
                continue
            if old:
                if r[6] !='':
                    continue
                d['reference'] = r[0]
                d['protein'] = r[2].replace("__","_").lower()
                d['mutation_pos'] = r[3]
                d['mutation_from'] = r[4]
                d['mutation_to'] = r[5]
                #r[6] is new double multi mutant group #FIXME FOR LATER
                d['ligand_name'] = r[7]
                d['ligand_type'] = r[8]
                d['ligand_id'] = r[9]
                d['ligand_class'] = r[10]
                #r[10] is new reference ligand #FIXME FOR LATER
                d['exp_type'] = r[12]
                d['exp_func'] = r[13]
                d['exp_wt_value'] = float(r[14]) if r[14] else 0
                d['exp_wt_unit'] = r[15]
                d['exp_mu_effect_sign'] = r[16]
                d['exp_mu_value_raw'] = float(r[17]) if r[17] else 0
                d['fold_effect'] = float(r[18]) if r[18] else 0
                d['exp_mu_effect_qual'] = r[19]
                d['exp_mu_effect_ligand_prop'] = '' #removed
                d['exp_mu_ligand_ref'] = r[11] #check if correct

                d['review'] =''
                d['submitting_group'] =''
                d['data_container'] =''
                d['data_container_number'] = ''

                d['opt_receptor_expression'] =  0
                d['opt_basal_activity'] =  0
                d['opt_gain_of_activity'] =  0
                d['opt_ligand_emax'] =  0
                d['opt_agonist'] =  0

            else:

                d['submitting_group'] = r[0]
                d['reference'] = r[1]
                d['data_container'] = r[2]
                d['data_container_number'] = r[3]
                d['review'] = r[4]
                d['protein'] = r[5].replace("__","_").lower()
                d['mutation_pos'] = r[6]
                d['mutation_from'] = r[7]
                d['mutation_to'] = r[8]
                #r[9] is new double multi mutant group #FIXME FOR LATER
                d['ligand_name'] = r[10]
                d['ligand_type'] = r[11]
                d['ligand_id'] = r[12]
                d['ligand_class'] = r[13]
                d['exp_mu_ligand_ref'] = r[14] #check if correct?
                #r[10] is new reference ligand #FIXME FOR LATER
                d['exp_type'] = r[15]
                d['exp_func'] = r[16]
                d['exp_wt_value'] = float(r[17]) if r[17] else 0
                d['exp_wt_unit'] = r[18]
                d['exp_mu_effect_sign'] = r[19]
                d['exp_mu_value_raw'] = float(r[20]) if r[20] else 0
                d['fold_effect'] = float(r[21]) if r[21] else 0
                d['exp_mu_effect_qual'] = r[22]
                d['exp_mu_effect_ligand_prop'] = '' #removed

                d['opt_receptor_expression'] = float(r[23]) if r[23] else 0
                d['opt_basal_activity'] = float(r[24]) if r[24] else 0
                d['opt_gain_of_activity'] = r[25]
                d['opt_ligand_emax'] = float(r[26]) if r[26] else 0
                d['opt_agonist'] = r[27]

            # d['opt_type'] = r[20]
            # d['opt_wt'] = float(r[21]) if r[21] else 0
            # d['opt_mu'] = float(r[23]) if r[23] else 0
            # d['opt_sign'] = r[22]
            # d['opt_percentage'] = float(r[24]) if r[24] else 0
            # d['opt_qual'] = r[25]
            # d['opt_agonist'] = r[26]
            d['source_file'] = source_file + "_" + str(i)

            if len(d['mutation_to'])>1 or len(d['mutation_from'])>1: #if something is off with amino acid
                continue



            if isinstance(d['ligand_id'], float): d['ligand_id'] = int(d['ligand_id'])
            if isinstance(d['mutation_pos'], float): d['mutation_pos'] = int(d['mutation_pos'])


            temp.append(d)
        return temp


    def insert_raw(self,r):
        obj = MutationRaw(
        reference=r['reference'],
        review=r['review'],
        submitting_group = r['submitting_group'],
        data_container = r['data_container'],
        data_container_number = r['data_container_number'],
        protein=r['protein'],
        mutation_pos=r['mutation_pos'],
        mutation_from=r['mutation_from'],
        mutation_to=r['mutation_to'],
        ligand_name=r['ligand_name'],
        ligand_idtype=r['ligand_type'],
        ligand_id=r['ligand_id'],
        ligand_class=r['ligand_class'],
        exp_type=r['exp_type'],
        exp_func=r['exp_func'],
        exp_wt_value=r['exp_wt_value'],
        exp_wt_unit=r['exp_wt_unit'],
        exp_fold_change=r['fold_effect'],
        exp_mu_effect_sign=r['exp_mu_effect_sign'],
        exp_mu_effect_value=r['exp_mu_value_raw'],
        exp_mu_effect_qual=r['exp_mu_effect_qual'],
        exp_mu_effect_ligand_prop=r['exp_mu_effect_ligand_prop'],
        exp_mu_ligand_ref=r['exp_mu_ligand_ref'],
        opt_receptor_expression=r['opt_receptor_expression'],
        opt_basal_activity=r['opt_basal_activity'],
        opt_gain_of_activity=r['opt_gain_of_activity'],
        opt_ligand_emax=r['opt_ligand_emax'],
        opt_agonist=r['opt_agonist'],
        # opt_type=r['opt_type'],
        # opt_wt=r['opt_wt'],
        # opt_mu=r['opt_mu'],
        # opt_sign=r['opt_sign'],
        # opt_percentage=r['opt_percentage'],
        # opt_qual=r['opt_qual'],
        # opt_agonist=r['opt_agonist'],
        added_by='munk',
        added_date=datetime.now()
        )

        raw_id = obj

        return raw_id


    def load_mutant_from_yaml(self,filename):
        pass

    def prepare_all_data(self, filenames):

        if not filenames:
            filenames = os.listdir(self.structure_data_dir)
        for source_file in filenames:
            source_file_path = os.sep.join([self.structure_data_dir, source_file])
            if os.path.isfile(source_file_path) and source_file[0] != '.':
                self.logger.info('Reading file {}'.format(source_file_path))
                print('Reading file {}'.format(source_file_path))
                # read the yaml file
                rows = []
                if source_file[-4:]=='xlsx' or source_file[-3:]=='xls':
                    if "~$" in source_file:
                        # ignore open excel files
                        continue
                    rows, old = self.loaddatafromexcel(source_file_path)
                    rows = self.analyse_rows(rows,source_file, old)
                elif source_file[-4:]=='yaml':
                    rows = yaml.load(open(source_file_path, 'r'), Loader=yaml.FullLoader)
                    temp = []
                    for i,r in enumerate(rows):
                        d = {}
                        d['reference'] = r['pubmed']
                        d['submitting_group'] = ''
                        d['data_container'] = ''
                        d['data_container_number'] = ''
                        d['review'] = ''
                        d['protein'] = r['entry_name'].replace("__","_").lower()
                        d['mutation_pos'] = r['seq']
                        d['mutation_from'] = r['from_res']
                        d['mutation_to'] = r['to_res']
                        d['ligand_name'] = ''
                        d['ligand_type'] = ''
                        d['ligand_id'] = ''
                        d['ligand_class'] = ''
                        d['exp_type'] = ''
                        d['exp_func'] = ''
                        d['exp_wt_value'] = 0
                        d['exp_wt_unit'] = ''
                        d['exp_mu_effect_sign'] = ''
                        d['exp_mu_value_raw'] = 0
                        d['fold_effect'] = 0
                        d['exp_mu_effect_qual'] = ''
                        d['exp_mu_effect_ligand_prop'] = ''
                        d['exp_mu_ligand_ref'] = ''
                        d['opt_receptor_expression'] = 0
                        d['opt_basal_activity'] = 0
                        d['opt_gain_of_activity'] = ''
                        d['opt_ligand_emax'] = 0
                        d['opt_agonist'] = ''
                        d['source_file'] = source_file + "_" + str(i)
                        if len(d['mutation_to'])>1 or len(d['mutation_from'])>1: #if something is off with amino acid
                            continue
                        temp.append(d)
                    rows = temp
                else:
                    self.logger.info('unknown format'.source_file)
                    continue

                self.data_all += rows
        print(len(self.data_all)," total data points")

    #def create_mutant_data(self, filenames):
    def main_func(self, positions, iteration,count,lock):
        # filenames
        # if not positions[1]:
        #     rows = self.data[positions[0]:]
        # else:
        #     rows = self.data[positions[0]:positions[1]]


        missing_proteins = {}
        mutants_for_proteins = {}
        wrong_uniport_ids = {}

        c = 0
        skipped = 0
        inserted = 0
        bulk_m = []
        bulk_r = []
        current_sheet = time.time()

        rows = self.data_all
        while count.value<len(rows):
            with lock:
                r = rows[count.value]
                count.value +=1
        # for r in rows:
            # print(r['source_file'],c)
            # PRINT IF ERRORS OCCUR
            #self.logger.info('File '+str(r['source_file'])+' number '+str(c))
            current = time.time()
            c += 1
            # if c%100==0:
            #     self.logger.info('Parsed '+str(c)+' mutant data entries')
            try:
                # publication
                try: #fix if it thinks it's float.
                    float(r['reference'])
                    r['reference'] = str(int(r['reference']))
                    float(r['review'])
                    r['review'] = str(int(r['review']))
                except ValueError:
                    pass

                if r['reference'].isdigit(): #assume pubmed
                    pub_type = 'pubmed'
                else: #assume doi
                    pub_type = 'doi'

                if r['reference'] not in self.publication_cache and len(r['reference']) > 0:
                    if pub_type == 'doi':
                        pub = Publication.get_or_create_from_doi(r['reference'])
                    elif pub_type == 'pubmed':
                        pub = Publication.get_or_create_from_pubmed(r['reference'])

                    if not pub:
                        self.logger.error('error with reference ' + str(r['reference']) + ' ' + pub_type)
                        continue #if something off with publication, skip.

                    self.publication_cache[r['reference']] = pub
                else:
                    pub = self.publication_cache[r['reference']]

                # print(r['review'],r['reference'])
                if r['review'].isdigit(): #assume pubmed
                    pub_type = 'pubmed'
                elif r['review'].startswith('http'):
                    pub_type = 'raw_link'
                    if r['review'].startswith("https://doi.org/"):
                        r['review'] = r['review'][len("https://doi.org/"):]
                        pub_type = 'doi'
                else: #assume doi
                    pub_type = 'doi'

                # print(r['review'],pub_type)
                if r['review']:
                    if r['review'] not in self.publication_cache:
                        if pub_type == "doi":
                            pub_review = Publication.get_or_create_from_doi(r['review'])
                        elif pub_type == "pubmed":
                            pub_review = Publication.get_or_create_from_pubmed(r['review'])
                        elif pub_type == "raw_link":
                            wr = WebResource.objects.get(slug=pub_type)
                            try:
                                wl, created = WebLink.objects.get_or_create(defaults={"index": r['review']}, index__iexact=r['review'], web_resource=wr)
                            except IntegrityError:
                                # Try again (paralellization)
                                wl, created = WebLink.objects.get_or_create(defaults={"index": r['review']}, index__iexact=r['review'], web_resource=wr)

                            try:
                                pub_review = Publication.objects.get(web_link=wl)
                            except Publication.DoesNotExist:
                                pub_review = Publication()
                                try:
                                    pub_review.web_link = wl
                                    pub_review.save()
                                except IntegrityError:
                                    pub_review = Publication.objects.get(web_link=wl)
                                    pub_review.save()

                        if not pub_review:
                            self.logger.error('error with review ' + str(r['review']) + ' ' + pub_type)
                            continue #if something off with publication, skip.

                        self.publication_cache[r['review']] = pub_review
                    else:
                        pub_review = self.publication_cache[r['review']]
                else:
                    pub_review = None

                l = None
                if str(r['ligand_name']) in self.ligand_cache:
                    if r['ligand_id'] in self.ligand_cache[str(r['ligand_name'])]:
                        l = self.ligand_cache[str(r['ligand_name'])][r['ligand_id']]
                else:
                    self.ligand_cache[str(r['ligand_name'])] = {}

                if not l:
                    try:
                        ids = {}
                        lig_type = "small-molecule"
                        if r['ligand_type'] in self.mol_types:
                            ids[self.mol_types[r['ligand_type']]] = r['ligand_id']
                            if self.mol_types[r['ligand_type']] in ["uniprot", "sequence"]:
                                lig_type = "peptide"
                        l = get_or_create_ligand(r['ligand_name'], ids, lig_type, True, True)
                    except Exception as msg:
                        print('Something errored with ligand, aborting entry of mutation',r['ligand_name'],r['ligand_type'],r['ligand_id'],r['source_file'])
                        print(msg)
                        traceback.print_exc()
                        continue
                    self.ligand_cache[str(r['ligand_name'])][r['ligand_id']] = l


                l_ref = None
                if str(r['exp_mu_ligand_ref']) in self.ref_ligand_cache:
                    l_ref = get_or_create_ligand(r['exp_mu_ligand_ref'], {})
                    self.ref_ligand_cache[str(r['exp_mu_ligand_ref'])] = l_ref

                protein_id = 0
                residue_id = 0

                protein=Protein.objects.filter(entry_name=r['protein'])
                if protein.exists():
                    protein=protein.get()
                    if r['protein'] in mutants_for_proteins:
                        mutants_for_proteins[r['protein']] += 1
                    else:
                        mutants_for_proteins[r['protein']] = 1

                elif r['protein'] not in missing_proteins:

                    try:
                        r['protein'] = wrong_uniport_ids[r['protein']]
                        # real_uniprot = wrong_uniport_ids[r['protein']]
                        # protein=Protein.objects.get(entry_name=r['protein'])
                        protein=Protein.objects.get(entry_name=r['protein'])
                        # print('fetched with lookup table',r['protein'])
                    except:
                        # look for it as uniprot
                        protein=Protein.objects.filter(web_links__web_resource__slug='uniprot', web_links__index=r['protein'].upper())
                        if protein.exists():
                            protein=protein.get()
                            real_uniprot = protein.entry_name
                            if r['protein'] in mutants_for_proteins:
                                mutants_for_proteins[r['protein']] += 1
                            else:
                                mutants_for_proteins[r['protein']] = 1
                        else:
                            # Try to lookup in uniprot to catch typing errors / variants in entry_name
                            url = 'http://www.uniprot.org/uniprot/$index.xml'
                            cache_dir = ['uniprot', 'id']
                            uniprot_protein = fetch_from_web_api(url, r['protein'], cache_dir, xml = True)
                            try:
                                real_uniprot = uniprot_protein.find('.//{http://uniprot.org/uniprot}name').text.lower()
                                protein=Protein.objects.get(entry_name=real_uniprot)
                            except:
                                skipped += 1
                                if r['protein'] in missing_proteins:
                                    missing_proteins[r['protein']] += 1
                                else:
                                    missing_proteins[r['protein']] = 1
                                    # print('Skipped due to no protein '+ r['protein'])
                                    self.logger.error('Skipped due to no protein '+ r['protein'])
                                continue
                        wrong_uniport_ids[r['protein']] = protein.entry_name
                        r['protein'] = real_uniprot
                else:
                    missing_proteins[r['protein']] += 1
                    continue


                res=Residue.objects.filter(protein_conformation__protein=protein,amino_acid=r['mutation_from'], sequence_number=r['mutation_pos']) #FIXME MAKE AA CHECK
                if res.exists():
                    res=res.get()
                else:
                    # In case of ada2a try again with correction for historical numbering
                    failed = True
                    if protein.entry_name in ["ada2a_human", "ada2a_pig", "ada2a_mouse"]:
                        if isinstance(r['mutation_pos'], str):
                            r['mutation_pos'] = int(r['mutation_pos'])
                        res=Residue.objects.filter(protein_conformation__protein=protein,amino_acid=r['mutation_from'], sequence_number=r['mutation_pos'] + 15)
                        if res.exists():
                            res=res.get()
                            failed = False

                    if failed:
                        self.logger.error('Skipped due to no residue or mismatch AA ' + r['protein'] + ' pos:'+str(r['mutation_pos']) + ' AA:'+r['mutation_from'])
                        # print('Skipped due to no residue or mismatch AA ' + r['protein'] + ' pos:'+str(r['mutation_pos']) + ' AA:'+r['mutation_from'],r['source_file'])
                        skipped += 1
                        continue

                if r['ligand_class']:
                    try:
                        l_role, created = LigandRole.objects.get_or_create(name=r['ligand_class'],
                            defaults={'slug': slugify(r['ligand_class'])[:50]}) # FIXME this should not be needed
                    except Exception as e:
                        if LigandRole.objects.filter(slug=slugify(r['ligand_class'])[:50]).exists():
                            l_role = LigandRole.objects.get(slug=slugify(r['ligand_class'])[:50])
                            if l_role.name == slugify(r['ligand_class'])[:50]:
                                #if name of role is same as slug, then it was created by constructs script, replace it
                                l_role.name = r['ligand_class']
                                l_role.save()
                        else:
                            print(e)
                            print("Error with",r['ligand_class'],slugify(r['ligand_class'])[:50] )
                            l_role, created = LigandRole.objects.get_or_create(slug=slugify(r['ligand_class'])[:50]) # FIXME this should not be needed
                else:
                    l_role = None

                if r['exp_type']:
                    exp_type_id, created = MutationExperimentalType.objects.get_or_create(type=r['exp_type'])
                else:
                    exp_type_id = None

                if r['exp_func']:
                    exp_func_id, created = MutationFunc.objects.get_or_create(func=r['exp_func'])
                else:
                    exp_func_id = None

                if r['exp_mu_effect_ligand_prop'] or r['exp_mu_effect_qual']:
                    exp_qual_id, created = MutationQual.objects.get_or_create(qual=r['exp_mu_effect_qual'], prop=r['exp_mu_effect_ligand_prop'])
                else:
                    exp_qual_id = None

                # if r['opt_type'] or r['opt_wt'] or r['opt_mu'] or r['opt_sign'] or r['opt_percentage'] or r['opt_qual'] or r['opt_agonist']:
                #     exp_opt_id, created =  MutationOptional.objects.get_or_create(type=r['opt_type'], wt=r['opt_wt'], mu=r['opt_mu'], sign=r['opt_sign'], percentage=r['opt_percentage'], qual=r['opt_qual'], agonist=r['opt_agonist'])
                # else:
                #     exp_opt_id = None

                try:
                    mutation, created =  Mutation.objects.get_or_create(amino_acid=r['mutation_to'],protein=protein, residue=res)
                except IntegrityError:
                    mutation = Mutation.objects.get(amino_acid=r['mutation_to'],protein=protein, residue=res)
                # logtypes = ['pEC50','pIC50','pK']

# I swapped a few things around, since lower values are better when not talking about log values or percentages.
# If it's a log value (pKi, pEC50, pIC50, pKx), FC = 10^-wt/10^-mut or FC=10^(mut-wt)
# if it's a percentage, FC = mut/wt (higher percentage is usually better)
# else, FC = wt/mut

                foldchange = 0
                typefold = ''
                if r['exp_wt_value']!=0 and r['exp_mu_value_raw']!=0: #fix for new format
                    if r['exp_type'].startswith(("Log", "log", "p")):
                    # if re.match("(" + ")|(".join(logtypes) + ")", r['exp_type']):  #-log values!
                        try:
                            # Previous calculation
                            # foldchange = round(math.pow(10,-r['exp_mu_value_raw'])/pow(10,-r['exp_wt_value']),3);
                            # If it's a log value (pKi, pEC50, pIC50, pKx), FC = 10^-wt/10^-mut or FC=10^(mut-wt)
                            tmp = r['exp_mu_value_raw'] - r['exp_wt_value']
                            foldchange = round(math.pow(10, tmp),3)
                        except:
                            print(r)
                        typefold = r['exp_type']+"_log"
                    elif r['exp_wt_unit'] == '%':
                        # if % then it's a difference case, then lower value is bad. Otherwise it's conc and lower is better
                        # foldchange = round(r['exp_wt_value']/r['exp_mu_value_raw'],3);
                        # if it's a percentage, FC = mut/wt (higher percentage is usually better)
                        foldchange = round(r['exp_mu_value_raw']/r['exp_wt_value'],3);
                    else:
                        # Previous calculation
                        # foldchange = round(r['exp_mu_value_raw']/r['exp_wt_value'],3);
                        # else, FC = wt/mut
                        foldchange = round(r['exp_wt_value']/r['exp_mu_value_raw'],3);
                        typefold = r['exp_type']+"_not_log"
                    if foldchange < 1 and foldchange != 0:
                        foldchange = -round((1/foldchange),3)
                elif r['fold_effect']!=0:
                        try:
                            if r['fold_effect']<1:
                                foldchange = -round((1/r['fold_effect']),3);
                            else:
                                foldchange = round(r['fold_effect'],3);
                        except:
                            print('FOLD ERROR',r)
                r['fold_effect'] = foldchange

                raw_experiment = self.insert_raw(r)
                # raw_experiment.save()
                bulk = MutationExperiment(
                    refs=pub,
                    review=pub_review,
                    submitting_group = r['submitting_group'],
                    data_container = r['data_container'],
                    data_container_number = r['data_container_number'],
                    protein=protein,
                    residue=res,
                    ligand=l,
                    ligand_role=l_role,
                    ligand_ref = l_ref,
                    # raw = raw_experiment, #raw_experiment, OR None
                    # optional = exp_opt_id,
                    exp_type=exp_type_id,
                    exp_func=exp_func_id,
                    exp_qual = exp_qual_id,

                    mutation=mutation,
                    wt_value=r['exp_wt_value'], #
                    wt_unit=r['exp_wt_unit'],

                    mu_value = r['exp_mu_value_raw'],
                    mu_sign = r['exp_mu_effect_sign'],
                    foldchange = foldchange,
                    opt_receptor_expression = r['opt_receptor_expression'],
                    opt_basal_activity = r['opt_basal_activity'],
                    opt_gain_of_activity = r['opt_gain_of_activity'],
                    opt_ligand_emax = r['opt_ligand_emax'],
                    opt_agonist =  r['opt_agonist'],
                    )
                # for line,val in r.items():
                #     val = str(val)
                #     if len(val)>100:
                #         print(line,"too long",val)
                # mut_id = obj.id
                bulk_r.append(raw_experiment)
                bulk_m.append(bulk)
                # try:
                #     bulk.save()
                # except Exception as e:
                #     print(e)
                #     print(r)
                #     break
                #print('saved ',r['source_file'])
                inserted += 1
                end = time.time()
                diff = round(end - current,2)

            except:
                print('error with mutation record')
                traceback.print_exc()
                exit(0)

            #print(diff)

        self.logger.info('Parsed '+str(c)+' mutant data entries. Skipped '+str(skipped))

        current = time.time()

        raws = MutationRaw.objects.bulk_create(bulk_r)
        for i,me in enumerate(bulk_m):
            me.raw = raws[i]
        MutationExperiment.objects.bulk_create(bulk_m)
        end = time.time()
        diff = round(end - current,2)
        # current_sheet
        diff_2 = round(end - current_sheet,2)
        print("overall",diff_2,"bulk",diff,len(bulk_m),"skipped",str(skipped))
        sorted_missing_proteins = sorted(missing_proteins.items(), key=operator.itemgetter(1),reverse=True)
        # print(missing_proteins)
        # sorted_mutants_for_proteins = sorted(mutants_for_proteins.items(), key=operator.itemgetter(1),reverse=True)
