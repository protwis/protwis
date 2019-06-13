from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db import IntegrityError
from django.utils.text import slugify


from build.management.commands.base_build import Command as BaseBuild
from common.tools import fetch_from_cache, save_to_cache, fetch_from_web_api
from residue.models import Residue
from protein.models import Protein
from ligand.models import Ligand, LigandProperities, LigandRole, LigandType
from ligand.functions import get_or_make_ligand
from common.models import WebLink, WebResource, Publication

import logging
import os
from datetime import datetime
import xlrd
import operator
import traceback
import time



class Command(BaseBuild):
    help = 'Reads bias data and imports it'

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
        parser.add_argument('--test_run', action='store_true', help='Skip this during a test run', default=False)

    logger = logging.getLogger(__name__)

    # source file directory
    structure_data_dir = os.sep.join([settings.DATA_DIR, 'ligand_data', 'bias'])

    publication_cache = {}
    ligand_cache = {}
    data_all = []

    def handle(self, *args, **options):
        if options['test_run']:
            print('Skipping in test run')
            return
        # delete any existing structure data
        if options['purge']:
            try:
                self.purge_bias_data()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
        # import the structure data
        try:
            self.logger.info('CREATING BIAS DATA')
            self.prepare_all_data(options['filename'])
            import random
            random.shuffle(self.data_all)
            self.prepare_input(options['proc'], self.data_all)
            self.logger.info('COMPLETED CREATING BIAS')

        except Exception as msg:
            print(msg)
            traceback.print_exc()
            self.logger.error(msg)

    def purge_bias_data(self):
        #Code to remove bias data
        pass

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


            d['source_file'] = source_file + "_" + str(i)

            if len(d['mutation_to'])>1 or len(d['mutation_from'])>1: #if something is off with amino acid
                continue



            if isinstance(d['ligand_id'], float): d['ligand_id'] = int(d['ligand_id'])
            if isinstance(d['mutation_pos'], float): d['mutation_pos'] = int(d['mutation_pos'])


            temp.append(d)
        return temp



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
                else:
                    self.logger.info('unknown format'.source_file)
                    continue

                self.data_all += rows
        print(len(self.data_all)," total data points")

    def main_func(self, positions, iteration,count,lock):

        missing_proteins = {}

        c = 0
        skipped = 0
        rows = self.data_all

        while count.value<len(rows):
            with lock:
                r = rows[count.value]
                count.value +=1 

            current = time.time()
            c += 1

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

            if r['reference'] not in self.publication_cache:
                try:
                    wl = WebLink.objects.get(index=r['reference'], web_resource__slug=pub_type)
                except WebLink.DoesNotExist:
                    try:
                        wl = WebLink.objects.create(index=r['reference'],
                                web_resource = WebResource.objects.get(slug=pub_type))
                    except IntegrityError:
                        wl = WebLink.objects.get(index=r['reference'], web_resource__slug=pub_type)


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
                        pub.update_from_doi(doi=r['reference'])
                    elif pub_type == 'pubmed':
                        pub.update_from_pubmed_data(index=r['reference'])
                    try:
                        pub.save()
                    except:
                        self.logger.error('error with reference ' + str(r['reference']) + ' ' + pub_type)
                        continue #if something off with publication, skip.
                self.publication_cache[r['reference']] = pub
            else:
                pub = self.publication_cache[r['reference']]


            l = None
            if str(r['ligand_name']) in self.ligand_cache:
                if r['ligand_id'] in self.ligand_cache[str(r['ligand_name'])]:
                    l = self.ligand_cache[str(r['ligand_name'])][r['ligand_id']]
            else:
                self.ligand_cache[str(r['ligand_name'])] = {}

            if not l:
                try:
                    l = get_or_make_ligand(r['ligand_id'],r['ligand_type'],str(r['ligand_name']))
                except Exception as msg:
                    print('Something errored with ligand, aborting entry of mutation',r['ligand_name'],r['ligand_type'],r['ligand_id'],r['source_file'])
                    print(msg)
                    traceback.print_exc()
                    continue
                self.ligand_cache[str(r['ligand_name'])][r['ligand_id']] = l


            protein=Protein.objects.filter(entry_name=r['protein'])
            if protein.exists():
                protein=protein.get()

            elif r['protein'] not in missing_proteins:
                # Can contain code to try to figure out what protein it is.
                pass
            else:
                missing_proteins[r['protein']] += 1
                continue

            res=Residue.objects.filter(protein_conformation__protein=protein,amino_acid=r['mutation_from'],sequence_number=r['mutation_pos']) #FIXME MAKE AA CHECK
            if res.exists():
                res=res.get()
            else:
                self.logger.error('Skipped due to no residue or mismatch AA ' + r['protein'] + ' pos:'+str(r['mutation_pos']) + ' AA:'+r['mutation_from'])
                # print('Skipped due to no residue or mismatch AA ' + r['protein'] + ' pos:'+str(r['mutation_pos']) + ' AA:'+r['mutation_from'],r['source_file'])
                skipped += 1
                continue

        self.logger.info('Parsed '+str(c)+' bias data entries. Skipped '+str(skipped))

        sorted_missing_proteins = sorted(missing_proteins.items(), key=operator.itemgetter(1),reverse=True)
        print(missing_proteins)
