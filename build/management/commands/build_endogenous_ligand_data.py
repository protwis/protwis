from django.db import IntegrityError
from build.management.commands.base_build import Command as BaseBuild
from protein.models import Protein, Species
from ligand.models import Endogenous_GTP, Ligand, LigandProperities, LigandType, LigandRole
from common.models import WebLink, WebResource, Publication
from bs4 import BeautifulSoup #will need to be removed
import logging
import time
import requests
import statistics
import pandas as pd
import numpy as np

#defining globals and URLs
missing_info = []
not_commented = []

#This command will build the endogenous data starting
#from an excel file. First iteration will be based on Excel file
#obtained from Guide to Pharmacologia via scraped data (scraper implemented here)
#Second iteration will be designed upon data directly provided by
#Guide to Pharmacology in csv/xls sheet to be further processed.

class Command(BaseBuild):
    mylog = logging.getLogger(__name__)
    mylog.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    file_handler = logging.FileHandler('biasDataTest.log')
    file_handler.setLevel(logging.ERROR)
    file_handler.setFormatter(formatter)
    mylog.addHandler(file_handler)
    endogenous_data_dir = os.sep.join([settings.DATA_DIR, 'ligand_data', 'endogenous_data'])

    help = 'Updates GuideToPharma data and imports it'
    publication_cache = {}
    receptor_cache = {}
    ligand_cache = {}
    ligand_info_cache = {}
    gtp_url = "https://www.guidetopharmacology.org/services/targets/families"
    DRUG = 'https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?tab=biology&ligandId={}'
    Summary = 'https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?&ligandId={}'
    URL = 'https://www.guidetopharmacology.org/GRAC/ObjectDisplayForward?objectId={}'
    pub_link = "https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?tab=refs&ligandId={}"

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
        #### START BLOCK WITH SCRAPER
        # print("\n#1 Get GPCR ids from target families")
        # target_list = self.get_gpcrs() #updated
        # print("\n#2 Get Guide to Pharmacology endogenous ligand data")
        # endogenous_data = self.scrape_data(target_list)   #updated
        # print("\n#3 Adding drug data")
        # endogenous_data = self.adding_drug_info(endogenous_data)
        # print("\n#4 Assessing principal endogenous from comments")
        # endogenous_data, to_be_ranked = self.labeling_principals(endogenous_data)
        # print("\n#5 Adding potency ranking (pEC50, pKi)")
        # endogenous_data, problematic_data = self.adding_potency_rankings(self, endogenous_data, to_be_ranked)
        #### END BLOCK WITH SCRAPER
        print('\n#1 Reading and Parsing Excel')
        endogenous_data = Command.read_and_convert_excel(endogenous_data_dir, 'GtoP_Endogenous_Full_Data.xlsx')
        print("\n#2 Clean dataframe and upload")
        self.create_model(endogenous_data)
        print('\n\n---Finished---')

    @staticmethod
    def read_and_convert_excel(datadir, filename):
        source_file_path = os.sep.join([Command.datadir, filename]).replace('//', '/')
        xls = pd.ExcelFile(source_file_path)
        df = pd.read_excel(xls, 'Data')
        df = df.astype(str)
        for column in df:
            df[column] = df[column].replace({'nan':None})
        return_list_of_dicts = df.to_dict('records')
        return return_list_of_dicts

    @staticmethod
    def get_soup(URL, id):
        r = requests.get(URL.format(id))
        soup = BeautifulSoup(r.text, "html.parser")
        return soup

    @staticmethod
    def get_infos(drug_info_rows, INFO, search_class = False):
        for row in drug_info_rows:
            if search_class == False:
                for info in useful_info:
                    try:
                        if info in str(row.find('a')) and info not in INFO.keys():
                            data = row.findAll('a')
                            if len(data) > 1:
                                info_add = []
                                for piece in data:
                                    info_add.append(piece.text)
                                INFO[info] = (' | ').join(info_add)
                            else:
                                INFO[info] = row.find('a').text
                    except TypeError:
                        pass
            else:
                try:
                    if str(row.find('a').text).strip() in compound_classes.keys():
                        INFO['Compound Class'] = compound_classes[str(row.find('a').text).strip()]
                except (TypeError, AttributeError):
                    pass
        return INFO

    @staticmethod
    def get_pub_info(drug_id, pub_list):
        pub_ids = []
        pub_page = get_soup(self.pub_link, drug_id)
        pub_data = pub_page.find('table', {'class' : 'receptor_data_tables'})
        pub_rows = pub_data.findAll('tr')
        for row in pub_rows[1:]:
            if row.find('span').text.split('.')[0] in pub_list:
                try:
                    pmid = row.find('a').text
                    pub_ids.append(pmid)
                except AttributeError:
                    pmid = row.find('td').text.strip().split('\n')[-1].strip()
                    pub_ids.append(pmid)
        pub_ids = (' | ').join(pub_ids)
        return pub_ids

    @staticmethod
    def get_gpcrs(self):
        gpcr_gtp_ids = []
        response = ''
        while response == '':
            try:
                response = requests.get(self.gtp_url)
            except:
                print("Connection refused by the server..")
                time.sleep(1)
                response == ''

        for entry in response.json():
            if entry['parentFamilyIds']:
                if entry['parentFamilyIds'][0] == 694 or entry['parentFamilyIds'][0] == 115:
                    gpcr_gtp_ids.extend(entry['targetIds'])
        return gpcr_gtp_ids

    @staticmethod
    def scrape_data(self, gpcr_gtp_ids):
        final = {}
        for id in gpcr_gtp_ids:
            soup = get_soup(self.URL, id)
            title = str(soup.title).split(' receptor')[0].strip('<title>').split(' |')[0]
            clean_title = str(soup.title.text).split(' |')[0]
            final[id] = {"Receptor": clean_title}
            tables = soup.findAll('table', { 'class' : 'receptor_data_tables' })
            #we need to find the correct table among the ones we have fetched
            to_parse = ''
            for table in tables:
                row = ''
                rows = table.findAll('tr')
                for row in rows:
                    if(row.text.find("Natural/Endogenous") > -1):
                        to_parse = table
                        break
            #now we have the actual table with all the info.
            #we need to parse THIS table and get all the other info by fetching data via links
            if to_parse != '':
                drugtable = to_parse.findAll('tr')
                for i in range(1, len(drugtable)):
                    final[id]['Comment'] = ''
                    if 'Comments' not in drugtable[i].text:
                        try:
                            #this has to be fixed
                            drug = drugtable[i].find('a').text.lower()
                            drug_ids = [x['href'].split('=')[1] for x in drugtable[i].findAll('a')]
                            for drug_id in drug_ids:
                                final[id][drug_id] = {"Name": drug}
                                dsoup = get_soup(self.DRUG, drug_id)
                                ligand_specie = dsoup.findAll('div', {'class': 'textright_ligsum'})[-1].text.strip()
                                if len(ligand_specie) > 0:
                                    try:
                                        ligand_specie = ligand_specie.split(u'\xa0')[1]
                                        final[id][drug_id]['Ligand Specie'] =  ligand_specie
                                    except IndexError:
                                        final[id][drug_id]['Ligand Specie'] = 'Same as target'
                                else:
                                    final[id][drug_id]['Ligand Specie'] = 'Same as target'
                                try:
                                    drug_data = dsoup.find('table', {'id' : 'Selectivity at GPCRs'})
                                    drug_rows = drug_data.findAll('tr')
                                    for k in range(len(drug_rows)):
                                        if drug_rows[k].find('a') and (drug_rows[k].find('a')['href'].split('=')[1] == str(id)):
                                            target_specie = drug_rows[k].findAll('td')[3].find('a')['title']
                                            if not target_specie:
                                                target_specie = 'No Specie'
                                            if target_specie not in final[id][drug_id].keys():
                                                final[id][drug_id][target_specie] = {"Target Specie": target_specie}
                                            if drug_rows[k].findAll('td')[2].find('img'):
                                                if 'endogenous' in drug_rows[k].findAll('td')[2].find('img')['alt']:
                                                    final[id][drug_id][target_specie]['Endogenous'] = 'True'
                                            else:
                                                final[id][drug_id][target_specie]['Endogenous'] = 'False'
                                            pubs = drug_rows[k].findAll('td')[-2].text
                                            if pubs != '':
                                                pubs = pubs.replace('-',',').split(',')
                                                final[id][drug_id][target_specie]['PMIDs'] = get_pub_info(drug_id, pubs)
                                            final[id][drug_id][target_specie]['Type'] = drug_rows[k].findAll('td')[4].text
                                            final[id][drug_id][target_specie]['Action'] = drug_rows[k].findAll('td')[5].text
                                            parameter = drug_rows[k].findAll('td')[7].text
                                            if '–' in drug_rows[k].findAll('td')[6].text:
                                                first = float(drug_rows[k].findAll('td')[6].text.split(' – ')[0])
                                                second = float(drug_rows[k].findAll('td')[6].text.split(' – ')[1])
                                                final[id][drug_id][target_specie][parameter+'_min'] = first
                                                final[id][drug_id][target_specie][parameter+'_max'] = second
                                                final[id][drug_id][target_specie][parameter+'_avg'] = statistics.mean([first, second])
                                            else:
                                                final[id][drug_id][target_specie][parameter+'_max'] = drug_rows[k].findAll('td')[6].text
                                except AttributeError:
                                    # final[id][drug_id]['Human'] = {"Name": drug}
                                    print('Something went wrong on ligand: ' + str(drug) + ' , ' + str(drug_id) + ' , Receptor: ' + str(id))
                                    pass
                        except AttributeError:
                            drug = drugtable[i].text.lower().split(', ')
                            final[id]['Drugs'] = []
                            if len(drug) > 1:
                                for entry in drug:
                                    final[id]['Drugs'].append(entry)
                            else:
                                final[id]['Drugs'].append(drug[0])
                            final[id]['Comment'] = "No drug link available"
                    else:
                        final[id]['Comment'] = str(drugtable[i].text).split(': ')[1]
        return final

    @staticmethod
    def adding_drug_info(self, final):
        keys_to_skip = ['Receptor', 'Comment', 'Drugs']
        useful_info = ['PubChem SID','PubChem CID','InChIKey', 'UniProtKB', 'Name', 'Ligand Specie', 'Compound Class']
        for gpcr in final.keys():
            for drug in final[gpcr]:
                if drug not in keys_to_skip:
                    if drug not in self.ligand_info_cache.keys():
                        INFO = {}
                        SummarySoup = get_soup(self.Summary, drug)
                        drug_info = SummarySoup.findAll('table', {'class' : 'receptor_data_tables'})[-1]
                        drug_class = SummarySoup.findAll('table', {'class' : 'receptor_data_tables'})[0]
                        drug_info_rows = drug_info.findAll('tr')
                        drug_class_rows = drug_class.findAll('tr')
                        if len(drug_info_rows) <= 2:
                            drug_info = SummarySoup.findAll('table', {'class' : 'receptor_data_tables'})[-2]
                            drug_info_rows = drug_info.findAll('tr')
                        INFO = get_infos(drug_info_rows, INFO)
                        INFO = get_infos(drug_class_rows, INFO, True)
                        self.ligand_info_cache[drug] = INFO
                        final[gpcr][drug].update(INFO)
                    else:
                        final[gpcr][drug].update(self.ligand_info_cache[drug])
        return final

    @staticmethod
    def generating_dataframe(self, final):
        GtoP_endogenous = pd.DataFrame(columns=
                       ['Receptor ID', 'Receptor Name', 'Ligand ID', 'UniProtKB', 'Ligand Specie', 'Compound Class',
                        'PubChem CID', 'PubChem SID', 'InChIKey', 'Name', 'Target Specie', 'Type', 'Action',
                        'pKi_min', 'pKi_avg', 'pKi_max', 'pEC50_min', 'pEC50_avg', 'pEC50_max',
                        'pKd_min', 'pKd_avg', 'pKd_max', 'pIC50_min', 'pIC50_avg', 'pIC50_max',
                        'Endogenous', 'Comment', 'Potency Ranking', 'Principal / Secondary', 'PMIDs'])
        for ID in final:
            for drug in final[ID].keys():
                row = {}
                row['Receptor ID'] = ID
                row['Receptor Name'] = final[ID]['Receptor']
                try:
                    row['Comment'] = final[ID]['Comment']
                except KeyError:
                    row['Comment'] =  ''
                if drug == 'Receptor':
                    continue
                if drug == 'Drugs':
                    for value in final[ID][drug]:
                        row['Name'] = value
                        GtoP_endogenous = GtoP_endogenous.append(row, ignore_index=True)
                    continue
                if drug == 'Comment':
                    continue
                row['Ligand ID'] = drug
                for key in useful_info:
                    try:
                        row[key] = final[ID][drug][key]
                    except KeyError:
                        pass
                for specie in final[ID][drug].keys():
                    if specie not in row.keys():
                        row['Target Specie'] = specie
                        temp = {}
                        for value in final[ID][drug][specie].keys():
                            if value in GtoP_endogenous.keys():
                                temp[value] = final[ID][drug][specie][value]
                        row = {**row, **temp}
                        GtoP_endogenous = GtoP_endogenous.append(row, ignore_index=True)
                if (len(row) > 3) and (len(row) < 13):
                    GtoP_endogenous = GtoP_endogenous.append(row, ignore_index=True)
        return GtoP_endogenous

    @staticmethod
    def labeling_principals(self, dataframe):
        IDS = list(dataframe['Receptor ID'].unique())
        not_commented = []
        for id in IDS:
            slice = dataframe.loc[dataframe['Receptor ID'] == id]
            comment = slice['Comment'].unique()[0].split('.')[0].lower()
            if len(slice['Name'].unique()) == 1:
                dataframe.loc[dataframe['Receptor ID'] == id, 'Principal / Secondary'] = 'Principal'
            if 'principal' in comment:
                if 'agonists' in comment:
                    drugs = comment.replace(' and ', ', ').split(' are')[0].split(', ')
                    drugs = [x.strip(',') for x in drugs]
                    dataframe.loc[(dataframe['Receptor ID'] == id) & (dataframe.Name.isin(drugs)), 'Principal / Secondary'] = 'Principal'
                    dataframe.loc[(dataframe['Receptor ID'] == id) & (~dataframe.Name.isin(drugs)), 'Principal / Secondary'] = 'Secondary'
                else:
                    drugs = comment.split(' is')[0]
                    dataframe.loc[(dataframe['Receptor ID'] == id) & (dataframe['Name'] == drugs), 'Principal / Secondary'] = 'Principal'
                    dataframe.loc[(dataframe['Receptor ID'] == id) & (dataframe['Name'] != drugs), 'Principal / Secondary'] = 'Secondary'
            else:
                not_commented.append(id)
        return dataframe, not_commented

    @staticmethod
    def adding_potency_rankings(self, GtoP_endogenous, not_commented):
        #fix things, drop unused values
        GtoP_endogenous.pKi_avg.fillna(GtoP_endogenous.pKi_max, inplace=True)
        GtoP_endogenous.pEC50_avg.fillna(GtoP_endogenous.pEC50_max, inplace=True)
        GtoP_endogenous.pKd_avg.fillna(GtoP_endogenous.pKd_max, inplace=True)
        GtoP_endogenous.pIC50_avg.fillna(GtoP_endogenous.pIC50_max, inplace=True)

        #Adding Potency Ranking to ligands in receptors without Principal status information
        #while tracking problematic values (missing info, symbols in data etc)
        for id in not_commented:
          slice = GtoP_endogenous.loc[GtoP_endogenous['Receptor ID'] == id]
          if len(slice['Name'].unique()) != 1:
              if slice['pEC50_avg'].isna().any() == False:
                  try:
                      #we have all pEC50 values
                      sorted_list = sorted(list(set([float(x) for x in slice['pEC50_avg'].to_list()])), reverse=True)
                      counter = 1
                      for item in sorted_list:
                          if item in slice['pEC50_avg'].to_list():
                              GtoP_endogenous.loc[(GtoP_endogenous['Receptor ID'] == id) & (GtoP_endogenous['pEC50_avg'] == item), 'Potency Ranking'] = counter
                          else:
                              GtoP_endogenous.loc[(GtoP_endogenous['Receptor ID'] == id) & (GtoP_endogenous['pEC50_avg'] == str(item)), 'Potency Ranking'] = counter
                          counter += 1
                  except ValueError:
                      missing_info.append(id)
              elif slice['pKi_avg'].isna().any() == False:
                  try:
                      #we have all pEC50 values
                      sorted_list = sorted(list(set([float(x) for x in slice['pKi_avg'].to_list()])), reverse=True)
                      counter = 1
                      for item in sorted_list:
                          if item in slice['pKi_avg'].to_list():
                              GtoP_endogenous.loc[(GtoP_endogenous['Receptor ID'] == id) & (GtoP_endogenous['pKi_avg'] == item), 'Potency Ranking'] = counter
                          else:
                              GtoP_endogenous.loc[(GtoP_endogenous['Receptor ID'] == id) & (GtoP_endogenous['pKi_avg'] == str(item)), 'Potency Ranking'] = counter
                          counter += 1
                  except ValueError:
                      missing_info.append(id)
              else:
                  #we don't have full values, grab higher pEC50 or higher pKi?
                  values_pEC50 = slice['pEC50_avg'].dropna().to_list()
                  values_pKi = slice['pKi_avg'].dropna().to_list()
                  if len(values_pEC50) > 0:
                      try:
                          #we have all pEC50 values
                          sorted_list = sorted(list(set([float(x) for x in slice['pEC50_avg'].dropna().to_list()])), reverse=True)
                          counter = 1
                          for item in sorted_list:
                              if item in slice['pEC50_avg'].to_list():
                                  GtoP_endogenous.loc[(GtoP_endogenous['Receptor ID'] == id) & (GtoP_endogenous['pEC50_avg'] == item), 'Potency Ranking'] = counter
                              else:
                                  GtoP_endogenous.loc[(GtoP_endogenous['Receptor ID'] == id) & (GtoP_endogenous['pEC50_avg'] == str(item)), 'Potency Ranking'] = counter
                              counter += 1
                          GtoP_endogenous.loc[(GtoP_endogenous['Receptor ID'] == id) & (GtoP_endogenous['pEC50_avg'].isna()), 'Potency Ranking'] = counter
                      except ValueError:
                          missing_info.append(id)
                  elif len(values_pKi) > 0:
                      try:
                          #we have all pEC50 values
                          sorted_list = sorted(list(set([float(x) for x in slice['pKi_avg'].dropna().to_list()])), reverse=True)
                          counter = 1
                          for item in sorted_list:
                              if item in slice['pKi_avg'].to_list():
                                  GtoP_endogenous.loc[(GtoP_endogenous['Receptor ID'] == id) & (GtoP_endogenous['pKi_avg'] == item), 'Potency Ranking'] = counter
                              else:
                                  GtoP_endogenous.loc[(GtoP_endogenous['Receptor ID'] == id) & (GtoP_endogenous['pKi_avg'] == str(item)), 'Potency Ranking'] = counter
                              counter += 1
                          GtoP_endogenous.loc[(GtoP_endogenous['Receptor ID'] == id) & (GtoP_endogenous['pKi_avg'].isna()), 'Potency Ranking'] = counter
                      except ValueError:
                          missing_info.append(id)
                  else:
                      missing_info.append(id)
        return GtoP_endogenous, missing_info

    @staticmethod
    def create_model(self, GtoP_endogenous):]
        types_dict = {'Inorganic': 'small-molecule',
                      'Metabolite': 'small-molecule',
                      'Natural product': 'small-molecule',
                      'Peptide': 'peptide',
                      'Synthetic organic': 'small-molecule',
                      '': 'na'}
        values = ['pKi', 'pEC50', 'pKd', 'pIC50']
        for row in GtoP_endogenous:
            # row = row.to_dict(orient='records')
            # row = row[0]
            numeric_data = {}
            if row['Receptor ID'] not in self.receptor_cache.keys():
                receptor = fetch_protein(row['Receptor ID'])
                receptor_cache[row['Receptor ID']] = receptor
            if row['Ligand ID'] not in self.ligand_cache.keys():
                ligand_id = LigandType.objects.get(slug=types_dict[row['Compound Class']])
                ligand = fetch_ligand(row['Ligand ID'], ligand_id, row['Name'])
                self.ligand_cache[row['Ligand ID']] = ligand
            try:
                role = fetch_role(row['Type'].lower(), row['Action'].lower())
            except AttributeError:
                role = None
            for v in values:
                try:
                    numeric_data[v] = (' | ').join([str(row[v+'_min']), str(row[v+'_avg']), str(row[v+'_max'])])
                except KeyError: #missing info on pEC50
                    numeric_data[v] = ''

            #Adding publications from the PMIDs section
            pmids = row['PMIDs'].split(' | ')

            #last step because it requires multiple uploads in case we have multiple species
            if len(row['Ligand Specie'].split(', ')) == 1:
                ligand_specie = fetch_specie(row['Ligand Specie'], row['Target Specie'])
                gtp_data = Endogenous_GTP(
                            ligand = ligand,
                            ligand_specie = ligand_specie,
                            ligand_action = role,
                            endogenous_status = row['Principal / Secondary'], #principal/secondary
                            potency_ranking = row['Potency Ranking'], #potency ranking
                            receptor = receptor, #link to protein model
                            pec50 = numeric_data['pEC50'],
                            pKi = numeric_data['pKi'],
                            pic50 = numeric_data['pIC50'],
                            pKd = numeric_data['pKd'],
                            publication = pubs, #manytomany field
                            )
                gtp_data.save()
                try:
                    for pmid in pmids:
                        publication = self.fetch_publication(pmid)
                        gtp_data.publication.add(publication)
                except:
                    publication= None
            else:
                species = row['Ligand Specie'].split(', ')
                for s in species:
                    specie = fetch_specie(row['Ligand Specie'], row['Target Specie'])
                    gtp_data = Endogenous_GTP(
                                ligand = ligand,
                                ligand_specie = specie,
                                ligand_action = role,
                                endogenous_status = row['Principal / Secondary'], #principal/secondary
                                potency_ranking = row['Potency Ranking'], #potency ranking
                                receptor = receptor, #link to protein model
                                pec50 = numeric_data['pEC50'],
                                pKi = numeric_data['pKi'],
                                pic50 = numeric_data['pIC50'],
                                pKd = numeric_data['pKd'],
                                publication = pubs, #manytomany field
                                )
                    gtp_data.save()
                    try:
                        for pmid in pmids:
                            publication = self.fetch_publication(reference['pmid'])
                            gtp_data.publication.add(publication)
                    except:
                        publication= None

    @staticmethod
    def fetch_protein(self, target):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        try:
            test = None
            test = list(Protein.objects.filter(web_links__index=target, web_links__web_resource__slug='gtop'))
            return test
        except:
            return None

    @staticmethod
    def fetch_ligand(self, ligand_id, ligand_type, name):
            """
            fetch ligands with Ligand model
            requires: ligand id, ligand id type, ligand name
            requires: source_file name
            """

            l = Ligand.objects.filter(properities__web_links__index=ligand_id, properities__web_links__web_resource__slug='gtoplig')
            if l.count() > 0:
                 return l.first()
            else:
                lig = Ligand()
                l = lig.load_by_gtop_id(name, ligand_id, ligand_type)
                return l

    @staticmethod
    def fetch_role(self, drug_type, drug_action):
        try:
            query = drug_type.title()
            if drug_type == 'agonist':
                if 'partial' in drug_action:
                    query == 'Agonist (partial)'
                else:
                    query == 'Agonist'
            elif drug_type == 'antagonist':
                if 'reverse' in drug_action:
                    query == 'Inverse agonist'
            elif drug_type == 'allosteric modulator':
                if drug_action == 'negative':
                    query = 'NAM'
                else:
                    query = 'PAM'
            else:
                query = 'unknown'
            lr = LigandRole.objects.get(name=query)
            return lr
        except:
            return None

    @staticmethod
    def fetch_publication(self, publication_doi):
        """
        fetch publication with Publication model
        requires: publication doi or pmid

        """
        if ("ISBN" in publication_doi) or (int(publication_doi) == 0):
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

    @staticmethod
    def fetch_specie(self, ligand_specie, target_specie):
        try:
            if ligand_specie == 'Same as target':
                specie = Species.objects.get(common_name=target_specie)
            elif ligand_specie == None:
                specie = None
            else:
                specie = Species.objects.get(common_name=ligand_specie)
            return specie
        except:
            return None

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
