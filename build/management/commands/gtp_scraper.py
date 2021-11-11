import requests
import statistics
import time
import json
import pandas as pd
import numpy as np
from bs4 import BeautifulSoup
from common.models import WebLink, WebResource, Publication

def get_soup(URL, id):
    r = requests.get(URL.format(id))
    soup = BeautifulSoup(r.text, "html.parser")
    return soup

def get_infos(drug_info_rows, INFO, search_class = False):
    for row in drug_info_rows:
        if search_class == False:
            for info in useful_info:
                try:
                    if info in str(row.find('a')) and info not in INFO.keys():
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

def get_pub_info(drug_id, pub_list):
    pub_ids = []
    pub_page = get_soup(pub_link, drug_id)
    pub_data = pub_page.find('table', {'class' : 'receptor_data_tables'})
    pub_rows = pub_data.findAll('tr')
    for row in pub_rows[1:]:
        if row.find('span').text.split('.')[0] in pub_list:
            pmid = row.find('a').text
            pub_ids.append(pmid)
            # pub = fetch_publication(pmid)
    return pub_ids

def fetch_publication(publication_doi):
    try:
        float(publication_doi)
        publication_doi = str(int(publication_doi))
    except ValueError:
        pass
    if publication_doi.isdigit():  # assume pubmed
        pub_type = 'pubmed'
    else:  # assume doi
        pub_type = 'doi'
    if publication_doi not in publication_cache:
        try:
            wl = WebLink.objects.get(index=publication_doi, web_resource__slug=pub_type)
        except WebLink.DoesNotExist:
            try:
                wl = WebLink.objects.create(index=publication_doi, web_resource=WebResource.objects.get(slug=pub_type))
            except IntegrityError:
                wl = WebLink.objects.get(index=publication_doi, web_resource__slug=pub_type)
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
                pass
        publication_cache[publication_doi] = pub
    else:
        pub = publication_cache[publication_doi]
    return pub

#defining globals and URLs
gpcr_gtp_ids = []
final = {}
missing_info = []
not_commented = []
publication_cache = {}
useful_info = ['PubChem SID','PubChem CID','InChIKey']
compound_classes = {'Metabolite or derivative': 'Metabolite',
                    'Natural product or derivative': 'Natural product',
                    'Endogenous peptide in human, mouse or rat': 'Peptide',
                    'Inorganic': 'Inorganic',
                    'Synthetic organic': 'Synthetic organic',
                    'Peptide or derivative': 'Peptide'}
keys_to_skip = ['Receptor', 'Comment', 'Drugs']
GtoP_endogenous = pd.DataFrame(columns=['Receptor ID', 'Receptor Name', 'Ligand ID', 'Ligand Specie', 'Compound Class',
                'PubChem CID', 'PubChem SID', 'InChIKey', 'Name', 'Target Specie', 'Type', 'Action',
                'pKi_min', 'pKi_avg', 'pKi_max', 'pEC50_min', 'pEC50_avg', 'pEC50_max',
                'pKd_min', 'pKd_avg', 'pKd_max', 'pIC50_min', 'pIC50_avg', 'pIC50_max',
                'Endogenous', 'Comment', 'Potency Ranking', 'Principal / Secondary', 'PMIDs'])

gtp_url = "https://www.guidetopharmacology.org/services/targets/families"
DRUG = 'https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?tab=biology&ligandId={}'
Summary = 'https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?&ligandId={}'
URL = 'https://www.guidetopharmacology.org/GRAC/ObjectDisplayForward?objectId={}'
interactions = "https://www.guidetopharmacology.org/services/ligands/{}/interactions"
pub_link = "https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?tab=refs&ligandId={}"

#fetching the list of GPCRs from GtoP
response = ''
while response == '':
    try:
        response = requests.get(gtp_url)
    except:
        print("Connection refused by the server..")
        time.sleep(1)
        response == ''

for entry in response.json():
    if entry['parentFamilyIds']:
        if entry['parentFamilyIds'][0] == 694 or entry['parentFamilyIds'][0] == 115:
            gpcr_gtp_ids.extend(entry['targetIds'])

#Parsing each GPCR receptor
for id in gpcr_gtp_ids:
    soup = get_soup(URL, id)
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
                        dsoup = get_soup(DRUG, drug_id)
                        try:
                            drug_data = dsoup.find('table', {'id' : 'Selectivity at GPCRs'})
                            drug_rows = drug_data.findAll('tr')
                            for k in range(len(drug_rows)):
                                if drug_rows[k].find('a') and (drug_rows[k].find('a')['href'].split('=')[1] == str(id)):
                                    target_specie = drug_rows[k].findAll('td')[3].find('a')['title']
                                    if not target_specie:
                                        target_specie = 'No Data'
                                    if target_specie not in final[id][drug_id].keys():
                                        final[id][drug_id][target_specie] = {"Target Specie": target_specie}
                                    if drug_rows[k].findAll('td')[2].find('img'):
                                        if 'endogenous' in drug_rows[k].findAll('td')[2].find('img')['alt']:
                                            final[id][drug_id][target_specie]['Endogenous'] = 'True'
                                    else:
                                        final[id][drug_id][target_specie]['Endogenous'] = 'False'
                                    pubs = drug_rows[k].findAll('td')[-2].text.replace('-',',').split(',')
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

                            ligand_specie = dsoup.findAll('div', {'class': 'textright_ligsum'})[-1].text.strip()
                            if len(ligand_specie) > 0:
                                try:
                                    ligand_specie = ligand_specie.split(u'\xa0')[1]
                                    final[id][drug_id]['Ligand Specie'] =  ligand_specie
                                except IndexError:
                                    final[id][drug_id]['Ligand Specie'] = 'Same as target'
                            else:
                                final[id][drug_id]['Ligand Specie'] = 'Same as target'
                        except AttributeError:
                            # final[id][drug_id]['Human'] = {"Name": drug}
                            print('Something went wrong on ligand: ' + str(drug) + ' , ' + str(drug_id))
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

#Implementing ligand info for each ligand
for gpcr in final.keys():
    for drug in final[gpcr]:
        if drug not in keys_to_skip:
            INFO = {}
            SummarySoup = get_soup(Summary, drug)
            drug_info = SummarySoup.findAll('table', {'class' : 'receptor_data_tables'})[-1]
            drug_class = SummarySoup.findAll('table', {'class' : 'receptor_data_tables'})[0]
            drug_info_rows = drug_info.findAll('tr')
            drug_class_rows = drug_class.findAll('tr')
            if len(drug_info_rows) <= 2:
                drug_info = SummarySoup.findAll('table', {'class' : 'receptor_data_tables'})[-2]
                drug_info_rows = drug_info.findAll('tr')
            INFO = get_infos(drug_info_rows, INFO)
            INFO = get_infos(drug_class_rows, INFO, True)
            for specie in final[gpcr][drug].keys():
                if specie not in ['Name', 'Ligand Specie']:
                    final[gpcr][drug][specie].update(INFO)

#Populating the Pandas Dataframe with all the scraped info
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
        row['Name'] = final[ID][drug]['Name']
        try:
            row['Ligand Specie'] = final[ID][drug]['Ligand Specie']
        except KeyError:
            pass
        for specie in final[ID][drug].keys():
            if specie not in row.keys():
                row['Target Specie'] = specie
                temp = {}
                for value in final[ID][drug][specie].keys():
                    if value in GtoP_endogenous.keys():
                        temp[value] = final[ID][drug][specie][value]
                to_add = {**row, **temp}
                GtoP_endogenous = GtoP_endogenous.append(to_add, ignore_index=True)
        if (len(row) > 3) and (len(row) < 6):
            GtoP_endogenous = GtoP_endogenous.append(row, ignore_index=True)

#Adding the Principal / Secondary labels where comments explicitly states principal
#and to receptors with a single reported endogenous ligand, while tracking not commented ones
IDS = list(GtoP_endogenous['Receptor ID'].unique())
for id in IDS:
    slice = GtoP_endogenous.loc[GtoP_endogenous['Receptor ID'] == id]
    comment = slice['Comment'].unique()[0].split('.')[0].lower()
    if len(slice['Name'].unique()) == 1:
        GtoP_endogenous.loc[GtoP_endogenous['Receptor ID'] == id, 'Principal / Secondary'] = 'Principal'
    if 'principal' in comment:
        if 'agonists' in comment:
            drugs = comment.replace(' and ', ', ').split(' are')[0].split(', ')
            drugs = [x.strip(',') for x in drugs]
            GtoP_endogenous.loc[(GtoP_endogenous['Receptor ID'] == id) & (GtoP_endogenous.Name.isin(drugs)), 'Principal / Secondary'] = 'Principal'
            GtoP_endogenous.loc[(GtoP_endogenous['Receptor ID'] == id) & (~GtoP_endogenous.Name.isin(drugs)), 'Principal / Secondary'] = 'Secondary'
        else:
            drugs = comment.split(' is')[0]
            GtoP_endogenous.loc[(GtoP_endogenous['Receptor ID'] == id) & (GtoP_endogenous['Name'] == drugs), 'Principal / Secondary'] = 'Principal'
            GtoP_endogenous.loc[(GtoP_endogenous['Receptor ID'] == id) & (GtoP_endogenous['Name'] != drugs), 'Principal / Secondary'] = 'Secondary'
    else:
        not_commented.append(id)

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

GtoP_endogenous_human = GtoP_endogenous.loc[GtoP_endogenous['Specie'] == 'Human']
GtoP_endogenous_mouse = GtoP_endogenous.loc[GtoP_endogenous['Specie'] == 'Mouse']
GtoP_endogenous_rat = GtoP_endogenous.loc[GtoP_endogenous['Specie'] == 'Rat']
GtoP_endogenous_monkey = GtoP_endogenous.loc[GtoP_endogenous['Specie'] == 'Monkey']
GtoP_endogenous_guinea_pig = GtoP_endogenous.loc[GtoP_endogenous['Specie'] == 'Guinea pig']


GtoP_endogenous.to_excel("GtoP_Endogenous_Full_Data.xlsx", sheet_name='Data', index=False)
GtoP_endogenous_human.to_excel("GtoP_Endogenous_Data.xlsx", sheet_name='Human', index=False)
GtoP_endogenous_mouse.to_excel("GtoP_Endogenous_Data.xlsx", sheet_name='Mouse', index=False)
GtoP_endogenous_rat.to_excel("GtoP_Endogenous_Data.xlsx", sheet_name='Rat', index=False)
GtoP_endogenous_monkey.to_excel("GtoP_Endogenous_Data.xlsx", sheet_name='Monkey', index=False)
GtoP_endogenous_guinea_pig.to_excel("GtoP_Endogenous_Data.xlsx", sheet_name='Guinea Pig', index=False)

#TODOs (easier matching with GPCRdb):
# Getting Entry Names for each receptor
# Handling formatting of the PMIDs values
# Handling string to integer value shifts
