from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db import IntegrityError
from django.utils.text import slugify
from django.http import HttpResponse, JsonResponse
from build.management.commands.base_build import Command as BaseBuild
from protein.models import Protein, Species
from ligand.models import Endogenous_GTP, Ligand, LigandType, LigandRole
from common.models import WebLink, WebResource, Publication
from build.management.commands.build_ligand_functions import *
from build.management.commands.build_all_gtp_ligands import Command as GtPLigand
import logging
import time
import requests
import statistics
import os
import pandas as pd
import numpy as np


class Command(BaseBuild):
    mylog = logging.getLogger(__name__)
    mylog.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    file_handler = logging.FileHandler('biasDataTest.log')
    file_handler.setLevel(logging.ERROR)
    file_handler.setFormatter(formatter)
    mylog.addHandler(file_handler)
    publication_cache = {}

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
                self.purge_data()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
        self.model_assemble()

# pylint: disable=R0201
    def purge_data(self):
        print("\n# Purging data")
        endo_data = Endogenous_GTP.objects.all()
        # endo_pub_data = Endogenous_GTP_publication.objects.all()
        endo_data.delete()
        # endo_pub_data.delete()
        print("\n# Old data removed")

    def model_assemble(self):
        """
        Fetch data to models
        Saves to DB
        """
        print('*** Starting *** \n')
        print('\n#1 Fetching and setting up the GTP endogenous data')
        gtp_detailed_endogenous_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/detailed_endogenous_ligands.csv", 7 * 24 * 3600)
        gtp_data = pd.read_csv(gtp_detailed_endogenous_link, dtype=str)
        gtp_interactions_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/interactions.csv", 7 * 24 * 3600)
        gtp_interactions = pd.read_csv(gtp_interactions_link, dtype=str)
        gtp_uniprot_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/GtP_to_UniProt_mapping.csv", 7 * 24 * 3600)
        gtp_uniprot = pd.read_csv(gtp_uniprot_link, dtype=str)
        iuphar_ids = GtPLigand.compare_proteins(gtp_uniprot)

        processed_data = self.data_preparation(gtp_data, gtp_interactions, iuphar_ids)

        print('\n#2 Labeling Principal and Secondary endogenous ligands')
        endogenous_data, to_be_ranked = self.labeling_principals(processed_data)

        print('\n#3 Adding potency ranking where required')
        ranked_data = self.adding_potency_rankings(endogenous_data, to_be_ranked)

        print('\n#4 Creating and filling the Endogenous_GTP model')
        endogenous_dicts = self.convert_dataframe(ranked_data)
        self.create_model(endogenous_dicts)
        print('\n *** Finished! ***')

    @staticmethod
    def convert_dataframe(df):
        df = df.astype(str)
        for column in df:
            df[column] = df[column].replace({'nan':None})
        return_list_of_dicts = df.to_dict('records')
        return return_list_of_dicts

    @staticmethod
    def data_preparation(endogenous_data, interactions, iuphar_ids):
        #Remove parameter, value and PMID is columns to have the complete dataset
        filtered_data = endogenous_data.loc[endogenous_data['Target ID'].isin(iuphar_ids)]
        uniq_rows = filtered_data.drop(columns=['Parameter','Value','PubMed IDs']).drop_duplicates()
        uniq_rows = uniq_rows.dropna(subset=['Ligand Name'])
        association = filtered_data[['Ligand ID','Target ID']].drop_duplicates().groupby('Target ID')['Ligand ID'].apply(list).to_dict()

        info_we_want = ['pKi', 'pIC50', 'pKd', 'pEC50']
        new_columns = ['pKi_min', 'pKi_avg', 'pKi_max',
                       'pKd_min', 'pKd_avg', 'pKd_max',
                       'pIC50_min', 'pIC50_avg', 'pIC50_max',
                       'pEC50_min', 'pEC50_avg', 'pEC50_max',
                       'Ligand Species', 'Ligand Action', 'Ligand Role', 'PubMed IDs']

        columns = ['Ligand ID', 'Ligand Name', 'Ligand Type', 'Ligand UniProt IDs',
                   'Ligand Ensembl Gene ID', 'Ligand Subunit ID', 'Ligand Subunit Name',
                   'Ligand Subunit UniProt IDs', 'Ligand Subunit Ensembl IDs', 'Target ID',
                   'Target Name', 'Target UniProt ID', 'Target Ensembl Gene ID',
                   'Subunit ID', 'Subunit Name', 'Subunit UniProt IDs',
                   'Subunit Ensembl IDs', 'Natural/Endogenous Ligand Comments',
                   'RankPotency', 'Interaction Species']

        uniq_rows = uniq_rows.reindex(columns=columns+new_columns)
        #remove spaces in the column names
        uniq_rows.columns = [c.replace(' ', '_') for c in uniq_rows.columns]

        for target in association.keys():
            for ligand in association[target]:
                #adding species and role info
                try:
                    role = interactions.loc[(interactions['target_id'] == target) & (interactions['ligand_id'] == ligand), 'action'].values[0]
                except IndexError:
                    role = None
                uniq_rows.loc[(uniq_rows['Target_ID'] == target) & (uniq_rows['Ligand_ID'] == ligand), 'Ligand_Role'] = role
                try:
                    action = interactions.loc[(interactions['target_id'] == target) & (interactions['ligand_id'] == ligand), 'type'].values[0]
                except IndexError:
                    action = None
                uniq_rows.loc[(uniq_rows['Target_ID'] == target) & (uniq_rows['Ligand_ID'] == ligand), 'Ligand_Action'] = action
                try:
                    species = interactions.loc[(interactions['target_id'] == target) & (interactions['ligand_id'] == ligand), 'ligand_species'].values[0]
                except IndexError:
                    species = None
                uniq_rows.loc[(uniq_rows['Target_ID'] == target) & (uniq_rows['Ligand_ID'] == ligand), 'Ligand_Species'] = species
                #fetching the parameters of interaction between receptor and ligand
                params = endogenous_data.loc[(endogenous_data['Target ID'] == target) & (endogenous_data['Ligand ID'] == ligand), 'Parameter'].to_list()
                pmids = ''.join(endogenous_data.loc[(endogenous_data['Target ID'] == target) & (endogenous_data['Ligand ID'] == ligand), 'PubMed IDs'].dropna().to_list())
                uniq_rows.loc[(uniq_rows['Target_ID'] == target) & (uniq_rows['Ligand_ID'] == ligand), 'PubMed_IDs'] = pmids
                #now parsing the data based on parameter
                for par in params:
                    #we want only pKi, pKd, pEC50 and pIC50, not nans or other weird stuff
                    if par in info_we_want:
                        species = endogenous_data.loc[(endogenous_data['Target ID'] == target) & (endogenous_data['Ligand ID'] == ligand) & (endogenous_data['Parameter'] == par)]['Interaction Species'].tolist()
                        for org in species:
                            data = endogenous_data.loc[(endogenous_data['Target ID'] == target) & (endogenous_data['Ligand ID'] == ligand) & (endogenous_data['Parameter'] == par) & (endogenous_data['Interaction Species'] == org)]['Value'].tolist()
                            if len(data) == 1:
                                if '-' not in data[0]:
                                    uniq_rows.loc[(uniq_rows['Target_ID'] == target) & (uniq_rows['Ligand_ID'] == ligand), par+'_avg'] = data[0]
                                    uniq_rows.loc[(uniq_rows['Target_ID'] == target) & (uniq_rows['Ligand_ID'] == ligand), par+'_max'] = data[0]
                                elif data[0] == '-':
                                    continue
                                else:
                                    uniq_rows.loc[(uniq_rows['Target_ID'] == target) & (uniq_rows['Ligand_ID'] == ligand), par+'_min'] = data[0].split(' - ')[0]
                                    uniq_rows.loc[(uniq_rows['Target_ID'] == target) & (uniq_rows['Ligand_ID'] == ligand), par+'_avg'] = statistics.mean([float(data[0].split(' - ')[0]), float(data[0].split(' - ')[1])])
                                    uniq_rows.loc[(uniq_rows['Target_ID'] == target) & (uniq_rows['Ligand_ID'] == ligand), par+'_max'] = data[0].split(' - ')[1]
                            else:
                                try:
                                    vals = [float(x) for x in data]
                                except ValueError:
                                    vals = [float(y) for x in data for y in x.split(' - ')]
                                uniq_rows.loc[(uniq_rows['Target_ID'] == target) & (uniq_rows['Ligand_ID'] == ligand), par+'_min'] = min(vals)
                                uniq_rows.loc[(uniq_rows['Target_ID'] == target) & (uniq_rows['Ligand_ID'] == ligand), par+'_avg'] = statistics.mean(vals)
                                uniq_rows.loc[(uniq_rows['Target_ID'] == target) & (uniq_rows['Ligand_ID'] == ligand), par+'_max'] = max(vals)
        return uniq_rows

    @staticmethod
    def labeling_principals(dataframe):
        dataframe['Principal/Secondary'] = np.nan
        IDS = list(dataframe['Target_ID'].unique())
        not_commented = []
        for id in IDS:
            slice = dataframe.loc[dataframe['Target_ID'] == id]
            try:
                comment = slice['Natural/Endogenous_Ligand_Comments'].unique()[0].split('.')[0]
            except AttributeError: #the comment is nan
                comment = ''
            if len(slice['Ligand_Name'].unique()) == 1:
                dataframe.loc[dataframe['Target_ID'] == id, 'Principal/Secondary'] = 'Principal'
            if 'principal' in comment:
                if 'agonists' in comment:
                    drugs = comment.replace(' and ', ', ').split(' are')[0].split(', ')
                    drugs = [x.strip(',') for x in drugs]
                    dataframe.loc[(dataframe['Target_ID'] == id) & (dataframe.Ligand_Name.isin(drugs)), 'Principal/Secondary'] = 'Principal'
                    dataframe.loc[(dataframe['Target_ID'] == id) & (~dataframe.Ligand_Name.isin(drugs)), 'Principal/Secondary'] = 'Secondary'
                else:
                    drugs = comment.split(' is')[0]
                    dataframe.loc[(dataframe['Target_ID'] == id) & (dataframe['Ligand_Name'] == drugs), 'Principal/Secondary'] = 'Principal'
                    dataframe.loc[(dataframe['Target_ID'] == id) & (dataframe['Ligand_Name'] != drugs), 'Principal/Secondary'] = 'Secondary'
            else:
                not_commented.append(id)
        return dataframe, not_commented

    @staticmethod
    def adding_potency_rankings(GtoP_endogenous, not_commented):
        #fix things, drop unused values
        GtoP_endogenous['Ranking'] = np.nan
        missing_info = []
        GtoP_endogenous.pKi_avg.fillna(GtoP_endogenous.pKi_max, inplace=True)
        GtoP_endogenous.pEC50_avg.fillna(GtoP_endogenous.pEC50_max, inplace=True)
        GtoP_endogenous.pKd_avg.fillna(GtoP_endogenous.pKd_max, inplace=True)
        GtoP_endogenous.pIC50_avg.fillna(GtoP_endogenous.pIC50_max, inplace=True)

        #Adding Ranking to ligands in receptors without Principal status information
        #while tracking problematic values (missing info, symbols in data etc)
        for id in not_commented:
          slice = GtoP_endogenous.loc[GtoP_endogenous['Target_ID'] == id]
          if len(slice['Ligand_Name'].unique()) != 1:
              if slice['pEC50_avg'].isna().any() == False:
                  try:
                      #we have all pEC50 values
                      sorted_list = sorted(list(set([float(x) for x in slice['pEC50_avg'].to_list()])), reverse=True)
                      counter = 1
                      for item in sorted_list:
                          if item in slice['pEC50_avg'].to_list():
                              GtoP_endogenous.loc[(GtoP_endogenous['Target_ID'] == id) & (GtoP_endogenous['pEC50_avg'] == item), 'Ranking'] = counter
                          else:
                              GtoP_endogenous.loc[(GtoP_endogenous['Target_ID'] == id) & (GtoP_endogenous['pEC50_avg'] == str(item)), 'Ranking'] = counter
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
                              GtoP_endogenous.loc[(GtoP_endogenous['Target_ID'] == id) & (GtoP_endogenous['pKi_avg'] == item), 'Ranking'] = counter
                          else:
                              GtoP_endogenous.loc[(GtoP_endogenous['Target_ID'] == id) & (GtoP_endogenous['pKi_avg'] == str(item)), 'Ranking'] = counter
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
                                  GtoP_endogenous.loc[(GtoP_endogenous['Target_ID'] == id) & (GtoP_endogenous['pEC50_avg'] == item), 'Ranking'] = counter
                              else:
                                  GtoP_endogenous.loc[(GtoP_endogenous['Target_ID'] == id) & (GtoP_endogenous['pEC50_avg'] == str(item)), 'Ranking'] = counter
                              counter += 1
                          GtoP_endogenous.loc[(GtoP_endogenous['Target_ID'] == id) & (GtoP_endogenous['pEC50_avg'].isna()), 'Ranking'] = counter
                      except ValueError:
                          missing_info.append(id)
                  elif len(values_pKi) > 0:
                      try:
                          #we have all pEC50 values
                          sorted_list = sorted(list(set([float(x) for x in slice['pKi_avg'].dropna().to_list()])), reverse=True)
                          counter = 1
                          for item in sorted_list:
                              if item in slice['pKi_avg'].to_list():
                                  GtoP_endogenous.loc[(GtoP_endogenous['Target_ID'] == id) & (GtoP_endogenous['pKi_avg'] == item), 'Ranking'] = counter
                              else:
                                  GtoP_endogenous.loc[(GtoP_endogenous['Target_ID'] == id) & (GtoP_endogenous['pKi_avg'] == str(item)), 'Ranking'] = counter
                              counter += 1
                          GtoP_endogenous.loc[(GtoP_endogenous['Target_ID'] == id) & (GtoP_endogenous['pKi_avg'].isna()), 'Ranking'] = counter
                      except ValueError:
                          missing_info.append(id)
                  else:
                      missing_info.append(id)
        return GtoP_endogenous

    @staticmethod
    def create_model(GtoP_endogenous):
        types_dict = {'Inorganic': 'small-molecule',
                      'Metabolite': 'small-molecule',
                      'Natural product': 'small-molecule',
                      'Peptide': 'peptide',
                      'Synthetic organic': 'small-molecule',
                      None: 'na'}
        values = ['pKi', 'pEC50', 'pKd', 'pIC50']
        ligands = {}
        for row in GtoP_endogenous:
            numeric_data = {}
            receptor = Command.fetch_protein(row['Target_ID'], row['Interaction_Species'])

            ligand = get_ligand_by_id("gtoplig", row['Ligand_ID'])
            if ligand != None:
                ligand.endogenous = True
                ligand.save()
            else:
                print("Ligand ", row['Ligand_ID'], "not found", row['Ligand_Name'])

            try:
                role = Command.fetch_role(row['Ligand_Action'].lower(), row['Ligand_Role'].lower())
            except AttributeError:
                role = None
            for v in values:
                #adapting average values to 2 decimal numbers when available
                try:
                    avg = "{:.2f}".format(float(row[v+'_avg']))
                except (ValueError, TypeError):
                    avg = row[v+'_avg']
                try:
                    numeric_data[v] = (' | ').join([str(row[v+'_min']), str(avg), str(row[v+'_max'])])
                except KeyError: #missing info on pEC50
                    numeric_data[v] = ''

            #Adding publications from the PMIDs section
            try:
                pmids = row['PubMed_IDs'].split('|')
            except AttributeError:
                pmids = None

            #Species check
            try:
                species = row['Ligand_Species'].split('|')
            except AttributeError:
                species = None

            try:
                potency = int(float(row['Ranking']))
            except (TypeError, ValueError):
                potency = None

            try:
                endo_status = str(row['Principal/Secondary'])
            except:
                endo_status = None

            if receptor is not None:
            #last step because it requires multiple uploads in case we have multiple species
                if species is None:
                    gtp_data = Endogenous_GTP(
                                ligand = ligand,
                                ligand_species = species,
                                ligand_action = role,
                                endogenous_status = row['Principal/Secondary'], #principal/secondary
                                potency_ranking = potency, #Ranking
                                receptor = receptor, #link to protein model
                                pec50 = numeric_data['pEC50'],
                                pKi = numeric_data['pKi'],
                                pic50 = numeric_data['pIC50'],
                                pKd = numeric_data['pKd'],
                                )
                    gtp_data.save()
                    try:
                        for pmid in pmids:
                            publication = Command.fetch_publication(pmid)
                            gtp_data.publication.add(publication)
                    except:
                        publication= None
                elif len(species) == 1:
                    ligand_species = Command.fetch_species(species[0], row['Interaction_Species'])
                    gtp_data = Endogenous_GTP(
                                ligand = ligand,
                                ligand_species = ligand_species,
                                ligand_action = role,
                                endogenous_status = row['Principal/Secondary'], #principal/secondary
                                potency_ranking = potency, #Ranking
                                receptor = receptor, #link to protein model
                                pec50 = numeric_data['pEC50'],
                                pKi = numeric_data['pKi'],
                                pic50 = numeric_data['pIC50'],
                                pKd = numeric_data['pKd'],
                                )
                    gtp_data.save()
                    try:
                        for pmid in pmids:
                            publication = Command.fetch_publication(pmid)
                            gtp_data.publication.add(publication)
                    except:
                        publication= None
                else:
                    for s in species:
                        species = Command.fetch_species(s, row['Interaction_Species'])
                        gtp_data = Endogenous_GTP(
                                    ligand = ligand,
                                    ligand_species = species,
                                    ligand_action = role,
                                    endogenous_status = row['Principal/Secondary'], #principal/secondary
                                    potency_ranking = potency, #Ranking
                                    receptor = receptor, #link to protein model
                                    pec50 = numeric_data['pEC50'],
                                    pKi = numeric_data['pKi'],
                                    pic50 = numeric_data['pIC50'],
                                    pKd = numeric_data['pKd'],
                                    )
                        gtp_data.save()
                        try:
                            for pmid in pmids:
                                publication = Command.fetch_publication(pmid)
                                gtp_data.publication.add(publication)
                        except:
                            publication= None
            else:
                print("SKIPPING", row["Ligand_ID"], row['Target_ID'], row['Interaction_Species'])

    @staticmethod
    def fetch_role(drug_type, drug_action):
        #This still need to be addressed and fixed
        conversion = {'Activator Agonist': 'Agonist',
                      'Activator Full agonist': 'Agonist',
                      'Agonist Agonist': 'Agonist',
                      'Agonist Binding': 'Agonist',
                      'Agonist Full agonist': 'Agonist',
                      'Agonist Partial agonist': 'Agonist (partial)',
                      'Agonist Unknown': 'Agonist',
                      'Allosteric modulator Inhibition': 'NAM',
                      'Allosteric modulator Negative': 'NAM',
                      'Allosteric modulator Positive': 'PAM',
                      'Allosteric modulator Potentiation': 'PAM',
                      'Antagonist Inverse agonist': 'Inverse agonist'}
        lig_function = str(drug_type) + str(drug_action)
        lr = None
        if lig_function in conversion.keys():
            query = conversion[lig_function]
            lr = LigandRole.objects.get(name=query)
        return lr

    @staticmethod
    def fetch_protein(target, species):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        try:
            if species == None or species == "None":
                # Sorting by species => human first, otherwise next species in line
                # TODO => potentially capture all species with GtP ID
                return Protein.objects.filter(web_links__index=target, web_links__web_resource__slug='gtop').order_by("species_id").first()
            else:
                prots = Protein.objects.filter(web_links__index=target, web_links__web_resource__slug='gtop', species__common_name__iexact=species)
                if prots.count() > 0:
                    return prots.first()
                else:
                    receptor_fam = list(Protein.objects.filter(web_links__index=target, web_links__web_resource__slug='gtop').values_list("family_id", flat = True))[0]
                    return Protein.objects.get(family_id=receptor_fam, species__common_name__iexact=species)
        except:
            return None

    @staticmethod
    def fetch_species(ligand_species, target_species):
        try:
            if ligand_species == 'Same as target':
                species = Species.objects.get(common_name=target_species)
            elif ligand_species == None:
                species = None
            else:
                species = Species.objects.get(common_name=ligand_species)
            return species
        except:
            return None

    @staticmethod
    def fetch_publication(publication_doi):
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

        if publication_doi not in Command.publication_cache:
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
                    self.mylog.debug(
                        "publication fetching error | module: fetch_publication. Row # is : " + str(publication_doi) + ' ' + pub_type)
                    # if something off with publication, skip.
            Command.publication_cache[publication_doi] = pub
        else:
            pub = Command.publication_cache[publication_doi]

        return pub
