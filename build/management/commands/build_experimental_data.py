from build.management.commands.base_build import Command as BaseBuild
from build.management.commands.build_ligand_functions import get_ligand_by_id, match_id_via_unichem, get_or_create_ligand, is_float
from django.conf import settings
from django.utils.text import slugify
from django.db import IntegrityError

from common.tools import get_or_create_url_cache, fetch_from_web_api
from common.models import WebLink, WebResource, Publication, PublicationJournal
from ligand.models import Ligand, LigandID, LigandType, LigandVendors, LigandVendorLink, AssayExperiment, Endogenous_GTP, LigandRole
from protein.models import Protein, Species

import math
import os
import statistics
import datamol as dm
import datetime
import pandas as pd
import numpy as np
import urllib.parse
import urllib.request

from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

class Command(BaseBuild):
    help = "Build ChEMBL, PDSP and GtoP bioactivities data"
    bulk_size = 50000
    mapper_cache = {}
    publication_cache = {}

    def add_arguments(self, parser):
        parser.add_argument("--test_run",
                            action="store_true",
                            help="Skip this during a test run",
                            default=False)
        parser.add_argument('-u', '--purge',
                            action='store_true',
                            dest='purge',
                            default=False,
                            help='Purge existing ligand records')

    def handle(self, *args, **options):
        if options["test_run"]:
            print("Skipping in test run")
            return

        if options['purge']:
            # delete any existing ChEMBL bioactivity data
            # For deleting the ligands - purse all ligand data using the GtP ligand build
            print("Started purging bioactivity and ligand data")
            self.purge_data()
            print("Ended purging data")

        #Fetching all the Guide to Pharmacology data
        print("\n\nStarted parsing Guide to Pharmacology bioactivities data")
        gtp_uniprot_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/GtP_to_UniProt_mapping.csv", 7 * 24 * 3600)
        gtp_uniprot = pd.read_csv(gtp_uniprot_link, dtype=str, header=1)
        self.normalize_gtp_headers(gtp_uniprot)
        gtp_complete_ligands_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/ligands.csv", 7 * 24 * 3600)
        gtp_complete_ligands = pd.read_csv(gtp_complete_ligands_link, dtype=str, header=1)
        self.normalize_gtp_headers(gtp_complete_ligands)
        gtp_ligand_mapping_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/ligand_id_mapping.csv", 7 * 24 * 3600)
        gtp_ligand_mapping = pd.read_csv(gtp_ligand_mapping_link, dtype=str, header=1)
        self.normalize_gtp_headers(gtp_ligand_mapping)
        gtp_interactions_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/interactions.csv", 7 * 24 * 3600)
        gtp_interactions = pd.read_csv(gtp_interactions_link, dtype=str, header=1)
        self.normalize_gtp_headers(gtp_interactions)
        gtp_detailed_endogenous_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/detailed_endogenous_ligands.csv", 7 * 24 * 3600)
        gtp_detailed_endogenous = pd.read_csv(gtp_detailed_endogenous_link, dtype=str, header = 1)
        self.normalize_gtp_headers(gtp_detailed_endogenous)
        gtp_peptides_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/peptides.csv", 7 * 24 * 3600)
        gtp_peptides = pd.read_csv(gtp_peptides_link, dtype=str, header=1)
        self.normalize_gtp_headers(gtp_peptides)

        # This gets all the info of the ligand and the interaction with the target
        iuphar_ids = self.compare_proteins(gtp_uniprot)
        bioactivity_ligands_ids = self.obtain_ligands(gtp_interactions, iuphar_ids, ['target_id','ligand_id'])
        #Now I have all the data I need
        bioactivity_data_gtp = self.get_ligands_data(bioactivity_ligands_ids, gtp_complete_ligands, gtp_ligand_mapping, ligand_interactions=gtp_interactions, target_ids=iuphar_ids)
        #Assess the assay type given info from affinity units and assay comments
        bioactivity_data_gtp = self.classify_assay(bioactivity_data_gtp)
        bioactivity_data_gtp.fillna('None', inplace=True)
        print("Ended parsing Guide to Pharmacology bioactivities data")

        print("\n\nStarted building all Guide to Pharmacology ligands")
        print('\n\nRetrieving IUPHAR ids from UniProt ids')

        print('\n\nRetrieving ALL ligands from GTP associated to GPCRs')
        endogenous_ligands_ids = self.obtain_ligands(gtp_detailed_endogenous, iuphar_ids, ['target_id','ligand_id'])
        ligand_ids = list(set(bioactivity_ligands_ids + endogenous_ligands_ids))

        print('\n\nCollating all info from GPCR related ligands in the GTP')
        ligand_data = self.get_ligands_data(ligand_ids, gtp_complete_ligands, gtp_ligand_mapping)

        print('\n\nSaving the ligands in the models')
        self.save_the_ligands_save_the_world(ligand_data, gtp_peptides)
        # Building GtP bioactivity data
        print("\n\nStarted building Guide to Pharmacology bioactivities")
        self.build_gtp_bioactivities(bioactivity_data_gtp)
        print("Ended building Guide to Pharmacology bioactivities")

        print("\n\nStarted building the Endogenous data from Guide to Pharmacology")

        print('\n#1 Preprocessing the data')
        processed_data = self.data_preparation(gtp_detailed_endogenous, gtp_interactions, iuphar_ids)
        print('\n#2 Labeling Principal and Secondary endogenous ligands')
        endogenous_data, to_be_ranked = self.labeling_principals(processed_data)
        print('\n#3 Adding potency ranking where required')
        ranked_data = self.adding_potency_rankings(endogenous_data, to_be_ranked)
        print('\n#4 Creating and filling the Endogenous_GTP model')
        endogenous_dicts = self.convert_dataframe(ranked_data)
        self.create_model(endogenous_dicts)

        print("\n\Ended building endogenous data")

        print("\n\nStarted building ChEBML ligands")
        self.build_chembl_ligands()
        print("\n\nEnded building ChEMBL ligands")
        # Parse ChEMBL bioactivity data
        print("\n\nStarted building ChEMBL bioactivities")
        self.build_chembl_bioactivities()
        print("Ended building ChEMBL bioactivities")
        # Parse ChEMBL/PubChem vendor data
        print("\n\nStarted building PubChem vendor data")
        self.build_pubchem_vendor_links()
        print("Ended building PubChem vendor data")

        # Building PDSP KiDatabase bioactivity data
        print("\n\nStarted building PDSP KiDatabase bioactivities")
        self.build_kidatabase_bioactivities() #14,562
        print("Ended building PDSP KiDatabase bioactivities")


    @staticmethod
    def purge_data():
        delete_experimental = AssayExperiment.objects.all() #New Model Biased Data
        delete_experimental.delete()
        endo_data = Endogenous_GTP.objects.all()
        endo_data.delete()
        Ligand.objects.all().delete()


    @staticmethod
    def data_preparation(endogenous_data, interactions, iuphar_ids):
        #Remove parameter, value and PMID is columns to have the complete dataset
        filtered_data = endogenous_data.loc[endogenous_data['target_id'].isin(iuphar_ids)]
        uniq_rows = filtered_data.drop(columns=['parameter','value','pubmed_ids']).drop_duplicates()
        uniq_rows = uniq_rows.dropna(subset=['ligand_name'])
        association = filtered_data[['ligand_id','target_id']].drop_duplicates().groupby('target_id')['ligand_id'].apply(list).to_dict()

        info_we_want = ['pki', 'pic50', 'pkd', 'pec50']
        new_columns = ['pki_min', 'pki_avg', 'pki_max',
                       'pkd_min', 'pkd_avg', 'pkd_max',
                       'pic50_min', 'pic50_avg', 'pic50_max',
                       'pec50_min', 'pec50_avg', 'pec50_max',
                       'ligand_species', 'ligand_action', 'ligand_role', 'pubmed_ids']

        columns = ['ligand_id', 'ligand_name', 'ligand_type', 'ligand_uniprot_ids',
                   'ligand_ensembl_gene_id', 'ligand_subunit_id', 'ligand_subunit_name',
                   'ligand_subunit_uniprot_ids', 'ligand_subunit_ensembl_ids', 'target_id',
                   'target_name', 'target_uniprot_id', 'target_ensembl_gene_id',
                   'subunit_id', 'subunit_name', 'subunit_uniprot_ids',
                   'subunit_ensembl_ids', 'natural/endogenous_ligand_comments',
                   'rankpotency', 'interaction_species']

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
                uniq_rows.loc[(uniq_rows['target_id'] == target) & (uniq_rows['ligand_id'] == ligand), 'ligand_role'] = role
                try:
                    action = interactions.loc[(interactions['target_id'] == target) & (interactions['ligand_id'] == ligand), 'type'].values[0]
                except IndexError:
                    action = None
                uniq_rows.loc[(uniq_rows['target_id'] == target) & (uniq_rows['ligand_id'] == ligand), 'ligand_action'] = action
                try:
                    species = interactions.loc[(interactions['target_id'] == target) & (interactions['ligand_id'] == ligand), 'ligand_species'].values[0]
                except IndexError:
                    species = None
                uniq_rows.loc[(uniq_rows['target_id'] == target) & (uniq_rows['ligand_id'] == ligand), 'ligand_species'] = species
                #fetching the parameters of interaction between receptor and ligand
                params = endogenous_data.loc[(endogenous_data['target_id'] == target) & (endogenous_data['ligand_id'] == ligand), 'parameter'].to_list()
                pmids = '|'.join(endogenous_data.loc[(endogenous_data['target_id'] == target) & (endogenous_data['ligand_id'] == ligand), 'pubmed_ids'].dropna().to_list())
                uniq_rows.loc[(uniq_rows['target_id'] == target) & (uniq_rows['ligand_id'] == ligand), 'pubmed_ids'] = pmids
                #now parsing the data based on parameter
                for par in params:
                    #we want only pKi, pKd, pEC50 and pIC50, not nans or other weird stuff
                    if par in info_we_want:
                        species = endogenous_data.loc[(endogenous_data['target_id'] == target) & (endogenous_data['ligand_id'] == ligand) & (endogenous_data['parameter'] == par)]['interaction_species'].tolist()
                        for org in species:
                            data = endogenous_data.loc[(endogenous_data['target_id'] == target) & (endogenous_data['ligand_id'] == ligand) & (endogenous_data['parameter'] == par) & (endogenous_data['interaction_species'] == org)]['value'].tolist()
                            if len(data) == 1:
                                if '-' not in data[0]:
                                    uniq_rows.loc[(uniq_rows['target_id'] == target) & (uniq_rows['ligand_id'] == ligand), par+'_avg'] = data[0]
                                    uniq_rows.loc[(uniq_rows['target_id'] == target) & (uniq_rows['ligand_id'] == ligand), par+'_max'] = data[0]
                                elif data[0] == '-':
                                    continue
                                else:
                                    uniq_rows.loc[(uniq_rows['target_id'] == target) & (uniq_rows['ligand_id'] == ligand), par+'_min'] = data[0].split(' - ')[0]
                                    uniq_rows.loc[(uniq_rows['target_id'] == target) & (uniq_rows['ligand_id'] == ligand), par+'_avg'] = statistics.mean([float(data[0].split(' - ')[0]), float(data[0].split(' - ')[1])])
                                    uniq_rows.loc[(uniq_rows['target_id'] == target) & (uniq_rows['ligand_id'] == ligand), par+'_max'] = data[0].split(' - ')[1]
                            else:
                                try:
                                    vals = [float(x) for x in data]
                                except ValueError:
                                    vals = [float(y) for x in data for y in x.split(' - ')]
                                uniq_rows.loc[(uniq_rows['target_id'] == target) & (uniq_rows['ligand_id'] == ligand), par+'_min'] = min(vals)
                                uniq_rows.loc[(uniq_rows['target_id'] == target) & (uniq_rows['ligand_id'] == ligand), par+'_avg'] = statistics.mean(vals)
                                uniq_rows.loc[(uniq_rows['target_id'] == target) & (uniq_rows['ligand_id'] == ligand), par+'_max'] = max(vals)
        return uniq_rows

    @staticmethod
    def labeling_principals(dataframe):
        dataframe['principal/secondary'] = np.nan
        ids = list(dataframe['target_id'].unique())
        not_commented = []
        #adding a new labeling for drugs that have been
        #defined as proposed endogenous ligands
        for val_id in ids:
            data_slice = dataframe.loc[dataframe['target_id'] == val_id]
            try:
                comment = data_slice['natural/endogenous_ligand_comments'].unique()[0].split('.')[0]
            except AttributeError: #the comment is nan
                comment = ''
            if len(data_slice['ligand_name'].unique()) == 1:
                label = 'Principal'
                if 'Proposed' in comment:
                    label = 'Proposed'
                dataframe.loc[dataframe['target_id'] == id, 'principal/secondary'] = label
            if 'principal' in comment:
                if 'agonists' in comment:
                    drugs = comment.replace(' and ', ', ').split(' are')[0].split(', ')
                    drugs = [x.strip(',') for x in drugs]
                    dataframe.loc[(dataframe['target_id'] == val_id) & (dataframe.ligand_name.isin(drugs)), 'principal/secondary'] = 'Principal'
                    dataframe.loc[(dataframe['target_id'] == val_id) & (~dataframe.ligand_name.isin(drugs)), 'principal/secondary'] = 'Secondary'
                else:
                    drugs = comment.split(' is')[0]
                    dataframe.loc[(dataframe['target_id'] == val_id) & (dataframe['ligand_name'] == drugs), 'principal/secondary'] = 'Principal'
                    dataframe.loc[(dataframe['target_id'] == val_id) & (dataframe['ligand_name'] != drugs), 'principal/secondary'] = 'Secondary'
            elif ('Proposed' in comment) and (len(data_slice['ligand_name'].unique()) > 1):
                dataframe.loc[dataframe['target_id'] == id, 'principal/secondary'] = 'Proposed'
                not_commented.append(val_id)
            else:
                not_commented.append(val_id)

        return dataframe, not_commented

    @staticmethod
    def adding_potency_rankings(gtop_endogenous, not_commented):
        #fix things, drop unused values
        gtop_endogenous['ranking'] = np.nan
        missing_info = []
        gtop_endogenous.pki_avg.fillna(gtop_endogenous.pki_max, inplace=True)
        gtop_endogenous.pec50_avg.fillna(gtop_endogenous.pec50_max, inplace=True)
        gtop_endogenous.pkd_avg.fillna(gtop_endogenous.pkd_max, inplace=True)
        gtop_endogenous.pic50_avg.fillna(gtop_endogenous.pic50_max, inplace=True)

        #adding ranking to ligands in receptors without principal status information
        #while tracking problematic values (missing info, symbols in data etc)
        for target_id in not_commented:
          data_slice = gtop_endogenous.loc[gtop_endogenous['target_id'] == target_id]
          if len(data_slice['ligand_name'].unique()) != 1:
              if data_slice['pec50_avg'].isna().any() == False:
                  try:
                      #we have all pec50 values
                      sorted_list = sorted(list(set([float(x) for x in data_slice['pec50_avg'].to_list()])), reverse=True)
                      counter = 1
                      for item in sorted_list:
                          if item in data_slice['pec50_avg'].to_list():
                              gtop_endogenous.loc[(gtop_endogenous['target_id'] == target_id) & (gtop_endogenous['pec50_avg'] == item), 'ranking'] = counter
                          else:
                              gtop_endogenous.loc[(gtop_endogenous['target_id'] == target_id) & (gtop_endogenous['pec50_avg'] == str(item)), 'ranking'] = counter
                          counter += 1
                  except ValueError:
                      missing_info.append(target_id)
              elif data_slice['pki_avg'].isna().any() == False:
                  try:
                      #we have all pec50 values
                      sorted_list = sorted(list(set([float(x) for x in data_slice['pki_avg'].to_list()])), reverse=True)
                      counter = 1
                      for item in sorted_list:
                          if item in data_slice['pki_avg'].to_list():
                              gtop_endogenous.loc[(gtop_endogenous['target_id'] == target_id) & (gtop_endogenous['pki_avg'] == item), 'ranking'] = counter
                          else:
                              gtop_endogenous.loc[(gtop_endogenous['target_id'] == target_id) & (gtop_endogenous['pki_avg'] == str(item)), 'ranking'] = counter
                          counter += 1
                  except ValueError:
                      missing_info.append(target_id)
              else:
                  #we don't have full values, grab higher pec50 or higher pki?
                  values_pec50 = data_slice['pec50_avg'].dropna().to_list()
                  values_pki = data_slice['pki_avg'].dropna().to_list()
                  if len(values_pec50) > 0:
                      try:
                          #we have all pec50 values
                          sorted_list = sorted(list(set([float(x) for x in data_slice['pec50_avg'].dropna().to_list()])), reverse=True)
                          counter = 1
                          for item in sorted_list:
                              if item in data_slice['pec50_avg'].to_list():
                                  gtop_endogenous.loc[(gtop_endogenous['target_id'] == target_id) & (gtop_endogenous['pec50_avg'] == item), 'ranking'] = counter
                              else:
                                  gtop_endogenous.loc[(gtop_endogenous['target_id'] == target_id) & (gtop_endogenous['pec50_avg'] == str(item)), 'ranking'] = counter
                              counter += 1
                          gtop_endogenous.loc[(gtop_endogenous['target_id'] == target_id) & (gtop_endogenous['pec50_avg'].isna()), 'ranking'] = counter
                      except ValueError:
                          missing_info.append(target_id)
                  elif len(values_pki) > 0:
                      try:
                          #we have all pec50 values
                          sorted_list = sorted(list(set([float(x) for x in data_slice['pki_avg'].dropna().to_list()])), reverse=True)
                          counter = 1
                          for item in sorted_list:
                              if item in data_slice['pki_avg'].to_list():
                                  gtop_endogenous.loc[(gtop_endogenous['target_id'] == target_id) & (gtop_endogenous['pki_avg'] == item), 'ranking'] = counter
                              else:
                                  gtop_endogenous.loc[(gtop_endogenous['target_id'] == target_id) & (gtop_endogenous['pki_avg'] == str(item)), 'ranking'] = counter
                              counter += 1
                          gtop_endogenous.loc[(gtop_endogenous['target_id'] == target_id) & (gtop_endogenous['pki_avg'].isna()), 'ranking'] = counter
                      except ValueError:
                          missing_info.append(target_id)
                  else:
                      missing_info.append(id)
        return gtop_endogenous

    @staticmethod
    def convert_dataframe(df):
        df = df.astype(str)
        for column in df:
            df[column] = df[column].replace({'nan':None})
        return_list_of_dicts = df.to_dict('records')
        return return_list_of_dicts

    @staticmethod
    def normalize_gtp_headers(df):
        # Fixing GtP whitespace issues in the headers
        df.columns = df.columns.str.strip()
        # Fixing changing format of GtP headers
        df.columns = [c.lower().replace(' ', '_') for c in df.columns]

    @staticmethod
    def create_model(gtop_endogenous):
        values = ['pki', 'pec50', 'pkd', 'pic50']

        human_entries = [entry["target_id"]+"|"+entry["ligand_id"] for entry in gtop_endogenous if entry["interaction_species"] in [None, "None", "", "Human"]]

        stereo_ligs = {}
        for row in gtop_endogenous:
            numeric_data = {}
            receptor = Command.fetch_protein(row['target_id'], 'GtoP', row['interaction_species'])

            # TODO Handle multiple matches (uniprot filter?)
            ligand = get_ligand_by_id("gtoplig", row['ligand_id'])
            if ligand != None:
                # Process stereoisomers when not specified:
                if ligand != None and ligand.smiles != None and row['ligand_id'] not in stereo_ligs:
                    stereo_ligs[row['ligand_id']] = []
                    # Check if compound has undefined chiral centers
                    input_mol = dm.to_mol(ligand.smiles, sanitize=False)
                    chiral_centers = Chem.FindMolChiralCenters(input_mol, useLegacyImplementation=False, includeUnassigned=True)
                    undefined_centers = [center for center in chiral_centers if center[1] == "?"]
                    # If up to 2 chiral centers => generate and match them via UniChem
                    if len(chiral_centers) > 0 and len(chiral_centers) <= 2 and len(chiral_centers) == len(undefined_centers):
                        isomers = tuple(EnumerateStereoisomers(input_mol, options = StereoEnumerationOptions(unique=True)))
                        for isomer in isomers:
                            new_key = dm.to_inchikey(isomer)
                            matches = match_id_via_unichem("inchikey", new_key)
                            if len(matches) > 0:
                                ster_ids = {}
                                ster_ids["inchikey"] = new_key
                                for match in matches:
                                    if match["type"] in ["chembl_ligand", "pubchem"]:
                                        ster_ids[match["type"]] = match["id"]
                                        stereo_ligs[row['ligand_id']].append(get_or_create_ligand("", ster_ids, unichem = True, extended_matching = True))
                                        continue
                elif row['ligand_id'] not in stereo_ligs:
                    stereo_ligs[row['ligand_id']] = []
            else:
                print("Ligand ", row['ligand_id'], "not found", row['ligand_name'])
                continue #Commented, but needed for my machine (relaxin-3 was not found)

            try:
                role = Command.fetch_role(row['ligand_action'], row['ligand_role'])
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
                pmids = row['pubmed_ids'].split('|')
            except AttributeError:
                pmids = None

            #Species check
            try:
                species = row['ligand_species'].split('|')
            except AttributeError:
                species = None

            try:
                potency = int(float(row['ranking']))
            except (TypeError, ValueError):
                potency = None

            try:
                endo_status = str(row['principal/secondary'])
            except:
                endo_status = None

            if receptor is not None:
            #last step because it requires multiple uploads in case we have multiple species
                if species is None:
                    gtp_data = Endogenous_GTP(
                                ligand = ligand,
                                ligand_species = species,
                                ligand_action = role,
                                endogenous_status = endo_status, #principal/secondary
                                potency_ranking = potency, #Ranking
                                receptor = receptor, #link to protein model
                                pec50 = numeric_data['pec50'],
                                pKi = numeric_data['pki'],
                                pic50 = numeric_data['pic50'],
                                pKd = numeric_data['pkd'],
                                )
                    gtp_data.save()
                    try:
                        for pmid in pmids:
                            publication = Command.fetch_publication(pmid)
                            gtp_data.publication.add(publication)
                    except:
                        publication= None

                    for isomer in stereo_ligs[row['ligand_id']]:
                        gtp_data.pk = None
                        gtp_data.ligand = isomer

                        # Affinity/Activity varies between isomers/mixture => do not copy
                        gtp_data.pec50 = gtp_data.pKi = gtp_data.pic50 = gtp_data.pKd = "None | None | None"
                        gtp_data.save()

                    # If endogenous is not present for HUMANs => add it
                    key_id = row['target_id']+"|"+row['ligand_id']
                    if gtp_data and key_id not in human_entries:
                        human_receptor = Command.fetch_protein(row['target_id'], 'GtoP', "Human")
                        if human_receptor:
                            human_entries.append(key_id)
                            gtp_data.pk = None
                            gtp_data.receptor = human_receptor
                            # Affinity/Activity varies between species => do not copy
                            gtp_data.pec50 = gtp_data.pKi = gtp_data.pic50 = gtp_data.pKd = "None | None | None"
                            gtp_data.save()
                elif len(species) == 1:
                    ligand_species = Command.fetch_species(species[0], row['interaction_species'])
                    gtp_data = Endogenous_GTP(
                                ligand = ligand,
                                ligand_species = ligand_species,
                                ligand_action = role,
                                endogenous_status = row['principal/secondary'], #principal/secondary
                                potency_ranking = potency, #Ranking
                                receptor = receptor, #link to protein model
                                pec50 = numeric_data['pec50'],
                                pKi = numeric_data['pki'],
                                pic50 = numeric_data['pic50'],
                                pKd = numeric_data['pkd'],
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
                        species = Command.fetch_species(s, row['interaction_species'])
                        gtp_data = Endogenous_GTP(
                                    ligand = ligand,
                                    ligand_species = species,
                                    ligand_action = role,
                                    endogenous_status = row['principal/secondary'], #principal/secondary
                                    potency_ranking = potency, #Ranking
                                    receptor = receptor, #link to protein model
                                    pec50 = numeric_data['pec50'],
                                    pKi = numeric_data['pki'],
                                    pic50 = numeric_data['pic50'],
                                    pKd = numeric_data['pkd'],
                                    )
                        gtp_data.save()
                        try:
                            for pmid in pmids:
                                publication = Command.fetch_publication(pmid)
                                gtp_data.publication.add(publication)
                        except:
                            publication= None
            else:
                print("SKIPPING", row["ligand_id"], row['target_id'], row['interaction_species'], receptor, ligand, species)



    @staticmethod
    def build_chembl_publications(publication_data):

        pub_db = {}
        wr_doi = WebResource.objects.get(slug="doi")
        existing_ids = list(WebLink.objects.filter(web_resource=wr_doi).values_list("index", flat=True).distinct())

        for _, row in publication_data.iterrows():
            if row['doi'] in existing_ids:
                wl = WebLink.objects.get(index=row['doi'], web_resource=wr_doi)
            else:
                try:
                    wl = WebLink.objects.create(index=row['doi'], web_resource=WebResource.objects.get(slug='doi'))
                except IntegrityError:
                    wl = WebLink.objects.get(index=row['doi'], web_resource=wr_doi)
            try:
                pub_db[row['doi']] = Publication.objects.get(web_link=wl)
            except Publication.DoesNotExist:
                pub = Publication()
                try:
                    pub.web_link = wl
                    pub.save()
                except IntegrityError:
                    pub = Publication.objects.get(web_link=wl)
                if not pd.isnull(row['journal']):
                    journal_slug = slugify(row['journal'].replace('.',' '))
                    if not pd.isnull(row['journal_full_title']):
                        journal_fullname =  row['journal_full_title']
                        # Remove dot when present at the end
                        if journal_fullname[-1] == ".":
                            journal_fullname = journal_fullname[:-1]
                    else:
                        journal_fullname =  row['journal']
                    pub.journal, _ = PublicationJournal.objects.get_or_create(defaults={"name": journal_fullname, 'slug': journal_slug}, name__iexact=journal_fullname)

                pub.title = row['title']
                pub.authors = row['authors']
                pub.year = row['year'] if not pd.isnull(row['year']) else None
                if (not pd.isnull(row['volume'])) and (not pd.isnull(row['first_page'])) and (not pd.isnull(row['last_page'])):
                    pub.reference = str(row['volume']) + ':' + str(row['first_page']) + '-' + str(row['last_page'])
                pub.web_link = wl
                pub.save()

                pub_db[row['doi']] = pub

        return pub_db

    @staticmethod
    def build_chembl_ligands():
        print("\n===============\n#1 Reading ChEMBL ligand data", datetime.datetime.now())
        ligand_input_file = os.sep.join([settings.DATA_DIR, "ligand_data", "assay_data", "chembl_cpds.csv.gz"])
        ligand_data = pd.read_csv(ligand_input_file, keep_default_na=False)
        for column in ligand_data:
            ligand_data[column] = ligand_data[column].replace({"":None})
        print("Found", len(ligand_data), "ligands")

        # Collect ChEMBL IDs from existing ligands
        print("\n#2 Collecting ChEMBL IDs from existing ligands", datetime.datetime.now())
        wr_chembl = WebResource.objects.get(slug="chembl_ligand")
        wr_pubchem = WebResource.objects.get(slug="pubchem")

        # Removing existing ligands based on ChEMBL IDs => second filter in loop necessary because of concatenated non-parent IDs
        existing_ids = list(LigandID.objects.filter(web_resource=wr_chembl).values_list("index", flat=True).distinct())
        filtered_ligands = ligand_data.loc[~ligand_data["molecule_chembl_id"].isin(existing_ids)].copy()

        # Set ChEMBL ID as name for ligands without pref_name
        filtered_ligands.loc[filtered_ligands['pref_name'].isnull(), 'pref_name'] = filtered_ligands.loc[filtered_ligands['pref_name'].isnull(), 'molecule_chembl_id']

        # Parse all new small-molecule ChEMBL ligands
        print("\n#3 Building new small-molecule ChEMBL ligands", datetime.datetime.now())
        sm_data = filtered_ligands.loc[filtered_ligands["molecule_type"].isin(["Small molecule", "Oligosaccharide", None])].reset_index()
        lig_entries = len(sm_data)
        print("Found", lig_entries, "new small molecules")

        # Additional matching via PubChem and InchiKeys
        existing_cids = list(LigandID.objects.filter(web_resource=wr_pubchem).values_list("index", flat=True).distinct())
        existing_inchis = list(Ligand.objects.exclude(inchikey=None).values_list("inchikey", flat=True).distinct())

        smallmol = LigandType.objects.get(slug="small-molecule")
        ligands = []
        weblinks = []
        for index, (_, row) in enumerate(sm_data.iterrows()):
            insert = True
            ids = [row['molecule_chembl_id']]

            # Filtering on non-parent ChEMBL IDs
            if row['other_ids'] != None:
                extra_ids = row['other_ids'].split(";")
                existing = list(set.intersection(set(extra_ids), set(existing_ids)))
                if len(existing) > 0:
                    # Add missing link to parent ChEMBL ID
                    match = LigandID.objects.get(index=existing[0], web_resource=wr_chembl)
                    LigandID(index=row['molecule_chembl_id'], web_resource=wr_chembl, ligand_id = match.ligand_id).save()
                    print("Found existing non-parent ChEMBL", existing[0], "for parent", row['molecule_chembl_id'])
                    insert = False # use this switch in case this is the last one in the list skipping the bulk insert
                else:
                    ids = ids + extra_ids

            # Filtering on PubChem CIDs
            if row['pubchem_cid'] != None:
                cids = row['pubchem_cid'].split(";")
                existing = list(set.intersection(set(cids), set(existing_cids)))
                if len(existing) > 0:
                    # Add missing link to parent ChEMBL ID
                    match = LigandID.objects.get(index=existing[0], web_resource=wr_pubchem)
                    LigandID(index=row['molecule_chembl_id'], web_resource=wr_chembl, ligand_id = match.ligand_id).save()
                    insert = False # use this switch in case this is the last one in the list skipping the bulk insert

            # For those rare cases in which neither the ChEMBL ID nor the PubChem ID was matched, but the InchiKey was
            if insert and row['standard_inchi_key'] in existing_inchis:
                ligand = Ligand.objects.get(inchikey=row['standard_inchi_key'])
                LigandID(index=row['molecule_chembl_id'], web_resource=wr_chembl, ligand = ligand).save()
                if row['pubchem_cid'] != None:
                    cids = row['pubchem_cid'].split(";")
                    for cid in cids:
                        LigandID(index=cid, web_resource=wr_pubchem, ligand = ligand).save()
                insert = False

            if insert:
                # creating ligand
                ligands.append(Ligand(name = row['pref_name'], ambiguous_alias = False))
                ligands[-1].ligand_type = smallmol
                ligands[-1].smiles = row['smiles']
                ligands[-1].inchikey = row['standard_inchi_key']
                ligands[-1].sequence = row["sequence"]

                try:
                    input_mol = dm.to_mol(row['smiles'], sanitize=True)
                    if input_mol:
                        # Check if InChIKey has been set
                        if ligands[-1].inchikey == None:
                            ligands[-1].inchikey = dm.to_inchikey(input_mol)

                        # Cleaned InChIKey has been set
                        # ligands[-1].clean_inchikey = get_cleaned_inchikey(row['smiles'])

                        # Calculate RDkit properties
                        ligands[-1].mw = dm.descriptors.mw(input_mol)
                        ligands[-1].rotatable_bonds = dm.descriptors.n_rotatable_bonds(input_mol)
                        ligands[-1].hacc = dm.descriptors.n_hba(input_mol)
                        ligands[-1].hdon = dm.descriptors.n_hbd(input_mol)
                        ligands[-1].logp = dm.descriptors.clogp(input_mol)
                except:
                    pass

                # Adding ligand IDs
                for val_id in ids:
                    weblinks.append({"link" : LigandID(index=val_id, web_resource=wr_chembl), "lig_idx" : len(ligands)-1})
                if row['pubchem_cid'] != None:
                    cids = row['pubchem_cid'].split(";")
                    for cid in cids:
                        weblinks.append({"link" : LigandID(index=cid, web_resource=wr_pubchem), "lig_idx" : len(ligands)-1})

            # BULK insert every X entries or last entry
            if len(ligands) == Command.bulk_size or (index == lig_entries - 1):
                Ligand.objects.bulk_create(ligands)

                # Ligands have been inserted => update LigandIDs for pairing
                for pair in weblinks:
                    pair["link"].ligand = ligands[pair["lig_idx"]]
                LigandID.objects.bulk_create([pair["link"] for pair in weblinks])

                print("Inserted", index + 1, "out of", lig_entries, "ligands")
                ligands = []
                weblinks = []

        # Parse all new non-small-molecule ChEMBL ligands
        print("\n#4 Building new non-small-molecule ChEMBL ligands", datetime.datetime.now())
        nonsm_data = filtered_ligands.loc[~filtered_ligands["molecule_type"].isin(["Small molecule", "Oligosaccharide", None])]
        print("Found", len(nonsm_data), "new non-small-molecules")

        ligands = []
        ligand_types = {"Unknown": "na", "Protein": "protein"}
        weblinks = []
        for _, row in nonsm_data.iterrows():
            ids = {}
            ids["smiles"] = row['smiles']
            ids["sequence"] = row['sequence']
            ids["inchikey"] = row['standard_inchi_key']
            ids["chembl_ligand"] = row['molecule_chembl_id']

            # Filter types
            ligand = get_or_create_ligand(row['pref_name'], ids, ligand_types[row['molecule_type']], False, True)

            # Add LigandIDs
            if row['other_ids'] != None:
                extra_ids = row['other_ids'].split(";")
                existing = list(set.intersection(set(extra_ids), set(existing_ids)))
                if len(existing) > 0:
                    continue # skip rest of creation
                else:
                    for extra_id in extra_ids:
                        weblinks.append(LigandID(ligand = ligand, index=extra_id, web_resource=wr_chembl))
        # Bulk insert all new ligandIDs
        LigandID.objects.bulk_create(weblinks)

    @staticmethod
    def build_chembl_bioactivities():
        print("\n===============\n#1 Reading ChEMBL bioacitivity data")

        cache_dir = ['chembl', 'document_entry']
        bioactivity_input_file = os.sep.join([settings.DATA_DIR, "ligand_data", "assay_data", "chembl_bioactivity_data.csv.gz"])
        chembl_document_conversion_file = os.sep.join([settings.DATA_DIR, "ligand_data", "assay_data", "chembl_document_data.csv.gz"])
        chembl_document_data = pd.read_csv(chembl_document_conversion_file, dtype=str)

        publication_array = Command.build_chembl_publications(chembl_document_data)

        bioactivity_data = pd.read_csv(bioactivity_input_file, dtype=str)

        url_doc_template = 'https://www.ebi.ac.uk/chembl/api/data/document?document_chembl_id=$index'
        bioactivity_data.fillna('None', inplace=True)
        bio_entries = len(bioactivity_data)
        print("Found", bio_entries, "bioactivities", datetime.datetime.now())

        print("\n#2 Building ChEMBL ligands cache", datetime.datetime.now())
        #ids = list(bioactivity_data["parent_molecule_chembl_id"].unique())  # not filtering is way faster
        ligands = list(LigandID.objects.filter(index__startswith="CHEMBL").values_list("ligand_id", "index"))
        lig_dict = {entry[1]: entry[0] for entry in ligands}

        print("\n#3 Building ChEMBL proteins cache", datetime.datetime.now())
        # NOTE => might need to switch to Accession as the Entry name changes more frequently
        # If so, keep isoform notations in mind
        names = list(bioactivity_data["Entry name"].unique())
        proteins = list(Protein.objects.filter(entry_name__in=names).values_list("pk", "entry_name"))
        prot_dict = {prot_entry[1]: prot_entry[0] for prot_entry in proteins}

        print("\n#4 Building ChEMBL bioactivity entries", datetime.datetime.now())
        bioacts = []
        pub_links = []
        for index, (_, row) in enumerate(bioactivity_data.iterrows()):
            try:
                if (row["parent_molecule_chembl_id"] in lig_dict.keys()) and (row["Entry name"] in prot_dict.keys()):
                    bioacts.append(AssayExperiment())
                    bioacts[-1].ligand_id = lig_dict[row["parent_molecule_chembl_id"]]
                    bioacts[-1].protein_id = prot_dict[row["Entry name"]]
                    bioacts[-1].assay_type = row["assay_type"]
                    bioacts[-1].assay_description = row["assay_description"]
                    bioacts[-1].standard_activity_value = row["standard_value"]
                    bioacts[-1].p_activity_value = row["pchembl_value"]
                    bioacts[-1].p_activity_ranges = None
                    bioacts[-1].standard_relation = row["standard_relation"]
                    bioacts[-1].value_type = row["standard_type"]
                    bioacts[-1].document_chembl_id = row["document_chembl_id"]
                    bioacts[-1].source = 'ChEMBL'

                    try:
                        doi = chembl_document_data.loc[chembl_document_data['document_chembl_id'] == row["document_chembl_id"], 'doi'].iloc[0]
                    except IndexError:
                        response = fetch_from_web_api(url_doc_template, row["document_chembl_id"], cache_dir, xml=True)
                        doi = response[0][0][4].text

                    if doi is not None:
                        publication = publication_array[doi]
                        if publication is not None:
                            pub_links.append([len(bioacts) -1, publication.id])
                    # except:
                    #     pass

                    # BULK insert every X entries or last entry
                    if (len(bioacts) == Command.bulk_size) or (index == bio_entries - 1):
                        print('Inserting bulk data')
                        AssayExperiment.objects.bulk_create(bioacts)

                        for pair in pub_links:
                            pair[0] = bioacts[pair[0]]

                        data_pub = AssayExperiment.publication.through
                        data_pub.objects.bulk_create([
                                            data_pub(assayexperiment_id=experiment.pk, publication_id=pub_id)
                                            for (experiment, pub_id) in pub_links], ignore_conflicts=True)
                        print("Inserted", index, "out of", bio_entries, "bioactivities")
                        bioacts = []
                        pub_links = []
            except KeyError:
                continue

    @staticmethod
    def build_pubchem_vendor_links():
        LigandVendors.objects.all().delete()
        print("\n===============\n#1 Reading and creating Vendors")
        vendor_url = os.sep.join([settings.DATA_DIR, "ligand_data", "assay_data", "pubchem_vendor_list.csv.gz"])
        vendor_data = pd.read_csv(vendor_url, dtype=str)
        vendors = []
        for _, row in vendor_data.iterrows():
            vendors.append(LigandVendors(slug=slugify(row["SourceName"]), name = row["SourceName"], url = row["SourceURL"]))
        LigandVendors.objects.bulk_create(vendors)
        vendor_dict = {vendor.name : vendor.pk for vendor in vendors}

        print("\n#2 Building ChEMBL ligands cache", datetime.datetime.now())
        ligands = list(LigandID.objects.filter(index__startswith="CHEMBL").values_list("ligand_id", "index"))
        lig_dict = {entry[1]: entry[0] for entry in ligands}

        print("\n#3 Creating all vendor links", datetime.datetime.now())
        vendor_links_url = os.sep.join([settings.DATA_DIR, "ligand_data", "assay_data", "pubchem_vendor_links.csv.gz"])
        vendor_links_data = pd.read_csv(vendor_links_url, dtype=str)
        links = []
        for _, row in vendor_links_data.iterrows():
            if len(row["SourceRecordURL"]) < 300:
                links.append(LigandVendorLink(vendor_id=vendor_dict[row["SourceName"]], ligand_id = lig_dict[row["chembl_id"]], url = row["SourceRecordURL"], external_id = row["RegistryID"]))

        LigandVendorLink.objects.bulk_create(links)

    @staticmethod
    def uniprot_mapper(protein, organism):
        organism_dict = {'PIG': 'sus_scrofa','RAT': 'rattus_norvegicus','HUMAN': 'homo_sapiens','MOUSE': 'mus_musculus',
                         'CANINE': 'canis_lupus_familiaris','BOVINE': 'bos_taurus','CALF': 'bos_taurus','COW': 'bos_taurus',
                         'GUINEA PIG': 'cavia_porcellus', 'CAT': 'felis_catus','NEONATAL RAT': 'rattus_norvegicus',
                         '? HUMAN': 'homo_sapiens','OPOSSUM': 'didelphis_marsupialis', 'Rat 6B': 'rattus_norvegicus',
                         'HUMAN M3': 'homo_sapiens','HUMAN M4': 'homo_sapiens', 'Chick': 'gallus_gallus',
                         'Frog': 'pseudis_balbodactyla','Newborn rats': 'rattus_norvegicus', 'Beef': 'bos_taurus',
                         'Sheep': 'ovis_aries','OX': 'bos_taurus','Dog': 'canis_lupus_familiaris',
                         'Rhesus': 'macaca_mulatta','Monkey': 'macaca_mulatta','PIGLET': 'sus_scrofa',
                         'Rat Y861': 'rattus_norvegicus','Zebra Finch': 'taeniopygia_guttata','Chicken': 'gallus_gallus',
                         'MICE': 'mus_musculus','Rhesus Monkey': 'macaca_mulatta', 'Zebrafish': 'danio_rerio'}
        if organism in organism_dict.keys():
            query = 'gene_exact:{0} AND organism:{1}'.format(protein.lower(), organism_dict[organism])
        else:
            query = 'gene_exact:{}'.format(protein.lower())

        if query not in Command.mapper_cache.keys():
            url = 'https://legacy.uniprot.org/uniprot/'
            params = {
                'query': query,
                'format': 'tab',
                'columns': 'entry_name'
            }
            data = urllib.parse.urlencode(params)
            data = data.encode('utf-8')
            req = urllib.request.Request(url, data)
            try:
                converted = urllib.request.urlopen(req).read().decode('utf-8').split('\n')[1].lower()
            except IndexError:
                converted = None

            Command.mapper_cache[query] = converted

        return Command.mapper_cache[query]


    @staticmethod
    def get_ligands_data(ligands, complete_ligands, ligand_mapping, ligand_interactions=pd.DataFrame(), target_ids=[]):
        full_info = ['ligand_id', 'name','species','type','approved','withdrawn','labelled','radioactive', 'pubchem_sid', 'pubchem_cid',
                     'uniprot_id','iupac_name', 'inn', 'synonyms','smiles','inchikey','inchi','gtoimmupdb','gtompdb']
        bioactivity_info = ['ligand_id', 'action','target', 'target_id', 'target_species',
                            'target_gene_symbol', 'target_uniprot_id', 'original_affinity_relation',
                            'action_comment', 'selectivity', 'primary_target',
                            'concentration_range', 'affinity_units', 'affinity_high',
                            'affinity_median', 'affinity_low', 'assay_description', 'pubmed_id']

        weblinks = ['ligand_id', 'chembl_id','chebi_id','cas','drugbank_id','drug_central_id']

        ligand_data = complete_ligands.loc[complete_ligands['ligand_id'].isin(ligands), full_info]
        ligand_weblinks = ligand_mapping.loc[ligand_mapping['ligand_id'].isin(ligands), weblinks]
        ligand_complete = ligand_data.merge(ligand_weblinks, on="ligand_id")

        if not ligand_interactions.empty:
            bioactivity_data = ligand_interactions.loc[ligand_interactions['ligand_id'].isin(ligands), bioactivity_info]
            if len(target_ids) > 0:
                bioactivity_data = bioactivity_data.loc[bioactivity_data["target_id"].isin(target_ids)]
            ligand_complete = ligand_complete.merge(bioactivity_data, on="ligand_id")

        return ligand_complete

    @staticmethod
    def obtain_ligands(data, compare_set, labels):
        interactions_targets = list(data[labels[0]].unique())
        targets_with_ligands = set(interactions_targets).intersection(set(compare_set))
        ligands = list(data.loc[data[labels[0]].isin(targets_with_ligands), labels[1]].unique())
        ligands = [x for x in ligands if x == x] #remove nan
        return ligands

    @staticmethod
    def compare_proteins(gtp_data):
        gpcrdb_proteins = Protein.objects.filter(family__slug__startswith="00", sequence_type__slug="wt").values_list('entry_name','accession')
        entries = gtp_data.loc[gtp_data['uniprotkb_id'].isin([protein[1].split("-")[0] for protein in gpcrdb_proteins]), ['uniprotkb_id', 'gtopdb_iuphar_id']]
        return list(entries['gtopdb_iuphar_id'].unique())

    @staticmethod
    def classify_assay(biodata):
        #starting to assess assay type
        biodata['assay_type'] = 'U'
        biodata['assay_description'].fillna('Unclassified',inplace=True)
        # find Binding assays (displacement, affinity, ) KI and KD are always binding
        biodata.loc[biodata['affinity_units'].isin(['pKi', 'pKd']), 'assay_type'] = 'B'
        # EC50, pKB and pA2 are always Functional
        biodata.loc[biodata['affinity_units'].isin(['pKB', 'pEC50', 'pA2']), 'assay_type'] = 'F'
        # IC50 can be both B: displacement, affinity, binding, radioligand else is F
        binding_words = ['isplacement', 'ffinity', 'inding', 'adioligand']
        biodata.loc[(biodata['affinity_units'] == 'pIC50') & (biodata['assay_description'].str.contains('|'.join(binding_words))), 'assay_type'] = 'B'
        biodata.loc[(biodata['affinity_units'] == 'pIC50') & ~(biodata['assay_description'].str.contains('|'.join(binding_words))), 'assay_type'] = 'F'
        # find Unclassified assayas
        biodata.loc[biodata['assay_description'].str.contains('Unclassified'), 'assay_type'] = 'U'

        return biodata

    @staticmethod
    def save_the_ligands_save_the_world(lig_df, pep_df):
        types_dict = {'Inorganic': 'small-molecule',
                      'Metabolite': 'small-molecule',
                      'Natural product': 'small-molecule',
                      'Peptide': 'peptide',
                      'Synthetic organic': 'small-molecule',
                      'Antibody': 'protein',
                      None: 'na'}
        weblink_keys = {'smiles': 'smiles',
                        'inchikey': 'inchikey',
                        'pubchem': 'pubchem_cid',
                        'gtoplig': 'ligand_id',
                        'chembl_ligand': 'chembl_id',
                        'drugbank': 'drugbank_id',
                        'drug_central': 'drug_central_id'}

        #remove radioactive ligands
        lig_df = lig_df.loc[lig_df['radioactive'] != 'yes']

        #Process dataframes and replace nan with None
        lig_df = lig_df.astype(str)
        for column in lig_df:
            lig_df[column] = lig_df[column].replace({'nan':None})
        pep_df = pep_df.astype(str)
        for column in pep_df:
            pep_df[column] = pep_df[column].replace({'nan':None})

        #start parsing the GtP ligands
        issues = []
        for _, row in lig_df.iterrows():
            if row['name']:
                ids = {}
                for key, value in weblink_keys.items():
                    if row[value] != None:
                        ids[key] = str(row[value])
                        if is_float(ids[key]):
                            ids[key] = str(int(float(ids[key])))

                # When protein or peptide => get UNIPROT AND FASTA Here
                uniprot_ids = []
                if types_dict[row['type']] != "small-molecule":
                    if pep_df["ligand_id"].eq(row["ligand_id"]).any():
                        peptide_entry = pep_df.loc[pep_df["ligand_id"] == row["ligand_id"]]

                        sequence = peptide_entry["single_letter_amino_acid_sequence"].item()
                        if sequence != None and sequence != "" and len(sequence) < 1000:
                            ids["sequence"] = sequence

                        uniprot = peptide_entry["uniprot_id"].item()
                        if uniprot != None and uniprot != "":
                            # TODO - when multiple UniProt IDs => clone ligand add new UniProt IDs for other species
                            uniprot_ids = uniprot.split("|")
                            ids["uniprot"] = uniprot_ids[0].strip()

                # Create or get ligand
                ligand = get_or_create_ligand(row['name'], ids, types_dict[row['type']], True, False)
                if ligand == None:
                    print("Issue with", row['name'])
                    exit()

                # Process multiple UniProt IDs
                if len(uniprot_ids) > 1:
                    ligand.uniprot = ",".join(uniprot_ids)
                    ligand.save()

        print("DONE - Ligands with issues:")
        print(issues)

    @staticmethod
    def build_gtp_bioactivities(gtp_biodata):
        print("# Start parsing the GTP Dataframe")
        for _, row in gtp_biodata.iterrows():
            receptor = Command.fetch_protein(row['target_id'], 'GtoP', row['target_species'])
            # TODO Handle multiple matches (uniprot filter?)
            ligand = get_ligand_by_id("gtoplig", row['ligand_id'])

            try:
                low_value = "{:.2f}".format(float(row['affinity_low']))
            except ValueError:
                low_value = 'None'

            try:
                high_value = "{:.2f}".format(float(row['affinity_high']))
            except ValueError:
                high_value =  'None'

            if row['affinity_median'] != 'None':
                activity_value = "{:.2f}".format(float(row['affinity_median']))
            elif row['affinity_high'] != 'None':
                if row['affinity_low'] != 'None':
                    activity_value = "{:.2f}".format(statistics.mean([float(row['affinity_high']),float(row['affinity_low'])]))
                else:
                    activity_value = "{:.2f}".format(float(row['affinity_high']))
            else:
                activity_value = 'None'

            ranges = '|'.join([low_value, activity_value, high_value])

            #Adding publications from the PMIDs section
            try:
                pmids = row['pubmed_id'].split('|')
            except AttributeError:
                pmids = None

            if (receptor is not None) and (ligand is not None):
            #last step because it requires multiple uploads in case we have multiple species
                    gtp_data = AssayExperiment(
                                ligand = ligand,
                                protein = receptor,
                                assay_type = row['assay_type'],
                                assay_description = row['assay_description'],
                                standard_activity_value = None,
                                p_activity_value = activity_value,
                                p_activity_ranges = ranges,
                                standard_relation = row['original_affinity_relation'],
                                value_type = row['affinity_units'],
                                source = 'Guide to Pharmacology',
                                document_chembl_id = None,
                                )
                    gtp_data.save()
                    try:
                        for pmid in pmids:
                            publication = Command.fetch_publication(pmid)
                            gtp_data.publication.add(publication)
                    except:
                        publication= None
            else:
                print("SKIPPING", ligand, row["ligand_id"], "|", receptor, row['target_id'], row['target_species'])

    @staticmethod
    def fetch_publication(publication_doi):
        """
        fetch publication with Publication model
        requires: publication doi or pmid

        """
        if pd.isna(publication_doi) == True:
            return None

        if ("ISBN" in publication_doi) or (publication_doi == '0'):
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
                    # if something off with publication, skip.
                    print("Publication fetching error | module: fetch_publication. Row # is : " + str(publication_doi) + ' ' + pub_type)

            Command.publication_cache[publication_doi] = pub
        else:
            pub = Command.publication_cache[publication_doi]

        return pub

    @staticmethod
    def fetch_protein(target, database, species=None):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        if database == 'GtoP':
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
        elif database == 'PDSP':
            try:
                test = None
                if Protein.objects.filter(entry_name=target):
                    protein = Protein.objects.filter(entry_name=target)
                    test = protein.get()
                elif Protein.objects.filter(web_links__index=target, web_links__web_resource__slug='uniprot'):
                    protein1 = Protein.objects.filter(web_links__index=target, web_links__web_resource__slug='uniprot')
                    test = protein1[0]
                return test
            except:
                return None

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
        lig_function = ' '.join([str(drug_type), str(drug_action)])
        lr = None
        if lig_function in conversion.keys():
            query = conversion[lig_function]
            role_slug = slugify(query)
            lr, _ = LigandRole.objects.get_or_create(slug=role_slug, defaults={'name': query})
        return lr

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
    def build_kidatabase_bioactivities():
        protein_names = {}
        ligand_cache = {}
        print("\n===============\n#1 Reading PDSP bioacitivity data")
        pdsp_link = get_or_create_url_cache("https://pdsp.unc.edu/databases/kiDownload/download.php", 7 * 24 * 3600)
        bioactivity_kidata = pd.read_csv(pdsp_link, dtype=str, encoding='mac_roman')
        #Keeping data that has either SMILES info OR CAS info
        #CAS number can be translated into pubchem CID
        bioactivity_data_filtered = bioactivity_kidata.loc[(~bioactivity_kidata['SMILES'].isnull()) | (~bioactivity_kidata['CAS'].isnull())]
        bioactivity_data_filtered = bioactivity_data_filtered.loc[(~bioactivity_data_filtered['Unigene'].isnull())]
        bioactivity_data_filtered.fillna('None', inplace=True)
        bio_entries = len(bioactivity_data_filtered)
        print("\n===============\n#2 Start parsing PDSP data")
        bioacts = []
        for index, (_, row) in enumerate(bioactivity_data_filtered.iterrows()):
            ids = {}
            receptor = None
            label = '_'.join([row['Unigene'], row['species']])
            if label not in protein_names.keys():
                protein = Command.uniprot_mapper(row['Unigene'], row['species'])
                if protein is not None:
                    protein_names[label] = protein

            if row['SMILES'] != 'None':
                ids['smiles'] = row['SMILES']
            if row['CAS'] != 'None':
                ids['CAS'] = row['CAS']
            if row[' Ligand Name'] not in ligand_cache.keys():
                ligand = get_or_create_ligand(row[' Ligand Name'], ids)
                ligand_cache[row[' Ligand Name']] = ligand
            if label in protein_names.keys():
                receptor = Command.fetch_protein(protein_names[label], 'PDSP')
            if (receptor is not None) and (ligand_cache[row[' Ligand Name']] is not None):
                bioacts.append(AssayExperiment())
                bioacts[-1].ligand_id = ligand_cache[row[' Ligand Name']].id
                bioacts[-1].protein_id = receptor.id
                bioacts[-1].assay_type = 'B'
                bioacts[-1].assay_description = None
                bioacts[-1].standard_activity_value  = round(float(row['ki Val']),2)
                bioacts[-1].p_activity_value = round(-math.log10(float(row['ki Val'])*10e-9),2)
                bioacts[-1].p_activity_ranges = None
                bioacts[-1].standard_relation = '='
                bioacts[-1].value_type = 'pKi'
                bioacts[-1].source = 'PDSP KiDatabase'
                bioacts[-1].document_chembl_id = None
                # BULK insert every X entries or last entry
            if (len(bioacts) == Command.bulk_size) or (index == bio_entries - 1):
                # AssayExperiment.objects.bulk_create(bioacts)
                print("Inserted", index, "out of", bio_entries, "bioactivities")
                bioacts = []

