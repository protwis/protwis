from build.management.commands.base_build import Command as BaseBuild
from build.management.commands.build_ligand_functions import get_ligand_by_id, match_id_via_unichem, get_or_create_ligand, is_float, standardize_smiles, generate_parent, apply_canonical_ligand_types
from django.conf import settings
from django.utils.text import slugify
from django.db import IntegrityError, transaction, connection
from django.db.models import Count, Q

from common.tools import get_or_create_url_cache, fetch_from_web_api, test_model_updates, find_role
from common.models import WebLink, WebResource, Publication, PublicationJournal
from ligand.models import Ligand, LigandID, LigandType, LigandVendors, LigandVendorLink, AssayExperiment, Endogenous_GTP, LigandRole, LigandEffect, LigandTargetPairing
from protein.models import Protein, Species

import requests
import math
import os
import statistics
import datamol as dm
import datetime
import pandas as pd
import numpy as np
import urllib.parse
import urllib.request
import django.apps

from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

class Command(BaseBuild):
    help = "Build ChEMBL, PDSP and GtoP bioactivities data"
    bulk_size = 10000
    mapper_cache = {}
    publication_cache = {}
    #Setting the variables for the test tracking of the model upadates
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    data_dir = os.sep.join([settings.DATA_DIR, 'ligand_data'])
    helm_chembl_filepath = os.sep.join([data_dir, 'HELM_CHEMBL.csv'])
    helm_cid_filepath = os.sep.join([data_dir, 'HELM_CID.csv'])
    helm_chembl = pd.read_csv(helm_chembl_filepath, index_col=0)
    helm_cid = pd.read_csv(helm_cid_filepath, index_col=0)
    test_model_updates(all_models, tracker, initialize=True)

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
            self.tracker = {}
            test_model_updates(self.all_models, self.tracker, initialize=True)
            print("Ended purging data")

        # Fetching all the Guide to Pharmacology data
        print("\n\nStarted parsing Guide to Pharmacology bioactivities data")
        gtp_uniprot_link = get_or_create_url_cache(
            "https://www.guidetopharmacology.org/DATA/GtP_to_UniProt_mapping.csv", 7 * 24 * 3600)
        gtp_uniprot = pd.read_csv(gtp_uniprot_link, dtype=str, header=1)
        self.normalize_gtp_headers(gtp_uniprot)
        gtp_complete_ligands_link = get_or_create_url_cache(
            "https://www.guidetopharmacology.org/DATA/ligands.csv", 7 * 24 * 3600)
        gtp_complete_ligands = pd.read_csv(
            gtp_complete_ligands_link, dtype=str, header=1)
        self.normalize_gtp_headers(gtp_complete_ligands)
        gtp_ligand_mapping_link = get_or_create_url_cache(
            "https://www.guidetopharmacology.org/DATA/ligand_id_mapping.csv", 7 * 24 * 3600)
        gtp_ligand_mapping = pd.read_csv(
            gtp_ligand_mapping_link, dtype=str, header=1)
        self.normalize_gtp_headers(gtp_ligand_mapping)
        gtp_interactions_link = get_or_create_url_cache(
            "https://www.guidetopharmacology.org/DATA/interactions.csv", 7 * 24 * 3600)
        gtp_interactions = pd.read_csv(
            gtp_interactions_link, dtype=str, header=1)
        self.normalize_gtp_headers(gtp_interactions)
        gtp_detailed_endogenous_link = get_or_create_url_cache(
            "https://www.guidetopharmacology.org/DATA/endogenous_ligand_detailed.csv", 7 * 24 * 3600)
        gtp_detailed_endogenous = pd.read_csv(
            gtp_detailed_endogenous_link, dtype=str, header=1)
        self.normalize_gtp_headers(gtp_detailed_endogenous)
        gtp_peptides_link = get_or_create_url_cache(
            "https://www.guidetopharmacology.org/DATA/peptides.csv", 7 * 24 * 3600)
        gtp_peptides = pd.read_csv(gtp_peptides_link, dtype=str, header=1)
        self.normalize_gtp_headers(gtp_peptides)

        # This gets all the info of the ligand and the interaction with the target
        iuphar_ids = self.compare_proteins(gtp_uniprot)
        bioactivity_ligands_ids = self.obtain_ligands(gtp_interactions, iuphar_ids, ['target_id', 'ligand_id'])
        # Now I have all the data I need
        bioactivity_data_gtp = self.get_ligands_data(
            bioactivity_ligands_ids, gtp_complete_ligands, gtp_ligand_mapping, ligand_interactions=gtp_interactions, target_ids=iuphar_ids)
        # Assess the assay type given info from affinity units and assay comments
        bioactivity_data_gtp = self.classify_assay(bioactivity_data_gtp, 'affinity_units', 'assay_description')
        bioactivity_data_gtp.fillna('None', inplace=True)
        print("Ended parsing Guide to Pharmacology bioactivities data")

        print("\n\nStarted building all Guide to Pharmacology ligands")
        print('\n\nRetrieving IUPHAR ids from UniProt ids')

        print('\n\nRetrieving ALL ligands from GTP associated to GPCRs')
        endogenous_ligands_ids = self.obtain_ligands(gtp_detailed_endogenous, iuphar_ids, ['target_id', 'ligand_id'])
        ligand_ids = list(set(bioactivity_ligands_ids + endogenous_ligands_ids))

        print('\n\nCollating all info from GPCR related ligands in the GTP')
        ligand_data = self.get_ligands_data(ligand_ids, gtp_complete_ligands, gtp_ligand_mapping)

        #Solving duplicated issues
        ligand_data["species"] = ligand_data["species"].fillna("")
        ligand_data["uniprot_id"] = ligand_data["uniprot_id"].fillna("")

        ligand_data["species_list"] = ligand_data["species"].apply(lambda x: [s.strip() for s in x.split(",") if s.strip()])
        ligand_data["uniprot_list"] = ligand_data["uniprot_id"].apply(lambda x: [u.strip() for u in x.split("|") if u.strip()])

        #Copy the same for gtp_peptides
        gtp_peptides["species"] = gtp_peptides["species"].fillna("")
        gtp_peptides["uniprot_id"] = gtp_peptides["uniprot_id"].fillna("")

        gtp_peptides["species_list"] = gtp_peptides["species"].apply(lambda x: [s.strip() for s in x.split(",") if s.strip()])
        gtp_peptides["uniprot_list"] = gtp_peptides["uniprot_id"].apply(lambda x: [u.strip() for u in x.split("|") if u.strip()])

        # 2) Apply the helper to each row
        ligand_expanded = []
        peptides_expanded = []
        for _, row in ligand_data.iterrows():
            ligand_expanded.extend(self.expand_row(row))
        for _, row in gtp_peptides.iterrows():
            peptides_expanded.extend(self.expand_row(row))

        ligand_data = pd.DataFrame(ligand_expanded)
        gtp_peptides = pd.DataFrame(peptides_expanded)

        # Build mask: not NaN, and not empty/whitespace-only
        mask_chembl = (
            Command.helm_chembl['helm_notation'].notna()
            & Command.helm_chembl['helm_notation'].astype(str).str.strip().ne('')
        )
        mask_cid = (
            Command.helm_cid['helm_notation'].notna()
        )

        # Apply mask
        helm_chembl_clean = Command.helm_chembl.loc[mask_chembl]
        helm_cid_clean = Command.helm_cid.loc[mask_cid]
        # And clean CID data
        helm_cid_clean['pubchem_cid'] = helm_cid_clean['pubchem_cid'].apply(
            lambda x: str(int(x)) if pd.notna(x) else np.nan
        )

        # First Merge
        ligand_data_merged = pd.merge(
            ligand_data,
            helm_chembl_clean,
            left_on='chembl_id',
            right_on='molecule_chembl_id',
            how='left')

        # Second Merge
        ligand_data_merged = pd.merge(
            ligand_data_merged,
            helm_cid_clean,
            on='pubchem_cid',
            how='left')

        ligand_data_merged['helm_notation'] = ligand_data_merged['helm_notation_x'].fillna(ligand_data_merged['helm_notation_y'])
        ligand_data_merged = ligand_data_merged.drop(columns=['helm_notation_x', 'helm_notation_y'])

        # First Merge
        gtp_peptides_merged = pd.merge(
            gtp_peptides,
            helm_chembl_clean,
            left_on='chembl_id',
            right_on='molecule_chembl_id',
            how='left')

        # Second Merge
        gtp_peptides_merged = pd.merge(
            gtp_peptides_merged,
            helm_cid_clean,
            on='pubchem_cid',
            how='left')

        gtp_peptides_merged['helm_notation'] = gtp_peptides_merged['helm_notation_x'].fillna(gtp_peptides_merged['helm_notation_y'])
        gtp_peptides_merged = gtp_peptides_merged.drop(columns=['helm_notation_x', 'helm_notation_y'])


        print('\n\nSaving the ligands in the models')
        self.save_the_ligands_save_the_world(ligand_data_merged, gtp_peptides_merged)
        # Assign the species to duplicated peptides
        # print('\n\nAssign the species to duplicated peptides')
        # self.assign_species_to_peptide()

        print('Performing checks')
        test_model_updates(self.all_models, self.tracker, check=True)
        print('\n\nFetching Drug Bank ligands and saving to model')
        self.build_drugbank_ligands()
        print("\n\nStarted building Guide to Pharmacology bioactivities")
        self.build_gtp_bioactivities(bioactivity_data_gtp)
        print("Ended building Guide to Pharmacology bioactivities")

        print('Performing checks')
        test_model_updates(self.all_models, self.tracker, check=True)

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
        print("\n\nEnded building endogenous data")

        print('Performing checks')
        test_model_updates(self.all_models, self.tracker, check=True)

        print("\n\nStarted building ChEBML ligands")
        self.build_chembl_ligands()
        print("\n\nEnded building ChEMBL ligands")
        print('Performing checks')
        test_model_updates(self.all_models, self.tracker, check=True)

        # Parse ChEMBL bioactivity data
        print("\n\nStarted building ChEMBL bioactivities")
        self.build_chembl_bioactivities()
        print("Ended building ChEMBL bioactivities")
        print('Performing checks')
        test_model_updates(self.all_models, self.tracker, check=True)

        # Parse ChEMBL/PubChem vendor data
        print("\n\nStarted building PubChem vendor data")
        self.build_pubchem_vendor_links()
        print("Ended building PubChem vendor data")
        print('Performing checks')
        test_model_updates(self.all_models, self.tracker, check=True)

        # Building PDSP KiDatabase bioactivity data
        print("\n\nStarted building PDSP KiDatabase bioactivities")
        self.build_kidatabase_bioactivities()  # 14,562
        print("Ended building PDSP KiDatabase bioactivities")
        print('Performing checks')
        test_model_updates(self.all_models, self.tracker, check=True)

        # Building Drug Central bioactivity data
        print("\n\nStarted building Drug Central bioactivities")
        self.build_drugcentral_bioactivities()  # 5,844
        print("Ended building Drug Central bioactivities")
        print('Performing checks')
        test_model_updates(self.all_models, self.tracker, check=True)

        # Building Drug Central bioactivity data
        print("\n\nStarted calculating potency and affinity indexes")
        self.calculate_potency_and_affinity()
        print("Potency and affinity indexes have been added to the model")

        print("\n\nFixing mismatched LigandType definition")
        n  = apply_canonical_ligand_types()
        print("\n\nUpdated LigandType on {} records to their canonical type".format(n))

    @staticmethod
    def reset_pk_sequence(model):
        """
        Resets the primary key sequence for a given model.
        """
        table_name = model._meta.db_table
        sequence_name = f"{table_name}_id_seq"

        with connection.cursor() as cursor:
            cursor.execute(f'ALTER SEQUENCE "{sequence_name}" RESTART WITH 1;')


    @staticmethod
    def purge_data():
        delete_experimental = AssayExperiment.objects.all()  # New Model Biased Data
        delete_experimental.delete()
        endo_data = Endogenous_GTP.objects.all()
        endo_data.delete()
        Ligand.objects.all().delete()
        LigandID.objects.all().delete()

        # Reset primary key ID sequences
        Command.reset_pk_sequence(AssayExperiment)
        Command.reset_pk_sequence(Endogenous_GTP)
        Command.reset_pk_sequence(Ligand)
        Command.reset_pk_sequence(LigandID)

    @staticmethod
    def data_preparation(endogenous_data, interactions, iuphar_ids):
        # Remove parameter, value and PMID is columns to have the complete dataset
        filtered_data = endogenous_data.loc[endogenous_data['target_id'].isin(
            iuphar_ids)]
        uniq_rows = filtered_data.drop(columns=['interaction_parameter','interaction_value','interaction_pubmed_ids']).drop_duplicates()
        uniq_rows = uniq_rows.dropna(subset=['ligand_name'])
        association = filtered_data[['ligand_id', 'target_id']].drop_duplicates(
        ).groupby('target_id')['ligand_id'].apply(list).to_dict()

        info_we_want = ['pKi', 'pIC50', 'pKd', 'pEC50']
        new_columns = ['pki_min', 'pki_avg', 'pki_max',
                       'pkd_min', 'pkd_avg', 'pkd_max',
                       'pic50_min', 'pic50_avg', 'pic50_max',
                       'pec50_min', 'pec50_avg', 'pec50_max',
                       'ligand_species', 'ligand_action', 'ligand_role', 'interaction_pubmed_ids']

        columns = ['ligand_id', 'ligand_name', 'ligand_type', 'ligand_uniprot_ids',
                   'ligand_ensembl_gene_id', 'ligand_subunit_id', 'ligand_subunit_name',
                   'ligand_subunit_uniprot_ids', 'ligand_subunit_ensembl_ids', 'target_id',
                   'target_name', 'target_uniprot_id', 'target_ensembl_gene_id',
                   'subunit_id', 'subunit_name', 'subunit_uniprot_ids',
                   'subunit_ensembl_ids', 'natural/endogenous_ligand_comments',
                   'rankpotency', 'interaction_species']

        uniq_rows = uniq_rows.reindex(columns=columns + new_columns)
        # remove spaces in the column names
        uniq_rows.columns = [c.replace(' ', '_') for c in uniq_rows.columns]
        for target in association.keys():
            for ligand in association[target]:
                # adding species and role info
                try:
                    role = interactions.loc[(interactions['target_id'] == target) & (
                        interactions['ligand_id'] == ligand), 'action'].values[0]
                except IndexError:
                    role = None
                uniq_rows.loc[(uniq_rows['target_id'] == target) & (
                    uniq_rows['ligand_id'] == ligand), 'ligand_role'] = role
                try:
                    action = interactions.loc[(interactions['target_id'] == target) & (
                        interactions['ligand_id'] == ligand), 'type'].values[0]
                except IndexError:
                    action = None
                uniq_rows.loc[(uniq_rows['target_id'] == target) & (
                    uniq_rows['ligand_id'] == ligand), 'ligand_action'] = action
                try:
                    species = interactions.loc[(interactions['target_id'] == target) & (
                        interactions['ligand_id'] == ligand), 'ligand_species'].values[0]
                except IndexError:
                    species = None
                uniq_rows.loc[(uniq_rows['target_id'] == target) & (
                    uniq_rows['ligand_id'] == ligand), 'ligand_species'] = species
                # fetching the parameters of interaction between receptor and ligand
                params = endogenous_data.loc[(endogenous_data['target_id'] == target) & (
                    endogenous_data['ligand_id'] == ligand), 'interaction_parameter'].to_list()
                pmids = '|'.join(endogenous_data.loc[(endogenous_data['target_id'] == target) & (
                    endogenous_data['ligand_id'] == ligand), 'interaction_pubmed_ids'].dropna().to_list())
                uniq_rows.loc[(uniq_rows['target_id'] == target) & (
                    uniq_rows['ligand_id'] == ligand), 'interaction_pubmed_ids'] = pmids
                # now parsing the data based on parameter
                for par in params:
                    # we want only pKi, pKd, pEC50 and pIC50, not nans or other weird stuff
                    if par in info_we_want:
                        par_normalized = par.lower()
                        species = endogenous_data.loc[(endogenous_data['target_id'] == target) & (
                            endogenous_data['ligand_id'] == ligand) & (endogenous_data['interaction_parameter'] == par)]['interaction_species'].tolist()
                        for org in species:
                            data = endogenous_data.loc[(endogenous_data['target_id'] == target) & (endogenous_data['ligand_id'] == ligand) & (
                                endogenous_data['interaction_parameter'] == par) & (endogenous_data['interaction_species'] == org)]['interaction_value'].tolist()
                            if len(data) == 1:
                                if '-' not in data[0]:
                                    uniq_rows.loc[(uniq_rows['target_id'] == target) & (
                                        uniq_rows['ligand_id'] == ligand), par_normalized + '_avg'] = data[0]
                                    uniq_rows.loc[(uniq_rows['target_id'] == target) & (
                                        uniq_rows['ligand_id'] == ligand), par_normalized + '_max'] = data[0]
                                elif data[0] == '-':
                                    continue
                                else:
                                    uniq_rows.loc[(uniq_rows['target_id'] == target) & (
                                        uniq_rows['ligand_id'] == ligand), par_normalized + '_min'] = data[0].split(' - ')[0]
                                    uniq_rows.loc[(uniq_rows['target_id'] == target) & (uniq_rows['ligand_id'] == ligand), par_normalized +
                                                  '_avg'] = statistics.mean([float(data[0].split(' - ')[0]), float(data[0].split(' - ')[1])])
                                    uniq_rows.loc[(uniq_rows['target_id'] == target) & (
                                        uniq_rows['ligand_id'] == ligand), par_normalized + '_max'] = data[0].split(' - ')[1]
                            else:
                                try:
                                    vals = [float(x) for x in data]
                                except ValueError:
                                    vals = [float(y)
                                            for x in data for y in x.split(' - ')]
                                uniq_rows.loc[(uniq_rows['target_id'] == target) & (
                                    uniq_rows['ligand_id'] == ligand), par_normalized + '_min'] = min(vals)
                                uniq_rows.loc[(uniq_rows['target_id'] == target) & (
                                    uniq_rows['ligand_id'] == ligand), par_normalized + '_avg'] = statistics.mean(vals)
                                uniq_rows.loc[(uniq_rows['target_id'] == target) & (
                                    uniq_rows['ligand_id'] == ligand), par_normalized + '_max'] = max(vals)
        return uniq_rows

    @staticmethod
    def labeling_principals(dataframe):
        dataframe['principal/secondary'] = np.nan
        ids = list(dataframe['target_id'].unique())
        not_commented = []
        # adding a new labeling for drugs that have been
        # defined as proposed endogenous ligands
        for val_id in ids:
            data_slice = dataframe.loc[dataframe['target_id'] == val_id]
            try:
                comment = data_slice['natural/endogenous_ligand_comments'].unique()[
                    0].split('.')[0]
            except AttributeError:  # the comment is nan
                comment = ''
            if len(data_slice['ligand_name'].unique()) == 1:
                label = 'Principal'
                if 'Proposed' in comment:
                    label = 'Proposed'
                dataframe.loc[dataframe['target_id'] == val_id, 'principal/secondary'] = label
            elif 'principal' in comment:
                if 'agonists' in comment:
                    drugs = comment.replace(' and ', ', ').split(' are')[
                        0].split(', ')
                    drugs = [x.strip(',') for x in drugs]
                    dataframe.loc[(dataframe['target_id'] == val_id) & (
                        dataframe.ligand_name.isin(drugs)), 'principal/secondary'] = 'Principal'
                    dataframe.loc[(dataframe['target_id'] == val_id) & (
                        ~dataframe.ligand_name.isin(drugs)), 'principal/secondary'] = 'Secondary'
                else:
                    drugs = comment.split(' is')[0]
                    dataframe.loc[(dataframe['target_id'] == val_id) & (
                        dataframe['ligand_name'] == drugs), 'principal/secondary'] = 'Principal'
                    dataframe.loc[(dataframe['target_id'] == val_id) & (
                        dataframe['ligand_name'] != drugs), 'principal/secondary'] = 'Secondary'
            elif ('Proposed' in comment) and (len(data_slice['ligand_name'].unique()) > 1):
                dataframe.loc[dataframe['target_id'] ==
                              id, 'principal/secondary'] = 'Proposed'
                not_commented.append(val_id)
            else:
                not_commented.append(val_id)

        return dataframe, not_commented

    @staticmethod
    def adding_potency_rankings(endogenous_data, to_be_ranked):
        # fix things, drop unused values
        endogenous_data['ranking'] = np.nan
        missing_info = []
        endogenous_data.pki_avg.fillna(endogenous_data.pki_max, inplace=True)
        endogenous_data.pec50_avg.fillna(endogenous_data.pec50_max, inplace=True)
        endogenous_data.pkd_avg.fillna(endogenous_data.pkd_max, inplace=True)
        endogenous_data.pic50_avg.fillna(endogenous_data.pic50_max, inplace=True)
        endogenous_data.loc[endogenous_data['principal/secondary'] == 'Principal', 'ranking'] = 1
        endogenous_data.loc[endogenous_data['principal/secondary'] == 'Proposed', 'ranking'] = 1
        endogenous_data.loc[endogenous_data['principal/secondary'] == 'Secondary', 'ranking'] = 2
        # adding ranking to ligands in receptors without principal status information
        # while tracking problematic values (missing info, symbols in data etc)
        for target_id in to_be_ranked:
            data_slice = endogenous_data.loc[endogenous_data['target_id'] == target_id]
            if len(data_slice['ligand_name'].unique()) != 1:
                if data_slice['pec50_avg'].isna().any() is False:
                    try:
                        # we have all pec50 values
                        sorted_list = sorted(
                            list(set([float(x) for x in data_slice['pec50_avg'].to_list()])), reverse=True)
                        counter = 1
                        for item in sorted_list:
                            endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                endogenous_data['pec50_avg'] == item), 'ranking'] = counter
                            endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                endogenous_data['pec50_avg'] == str(item)), 'ranking'] = counter
                            counter += 1
                    except ValueError:
                        missing_info.append(target_id)
                elif data_slice['pki_avg'].isna().any() is False:
                    try:
                        # we have all pec50 values
                        sorted_list = sorted(
                            list(set([float(x) for x in data_slice['pki_avg'].to_list()])), reverse=True)
                        counter = 1
                        for item in sorted_list:
                            endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                endogenous_data['pki_avg'] == item), 'ranking'] = counter
                            endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                endogenous_data['pki_avg'] == str(item)), 'ranking'] = counter
                            counter += 1
                    except ValueError:
                        missing_info.append(target_id)
                elif data_slice['pkd_avg'].isna().any() is False:
                    try:
                        # we have all pec50 values
                        sorted_list = sorted(
                            list(set([float(x) for x in data_slice['pkd_avg'].to_list()])), reverse=True)
                        counter = 1
                        for item in sorted_list:
                            endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                endogenous_data['pkd_avg'] == item), 'ranking'] = counter
                            endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                endogenous_data['pkd_avg'] == str(item)), 'ranking'] = counter
                            counter += 1
                    except ValueError:
                        missing_info.append(target_id)
                if data_slice['pic50_avg'].isna().any() is False:
                    try:
                        # we have all pec50 values
                        sorted_list = sorted(
                            list(set([float(x) for x in data_slice['pic50_avg'].to_list()])), reverse=True)
                        counter = 1
                        for item in sorted_list:
                            endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                endogenous_data['pic50_avg'] == item), 'ranking'] = counter
                            endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                endogenous_data['pic50_avg'] == str(item)), 'ranking'] = counter
                            counter += 1
                    except ValueError:
                        missing_info.append(target_id)
                else:
                    # we don't have full values, grab higher pec50 or higher pki?
                    values_pec50 = data_slice['pec50_avg'].dropna().to_list()
                    values_pki = data_slice['pki_avg'].dropna().to_list()
                    values_pic50 = data_slice['pic50_avg'].dropna().to_list()
                    values_pkd = data_slice['pkd_avg'].dropna().to_list()
                    if len(values_pec50) > 0:
                        try:
                            # we have all pec50 values
                            sorted_list = sorted(list(set(
                                [float(x) for x in data_slice['pec50_avg'].dropna().to_list()])), reverse=True)
                            counter = 1
                            for item in sorted_list:
                                endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                    endogenous_data['pec50_avg'] == item), 'ranking'] = counter
                                endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                    endogenous_data['pec50_avg'] == str(item)), 'ranking'] = counter
                                counter += 1
                            endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                endogenous_data['pec50_avg'].isna()), 'ranking'] = counter
                        except ValueError:
                            missing_info.append(target_id)
                    elif len(values_pki) > 0:
                        try:
                            # we have all pec50 values
                            sorted_list = sorted(list(
                                set([float(x) for x in data_slice['pki_avg'].dropna().to_list()])), reverse=True)
                            counter = 1
                            for item in sorted_list:
                                endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                    endogenous_data['pki_avg'] == item), 'ranking'] = counter
                                endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                    endogenous_data['pki_avg'] == str(item)), 'ranking'] = counter
                                counter += 1
                            endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                endogenous_data['pki_avg'].isna()), 'ranking'] = counter
                        except ValueError:
                            missing_info.append(target_id)
                    elif len(values_pkd) > 0:
                        try:
                            # we have all pec50 values
                            sorted_list = sorted(list(
                                set([float(x) for x in data_slice['pkd_avg'].dropna().to_list()])), reverse=True)
                            counter = 1
                            for item in sorted_list:
                                endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                    endogenous_data['pkd_avg'] == item), 'ranking'] = counter
                                endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                    endogenous_data['pkd_avg'] == str(item)), 'ranking'] = counter
                                counter += 1
                            endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                endogenous_data['pkd_avg'].isna()), 'ranking'] = counter
                        except ValueError:
                            missing_info.append(target_id)
                    elif len(values_pic50) > 0:
                        try:
                            # we have all pec50 values
                            sorted_list = sorted(list(
                                set([float(x) for x in data_slice['pic50_avg'].dropna().to_list()])), reverse=True)
                            counter = 1
                            for item in sorted_list:
                                endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                    endogenous_data['pic50_avg'] == item), 'ranking'] = counter
                                endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                    endogenous_data['pic50_avg'] == str(item)), 'ranking'] = counter
                                counter += 1
                            endogenous_data.loc[(endogenous_data['target_id'] == target_id) & (
                                endogenous_data['pic50_avg'].isna()), 'ranking'] = counter
                        except ValueError:
                            missing_info.append(target_id)
                    else:
                        missing_info.append(target_id)
        return endogenous_data

    @staticmethod
    def convert_dataframe(df):
        df = df.astype(str)
        for column in df:
            df[column] = df[column].replace({'nan': None})
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

        human_entries = [entry["target_id"] + "|" + entry["ligand_id"]
                         for entry in gtop_endogenous if entry["interaction_species"] in [None, "None", "", "Human"]]

        stereo_ligs = {}
        for row in gtop_endogenous:
            numeric_data = {}
            receptor = Command.fetch_protein(
                row['target_id'], 'GtoP', row['interaction_species'])

            # TODO Handle multiple matches (uniprot filter?)
            ligand = get_ligand_by_id("gtoplig", row['ligand_id'])
            if ligand is not None:
                # Process stereoisomers when not specified:
                if ligand is not None and ligand.smiles is not None and row['ligand_id'] not in stereo_ligs:
                    stereo_ligs[row['ligand_id']] = []
                    # Check if compound has undefined chiral centers
                    input_mol = dm.to_mol(ligand.smiles, sanitize=False)
                    chiral_centers = Chem.FindMolChiralCenters(
                        input_mol, useLegacyImplementation=False, includeUnassigned=True)
                    undefined_centers = [
                        center for center in chiral_centers if center[1] == "?"]
                    # If up to 2 chiral centers => generate and match them via UniChem
                    if len(chiral_centers) > 0 and len(chiral_centers) <= 2 and len(chiral_centers) == len(undefined_centers):
                        isomers = tuple(EnumerateStereoisomers(
                            input_mol, options=StereoEnumerationOptions(unique=True)))
                        for isomer in isomers:
                            new_key = dm.to_inchikey(isomer)
                            matches = match_id_via_unichem("inchikey", new_key)
                            if len(matches) > 0:
                                ster_ids = {}
                                ster_ids["inchikey"] = new_key
                                for match in matches:
                                    if match["type"] in ["chembl_ligand", "pubchem"]:
                                        ster_ids[match["type"]] = match["id"]
                                        stereo_ligs[row['ligand_id']].append(get_or_create_ligand("", ster_ids, unichem=True, extended_matching=True))
                                        continue
                elif row['ligand_id'] not in stereo_ligs:
                    stereo_ligs[row['ligand_id']] = []
            else:
                print("Ligand ", row['ligand_id'], "not found", row['ligand_name'])
                # Commented, but needed for my machine (relaxin-3 was not found)
                continue

            try:
                role = Command.fetch_role(row['ligand_action'], row['ligand_role'])
            except AttributeError:
                role = None
            for v in values:
                # adapting average values to 2 decimal numbers when available
                try:
                    avg = "{:.2f}".format(float(row[v + '_avg']))
                except (ValueError, TypeError):
                    avg = row[v + '_avg']
                try:
                    numeric_data[v] = (' | ').join(
                        [str(row[v + '_min']), str(avg), str(row[v + '_max'])])
                except KeyError:  # missing info on pEC50
                    numeric_data[v] = ''

            # Adding publications from the PMIDs section
            try:
                pmids = row['interaction_pubmed_ids'].split('|')
            except AttributeError:
                pmids = None

            # Species check
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
                # last step because it requires multiple uploads in case we have multiple species
                if species is None:
                    gtp_data = Endogenous_GTP(
                        ligand=ligand,
                        ligand_species=species,
                        ligand_action=role,
                        endogenous_status=endo_status,  # principal/secondary
                        potency_ranking=potency,  # Ranking
                        receptor=receptor,  # link to protein model
                        pec50=numeric_data['pec50'],
                        pKi=numeric_data['pki'],
                        pic50=numeric_data['pic50'],
                        pKd=numeric_data['pkd'],
                    )
                    gtp_data.save()
                    try:
                        for pmid in pmids:
                            publication = Command.fetch_publication(pmid)
                            gtp_data.publication.add(publication)
                    except:
                        publication = None

                    for isomer in stereo_ligs[row['ligand_id']]:
                        gtp_data.pk = None
                        gtp_data.ligand = isomer

                        # Affinity/Activity varies between isomers/mixture => do not copy
                        gtp_data.pec50 = gtp_data.pKi = gtp_data.pic50 = gtp_data.pKd = "None | None | None"
                        gtp_data.save()

                    # If endogenous is not present for HUMANs => add it
                    key_id = row['target_id'] + "|" + row['ligand_id']
                    if gtp_data and key_id not in human_entries:
                        human_receptor = Command.fetch_protein(
                            row['target_id'], 'GtoP', "Human")
                        if human_receptor:
                            human_entries.append(key_id)
                            gtp_data.pk = None
                            gtp_data.receptor = human_receptor
                            # Affinity/Activity varies between species => do not copy
                            gtp_data.pec50 = gtp_data.pKi = gtp_data.pic50 = gtp_data.pKd = "None | None | None"
                            gtp_data.save()
                elif len(species) == 1:
                    ligand_species = Command.fetch_species(
                        species[0], row['interaction_species'])
                    gtp_data = Endogenous_GTP(
                        ligand=ligand,
                        ligand_species=ligand_species,
                        ligand_action=role,
                        # principal/secondary
                        endogenous_status=row['principal/secondary'],
                        potency_ranking=potency,  # Ranking
                        receptor=receptor,  # link to protein model
                        pec50=numeric_data['pec50'],
                        pKi=numeric_data['pki'],
                        pic50=numeric_data['pic50'],
                        pKd=numeric_data['pkd'],
                    )
                    gtp_data.save()

                    try:
                        for pmid in pmids:
                            publication = Command.fetch_publication(pmid)
                            gtp_data.publication.add(publication)
                    except:
                        publication = None
                else:
                    for s in species:
                        species = Command.fetch_species(
                            s, row['interaction_species'])
                        gtp_data = Endogenous_GTP(
                            ligand=ligand,
                            ligand_species=species,
                            ligand_action=role,
                            # principal/secondary
                            endogenous_status=row['principal/secondary'],
                            potency_ranking=potency,  # Ranking
                            receptor=receptor,  # link to protein model
                            pec50=numeric_data['pec50'],
                            pKi=numeric_data['pki'],
                            pic50=numeric_data['pic50'],
                            pKd=numeric_data['pkd'],
                        )
                        gtp_data.save()
                        try:
                            for pmid in pmids:
                                publication = Command.fetch_publication(pmid)
                                gtp_data.publication.add(publication)
                        except:
                            publication = None
            else:
                print("SKIPPING", row["ligand_id"], row['target_id'],
                      row['interaction_species'], receptor, ligand, species)

    @staticmethod
    def build_chembl_publications(publication_data):

        pub_db = {}
        wr_doi = WebResource.objects.get(slug="doi")
        existing_ids = list(WebLink.objects.filter(
            web_resource=wr_doi).values_list("index", flat=True).distinct())

        for _, row in publication_data.iterrows():
            if row['doi'] in existing_ids:
                wl = WebLink.objects.get(index=row['doi'], web_resource=wr_doi)
            else:
                try:
                    wl = WebLink.objects.create(
                        index=row['doi'], web_resource=WebResource.objects.get(slug='doi'))
                except IntegrityError:
                    wl = WebLink.objects.get(
                        index=row['doi'], web_resource=wr_doi)
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
                    journal_slug = slugify(row['journal'].replace('.', ' '))
                    if not pd.isnull(row['journal_full_title']):
                        journal_fullname = row['journal_full_title']
                        # Remove dot when present at the end
                        if journal_fullname[-1] == ".":
                            journal_fullname = journal_fullname[:-1]
                    else:
                        journal_fullname = row['journal']
                    pub.journal, _ = PublicationJournal.objects.get_or_create(
                        defaults={"name": journal_fullname, 'slug': journal_slug}, name__iexact=journal_fullname)

                pub.title = row['title']
                pub.authors = row['authors']
                pub.year = row['year'] if not pd.isnull(row['year']) else None
                if (not pd.isnull(row['volume'])) and (not pd.isnull(row['first_page'])) and (not pd.isnull(row['last_page'])):
                    pub.reference = str(
                        row['volume']) + ':' + str(row['first_page']) + '-' + str(row['last_page'])
                pub.web_link = wl
                pub.save()

                pub_db[row['doi']] = pub

        return pub_db

    @staticmethod
    def build_chembl_ligands():
        print("\n===============\n#1 Reading ChEMBL ligand data", datetime.datetime.now())

        ligand_input_file = os.path.join(settings.DATA_DIR, "ligand_data", "assay_data", "chembl_cpds.csv.gz")
        ligand_data = pd.read_csv(ligand_input_file, keep_default_na=False)
        ligand_data.replace(["", "None", "null", "NaN"], np.nan, inplace=True)
        print(f"Found {len(ligand_data)} ligands")

        # Build mask: not NaN, and not empty/whitespace-only
        mask_chembl = (
            Command.helm_chembl['helm_notation'].notna()
            & Command.helm_chembl['helm_notation'].astype(str).str.strip().ne('')
        )
        mask_cid = (
            Command.helm_cid['helm_notation'].notna()
        )
        # Apply mask
        helm_chembl_clean = Command.helm_chembl.loc[mask_chembl]
        helm_cid_clean = Command.helm_cid.loc[mask_cid]
        # And clean CID data
        helm_cid_clean['pubchem_cid'] = helm_cid_clean['pubchem_cid'].apply(
            lambda x: str(int(x)) if pd.notna(x) else np.nan
        )

        # Get WebResource objects
        wr_chembl = WebResource.objects.get(slug="chembl_ligand")
        wr_pubchem = WebResource.objects.get(slug="pubchem")

        # Fetch existing ligand IDs
        print("\n#2 Collecting ChEMBL IDs from existing ligands", datetime.datetime.now())

        existing_ids = set(LigandID.objects.filter(web_resource=wr_chembl).values_list("index", flat=True))
        filtered_ligands = ligand_data[~ligand_data["molecule_chembl_id"].isin(existing_ids)].copy()

        # Set default name
        filtered_ligands["pref_name"].fillna(filtered_ligands["molecule_chembl_id"], inplace=True)

        # First Merge
        merged = pd.merge(
            filtered_ligands,
            helm_chembl_clean,
            on='molecule_chembl_id',
            how='left')

        # Fix the double CID values in the merged df
        merged['pubchem_cid'] = merged['pubchem_cid'].astype(str).str.split(';')
        merged = merged.explode('pubchem_cid').reset_index(drop=True)

        # Second Merge
        merged = pd.merge(
            merged,
            helm_cid_clean,
            on='pubchem_cid',
            how='left')

        # Remove double helm notation columns and collate them
        merged['helm_notation'] = merged['helm_notation_x'].fillna(merged['helm_notation_y'])
        filtered_ligands = merged.drop(columns=['helm_notation_x', 'helm_notation_y'])

        # Fetch additional existing ligand data
        existing_cids = set(LigandID.objects.filter(web_resource=wr_pubchem).values_list("index", flat=True))
        existing_inchis = set(Ligand.objects.exclude(inchikey=None).values_list("inchikey", flat=True))

        smallmol = LigandType.objects.get(slug="small-molecule")
        ligands, weblinks = [], []

        sm_data = filtered_ligands[
            (filtered_ligands["molecule_type"].isin(["Small molecule", "Oligosaccharide"])) |
            (pd.isna(filtered_ligands["molecule_type"]))
        ].reset_index()

        #check that we substitute all the 'nan' with actual np.nan values
        sm_data.replace('nan', np.nan, inplace=True)

        lig_entries = len(sm_data)

        print(f"\n#3 Building {lig_entries} new small-molecule ChEMBL ligands", datetime.datetime.now())

        parent_cache = {
            "smiles": {},
            "inchikey": {},
            "sequence": {},
            "name": {}
        }
        for index, row in sm_data.iterrows():
            insert = True
            chembl_id = row["molecule_chembl_id"]
            ids = [chembl_id]

            # Check other ChEMBL IDs
            if pd.notna(row["other_ids"]) and row["other_ids"]:
                extra_ids = set(row["other_ids"].split(";"))
                existing_matches = extra_ids & existing_ids
                if existing_matches:
                    try:
                        match = LigandID.objects.get(index=next(iter(existing_matches)))
                    except LigandID.MultipleObjectsReturned:
                        match = LigandID.objects.filter(index=next(iter(existing_matches))).first()

                    if match and (match.web_resource == wr_pubchem):
                        # Before saving, only save if it does not already exist:
                        if not LigandID.objects.filter(
                                ligand_id=match.ligand_id,
                                index=chembl_id,
                                web_resource=wr_chembl
                        ).exists():
                            LigandID(
                                index=chembl_id,
                                web_resource=wr_chembl,
                                ligand_id=match.ligand_id
                            ).save()
                            print(f"Found existing non-parent ChEMBL {next(iter(existing_matches))} for parent {chembl_id}")
                        insert = False
                else:
                    ids.extend(extra_ids)

            # Check PubChem CIDs
            if pd.notna(row["pubchem_cid"]) and row["pubchem_cid"]:
                cids = set(row["pubchem_cid"].split(";"))
                existing_matches = cids & existing_cids
                if existing_matches:
                    try:
                        match = LigandID.objects.get(index=next(iter(existing_matches)))
                    except LigandID.MultipleObjectsReturned:
                        match = LigandID.objects.filter(index=next(iter(existing_matches))).first()

                    if match and (match.web_resource == wr_pubchem):
                        # Only create if not already in place:
                        if not LigandID.objects.filter(
                                ligand_id=match.ligand_id,
                                index=chembl_id,
                                web_resource=wr_chembl
                        ).exists():
                            LigandID(
                                index=chembl_id,
                                web_resource=wr_chembl,
                                ligand_id=match.ligand_id
                            ).save()
                        insert = False

            # Check InChIKey
            if insert and row["standard_inchi_key"] in existing_inchis:
                try:
                    ligand = Ligand.objects.get(
                        inchikey=row["standard_inchi_key"],
                        parent__isnull=False
                    )
                except Ligand.DoesNotExist:
                    # If no matching Ligand, skip to next row
                    continue

                # Before saving the new LigandID for the ChEMBL ID:
                if not LigandID.objects.filter(
                        ligand=ligand,
                        index=chembl_id,
                        web_resource=wr_chembl
                    ).exists():
                    LigandID(
                        index=chembl_id,
                        web_resource=wr_chembl,
                        ligand=ligand
                    ).save()

                # Then link any PubChem CIDs, but only if they don’t already exist:
                if pd.notna(row["pubchem_cid"]) and row["pubchem_cid"]:
                    for cid in row["pubchem_cid"].split(";"):
                        if not LigandID.objects.filter(
                                ligand=ligand,
                                index=cid,
                                web_resource=wr_pubchem
                        ).exists():
                            LigandID(
                                index=cid,
                                web_resource=wr_pubchem,
                                ligand=ligand
                            ).save()

                insert = False

            if insert:
                parent = None
                keys = {}
                query = Q(parent__isnull=True)  # This must always be true
                optional_conditions = Q()  # Will hold the OR conditions

                match_fields = []  # Store which field(s) caused a match

                smiles = row.get("smiles")
                if pd.notna(smiles) and smiles:
                    std_smiles = standardize_smiles(smiles)
                    if std_smiles:
                        keys["smiles"] = std_smiles
                        optional_conditions |= Q(smiles=std_smiles)  # OR condition

                inchi = row.get("standard_inchi_key")
                if pd.notna(inchi) and inchi:
                    head_inchi = inchi.split("-")[0]
                    if head_inchi:
                        keys["inchikey"] = head_inchi
                        optional_conditions |= Q(clean_inchikey=head_inchi)  # OR condition

                sequence = row.get("sequence")
                if pd.notna(sequence) and sequence:
                    keys["sequence"] = sequence
                    optional_conditions |= Q(sequence=sequence)  # OR condition

                name = row.get("pref_name")
                if pd.notna(name) and name:
                    keys["name"] = name
                    optional_conditions |= Q(name=name)  # OR condition

                # Final query: (parent__isnull=True) AND (at least one optional condition)
                if optional_conditions:
                    query &= optional_conditions

                # Check cache first
                for key_type, key_value in keys.items():
                    if key_value in parent_cache[key_type]:
                        parent = parent_cache[key_type][key_value]
                        match_fields.append(key_type)
                        break  # Use the first cached match

                # If parent is not found, query the database
                if not parent:
                    parent = Ligand.objects.filter(query).first()
                    if parent:
                        # Store the parent in the cache
                        for key_type, key_value in keys.items():
                            parent_cache[key_type][key_value] = parent

                # Generate parent if still not found
                if not parent:
                    input_ids = {k: v for k, v in keys.items() if v}
                    parent = generate_parent(name, input_ids, "small-molecule")

                    # Store the generated parent in the cache
                    for key_type, key_value in keys.items():
                        parent_cache[key_type][key_value] = parent

                # Create the new ligand and associate it with the parent
                ligand = Ligand(
                    name=row['pref_name'],
                    ambiguous_alias=False,
                    ligand_type=smallmol,  # assuming smallmol is already defined
                    smiles=row.get('smiles'),
                    inchikey=row.get('standard_inchi_key'),
                    sequence=row.get("sequence"),
                    source="ChEMBL_sm",
                    helm=row.get('helm_notation'),
                    parent=parent  # directly assign the parent (which was already created)
                )

                try:
                    input_mol = dm.to_mol(row['smiles'], sanitize=True)
                    if input_mol:
                        # If the ligand's InChIKey wasn't set from the row, compute it.
                        if not ligand.inchikey:
                            ligand.inchikey = dm.to_inchikey(input_mol)
                        ligand.mw = dm.descriptors.mw(input_mol)
                        ligand.rotatable_bonds = dm.descriptors.n_rotatable_bonds(input_mol)
                        ligand.hacc = dm.descriptors.n_hba(input_mol)
                        ligand.hdon = dm.descriptors.n_hbd(input_mol)
                        ligand.logp = dm.descriptors.clogp(input_mol)
                except Exception:
                    pass

                ligands.append(ligand)

                # Add LigandIDs as before, e.g.:
                for val_id in ids:  # or iterate over appropriate keys
                    weblinks.append({
                        "link": LigandID(index=val_id, web_resource=wr_chembl),
                        "lig_idx": len(ligands) - 1
                    })
                if pd.notna(row["pubchem_cid"]) and row["pubchem_cid"]:
                    for cid in row['pubchem_cid'].split(";"):
                        weblinks.append({
                            "link": LigandID(index=cid, web_resource=wr_pubchem),
                            "lig_idx": len(ligands) - 1
                        })

                # Bulk insert every X entries or on the last row
                if len(ligands) == Command.bulk_size or (index == lig_entries - 1):
                    # 1) Insert all pending Ligand objects
                    Ligand.objects.bulk_create(ligands)

                    # 2) Now assign each weblink to its newly‐saved Ligand and check "exists"
                    to_create = []
                    seen_keys = set()  # will hold (ligand_id, index, web_resource_id) tuples

                    for pair in weblinks:
                        # pair["lig_idx"] points into `ligands` list. Because we just bulk‐created,
                        # each `ligands[...]` now has a primary key.
                        ligand_instance = ligands[pair["lig_idx"]]
                        link_obj = pair["link"]

                        # build a “deduplication key” based on the FK IDs and index
                        key = (
                            ligand_instance.id,
                            link_obj.index,
                            # if web_resource is a FK, use its .id; otherwise, use link_obj.web_resource
                            getattr(link_obj.web_resource, "id", link_obj.web_resource),
                        )

                        # Skip immediately if we’ve already queued this exact combination
                        if key in seen_keys:
                            continue

                        # Check if a record with (ligand, index, web_resource) already exists in the DB:
                        already_exists = LigandID.objects.filter(
                            ligand=ligand_instance,
                            index=link_obj.index,
                            web_resource=link_obj.web_resource,
                        ).exists()

                        if not already_exists:
                            # Assign the real Ligand instance and queue for bulk_create
                            link_obj.ligand = ligand_instance
                            to_create.append(link_obj)
                            seen_keys.add(key)

                    # 3) Bulk‐insert only those new (and now‐deduplicated) links
                    if to_create:
                        LigandID.objects.bulk_create(to_create)

                    # -- PROGRESS REPORTING ADDED HERE --
                    completed = index + 1
                    percent = completed / lig_entries * 100
                    print(f"Inserted {completed} of {lig_entries} ligands — {percent:.1f}% complete")

                    # 4) Clear the lists for the next batch
                    ligands = []
                    weblinks = []

        # Parse all new non-small-molecule ChEMBL ligands
        print("\n#4 Building new non-small-molecule ChEMBL ligands", datetime.datetime.now())

        nonsm_data = filtered_ligands[
                        ~((filtered_ligands["molecule_type"].isin(["Small molecule", "Oligosaccharide"])) |
                        (pd.isna(filtered_ligands["molecule_type"])))
                     ].reset_index()

        nonsm_entries = len(nonsm_data)

        print("Found", nonsm_entries, "new non-small-molecules")

        ligands = []
        ligand_types = {"Unknown": "na", "Protein": "protein"}
        weblinks = []
        for _, row in nonsm_data.iterrows():
            nonsm_ids = {}
            if pd.notna(row["smiles"]):
                nonsm_ids["smiles"] = row['smiles']
            if pd.notna(row["sequence"]):
                nonsm_ids["sequence"] = row['sequence']
            if pd.notna(row["standard_inchi_key"]):
                nonsm_ids["inchikey"] = row['standard_inchi_key']
            if pd.notna(row["molecule_chembl_id"]):
                nonsm_ids["chembl_ligand"] = row['molecule_chembl_id']

            # Filter types
            ligand = get_or_create_ligand(row['pref_name'], nonsm_ids, ligand_types[row['molecule_type']], False, True, "ChEMBL_peptide", row.get('helm_notation'))
            # Add LigandIDs
            if pd.notna(row["other_ids"]):
                # 1) Strip and dedupe
                raw_ids = [s.strip() for s in row["other_ids"].split(";")]
                extra_ids = set(raw_ids)
                # 2) Skip if any of these already match existing_ids
                if extra_ids & existing_ids:
                    continue
                # 3) For each trimmed, deduped extra_id, only queue it if
                #    a) Not already in DB
                #    b) Not already in our weblinks list
                for link_index in extra_ids:
                    wr = wr_chembl
                    already_in_db = LigandID.objects.filter(
                        ligand=ligand,
                        index=link_index,
                        web_resource=wr
                    ).exists()
                    already_queued = any(
                        (w.ligand_id == ligand.id and
                         w.index == link_index and
                         w.web_resource_id == wr.id)
                        for w in weblinks
                    )
                    if not (already_in_db or already_queued):
                        weblinks.append(
                            LigandID(
                                ligand=ligand,
                                index=link_index,
                                web_resource=wr
                            )
                        )

            # -- PROGRESS REPORTING ADDED HERE --
            completed = index + 1
            if completed % 200 == 0 or completed == nonsm_entries:
                percent = completed / nonsm_entries * 100
                print(f"Processed {completed} of {nonsm_entries} ligands — {percent:.1f}% complete")

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
        # ids = list(bioactivity_data["parent_molecule_chembl_id"].unique())  # not filtering is way faster
        # Why this was based on LigandID and not on Ligand?
        ligands = list(Ligand.objects.filter(name__startswith="CHEMBL", parent__isnull=False).values_list("id", "name").distinct())
        # ligands = list(Ligand.objects.filter(name__startswith="CHEMBL").values_list("id", "name").distinct())
        # ligands = list(LigandID.objects.filter(index__startswith="CHEMBL").values_list("ligand_id", "index"))
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
        count = 0
        missed = 0
        for index, (_, row) in enumerate(bioactivity_data.iterrows()):
            try:
                if (row["parent_molecule_chembl_id"] in lig_dict.keys()) and (row["Entry name"] in prot_dict.keys()):
                    count += 1
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

                    input_ligand = Ligand.objects.get(id = lig_dict[row["parent_molecule_chembl_id"]])
                    input_prot = Protein.objects.get(id = prot_dict[row["Entry name"]])
                    Command.assign_ligand_target_pairing(input_ligand, input_prot, row["assay_description"], row["standard_type"])

                    try:
                        doi = chembl_document_data.loc[chembl_document_data['document_chembl_id'] == row["document_chembl_id"], 'doi'].iloc[0]
                    except IndexError:
                        response = fetch_from_web_api(
                            url_doc_template, row["document_chembl_id"], cache_dir, xml=True)
                        doi = response[0][0][4].text

                    if doi is not None:
                        publication = publication_array[doi]
                        if publication is not None:
                            pub_links.append(
                                [len(bioacts) - 1, publication.id])
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
                            data_pub(assayexperiment_id=experiment.pk,
                                     publication_id=pub_id)
                            for (experiment, pub_id) in pub_links], ignore_conflicts=True)
                        print("Inserted", index, "out of",
                              bio_entries, "bioactivities")

                        # —— PROGRESS REPORTING ——
                        completed = index + 1
                        percent = completed / bio_entries * 100
                        print(f"Inserted {completed} of {bio_entries} bioactivities — {percent:.1f}% complete")

                        # reset accumulators
                        bioacts = []
                        pub_links = []
                else:
                    missed += 1
                    if row["parent_molecule_chembl_id"] not in lig_dict.keys():
                        print(f'{row["parent_molecule_chembl_id"]} not in lig_dict')
                    elif row["Entry name"] not in prot_dict.keys():
                        print(f'{row["Entry name"]} not in prot_dict')
            except KeyError:
                continue

    @staticmethod
    def build_pubchem_vendor_links():
        LigandVendors.objects.all().delete()
        print("\n===============\n#1 Reading and creating Vendors")
        vendor_url = os.sep.join(
            [settings.DATA_DIR, "ligand_data", "assay_data", "pubchem_vendor_list.csv.gz"])
        vendor_data = pd.read_csv(vendor_url, dtype=str)
        vendors = []
        for _, row in vendor_data.iterrows():
            vendors.append(LigandVendors(slug=slugify(
                row["SourceName"]), name=row["SourceName"], url=row["SourceURL"]))
        LigandVendors.objects.bulk_create(vendors)
        vendor_dict = {vendor.name: vendor.pk for vendor in vendors}

        print("\n#2 Building ChEMBL ligands cache", datetime.datetime.now())
        ligands = list(LigandID.objects.filter(
            index__startswith="CHEMBL").values_list("ligand_id", "index"))
        lig_dict = {entry[1]: entry[0] for entry in ligands}

        print("\n#3 Creating all vendor links", datetime.datetime.now())
        vendor_links_url = os.sep.join(
            [settings.DATA_DIR, "ligand_data", "assay_data", "pubchem_vendor_links.csv.gz"])
        vendor_links_data = pd.read_csv(vendor_links_url, dtype=str)
        links = []
        for _, row in vendor_links_data.iterrows():
            if len(row["SourceRecordURL"])<=400 and len(row["RegistryID"])<=500:
                links.append(LigandVendorLink(
                    vendor_id=vendor_dict[row["SourceName"]], ligand_id=lig_dict[row["chembl_id"]], url=row["SourceRecordURL"], external_id=row["RegistryID"]))

        LigandVendorLink.objects.bulk_create(links)

    @staticmethod
    def uniprot_mapper(protein, organism):
        organism_dict = {'PIG': 'sus_scrofa', 'RAT': 'rattus_norvegicus', 'HUMAN': 'homo_sapiens', 'MOUSE': 'mus_musculus',
                         'CANINE': 'canis_lupus_familiaris', 'BOVINE': 'bos_taurus', 'CALF': 'bos_taurus', 'COW': 'bos_taurus',
                         'GUINEA PIG': 'cavia_porcellus', 'CAT': 'felis_catus', 'NEONATAL RAT': 'rattus_norvegicus',
                         '? HUMAN': 'homo_sapiens', 'OPOSSUM': 'didelphis_marsupialis', 'Rat 6B': 'rattus_norvegicus',
                         'HUMAN M3': 'homo_sapiens', 'HUMAN M4': 'homo_sapiens', 'Chick': 'gallus_gallus',
                         'Frog': 'pseudis_balbodactyla', 'Newborn rats': 'rattus_norvegicus', 'Beef': 'bos_taurus',
                         'Sheep': 'ovis_aries', 'OX': 'bos_taurus', 'Dog': 'canis_lupus_familiaris',
                         'Rhesus': 'macaca_mulatta', 'Monkey': 'macaca_mulatta', 'PIGLET': 'sus_scrofa',
                         'Rat Y861': 'rattus_norvegicus', 'Zebra Finch': 'taeniopygia_guttata', 'Chicken': 'gallus_gallus',
                         'MICE': 'mus_musculus', 'Rhesus Monkey': 'macaca_mulatta', 'Zebrafish': 'danio_rerio'}
        if organism in organism_dict.keys():
            query = 'gene_exact:{0}+AND+organism_name:{1}'.format(
                urllib.parse.quote(protein.lower()), organism_dict[organism])
        else:
            query = 'gene_exact:{}'.format(urllib.parse.quote(protein.lower()))
        if query not in Command.mapper_cache.keys():
            url = 'https://rest.uniprot.org/uniprotkb/search?query={}&fields=id&format=tsv'.format(query)
            req = urllib.request.Request(url)
            try:
                converted = urllib.request.urlopen(
                    req).read().decode('utf-8').split('\n')[1].lower()
            except IndexError:
                converted = None

            Command.mapper_cache[query] = converted

        return Command.mapper_cache[query]

    @staticmethod
    def get_ligands_data(ligands, complete_ligands, ligand_mapping, ligand_interactions=pd.DataFrame(), target_ids=False):
        full_info = ['ligand_id', 'name', 'species', 'type', 'approved', 'withdrawn', 'labelled', 'radioactive', 'pubchem_sid', 'pubchem_cid',
                     'uniprot_id', 'iupac_name', 'inn', 'synonyms', 'smiles', 'inchikey', 'inchi', 'gtoimmupdb', 'gtompdb']
        bioactivity_info = ['ligand_id', 'action', 'target', 'target_id', 'target_species',
                            'target_gene_symbol', 'target_uniprot_id', 'original_affinity_relation',
                            'action_comment', 'selectivity', 'primary_target',
                            'concentration_range', 'affinity_units', 'affinity_high',
                            'affinity_median', 'affinity_low', 'assay_description', 'pubmed_id']

        weblinks = ['ligand_id', 'chembl_id', 'chebi_id',
                    'cas', 'drugbank_id', 'drug_central_id']

        ligand_data = complete_ligands.loc[complete_ligands['ligand_id'].isin(
            ligands), full_info]
        ligand_weblinks = ligand_mapping.loc[ligand_mapping['ligand_id'].isin(
            ligands), weblinks]
        ligand_complete = ligand_data.merge(ligand_weblinks, on="ligand_id")

        if not ligand_interactions.empty:
            bioactivity_data = ligand_interactions.loc[ligand_interactions['ligand_id'].isin(
                ligands), bioactivity_info]
            if target_ids is not False and len(target_ids) > 0:
                bioactivity_data = bioactivity_data.loc[bioactivity_data["target_id"].isin(
                    target_ids)]
            ligand_complete = ligand_complete.merge(
                bioactivity_data, on="ligand_id")

        return ligand_complete

    @staticmethod
    def obtain_ligands(data, compare_set, labels):
        interactions_targets = list(data[labels[0]].unique())
        targets_with_ligands = set(
            interactions_targets).intersection(set(compare_set))
        ligands = list(data.loc[data[labels[0]].isin(
            targets_with_ligands), labels[1]].unique())
        ligands = [x for x in ligands if x == x]  # remove nan
        return ligands

    @staticmethod
    def compare_proteins(gtp_data):
        gpcrdb_proteins = Protein.objects.filter(
            family__slug__startswith="00", sequence_type__slug="wt").values_list('entry_name', 'accession')
        entries = gtp_data.loc[gtp_data['uniprotkb_id'].isin([protein[1].split(
            "-")[0] for protein in gpcrdb_proteins]), ['uniprotkb_id', 'gtopdb_iuphar_id']]
        return list(entries['gtopdb_iuphar_id'].unique())

    @staticmethod
    def classify_assay(biodata, type_column, description_column):
        # starting to assess assay type
        biodata['assay_type'] = 'U'
        biodata[description_column].fillna('Unclassified', inplace=True)
        # find Binding assays (displacement, affinity, ) KI and KD are always binding
        biodata.loc[biodata[type_column].isin(
            ['pKi', 'pKd', 'Ki', 'Kd']), 'assay_type'] = 'B'
        # EC50, pKB and pA2 are always Functional
        biodata.loc[biodata[type_column].isin(
            ['pKB', 'pEC50', 'pA2', 'A2', 'Kb', 'KB', 'EC50']), 'assay_type'] = 'F'
        # IC50 can be both B: displacement, affinity, binding, radioligand else is F
        binding_words = ['isplacement', 'ffinity', 'inding', 'adioligand']
        biodata.loc[(biodata[type_column] == 'pIC50') & (
            biodata[description_column].str.contains('|'.join(binding_words))), 'assay_type'] = 'B'
        biodata.loc[(biodata[type_column] == 'IC50') & (
            biodata[description_column].str.contains('|'.join(binding_words))), 'assay_type'] = 'B'
        biodata.loc[(biodata[type_column] == 'pIC50') & ~(
            biodata[description_column].str.contains('|'.join(binding_words))), 'assay_type'] = 'F'
        biodata.loc[(biodata[type_column] == 'IC50') & ~(
            biodata[description_column].str.contains('|'.join(binding_words))), 'assay_type'] = 'F'
        # find Unclassified assayas
        # biodata.loc[biodata[description_column].str.contains(
        #     'Unclassified'), 'assay_type'] = 'U'

        return biodata

    @staticmethod
    def save_the_ligands_save_the_world(lig_df, pep_df):
        types_dict = {
            'Inorganic': 'small-molecule',
            'Metabolite': 'small-molecule',
            'Natural product': 'small-molecule',
            'Peptide': 'peptide',
            'Synthetic organic': 'small-molecule',
            'Antibody': 'protein',
            None: 'na'
        }
        weblink_keys = {
            'smiles': 'smiles',
            'inchikey': 'inchikey',
            'pubchem': 'pubchem_cid',
            'gtoplig': 'ligand_id',
            'chembl_ligand': 'chembl_id',
            'drugbank': 'drugbank_id',
            'drug_central': 'drug_central_id'
        }

        # Remove radioactive ligands
        lig_df = lig_df[lig_df['radioactive'] != 'yes']

        # Convert all columns to string, then replace 'nan' with None
        lig_df = lig_df.astype(str).replace({'nan': None})
        pep_df = pep_df.astype(str).replace({'nan': None})

        total = len(lig_df)
        if total == 0:
            print("No ligands to process.")
            return

        print(f"Starting processing {total} ligands...")
        for idx, (_, row) in enumerate(lig_df.iterrows(), start=1):
            # print progress at start, every 100 rows, and at the end
            if idx == 1 or idx % 100 == 0 or idx == total:
                pct = idx / total * 100
                print(f"Processed {idx}/{total} ({pct:.1f}%)")

            # Skip if there's no name
            if not row['name']:
                continue

            # Build dictionary of IDs
            ids = {}
            for key, col_name in weblink_keys.items():
                raw_val = row.get(col_name)
                if raw_val is not None:
                    raw_val = str(raw_val)
                    if is_float(raw_val):
                        raw_val = str(int(float(raw_val)))
                    ids[key] = raw_val

            ligand_type = types_dict.get(row['type'], 'na')

            if ligand_type != "small-molecule":
                mask = (
                    (pep_df["ligand_id"] == row["ligand_id"]) &
                    (pep_df["species"] == row["species"])
                )
                if mask.any():
                    peptide_entries = pep_df[mask]
                    if len(peptide_entries) > 1:
                        for _, pep_row in peptide_entries.iterrows():
                            seq = pep_row["single_letter_amino_acid_sequence"]
                            if seq and len(seq) < 1000:
                                ids["sequence"] = seq
                            uniprot = pep_row["uniprot_id"]
                            if uniprot:
                                ids["uniprot"] = uniprot
                            new_name = row['name']
                            if row['species']:
                                new_name = f"{row['name']} ({row['species']})"
                            ligand = get_or_create_ligand(new_name, ids, ligand_type, True, False, 'GuideToPharma', pep_row['helm_notation'])
                            if ligand is None:
                                print("Issue with", row['name'])
                                print(row)
                                break
                    else:
                        pep_row = peptide_entries.iloc[0]
                        seq = pep_row["single_letter_amino_acid_sequence"]
                        if seq and len(seq) < 1000:
                            ids["sequence"] = seq
                        uniprot = pep_row["uniprot_id"]
                        if uniprot:
                            ids["uniprot"] = uniprot
                        new_name = row['name']
                        if row['species']:
                            new_name = f"{row['name']} ({row['species']})"
                        ligand = get_or_create_ligand(new_name, ids, ligand_type, True, False, 'GuideToPharma', pep_row['helm_notation'])
                        if ligand is None:
                            print("Issue with", row['name'])
                            print(row)
                else:
                    new_name = row['name']
                    ligand = get_or_create_ligand(new_name, ids, ligand_type, True, False, 'GuideToPharma')
                    if ligand is None:
                        print("Issue with", row['name'])
                        print(row)
            else:
                ligand = get_or_create_ligand(row['name'], ids, ligand_type, True, False, 'GuideToPharma', row['helm_notation'])
                if ligand is None:
                    print("Issue with", row['name'])
                    print(row)

    @staticmethod
    def build_gtp_bioactivities(gtp_biodata):
        print("# Start parsing the GTP Dataframe")
        for _, row in gtp_biodata.iterrows():
            receptor = Command.fetch_protein(
                row['target_id'], 'GtoP', row['target_species'])
            # TODO Handle multiple matches (uniprot filter?)
            ligand = get_ligand_by_id("gtoplig", row['ligand_id'])

            try:
                low_value = "{:.2f}".format(float(row['affinity_low']))
            except ValueError:
                low_value = 'None'

            try:
                high_value = "{:.2f}".format(float(row['affinity_high']))
            except ValueError:
                high_value = 'None'

            if row['affinity_median'] != 'None':
                activity_value = "{:.2f}".format(float(row['affinity_median']))
            elif row['affinity_high'] != 'None':
                if row['affinity_low'] != 'None':
                    activity_value = "{:.2f}".format(statistics.mean(
                        [float(row['affinity_high']), float(row['affinity_low'])]))
                else:
                    activity_value = "{:.2f}".format(
                        float(row['affinity_high']))
            else:
                activity_value = 'None'

            ranges = '|'.join([low_value, activity_value, high_value])

            # Adding publications from the PMIDs section
            try:
                pmids = row['pubmed_id'].split('|')
            except AttributeError:
                pmids = None

            if (receptor is not None) and (ligand is not None):
                # last step because it requires multiple uploads in case we have multiple species
                gtp_data = AssayExperiment(
                    ligand=ligand,
                    protein=receptor,
                    assay_type=row['assay_type'],
                    assay_description=row['assay_description'],
                    standard_activity_value=None,
                    p_activity_value=activity_value,
                    p_activity_ranges=ranges,
                    standard_relation=row['original_affinity_relation'],
                    value_type=row['affinity_units'],
                    source='Guide to Pharmacology',
                    document_chembl_id=None,
                )
                gtp_data.save()
                try:
                    for pmid in pmids:
                        publication = Command.fetch_publication(pmid)
                        gtp_data.publication.add(publication)
                except:
                    publication = None

                Command.assign_ligand_target_pairing(ligand, receptor, row['assay_description'], row['affinity_units'])
            else:
                print("SKIPPING", ligand, row["ligand_id"], "|",
                      receptor, row['target_id'], row['target_species'])

    @staticmethod
    def fetch_publication(publication_doi):
        """
        fetch publication with Publication model
        requires: publication doi or pmid

        """
        if pd.isna(publication_doi) is True:
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
                wl = WebLink.objects.get(
                    index=publication_doi, web_resource__slug=pub_type)
            except WebLink.DoesNotExist:
                try:
                    wl = WebLink.objects.create(
                        index=publication_doi, web_resource=WebResource.objects.get(slug=pub_type))
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
                    # if something off with publication, skip.
                    print("Publication fetching error | module: fetch_publication. Row # is : " +
                          str(publication_doi) + ' ' + pub_type)

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
                if species is None or species == "None":
                    # Sorting by species => human first, otherwise next species in line
                    # TODO => potentially capture all species with GtP ID
                    return Protein.objects.filter(web_links__index=target, web_links__web_resource__slug='gtop').order_by("species_id").first()
                else:
                    prots = Protein.objects.filter(
                        web_links__index=target, web_links__web_resource__slug='gtop', species__common_name__iexact=species)
                    if prots.count() > 0:
                        return prots.first()
                    else:
                        receptor_fam = list(Protein.objects.filter(
                            web_links__index=target, web_links__web_resource__slug='gtop').values_list("family_id", flat=True))[0]
                        return Protein.objects.get(family_id=receptor_fam, species__common_name__iexact=species)
            except:
                return None
        elif database == 'PDSP':
            try:
                test = None
                if len(Protein.objects.filter(entry_name=target))>0:
                    test = Protein.objects.get(entry_name=target)
                elif len(Protein.objects.filter(web_links__index=target, web_links__web_resource__slug='uniprot'))>0:
                    test = Protein.objects.get(web_links__index=target, web_links__web_resource__slug='uniprot')
                return test
            except:
                return None
        elif database == 'DrugCentral':
            try:
                protein = Protein.objects.get(accession=target)
                return protein
            except:
                return None

    @staticmethod
    def fetch_role(drug_type, drug_action):
        # This still need to be addressed and fixed
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
            lr = find_role(query)
            # role_slug = slugify(query)
            # lr, _ = LigandRole.objects.get_or_create(slug=role_slug, defaults={'name': query})
        if lr:
            return lr[0]
        else:
            return lr

    @staticmethod
    def fetch_species(ligand_species, target_species):
        try:
            if ligand_species == 'Same as target':
                species = Species.objects.get(common_name=target_species)
            elif ligand_species is None:
                species = None
            else:
                species = Species.objects.get(common_name=ligand_species)
            return species
        except:
            return None

    @staticmethod
    def fetch_from_accession(code):
        try:
            protein = Protein.objects.filter(accession=code)
            test = protein.get()
        except:
            test = None
        return test

    @staticmethod
    def build_kidatabase_bioactivities():
        protein_names = {}
        ligand_cache = {}
        print("\n===============\n#1 Reading PDSP bioacitivity data")
        pdsp_link = get_or_create_url_cache(
            "https://pdsp.unc.edu/databases/kiDownload/download.php", 7 * 24 * 3600)
        bioactivity_kidata = pd.read_csv(pdsp_link, dtype=str, encoding='mac_roman')
        # Keeping data that has either SMILES info OR CAS info
        # CAS number can be translated into pubchem CID
        bioactivity_data_filtered = bioactivity_kidata.loc[(
            ~bioactivity_kidata['SMILES'].isnull()) | (~bioactivity_kidata['CAS'].isnull())]
        bioactivity_data_filtered = bioactivity_data_filtered.loc[(
            ~bioactivity_data_filtered['Unigene'].isnull())]
        bioactivity_data_filtered.fillna('None', inplace=True)
        bio_entries = len(bioactivity_data_filtered)
        print("\n===============\n#2 Start parsing PDSP data")
        bioacts = []
        for index, (_, row) in enumerate(bioactivity_data_filtered.iterrows()):
            receptor = None
            label = '_'.join([row['Unigene'], row['species']])
            if label not in protein_names.keys():
                protein = Command.uniprot_mapper(row['Unigene'], row['species'])
                if protein is not None:
                    protein_names[label] = Command.fetch_protein(protein, 'PDSP')

            ids = {}
            if row['SMILES'] != 'None':
                ids['smiles'] = row['SMILES']
                input_mol = dm.to_mol(row['SMILES'], sanitize=True)
                if input_mol:
                # Check if InChIKey has been set
                    ids['inchikey'] = dm.to_inchikey(input_mol)
            if row['CAS'] != 'None':
                ids['CAS'] = row['CAS']
            if row[' Ligand Name'] not in ligand_cache.keys():
                print(ids)
                ligand = get_or_create_ligand(row[' Ligand Name'], ids, source='PDSP KiDatabase')
                ligand_cache[row[' Ligand Name']] = ligand
            if label in protein_names.keys():
                receptor = protein_names[label]
            if (receptor is not None) and (ligand_cache[row[' Ligand Name']] is not None):
                bioacts.append(AssayExperiment())
                bioacts[-1].ligand_id = ligand_cache[row[' Ligand Name']].id
                bioacts[-1].protein_id = receptor.id
                bioacts[-1].assay_type = 'B'
                bioacts[-1].assay_description = None
                bioacts[-1].standard_activity_value = round(float(row['ki Val']), 2)
                bioacts[-1].p_activity_value = round(-math.log10(float(row['ki Val']) * 1e-9), 2)
                bioacts[-1].p_activity_ranges = None
                bioacts[-1].standard_relation = '='
                bioacts[-1].value_type = 'pKi'
                bioacts[-1].source = 'PDSP KiDatabase'
                bioacts[-1].document_chembl_id = None
                bioacts[-1].reference_ligand = row['Hotligand']
                # BULK insert every X entries or last entry
            if (len(bioacts) == Command.bulk_size) or (index == bio_entries - 1):
                AssayExperiment.objects.bulk_create(bioacts)
                print("Inserted", index, "out of",
                      bio_entries, "bioactivities")
                bioacts = []

            Command.assign_ligand_target_pairing(ligand_cache[row[' Ligand Name']], receptor, None, 'pKi')

    @staticmethod
    def build_drugcentral_bioactivities():
        print("# Collecting DrugCentral data")
        accession_numbers = {}
        ligand_cache = {}
        drugcentral_ligands_link = get_or_create_url_cache("https://unmtid-shinyapps.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz", 7 * 24 * 3600)
        # drugcentral_ligands_link = "https://unmtid-shinyapps.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz"
        drugcentral_smiles_link = get_or_create_url_cache("https://unmtid-shinyapps.net/download/DrugCentral/2021_09_01/structures.smiles.tsv", 7 * 24 * 3600)
        # drugcentral_smiles_link = "https://unmtid-shinyapps.net/download/DrugCentral/2021_09_01/structures.smiles.tsv"
        drugcentral_ligands = pd.read_csv(drugcentral_ligands_link, sep='\t', header=0, compression='gzip')
        drugcentral_smiles = pd.read_csv(drugcentral_smiles_link, sep='\t', header=0)
        #Adjusting the data and filter
        drugcentral_ligands['ID'] =  drugcentral_ligands['STRUCT_ID']
        drugcentral_ligands.drop('STRUCT_ID', axis=1, inplace=True)
        drugcentral_ligands = drugcentral_ligands.loc[drugcentral_ligands['TARGET_CLASS'] == 'GPCR']
        merged_data = pd.merge(drugcentral_ligands, drugcentral_smiles, on='ID', how='left')
        # Keeping data that has either SMILES info OR CAS info
        # CAS number can be translated into pubchem CID
        merged_data_filtered = merged_data.loc[(~merged_data['SMILES'].isnull()) | (~merged_data['CAS_RN'].isnull())]
        merged_data_filtered = merged_data_filtered.loc[(~merged_data_filtered['GENE'].isnull())]
        merged_data_filtered = merged_data_filtered.loc[(~merged_data_filtered['ACT_VALUE'].isnull())]
        merged_data_filtered.fillna('None', inplace=True)
        merged_data_filtered = Command.classify_assay(merged_data_filtered, 'ACT_TYPE', 'ACT_COMMENT')
        bio_entries = len(merged_data_filtered)

        print("# Parsing DrugCentral data")
        bioacts = []
        for index, (_, row) in enumerate(merged_data_filtered.iterrows()):
            receptor = None
            code = row['ACCESSION']
            if code not in accession_numbers.keys():
                protein = Command.fetch_protein(code, 'DrugCentral')
                if protein is not None:
                    accession_numbers[code] = protein

            ids = {}
            if row['SMILES'] != 'None':
                ids['smiles'] = row['SMILES']
            # if row['CAS_RN'] != 'None':
            #     ids['CAS'] = row['CAS_RN']
            if row['InChIKey'] != 'None':
                ids['inchikey'] = row['InChIKey']
            if row['DRUG_NAME'] not in ligand_cache.keys():
                ligand = get_or_create_ligand(row['DRUG_NAME'], ids, source='Drug Central')
                ligand_cache[row['DRUG_NAME']] = ligand
            if code in accession_numbers.keys():
                receptor = accession_numbers[code]
            if (receptor is not None) and (ligand_cache[row['DRUG_NAME']] is not None):
                calc_val = round(-math.log10(float(row['ACT_VALUE']) * 1e-9), 2) if row['ACT_TYPE'] != 'pA2' else round(float(row['ACT_VALUE']), 2)
                bioacts.append(AssayExperiment())
                bioacts[-1].ligand_id = ligand_cache[row['DRUG_NAME']].id
                bioacts[-1].protein_id = receptor.id
                bioacts[-1].assay_type = row['assay_type']
                bioacts[-1].assay_description = row['ACT_COMMENT']
                bioacts[-1].standard_activity_value = round(float(row['ACT_VALUE']), 2) if row['ACT_TYPE'] != 'pA2' else None
                bioacts[-1].p_activity_value = calc_val
                bioacts[-1].p_activity_ranges = None
                bioacts[-1].standard_relation = row['RELATION']
                bioacts[-1].value_type = 'p'+row['ACT_TYPE'] if row['ACT_TYPE'] != 'pA2' else row['ACT_TYPE']
                bioacts[-1].source = 'Drug Central'
                bioacts[-1].document_chembl_id = None
                # BULK insert every X entries or last entry
            if (len(bioacts) == Command.bulk_size) or (index == bio_entries - 1):
                AssayExperiment.objects.bulk_create(bioacts)
                print("Inserted", index, "out of",
                      bio_entries, "bioactivities")
                bioacts = []
            value_type = 'p'+row['ACT_TYPE'] if row['ACT_TYPE'] != 'pA2' else row['ACT_TYPE']

            Command.assign_ligand_target_pairing(ligand_cache[row['DRUG_NAME']], receptor, row['ACT_COMMENT'], value_type)

    @staticmethod
    def build_drugbank_ligands():
        print("# Collecting Drug Bank data")
        data_link = os.sep.join([settings.DATA_DIR, 'ligand_data', 'assay_data', 'structure_links.csv'])

        data = pd.read_csv(data_link)
        ligand_cache = {}

        data.fillna('None', inplace=True)
        # Keep only rows with at least one relevant ID
        filtered = data.loc[
            (data['SMILES'] != 'None') |
            (data['CAS Number'] != 'None') |
            (data['InChIKey'] != 'None') |
            (data['PubChem Compound ID'] != 'None')
        ]

        total = len(filtered)
        if total == 0:
            print("No DrugBank ligands to process.")
            return

        print(f"# Parsing Drug Bank data ({total} ligands)…")
        for idx, row in enumerate(filtered.itertuples(index=False), start=1):
            # progress update
            if idx == 1 or idx % 100 == 0 or idx == total:
                pct = idx / total * 100
                print(f"Processed {idx}/{total} ({pct:.1f}%)")

            ids = {}
            if row.SMILES != 'None':
                ids['smiles'] = row.SMILES
            # if row['CAS Number'] != 'None':
            #     ids['CAS'] = row['CAS Number']
            if row.InChIKey != 'None':
                ids['inchikey'] = row.InChIKey
            #PubChem ID
            if row._10 != 'None':
                ids['pubchem'] = int(row._10)

            name = row.Name
            if name not in ligand_cache:
                ligand = get_or_create_ligand(name, ids, source='DrugBank')
                ligand_cache[name] = ligand

    @staticmethod
    def calculate_potency_and_affinity():
        ligand_target_couples = AssayExperiment.objects.exclude(p_activity_value='None').values_list('ligand_id',
                                                                                                     'protein_id',
                                                                                                     'value_type',
                                                                                                     'p_activity_value').distinct()
        connections = {}
        SI_dict = {}
        B_values = ['pKi', 'pKd']
        F_values = ['pEC50', 'pIC50', 'pA2', 'pKB', 'pKb', 'Potency', 'pAC50']
        for pair in ligand_target_couples:
            value_type = pair[2]
            if not pair[2].startswith(('p','P')):
                value_type = 'p'+pair[2]
            if value_type in B_values:
                sample = 'Affinity'
            if value_type in F_values:
                sample = 'Potency'
            if pair[0] not in connections.keys():
                connections[pair[0]] = {}
            if pair[1] not in connections[pair[0]].keys():
                connections[pair[0]][pair[1]] = {}
            if sample not in connections[pair[0]][pair[1]].keys():
                connections[pair[0]][pair[1]][sample] = []
            connections[pair[0]][pair[1]][sample].append(float(pair[3]))

        for ligand in connections:
            for target in connections[ligand]:
                for value in connections[ligand][target]:
                    connections[ligand][target][value] = round(statistics.mean(connections[ligand][target][value]),2)

        #Expected: ligand[target] = SI
        #ligand = 212224
        for ligand in connections:
            SI_dict[ligand] = {}
            Max_Affinity = []
            Max_Potency = []
            # target_count = len(connections[ligand].keys())
            #we have a list of targets for the ligand now
            #need to assess target with highest B and highest F
            for target in connections[ligand].keys():
                for measurement in connections[ligand][target].keys():
                    if measurement == 'Affinity':
                        Max_Affinity.append((target, connections[ligand][target][measurement]))
                    elif measurement == 'Potency':
                        Max_Potency.append((target, connections[ligand][target][measurement]))
            #Generate sorted lists
            affinity_sort = sorted(Max_Affinity, key=lambda x: x[1], reverse=True)
            potency_sort = sorted(Max_Potency, key=lambda x: x[1], reverse=True)
            SI_dict[ligand] = {"Affinity Count": len(affinity_sort),
                               "Potency Count": len(potency_sort)}
            #Max_B and Max_F sets the reference GPCR and
            #need to assess target with highest B and highest F
            for target in connections[ligand].keys():
                SI_dict[ligand][target] = {}
                if (len(affinity_sort) > 1) and ('Affinity' in connections[ligand][target].keys()):
                    if target == affinity_sort[0][0]:
                        va = connections[ligand][target]['Affinity'] - affinity_sort[1][1]
                        BAssay = int(10**(abs(va)))
                        SI_dict[ligand][target]['Affinity'] = BAssay
                    else:
                        va = connections[ligand][target]['Affinity'] - affinity_sort[0][1]
                        BAssay = int(10**(abs(va)))
                        SI_dict[ligand][target]['Affinity'] = -BAssay
                if (len(potency_sort) > 1) and ('Potency' in connections[ligand][target].keys()):
                    if target == potency_sort[0][0]:
                        va = connections[ligand][target]['Potency'] - potency_sort[1][1]
                        FAssay = int(10**(abs(va)))
                        SI_dict[ligand][target]['Potency'] = FAssay
                    else:
                        va = connections[ligand][target]['Potency'] - potency_sort[0][1]
                        FAssay = int(10**(abs(va)))
                        SI_dict[ligand][target]['Potency'] = -FAssay
        #This step is kinda slow and may be sped up using bulk_update
        #but I don't know if you can apply bulk_update with filters
        for lig in SI_dict:
            for tg in list(SI_dict[lig].keys())[2:]:
                affinity = SI_dict[lig][tg]['Affinity'] if 'Affinity' in SI_dict[lig][tg].keys() else '-'
                potency = SI_dict[lig][tg]['Potency'] if 'Potency' in SI_dict[lig][tg].keys() else '-'
                AssayExperiment.objects.filter(ligand_id=lig, protein_id=tg).update(
                    affinity=affinity,
                    potency=potency,
                    count_potency_test=SI_dict[lig]['Potency Count'],
                    count_affinity_test=SI_dict[lig]['Affinity Count']
                )

    @staticmethod
    def expand_row(row):
        """
        Given one row with species_list and uniprot_list:
          - If lengths match, explode 1:1.
          - If one list is empty, replicate the other list with "".
          - If both lists are non-empty but lengths differ, pad the smaller one with "".
        Returns a list of dictionaries with the expanded data.
        """
        row_dict = row.to_dict()
        sp_list = row_dict["species_list"]
        up_list = row_dict["uniprot_list"]

        # If both are empty, just return the single row
        if not sp_list and not up_list:
            return [row_dict]

        # If one is empty, replicate with ""
        if not sp_list and up_list:
            sp_list = [""] * len(up_list)
        elif not up_list and sp_list:
            up_list = [""] * len(sp_list)

        # Pad if lengths differ
        if len(sp_list) != len(up_list):
            if len(sp_list) < len(up_list):
                sp_list += [""] * (len(up_list) - len(sp_list))
            else:
                up_list += [""] * (len(sp_list) - len(up_list))

        expanded = []
        for s_val, u_val in zip(sp_list, up_list):
            new_row = dict(row_dict)
            # Store single values (not as lists)
            new_row["species"] = s_val
            new_row["uniprot_id"] = u_val
            expanded.append(new_row)
        return expanded

    @staticmethod
    def get_species_from_uniprot(uniprot_id):
        """
        Returns the scientific name for the organism from a UniProt entry, or None if not found.
        """
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        response = requests.get(url)

        if response.status_code == 200:
            data = response.json()
            # 'organism' entry usually contains at least 'scientificName'
            organism_info = data.get("organism", {})
            return organism_info.get("commonName", None)
        else:
            # Could not retrieve, or ID not found
            return None

    @staticmethod
    def assign_species_to_peptide():
        duplicates = (
            Ligand.objects
            .values('name')
            .annotate(name_count=Count('name'))
            .filter(name_count__gt=1, ligand_type=2)
        )

        for duplicate in duplicates:
            # Get actual model instances, not dictionary values
            ligands = Ligand.objects.filter(name=duplicate['name'])

            for lig in ligands:
                species = Command.get_species_from_uniprot(lig.uniprot)
                if species:
                    lig.name = lig.name + ' (' + species + ')'
                    lig.save()

    @staticmethod
    def assign_ligand_target_pairing(ligand, receptor, assay_description, unit):
        """
        Create a LigandTargetPairing linking `ligand` and `receptor` with:
          - effect determined by `unit` (and, for AC50, by keywords in `assay_description`)
          - role determined by `assay_type` slug lookup in LigandType
        """
        # 1. Determine effect slug
        effect_slug = 'unknown'
        if unit in ['EC50', 'pEC50']:
            effect_slug = 'stimulatory'
        elif unit in ['pA2', 'pKB', 'IC50', 'pIC50', 'pKb']:
            effect_slug = 'inhibitory'
        elif unit in ['pKi', 'pKd', 'Ki', 'Kd']:
            effect_slug = 'binding'
        elif unit == 'AC50':
            desc = assay_description.lower()
            # stimulatory if “agonist” but not part of “antagonist”
            if 'agonist' in desc and 'antagonist' not in desc:
                effect_slug = 'stimulatory'
            # inhibitory if either “antagonist” or “inhibitor”
            elif 'antagonist' in desc or 'inhibitor' in desc:
                effect_slug = 'inhibitory'
        else:
            print(f"This value was not in the selection: {unit}")
        # 2. Lookup or fall back
        effect = None
        if effect_slug:
            try:
                effect = LigandEffect.objects.get(slug=effect_slug)
            except LigandEffect.DoesNotExist:
                logger.warning("No LigandEffect with slug=%r found; leaving effect NULL", effect_slug)
        # 3. Create + save
        with transaction.atomic():
            pairing = LigandTargetPairing(
                ligand=ligand,
                target=receptor,
                effect=effect,
            )
            pairing.save()

        return pairing
