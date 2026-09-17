from django.core.management.base import CommandError
from build.management.commands.base_build import Command as BaseBuild
from build.management.commands.build_ligand_functions import get_ligand_by_id, match_id_via_unichem, get_or_create_ligand, get_or_create_ligands_batch, is_float, standardize_smiles, generate_parent, apply_canonical_ligand_types, classify_lipids, allocate_next_gpcrdb_id, load_gtp_source_data
from django.conf import settings
from django.utils.text import slugify
from django.db import IntegrityError, transaction, connection
from django.db.models import Count, Q

from common.tools import get_or_create_url_cache, fetch_from_web_api, test_model_updates, find_role, dump_checker, generate_helm_universal, parse_cycl_pos
from common.models import WebLink, WebResource, Publication, PublicationJournal
from ligand.models import Ligand, LigandID, LigandType, LigandVendors, LigandVendorLink, AssayExperiment, Endogenous_GTP, LigandRole, LigandEffect, LigandTargetPairing, AssayClassification
from protein.models import Protein, Species

import django.apps
import requests
import math
import os
import csv
import re
import time
import gzip
import statistics
import datamol as dm
import datetime
import pandas as pd
import numpy as np
import urllib.parse
import urllib.request
import urllib.error
import socket
import ssl
from io import StringIO
from typing import List, Tuple, Optional, Dict
from rdkit import Chem
from chembl_structure_pipeline import standardizer
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
    dump_dir = os.sep.join([settings.DATA_DIR, 'ligand_data', 'model_snapshots'])
    helm_chembl_filepath = os.sep.join([data_dir, 'HELM_CHEMBL.csv'])
    helm_cid_filepath = os.sep.join([data_dir, 'HELM_CID.csv'])
    helm_chembl = pd.read_csv(helm_chembl_filepath, index_col=0)
    helm_cid = pd.read_csv(helm_cid_filepath, index_col=0)
    ### These are here as a safety valve in case you want to manually check
    #ligand_dump = os.sep.join([dump_dir, 'ligand_reload_dump.csv'])
    #id_dump = os.sep.join([dump_dir, 'ligandid_reload_dump.csv'])
    #######################################################################
    ligand_dump = dump_checker('ligand.Ligand')
    ligand_csv = pd.read_csv(ligand_dump['latest_dump'], sep=';', index_col=0)
    id_dump = dump_checker('ligand.LigandID')
    id_csv = pd.read_csv(id_dump['latest_dump'], sep=';', index_col=0)

    #Fetch assay classification dictionary - not to query the db for every bioactivity
    assay_classification = {a.unit:a for a in AssayClassification.objects.all()}

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
        parser.add_argument('--make_dump',
                            action='store_true',
                            dest='make_dump',
                            default=False,
                            help='Create dump')
        parser.add_argument('--only_make_dump',
                            action='store_true',
                            dest='only_make_dump',
                            default=False,
                            help='Create dump, then quit')
        parser.add_argument('--no_reload',
                            action='store_true',
                            dest='no_reload',
                            default=False,
                            help='Skip ligand dump reload scripts')

    def handle(self, *args, **options):
        if options["test_run"]:
            print("Skipping in test run")
            return

        if options["only_make_dump"]:
            return

        if options['purge']:
            # delete any existing ChEMBL bioactivity data
            # For deleting the ligands - purse all ligand data using the GtP ligand build
            print("Started purging bioactivity and ligand data")
            self.purge_data()
            self.tracker = {}
            test_model_updates(self.all_models, self.tracker, initialize=True)
            print("Ended purging data")

        if not options['no_reload']:
            print("\n\nRebuilding the Ligand Model based on latest dump")
            self.reload_dump()
            print("Ended reloading data from ligand dump")

        # Fetching all the Guide to Pharmacology data
        print("\n\nStarted parsing Guide to Pharmacology bioactivities data")
        gtp = load_gtp_source_data()
        gtp_uniprot             = gtp["gtp_uniprot"]
        gtp_complete_ligands    = gtp["gtp_complete_ligands"]
        gtp_ligand_mapping      = gtp["gtp_ligand_mapping"]
        gtp_interactions        = gtp["gtp_interactions"]
        gtp_detailed_endogenous = gtp["gtp_detailed_endogenous"]
        gtp_peptides            = gtp["gtp_peptides"]
        iuphar_ids              = gtp["iuphar_ids"]
        ligand_ids              = gtp["ligand_ids"]
        bioactivity_ligands_ids = [lid for lid in ligand_ids
                                   if lid in set(gtp_interactions["ligand_id"].unique())]

        # This gets all the info of the ligand and the interaction with the target
        # Now I have all the data I need
        bioactivity_data_gtp = self.get_ligands_data(
            bioactivity_ligands_ids, gtp_complete_ligands, gtp_ligand_mapping, ligand_interactions=gtp_interactions, target_ids=iuphar_ids)
        # Assess the assay type given info from affinity units and assay comments
        bioactivity_data_gtp = self.classify_assay(bioactivity_data_gtp, 'affinity_units', 'assay_description')
        bioactivity_data_gtp.fillna('None', inplace=True)
        print("Ended parsing Guide to Pharmacology bioactivities data")

        print("\n\nStarted building all Guide to Pharmacology ligands")
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

        if not options['no_reload']:
            #HERE WE HAVE THE ACTUAL DATA TO COMPARE, WE NEED TO COMPARE AND THEN CREATE MISSING LIGANDS
            #GTP
            print('\n\nStarted comparing GTP data to reloaded database')
            ligand_data_merged, gtp_peptides_merged = self.find_unmatched_gtp(Command.ligand_csv, ligand_data_merged, gtp_peptides_merged)
            print(f"Building {len(ligand_data_merged)} small molecules and {len(gtp_peptides_merged)} peptide from GTP that were missing from dump")

        with open('./gtp_sm.txt','w') as f:
            for l in ligand_data_merged:
                f.write(l)
        with open('./gtp_peptides.txt','w') as f2:
            for k in gtp_peptides_merged:
                f2.write(k)

        self.build_gtoplig_data(ligand_data_merged, gtp_peptides_merged)
        print('Performing checks')
        test_model_updates(self.all_models, self.tracker, check=True, rebuild=True)

        #CHEMBL
        print("\n\nStarted comparing ChEBML ligands")
        self.build_chembl_ligands()
        print("\n\nEnded building ChEMBL ligands")
        print('Performing checks')
        test_model_updates(self.all_models, self.tracker, check=True, rebuild=True)

        #DrugBank
        print('\n\nFetching Drug Bank ligands and saving to model')
        self.build_drugbank_ligands()
        test_model_updates(self.all_models, self.tracker, check=True, rebuild=True)

        #BUILDING BIOACTIVITIES
        #GTP bioactivity data
        print("\n\nStarted building Guide to Pharmacology bioactivities")
        # bioactivities_to_update = self.find_unmatched_bioactivities(Command.ligand_csv, bioactivity_data_gtp)
        self.build_gtp_bioactivities(bioactivity_data_gtp)
        print("Ended building Guide to Pharmacology bioactivities")
        print('Performing checks')
        test_model_updates(self.all_models, self.tracker, check=True, rebuild=True)

        #ChEMBL bioactivity data
        print("\n\nStarted building ChEMBL bioactivities")
        self.build_chembl_bioactivities()
        print("Ended building ChEMBL bioactivities")
        print('Performing checks')
        test_model_updates(self.all_models, self.tracker, check=True, rebuild=True)

        #ChEMBL/PubChem vendor data
        print("\n\nStarted building PubChem vendor data")
        self.build_pubchem_vendor_links()
        print("Ended building PubChem vendor data")
        print('Performing checks')
        test_model_updates(self.all_models, self.tracker, check=True, rebuild=True)

        #PDSP KiDatabase bioactivity data
        print("\n\nStarted building PDSP KiDatabase bioactivities")
        # to_update = self.comparePDSP(Command.ligand_csv)
        # print(f"Building {len(to_update)} PDSP KiDatabase bioactivities missing from dump")
        self.build_kidatabase_bioactivities()  # 14,562
        print("Ended building PDSP KiDatabase bioactivities")
        print('Performing checks')
        test_model_updates(self.all_models, self.tracker, check=True, rebuild=True)

        #Drug Central bioactivity data
        print("\n\nStarted building Drug Central bioactivities")
        # to_update = self.compareDrugCentral(Command.ligand_csv)
        # print(f"Building {len(to_update)}  Drug Central bioactivities missing from dump")
        self.build_drugcentral_bioactivities()  # 5,844
        print("Ended building Drug Central bioactivities")
        print('Performing checks')
        test_model_updates(self.all_models, self.tracker, check=True, rebuild=True)

        # #Evolvus bioactivity data
        # print("\n\nStarted building Evolvus bioactivities")
        # self.build_evolvus_bioactivities()
        # print("\n\nEnded building Evolvus bioactivities")
        # print('Performing checks')
        # test_model_updates(self.all_models, self.tracker, check=True, rebuild=True)

        #ENDOGENOUS LIGANDS
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
        test_model_updates(self.all_models, self.tracker, check=True, rebuild=True)

        #AFTERMATH FIXES
        print("\n\nStarted calculating potency and affinity indexes")
        self.calculate_potency_and_affinity()
        print("Potency and affinity indexes have been added to the model")

        print("\n\nFixing mismatched LigandType definition")
        n  = apply_canonical_ligand_types()
        print("\n\nUpdated LigandType on {} records to their canonical type".format(n))

        print("\n\nClassifying lipid ligands")
        n = classify_lipids()
        print("\n\nReclassified {} records as lipid".format(n))

        if options["make_dump"]:
            print("\n\nRunning the Dump Checker and saving new csvs")
            ligand_dump = dump_checker('ligand.Ligand')
            id_dump = dump_checker('ligand.LigandID')

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
    def sync_pk_sequence_to_max(model):
        """
        Set the model's PostgreSQL sequence to MAX(id), so the next INSERT gets MAX(id)+1.
        Works even if rows were inserted with explicit IDs.
        """
        table = model._meta.db_table
        with connection.cursor() as cursor:
            cursor.execute(f'''
                SELECT setval(
                    pg_get_serial_sequence('"{table}"','id'),
                    COALESCE((SELECT MAX(id) FROM "{table}"), 0),
                    true
                );
            ''')

    @staticmethod
    def sync_sequence_to_max(model, column):
        """
        Sync the PostgreSQL sequence for `column` on `model` to MAX(column).
        Use after imports that set the column explicitly.
        Works for serial or identity columns.
        """
        table = model._meta.db_table
        with connection.cursor() as cursor:
            cursor.execute(f'''
                SELECT setval(
                    pg_get_serial_sequence('"{table}"','{column}'),
                    COALESCE((SELECT MAX("{column}") FROM "{table}"), 0),
                    true
                );
            ''')

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
        # Command.reset_pk_sequence(Ligand)
        Command.reset_pk_sequence(LigandID)

    @staticmethod
    def try_std(smiles):
        """
        Safely standardize a SMILES string:
          - returns None for missing/invalid SMILES
          - handles pipelines that return either Mol or (Mol, status)
          - strips stereochemistry and returns non-isomeric SMILES
        """
        if not isinstance(smiles, str) or not smiles.strip():
            return None
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            out = standardizer.standardize_mol(mol)
            std_mol = out[0] if isinstance(out, tuple) else out
            Chem.RemoveStereochemistry(std_mol)
            return Chem.MolToSmiles(std_mol, isomericSmiles=False)
        except Exception:
            return None

    @staticmethod
    def comparePDSP(dump):
        pdsp_link = get_or_create_url_cache(
            "https://pdsp.unc.edu/databases/kiDownload/download.php", 7 * 24 * 3600)
        bioactivity_kidata = pd.read_csv(pdsp_link, dtype=str, encoding='mac_roman')
        # Keeping data that has either SMILES info OR CAS info
        # CAS number can be translated into pubchem CID
        pdsp = bioactivity_kidata.loc[(
            ~bioactivity_kidata['SMILES'].isnull()) | (~bioactivity_kidata['CAS'].isnull())]
        pdsp = pdsp.loc[(
            ~pdsp['Unigene'].isnull())]
        pdsp.fillna('None', inplace=True)
        # 1) Normalize column names
        up = pdsp.rename(columns={
            ' Ligand Name': 'name',
            'SMILES':        'smiles',
        })

        # 2) Deduplicate the key columns
        up_sub = up[['name','smiles']].drop_duplicates().copy()
        lc_sub = dump[['name','smiles']].drop_duplicates().copy()

        up_sub['smiles'] = up_sub['smiles'].apply(Command.try_std)
        # 3) Merge on both name & smiles → find those missing entirely
        merged = up_sub.merge(
            lc_sub,
            on=['name','smiles'],
            how='left',
            indicator=True
        )
        missing_both = merged[merged['_merge']=='left_only'][['name','smiles']]

        # 5) Try matching again on (name, std_smiles)
        std_merged = missing_both.merge(
            lc_sub,
            on=['name','smiles'],
            how='left',
            indicator=True
        )
        std_matched   = std_merged[std_merged['_merge']=='both'][['name','smiles']]
        still_missing = std_merged[std_merged['_merge']=='left_only'][['name','smiles']]

        #  → truly no match
        true_no_matches = still_missing[
            ~still_missing['name'].isin(lc_sub['name']) &
            ~still_missing['smiles'].isin(lc_sub['smiles'])
        ]

        # get all names from up that are NOT in true_no_matches
        reverse_no_matches = up.loc[~up['name'].isin(true_no_matches['name']), 'name'].unique()

        reverse_no_matches = list(reverse_no_matches)

        return true_no_matches

    @staticmethod
    def compareDrugCentral(dump):
        """
        Compare ligand_csv vs merged_df on DRUG_NAME⇄name, SMILES⇄smiles, InChIKey⇄inchikey
        Returns a dict of DataFrames: std_matched, still_missing, name_only, smiles_only,
        inchikey_only, none.
        """
        # 1) Get DrugCentral data
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

        # 2) Normalize column names
        md = merged_data_filtered.rename(columns={
            'DRUG_NAME':'name',
            'SMILES':   'smiles',
            'InChIKey': 'inchikey',
        })

        # 3) Pull just the three key cols and dedupe
        lc_sub = dump[['name','smiles','inchikey']].drop_duplicates().copy()
        md_sub = md[['name','smiles','inchikey']].drop_duplicates().copy()

        # 4) Standardize SMILES
        md_sub['smiles'] = md_sub['smiles'].apply(Command.try_std)

        # 5) Exact-match on all three → find missing in ligand_csv
        merged_exact = md_sub.merge(
            lc_sub,
            on=['name','smiles','inchikey'],
            how='left',
            indicator=True
        )
        missing_exact = merged_exact[merged_exact['_merge']=='left_only'][['name','smiles','inchikey']]

        merged_std = missing_exact.merge(
            lc_sub[['name','smiles','inchikey']],
            left_on=['name','smiles','inchikey'],
            right_on=['name','smiles','inchikey'],
            how='left',
            indicator=True
        )
        still_missing = merged_std[merged_std['_merge']=='left_only'][['name','smiles','inchikey']]

        # 6) And the ones matching none of the three
        none = still_missing[
            ~still_missing['name'].isin(lc_sub['name']) &
            ~still_missing['smiles'].isin(lc_sub['smiles']) &
            ~still_missing['inchikey'].isin(lc_sub['inchikey'])
        ]

        to_update = pd.merge(
              none,
              md,
              on='inchikey',
              how='inner',       # keeps only rows where C exists in both
              suffixes=('_none','_og')   # if A/B differ, you’ll need suffixes; here we assume they’re identical
        )

        # 2) Select the “preferred” columns
        keep = []
        for col in to_update.columns:
            # if it’s a “_none” column but we have a corresponding “_og”, skip it
            if col.endswith('_none'):
                base = col[:-5]
                if f'{base}_og' in to_update.columns:
                    continue
            keep.append(col)

        # Sub‑select
        unmatched = to_update[keep].copy()

        # 3) Strip any of the suffixes off the names
        unmatched.rename(
            columns=lambda c: re.sub(r'_(none|og)$','',c),
            inplace=True
        )

        unmatched = unmatched.loc[:, ~unmatched.columns.duplicated()]

        return unmatched

    @staticmethod
    def find_unmatched_gtp(reference_df, small_mols_df, peptide_df):

        # 1) normalize columns
        small = small_mols_df.rename(columns={
            'ligand_name':'name'
        })[['name','smiles','inchikey']].drop_duplicates().copy()

        small['smiles'] = small['smiles'].apply(Command.try_std)

        peptide = peptide_df.rename(columns={
            'single_letter_amino_acid_sequence':'sequence'
        })[['name','smiles','inchikey','helm','sequence']].drop_duplicates().copy()

        peptide['smiles'] = peptide['smiles'].apply(Command.try_std)

        # 2) grab just name, smiles, inchikey, sequence and helm from reference
        ref = reference_df[['name','smiles','inchikey','sequence','helm']].drop_duplicates().copy()

        # 3) exact‐match on all three → which small are missing in reference?
        small_cmp = small.merge(
            ref, on=['name','smiles','inchikey'],
            how='left', indicator=True
        )
        small_missing = small_cmp[small_cmp['_merge']=='left_only'][['name','smiles','inchikey']].copy()

        # 4) those matching none of the three
        small_none = small_missing[
            ~small_missing['name'].isin(ref['name']) &
            ~small_missing['smiles'].isin(ref['smiles']) &
            ~small_missing['inchikey'].isin(ref['inchikey'] )
        ]

        small_to_update = pd.merge(
              small_none,
              small_mols_df,
              on='inchikey',
              how='inner',       # keeps only rows where C exists in both
              suffixes=('_none','_og')   # if A/B differ, you’ll need suffixes; here we assume they’re identical
        )

        # 2) Select the “preferred” columns
        keep = []
        for col in small_to_update.columns:
            # if it’s a “_none” column but we have a corresponding “_og”, skip it
            if col.endswith('_none'):
                base = col[:-5]
                if f'{base}_og' in small_to_update.columns:
                    continue
            keep.append(col)

        # Sub‑select
        unmatched_small = small_to_update[keep].copy()

        # 3) Strip any of the suffixes off the names
        unmatched_small.rename(
            columns=lambda c: re.sub(r'_(none|og)$','',c),
            inplace=True
        )

        unmatched_small = unmatched_small.loc[:, ~unmatched_small.columns.duplicated()]

        #6) Let's redo everything for peptides

        peptide_cmp = peptide.merge(
            ref, on=['name','smiles','inchikey','sequence','helm'],
            how='left', indicator=True
        )
        peptide_missing = peptide_cmp[peptide_cmp['_merge']=='left_only'][['name','smiles','inchikey','sequence','helm']].copy()

        # 5) those matching none of the three
        peptide_none = peptide_missing[
            ~peptide_missing['name'].isin(ref['name']) &
            ~peptide_missing['smiles'].isin(ref['smiles']) &
            ~peptide_missing['sequence'].isin(ref['sequence']) &
            ~peptide_missing['helm'].isin(ref['helm']) &
            ~peptide_missing['inchikey'].isin(ref['inchikey'] )
        ]

        peptides_to_update = pd.merge(
              peptide_none,
              peptide_df,
              on='inchikey',
              how='inner',       # keeps only rows where C exists in both
              suffixes=('_none','_og')   # if A/B differ, you’ll need suffixes; here we assume they’re identical
        )

        # 2) Select the “preferred” columns
        keep = []
        for col in peptides_to_update.columns:
            # if it’s a “_none” column but we have a corresponding “_og”, skip it
            if col.endswith('_none'):
                base = col[:-5]
                if f'{base}_og' in peptides_to_update.columns:
                    continue
            keep.append(col)

        # Sub‑select
        unmatched_peptide = peptides_to_update[keep].copy()

        # 3) Strip any of the suffixes off the names
        unmatched_peptide.rename(
            columns=lambda c: re.sub(r'_(none|og)$','',c),
            inplace=True
        )

        unmatched_peptide = unmatched_peptide.loc[:, ~unmatched_peptide.columns.duplicated()]

        return unmatched_small, unmatched_peptide

    @staticmethod
    def find_unmatched_bioactivities(reference_df, bio_df):

        # 1) normalize columns
        bio = bio_df[['name','smiles','inchikey']].drop_duplicates().copy()

        bio['smiles'] = small['smiles'].apply(Command.try_std)

        # 2) grab just name, smiles, inchikey, sequence and helm from reference
        ref = reference_df[['name','smiles','inchikey']].drop_duplicates().copy()

        # 3) exact‐match on all three → which small are missing in reference?
        bio_cmp = bio.merge(
            ref, on=['name','smiles','inchikey'],
            how='left', indicator=True
        )
        bio_missing = bio_cmp[bio_cmp['_merge']=='left_only'][['name','smiles','inchikey']].copy()

        # 4) those matching none of the three
        bio_none = bio_missing[
            ~bio_missing['name'].isin(ref['name']) &
            ~bio_missing['smiles'].isin(ref['smiles']) &
            ~bio_missing['inchikey'].isin(ref['inchikey'] )
        ]

        bio_to_update = pd.merge(
              bio_none,
              bio_df,
              on='name',
              how='inner',       # keeps only rows where C exists in both
              suffixes=('_none','_og')   # if A/B differ, you’ll need suffixes; here we assume they’re identical
        )

        # 2) Select the “preferred” columns
        keep = []
        for col in bio_to_update.columns:
            # if it’s a “_none” column but we have a corresponding “_og”, skip it
            if col.endswith('_none'):
                base = col[:-5]
                if f'{base}_og' in bio_to_update.columns:
                    continue
            keep.append(col)

        # Sub‑select
        unmatched_bio = bio_to_update[keep].copy()

        # 3) Strip any of the suffixes off the names
        unmatched_bio.rename(
            columns=lambda c: re.sub(r'_(none|og)$','',c),
            inplace=True
        )

        unmatched_bio = unmatched_bio.loc[:, ~unmatched_bio.columns.duplicated()]

        return unmatched_bio

    @staticmethod
    def reload_dump():
        # --- Helper to print progress every 10% ---
        def _progress_printer(total, prefix):
            """Returns a closure that you can call with current index to print at 10% intervals."""
            percent = 10
            def tick(idx):
                nonlocal percent
                # when idx crosses the next threshold, print and bump
                if idx >= total * percent / 100:
                    print(f"{prefix}: {percent}%")
                    percent += 10
            return tick

        # --- Setting up the opened for gzip or not ---
        if Command.ligand_dump['latest_dump'].endswith(".gz"):
            opener_ligand = lambda path, **kw: gzip.open(path, mode="rt", **kw)
        else:
            opener_ligand = lambda path, **kw: open(path, mode="r", **kw)
        if Command.id_dump['latest_dump'].endswith(".gz"):
            opener_ligand_id = lambda path, **kw: gzip.open(path, mode="rt", **kw)
        else:
            opener_ligand_id = lambda path, **kw: open(path, mode="r", **kw)

        # opener_ligand = gzip.open if Command.ligand_dump['latest_dump'].endswith(".gz") else open
        # opener_ligand_id = gzip.open if Command.id_dump['latest_dump'].endswith(".gz") else open

        compounds = {}
        # --- PASS 1: create all compounds without parent ---
        # 1a) count rows
        with opener_ligand(Command.ligand_dump['latest_dump'], newline='', encoding='utf-8-sig') as f:
            total = sum(1 for _ in f) - 1

        print("Pass 1 (compounds): 0%")
        tick1 = _progress_printer(total, "Pass 1 (compounds)")

        # 1b) actual work
        with opener_ligand(Command.ligand_dump['latest_dump'], newline='', encoding='utf-8-sig') as csvfile:
            reader = csv.DictReader(csvfile, delimiter=';')
            ligand_objects = []
            ligand_types = {'none': LigandType.objects.get(slug='none')}
            for idx, row in enumerate(reader, start=1):
                # print progress if needed
                tick1(idx)
                if row.get('gpcrdb_id')=='':
                    continue
                if row['ligand_type_id__slug'] in ligand_types:
                    lig_type = ligand_types[row['ligand_type_id__slug']]
                else:
                    try:
                        lig_type = LigandType.objects.get(slug=row['ligand_type_id__slug'])
                        ligand_types[lig_type.slug] = lig_type
                    except LigandType.DoesNotExist:
                        lig_type = ligand_types['none']

                compound = Ligand(
                    # id=int(row['id']),
                    name=row.get('name'),
                    pdbe=row.get('pdbe') or None,
                    ambiguous_alias=row.get('ambiguous_alias') or None,
                    clean_inchikey=row.get('clean_inchikey') or None,
                    hacc=row.get('hacc') or None,
                    hdon=row.get('hdon') or None,
                    inchikey=row.get('inchikey') or None,
                    ligand_type=lig_type,
                    logp=row.get('logp') or None,
                    mw=row.get('mw') or None,
                    rotatable_bonds=row.get('rotatable_bonds') or None,
                    sequence=row.get('sequence') or None,
                    smiles=row.get('smiles') or None,
                    uniprot=row.get('uniprot') or None,
                    source=row.get('source') or None,
                    helm=row.get('helm') or None,
                    gpcrdb_id=int(row.get('gpcrdb_id')),
                    radioactive=row.get('radioactive') or None,
                    stereo_status=row.get('stereo_status') or None,
                    parent=None
                )
                compounds[compound.gpcrdb_id] = compound
                ligand_objects.append(compound)
            Ligand.objects.bulk_create(ligand_objects)
        print("Pass 1 (compounds): 100%")

        # --- PASS 2: set parent relationships ---
        with opener_ligand(Command.ligand_dump['latest_dump'], newline='', encoding='utf-8-sig') as f:
            total = sum(1 for _ in f) - 1

        print("Pass 2 (parents): 0%")
        tick2 = _progress_printer(total, "Pass 2 (parents)")

        with opener_ligand(Command.ligand_dump['latest_dump'], newline='', encoding='utf-8-sig') as csvfile:
            reader = csv.DictReader(csvfile, delimiter=';')
            for idx, row in enumerate(reader, start=1):
                tick2(idx)
                if row.get('gpcrdb_id')=='':
                    continue
                parent_id = row.get('parent_id__gpcrdb_id')
                if parent_id:
                    try:
                        compound = compounds[int(row.get('gpcrdb_id'))]
                        compound.parent = compounds.get(int(parent_id))
                        compound.save()
                    except KeyError:
                        print(f"Parent ID {parent_id} not found for compound ID {row['id']}")
        print("Pass 2 (parents): 100%")

        # --- PASS 3: create all the LigandIDs ---
        with opener_ligand_id(Command.id_dump['latest_dump'], newline='', encoding='utf-8-sig') as f:
            total = sum(1 for _ in f) - 1

        print("Pass 3 (IDs): 0%")
        tick3 = _progress_printer(total, "Pass 3 (IDs)")

        with opener_ligand_id(Command.id_dump['latest_dump'], newline='', encoding='utf-8-sig') as idsfile:
            reader = csv.DictReader(idsfile, delimiter=';')
            web_resources = {}
            ligandid_objects = []
            for idx, row in enumerate(reader, start=1):
                tick3(idx)
                if row['web_resource_id__slug'] in web_resources:
                    wr = web_resources[row['web_resource_id__slug']]
                else:
                    wr = WebResource.objects.get(slug=row['web_resource_id__slug'])
                    web_resources[row['web_resource_id__slug']] = wr
                try:
                    record = LigandID(
                        id=int(row['id']),
                        index=row['index'],
                        ligand=Ligand.objects.get(gpcrdb_id=int(row['ligand_id__gpcrdb_id'])),
                        web_resource=wr
                    )
                    ligandid_objects.append(record)
                except Exception as e:
                    print(f"Impossible to import LigandID {row['index']}: {e!r}")
            LigandID.objects.bulk_create(ligandid_objects)
        print("Pass 3 (IDs): 100%")

        print("Syncing the MAX id values")
        Command.sync_sequence_to_max(Ligand, 'gpcrdb_id')
        # Command.sync_pk_sequence_to_max(Ligand)
        Command.sync_pk_sequence_to_max(LigandID)

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
            ligand = get_ligand_by_id("gtoplig", row['ligand_id'], forced=False)
            if ligand is not None:
                print(ligand)
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
        ligand_data.replace(["", "None", "null", "NaN", "nan"], np.nan, inplace=True)
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

        sm_data = filtered_ligands[
            (filtered_ligands["molecule_type"].isin(["Small molecule", "Oligosaccharide"])) |
            (pd.isna(filtered_ligands["molecule_type"]))
        ].reset_index()

        #check that we substitute all the 'nan' with actual np.nan values
        sm_data.replace('nan', np.nan, inplace=True)


        lig_entries = len(sm_data)

        print(f"\n#3 Building {lig_entries} new small-molecule ChEMBL ligands", datetime.datetime.now())

        entries_batch = []
        for index, row in sm_data.iterrows():
            _name = str(row.get('pref_name') or row.get('molecule_chembl_id') or "")
            _smiles = row.get('smiles') if pd.notna(row.get('smiles')) else None
            _ik = row.get('standard_inchi_key') if pd.notna(row.get('standard_inchi_key')) else None
            _chembl_id = row.get('molecule_chembl_id')
            _pubchem_cid = None
            if pd.notna(row.get('pubchem_cid')) and row.get('pubchem_cid'):
                try:
                    _pubchem_cid = str(int(float(str(row['pubchem_cid']))))
                except (ValueError, TypeError):
                    _pubchem_cid = str(row['pubchem_cid']).strip() or None

            source_ids = {}
            if _chembl_id:
                source_ids['chembl_ligand'] = str(_chembl_id)
            if _pubchem_cid:
                source_ids['pubchem'] = _pubchem_cid

            entries_batch.append({
                "name": _name,
                "raw_smiles": _smiles,
                "raw_inchikey": _ik,
                "raw_helm": None,
                "raw_sequence": None,
                "raw_uniprot": None,
                "source_ids": source_ids,
                "radioactive": False,
                "synonyms": False,
            })

            if len(entries_batch) >= Command.bulk_size:
                get_or_create_ligands_batch(
                    entries_batch, source="ChEMBL_sm", lig_type="small-molecule"
                )
                entries_batch = []

        if entries_batch:
            get_or_create_ligands_batch(
                entries_batch, source="ChEMBL_sm", lig_type="small-molecule"
            )

        # Parse all new non-small-molecule ChEMBL ligands
        print("\n#4 Building new non-small-molecule ChEMBL ligands", datetime.datetime.now())

        nonsm_data = filtered_ligands[
                        ~((filtered_ligands["molecule_type"].isin(["Small molecule", "Oligosaccharide"])) |
                        (pd.isna(filtered_ligands["molecule_type"])))
                     ].reset_index()
        print("Found", len(nonsm_data), "new non-small-molecules")

        ligand_types = {"Unknown": "na", "Protein": "protein"}

        alias_rows = []  # (_lig_type, row, _smiles, _ik, _sequence, _helm) for Phase 2

        for _lig_type in ("na", "protein"):
            typed_rows = [row for _, row in nonsm_data.iterrows()
                          if ligand_types.get(row.get('molecule_type'), 'na') == _lig_type]
            if not typed_rows:
                continue
            batch = []
            for row in typed_rows:
                _name = str(row.get('pref_name') or row.get('molecule_chembl_id') or "")
                _smiles = row.get('smiles') if pd.notna(row.get('smiles')) else None
                _ik = row.get('standard_inchi_key') if pd.notna(row.get('standard_inchi_key')) else None
                _sequence = row.get('sequence') if pd.notna(row.get('sequence')) else None
                _helm = row.get('helm_notation') if pd.notna(row.get('helm_notation')) else None
                _chembl_id = row.get('molecule_chembl_id')
                batch.append({
                    "name": _name,
                    "raw_smiles": _smiles,
                    "raw_inchikey": _ik,
                    "raw_helm": _helm,
                    "raw_sequence": _sequence,
                    "raw_uniprot": None,
                    "source_ids": {"chembl_ligand": str(_chembl_id)} if _chembl_id else {},
                    "radioactive": False,
                    "synonyms": False,
                    "skip_normalization": True,
                })
                if pd.notna(row.get("other_ids")) and row.get("other_ids"):
                    alias_rows.append((_lig_type, row, _smiles, _ik, _sequence, _helm))
            get_or_create_ligands_batch(batch, source="ChEMBL_peptide", lig_type=_lig_type)

        # Phase 2: alias IDs — children already in DB after Phase 1, Step E finds by
        # InChIKey and queues a new LigandID without creating a duplicate child
        alias_batch = []
        for _lig_type, row, _smiles, _ik, _sequence, _helm in alias_rows:
            for alias_id in str(row["other_ids"]).split(";"):
                alias_id = alias_id.strip()
                if not alias_id:
                    continue
                alias_batch.append({
                    "name": str(row.get('pref_name') or alias_id),
                    "raw_smiles": _smiles,
                    "raw_inchikey": _ik,
                    "raw_helm": _helm,
                    "raw_sequence": _sequence,
                    "raw_uniprot": None,
                    "source_ids": {"chembl_ligand": alias_id},
                    "radioactive": False,
                    "synonyms": False,
                    "skip_normalization": True,
                })
        if alias_batch:
            get_or_create_ligands_batch(alias_batch, source="ChEMBL_peptide", lig_type="na")

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
        # ligands = list(Ligand.objects.filter(name__startswith="CHEMBL", parent__isnull=False).values_list("id", "name").distinct())
        # ligands = list(Ligand.objects.filter(name__startswith="CHEMBL").values_list("id", "name").distinct())
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
                    bioacts[-1].assay_classification = Command.fetch_assay_classification(row["standard_type"], row["assay_description"])
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
        missed = 0
        for _, row in vendor_links_data.iterrows():
            if len(row["SourceRecordURL"])<=400 and len(row["RegistryID"])<=500:
                try:
                    links.append(LigandVendorLink(
                        vendor_id=vendor_dict[row["SourceName"]], ligand_id=lig_dict[row["chembl_id"]], url=row["SourceRecordURL"], external_id=row["RegistryID"]))
                except KeyError:
                    missed += 1
                    continue

        LigandVendorLink.objects.bulk_create(links)
        if missed:
            print(f"Skipped {missed} vendor links with an unknown vendor or ligand")

    @staticmethod
    def uniprot_mapper_update(protein, organism):
        organism_dict = {
            'PIG': 'sus_scrofa', 'RAT': 'rattus_norvegicus', 'HUMAN': 'homo_sapiens', 'MOUSE': 'mus_musculus',
            'CANINE': 'canis_lupus_familiaris', 'BOVINE': 'bos_taurus', 'CALF': 'bos_taurus', 'COW': 'bos_taurus',
            'GUINEA PIG': 'cavia_porcellus', 'CAT': 'felis_catus', 'NEONATAL RAT': 'rattus_norvegicus',
            '? HUMAN': 'homo_sapiens', 'OPOSSUM': 'didelphis_marsupialis', 'Rat 6B': 'rattus_norvegicus',
            'HUMAN M3': 'homo_sapiens', 'HUMAN M4': 'homo_sapiens', 'Chick': 'gallus_gallus',
            'Frog': 'pseudis_balbodactyla', 'Newborn rats': 'rattus_norvegicus', 'Beef': 'bos_taurus',
            'Sheep': 'ovis_aries', 'OX': 'bos_taurus', 'Dog': 'canis_lupus_familiaris',
            'Rhesus': 'macaca_mulatta', 'Monkey': 'macaca_mulatta', 'PIGLET': 'sus_scrofa',
            'Rat Y861': 'rattus_norvegicus', 'Zebra Finch': 'taeniopygia_guttata', 'Chicken': 'gallus_gallus',
            'MICE': 'mus_musculus', 'Rhesus Monkey': 'macaca_mulatta', 'Zebrafish': 'danio_rerio'
        }

        # Build the UniProt search query
        gene_term = f"gene_exact:{protein.lower()}"
        if organism in organism_dict:
            org_term = f'organism_name:"{organism_dict[organism]}"'
            query = f"{gene_term} AND {org_term}"
        else:
            query = gene_term

        if query in Command.mapper_cache:
            return Command.mapper_cache[query]

        params = {
            "query": query,
            "fields": "id",
            "format": "tsv",
            "size": 1  # only first hit
        }
        url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(params)

        headers = {"User-Agent": "Mozilla/5.0 (compatible; pdsp-mapper/1.0)"}

        # Retry with exponential backoff for transient network issues
        retries = 5
        backoff = 0.5
        last_err = None
        text = None

        for attempt in range(1, retries + 1):
            try:
                req = urllib.request.Request(url, headers=headers)
                with urllib.request.urlopen(req, timeout=20) as resp:
                    text = resp.read().decode("utf-8", errors="replace")
                break
            except (urllib.error.URLError, socket.timeout, ConnectionResetError, ssl.SSLError) as e:
                last_err = e
                if attempt < retries:
                    time.sleep(backoff * (2 ** (attempt - 1)))
                else:
                    # do not cache a transient failure
                    print(f"UniProt request failed for query={query!r}: {e}")
                    return None

        # Robust TSV parsing
        reader = csv.DictReader(StringIO(text), delimiter="\t")
        first_row = next(reader, None)
        converted = (first_row["Entry"] if first_row and "Entry" in first_row else None)

        # Optional: keep original case (recommended)
        # If you *must* lower for keys, do it at the call site, not the accession.
        Command.mapper_cache[query] = converted
        return converted

    # @staticmethod
    # def uniprot_mapper(protein, organism):
    #     organism_dict = {'PIG': 'sus_scrofa', 'RAT': 'rattus_norvegicus', 'HUMAN': 'homo_sapiens', 'MOUSE': 'mus_musculus',
    #                      'CANINE': 'canis_lupus_familiaris', 'BOVINE': 'bos_taurus', 'CALF': 'bos_taurus', 'COW': 'bos_taurus',
    #                      'GUINEA PIG': 'cavia_porcellus', 'CAT': 'felis_catus', 'NEONATAL RAT': 'rattus_norvegicus',
    #                      '? HUMAN': 'homo_sapiens', 'OPOSSUM': 'didelphis_marsupialis', 'Rat 6B': 'rattus_norvegicus',
    #                      'HUMAN M3': 'homo_sapiens', 'HUMAN M4': 'homo_sapiens', 'Chick': 'gallus_gallus',
    #                      'Frog': 'pseudis_balbodactyla', 'Newborn rats': 'rattus_norvegicus', 'Beef': 'bos_taurus',
    #                      'Sheep': 'ovis_aries', 'OX': 'bos_taurus', 'Dog': 'canis_lupus_familiaris',
    #                      'Rhesus': 'macaca_mulatta', 'Monkey': 'macaca_mulatta', 'PIGLET': 'sus_scrofa',
    #                      'Rat Y861': 'rattus_norvegicus', 'Zebra Finch': 'taeniopygia_guttata', 'Chicken': 'gallus_gallus',
    #                      'MICE': 'mus_musculus', 'Rhesus Monkey': 'macaca_mulatta', 'Zebrafish': 'danio_rerio'}
    #     if organism in organism_dict.keys():
    #         query = 'gene_exact:{0}+AND+organism_name:{1}'.format(
    #             urllib.parse.quote(protein.lower()), organism_dict[organism])
    #     else:
    #         query = 'gene_exact:{}'.format(urllib.parse.quote(protein.lower()))
    #     if query not in Command.mapper_cache.keys():
    #         url = 'https://rest.uniprot.org/uniprotkb/search?query={}&fields=id&format=tsv'.format(query)
    #         req = urllib.request.Request(url)
    #         try:
    #             converted = urllib.request.urlopen(
    #                 req).read().decode('utf-8').split('\n')[1].lower()
    #         except IndexError:
    #             converted = None
    #
    #         Command.mapper_cache[query] = converted
    #
    #     return Command.mapper_cache[query]

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
            family__slug__startswith="0", sequence_type__slug="wt").values_list('entry_name', 'accession')
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
    def build_gtoplig_data(lig_df, pep_df):
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

            radioactive = False
            if row['radioactive']=='yes':
                radioactive = True

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
                            ligand = get_or_create_ligand(new_name, ids, ligand_type, True, False, 'GuideToPharma', pep_row['helm_notation'], radioactive, row['synonyms'])
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
                        ligand = get_or_create_ligand(new_name, ids, ligand_type, True, False, 'GuideToPharma', pep_row['helm_notation'], radioactive, row['synonyms'])
                        if ligand is None:
                            print("Issue with", row['name'])
                            print(row)
                else:
                    new_name = row['name']
                    ligand = get_or_create_ligand(new_name, ids, ligand_type, True, False, 'GuideToPharma', radioactive=radioactive, synonyms=row['synonyms'])
                    if ligand is None:
                        print("Issue with", row['name'])
                        print(row)
            else:
                ligand = get_or_create_ligand(row['name'], ids, ligand_type, True, False, 'GuideToPharma', row['helm_notation'], radioactive, row['synonyms'])
                if ligand is None:
                    print("Issue with", row['name'])
                    print(row)

    @staticmethod
    def fetch_assay_classification(unit, assay_description=""):
        # Handling edge cases
        if unit in ['AC50','ED50']:
            desc = assay_description.lower()
            effect_slug = "Unknown"
            # stimulatory if “agonist” but not part of “antagonist”
            if 'agonist' in desc and 'antagonist' not in desc:
                effect_slug = 'Stimulatory'
            # inhibitory if either “antagonist” or “inhibitor”
            elif 'antagonist' in desc or 'inhibitor' in desc:
                effect_slug = 'Inhibitory'
            return AssayClassification.objects.get(unit=unit, modality=effect_slug)
        if unit in Command.assay_classification:
            return Command.assay_classification[unit]
        else:
            if unit and unit!='-':
                print(f'WARNING: unit type {unit} not in assay classification table')
            return None

    @staticmethod
    def build_gtp_bioactivities(gtp_biodata):
        print("# Start parsing the GTP Dataframe")
        for _, row in gtp_biodata.iterrows():
            receptor = Command.fetch_protein(
                row['target_id'], 'GtoP', row['target_species'])
            ligand = get_ligand_by_id("gtoplig", row['ligand_id'], forced=False, name=row.get('ligand_name'))

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
                    assay_classification=Command.fetch_assay_classification(row['affinity_units'], row['assay_description']),
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
        src = 'Remote'
        conversion = {'Gene': {'Remote': 'Unigene',
                              'Local': 'UniGene Code'},
                      'Species': {'Remote': 'species',
                                  'Local': 'Species'},
                      'Ligand': {'Remote': ' Ligand Name',
                                 'Local': 'Ligand name'},
                      'Ki': {'Remote': 'ki Val',
                             'Local': 'Ki value'},
                      'Hot': {'Remote': 'Hotligand',
                              'Local': 'Hot Ligand'}
                      }
        pdsp_link = get_or_create_url_cache("https://pdspdb.unc.edu/databases/kiDownload/download.php", 7 * 24 * 3600)
        pdsp_file = os.sep.join([Command.data_dir, 'pdsp_ki_backup_august25.csv'])
        bioactivity_kidata = pd.read_csv(pdsp_link, dtype=str, encoding='mac_roman')
        if 'Could not connect: ' in bioactivity_kidata.columns:
            print("\n===============\Prep: Error is fetching data from remote, using local backup file")
            bioactivity_kidata = pd.read_csv(pdsp_file, dtype=str, encoding='mac_roman')
            bioactivity_kidata["Ki value"].replace(">10000", 10000.0, inplace=True)
            src = 'Local'
        # Keeping data that has either SMILES info OR CAS info
        # CAS number can be translated into pubchem CID
        pdsp = bioactivity_kidata.loc[(
            ~bioactivity_kidata['SMILES'].isnull()) | (~bioactivity_kidata['CAS'].isnull())]
        pdsp = pdsp.loc[(
            ~pdsp[conversion['Gene'][src]].isnull())]
        pdsp.fillna('None', inplace=True)
        bio_entries = len(pdsp)
        print("\n===============\n#1 Start parsing PDSP data")
        bioacts = []
        missed_no_name = 0
        for index, (_, row) in enumerate(pdsp.iterrows()):
            receptor = None
            label = '_'.join([row[conversion['Gene'][src]], row[conversion['Species'][src]]])
            if label not in protein_names.keys():
                protein = Command.uniprot_mapper_update(row[conversion['Gene'][src]], row[conversion['Species'][src]])
                if protein is not None:
                    protein_names[label] = Command.fetch_protein(protein, 'PDSP')

            ligand_name = row[conversion['Ligand'][src]]
            if not ligand_name or not ligand_name.strip() or ligand_name == 'None':
                missed_no_name += 1
                continue

            ids = {}
            if row['SMILES'] != 'None':
                ids['smiles'] = row['SMILES']
                input_mol = dm.to_mol(row['SMILES'], sanitize=True)
                if input_mol:
                # Check if InChIKey has been set
                    ids['inchikey'] = dm.to_inchikey(input_mol)
            if row['CAS'] != 'None':
                ids['CAS'] = row['CAS']
            # Keyed on name + structural identifiers (not name alone) so two different
            # compounds that happen to share a display name aren't merged into one ligand.
            cache_key = (ligand_name, row['SMILES'], row['CAS'])
            if cache_key not in ligand_cache.keys():
                ligand = get_or_create_ligand(ligand_name, ids, source='PDSP KiDatabase')
                ligand_cache[cache_key] = ligand
            if label in protein_names.keys():
                receptor = protein_names[label]
            if (receptor is not None) and (ligand_cache[cache_key] is not None):
                bioacts.append(AssayExperiment())
                bioacts[-1].ligand_id = ligand_cache[cache_key].id
                bioacts[-1].protein_id = receptor.id
                bioacts[-1].assay_type = 'B'
                bioacts[-1].assay_description = None
                bioacts[-1].standard_activity_value = round(float(row[conversion['Ki'][src]]), 2)
                bioacts[-1].p_activity_value = round(-math.log10(float(row[conversion['Ki'][src]]) * 1e-9), 2)
                bioacts[-1].p_activity_ranges = None
                bioacts[-1].standard_relation = '='
                bioacts[-1].value_type = 'pKi'
                bioacts[-1].assay_classification = Command.fetch_assay_classification('pKi')
                bioacts[-1].source = 'PDSP KiDatabase'
                bioacts[-1].document_chembl_id = None
                bioacts[-1].reference_ligand = row[conversion['Hot'][src]]
                # BULK insert every X entries or last entry
            if (len(bioacts) == Command.bulk_size) or (index == bio_entries - 1):
                AssayExperiment.objects.bulk_create(bioacts)
                print("Inserted", index, "out of",
                      bio_entries, "bioactivities")
                bioacts = []

            Command.assign_ligand_target_pairing(ligand_cache[cache_key], receptor, None, 'pKi')
        if missed_no_name:
            print(f"Skipped {missed_no_name} PDSP KiDatabase rows with a missing ligand name")

    @staticmethod
    def build_drugcentral_bioactivities():
        print("# Collecting DrugCentral data")
        accession_numbers = {}
        ligand_cache = {}
        # 1) Get DrugCentral data
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
        missed_no_name = 0
        for index, (_, row) in enumerate(merged_data_filtered.iterrows()):
            receptor = None
            code = row['ACCESSION']
            if code not in accession_numbers.keys():
                protein = Command.fetch_protein(code, 'DrugCentral')
                if protein is not None:
                    accession_numbers[code] = protein

            if not row['DRUG_NAME'] or not row['DRUG_NAME'].strip() or row['DRUG_NAME'] == 'None':
                missed_no_name += 1
                continue

            ids = {'drug_central':row['ID']}
            if row['SMILES'] != 'None':
                ids['smiles'] = row['SMILES']
            # if row['CAS_RN'] != 'None':
            #     ids['CAS'] = row['CAS_RN']
            if row['InChIKey'] != 'None':
                ids['inchikey'] = row['InChIKey']
            # Keyed on the DrugCentral STRUCT_ID (row['ID']), not the display name, so two
            # different structures that happen to share a DRUG_NAME aren't merged into one
            # ligand and each structure's own drug_central id reaches get_or_create_ligand.
            cache_key = row['ID']
            if cache_key not in ligand_cache.keys():
                ligand = get_or_create_ligand(row['DRUG_NAME'], ids, source='Drug Central')
                ligand_cache[cache_key] = ligand
            if code in accession_numbers.keys():
                receptor = accession_numbers[code]

            # Case sensitive fix for pKB
            if row['ACT_TYPE']=='Kb':
                act_type = 'KB'
            elif row['ACT_TYPE']=='pKb':
                act_type = 'pKB'
            else:
                act_type = row['ACT_TYPE']
            if (receptor is not None) and (ligand_cache[cache_key] is not None):
                calc_val = round(-math.log10(float(row['ACT_VALUE']) * 1e-9), 2) if act_type != 'pA2' else round(float(row['ACT_VALUE']), 2)
                bioacts.append(AssayExperiment())
                bioacts[-1].ligand_id = ligand_cache[cache_key].id
                bioacts[-1].protein_id = receptor.id
                bioacts[-1].assay_type = row['assay_type']
                bioacts[-1].assay_description = row['ACT_COMMENT']
                bioacts[-1].standard_activity_value = round(float(row['ACT_VALUE']), 2) if act_type != 'pA2' else None
                bioacts[-1].p_activity_value = calc_val
                bioacts[-1].p_activity_ranges = None
                bioacts[-1].standard_relation = row['RELATION']
                bioacts[-1].value_type = 'p'+act_type if act_type != 'pA2' else act_type
                bioacts[-1].assay_classification = Command.fetch_assay_classification(bioacts[-1].value_type, row['ACT_COMMENT'])
                bioacts[-1].source = 'Drug Central'
                bioacts[-1].document_chembl_id = None
                # BULK insert every X entries or last entry
            if (len(bioacts) == Command.bulk_size) or (index == bio_entries - 1):
                AssayExperiment.objects.bulk_create(bioacts)
                print("Inserted", index, "out of",
                      bio_entries, "bioactivities")
                bioacts = []
            value_type = 'p'+act_type if act_type != 'pA2' else act_type

            Command.assign_ligand_target_pairing(ligand_cache[cache_key], receptor, row['ACT_COMMENT'], value_type)
        if missed_no_name:
            print(f"Skipped {missed_no_name} Drug Central rows with a missing ligand name")

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
            if row.SMILES not in ['None','nan']:
                ids['smiles'] = row.SMILES
            # if row['CAS Number'] != 'None':
            #     ids['CAS'] = row['CAS Number']
            if row.InChIKey not in ['None','nan']:
                ids['inchikey'] = row.InChIKey
            #PubChem ID
            if row._10 not in ['None','nan']:
                ids['pubchem'] = int(row._10)
            #DrugBank ID
            if row._0 not in ['None','nan']:
                ids['drugbank'] = row._0
            

            name = row.Name
            if name not in ligand_cache:
                ligand = get_or_create_ligand(name, ids, source='DrugBank')
                ligand_cache[name] = ligand

    # @staticmethod
    def build_evolvus_bioactivities(self):
        self.protein_names = {}
        self.ligand_cache = {}
        self.errors = {}
        new_columns = ['Reference',
                       'Fig/Table with data',
                       'Fig/Table no.',
                       'GPCR Name',
                       'UniProt',
                       'RAMP UniProt',
                       'Ligand Name',
                       'N-term mod',
                       'Sequence',
                       'C-term mod',
                       'Cyclization positions',
                       'Cyclization type',
                       'Label',
                       'Activity Type',
                       'Sign',
                       'Value',
                       'Unit',
                       'Emax (%)',
                       'Qualitative activity',
                       'Assay Type',
                       'Reference ligand name',
                       'Reference ligand PubChem CID',
                       'Curator',
                       'Remarks',
                       'Errors']

        evolvus_data = os.path.join(settings.DATA_DIR, "ligand_data", "evolvus_data")
        class_a_data = Command.read_data(evolvus_data, 'Annotation_AB1.xlsx', 'cl_A')
        class_b1_data = Command.read_data(evolvus_data, 'Annotation_AB1.xlsx', 'cl_B1')
        class_a_data = class_a_data.loc[:, ~class_a_data.columns.str.lower().str.startswith("unnamed")]
        class_b1_data = class_b1_data.loc[:, ~class_b1_data.columns.str.lower().str.startswith("unnamed")]
        class_a_data.columns = new_columns
        class_b1_data.columns = new_columns
        print("\n===============\n#1 Start parsing Evolvus data")
        self.evolvus_main(class_a_data)
        self.evolvus_main(class_b1_data)

        c = 0
        for p, ligs in self.errors.items():
            for l in ligs:
                c+=1
        print('Evolvus error count: ', c)

    # @staticmethod
    def evolvus_main(self, data):
        bioacts = []
        pub_links = []
        bio_entries = len(data)
        unit_to_factor = {
            'pM': 1e-12,
            'nM': 1e-9,
            'µM': 1e-6,
            'uM': 1e-6,
            'mM': 1e-3,
            'M':  1.0,
        }

        for index, (_, row) in enumerate(data.iterrows()):
            assay_type = 'U'
            receptor = None
            helm = None
            if pd.notna(row['Assay Type']):
                assay_type = row['Assay Type'][0]
            protein = str(row['UniProt']).lower() if pd.notna(row['UniProt']) else None

            if protein:
                if protein not in self.protein_names:
                    self.protein_names[protein] = Command.fetch_protein(protein, 'PDSP')  # may be None
                receptor = self.protein_names[protein]

            if pd.notna(row['Sequence']):
                helm = generate_helm_universal(sequence = row['Sequence'],
                                               n_term = row['N-term mod'],
                                               c_term = row['C-term mod'],
                                               cyclizations = parse_cycl_pos(row['Cyclization positions']),
                                               cyclization_type = row['Cyclization type'])

            ligand_label = f"{row['Ligand Name']}_{hash(helm)}"  # shorter/stable key

            # Handling sequence data - setting seq to None if special characters are present
            if pd.notna(row['Sequence']):
                seq = row['Sequence']
                char_set_not_wanted = set("[]-()=d")
                has_any = not char_set_not_wanted.isdisjoint(seq)
                if has_any:
                    seq = None
            else:
                seq = None

            if ligand_label not in self.ligand_cache.keys():
                ids = {}
                ids['sequence'] = seq
                ligand = get_or_create_ligand(row['Ligand Name'], ids, lig_type='peptide', source='Evolvus', helm=helm)
                self.ligand_cache[ligand_label] = ligand

            if (receptor is not None) and (self.ligand_cache[ligand_label] is not None):
                exp = AssayExperiment()
                exp.ligand_id = self.ligand_cache[ligand_label].id
                exp.protein_id = receptor.id
                exp.assay_type = assay_type
                exp.assay_description = row['Assay Type'] if pd.notna(row['Assay Type']) else None

                act_type = row['Activity Type']

                if pd.notna(act_type) and act_type.startswith('p'):
                    # Already a p-log value — both fields are the same, no conversion needed
                    if pd.notna(row['Value']):
                        val = round(float(row['Value']), 2)
                        exp.standard_activity_value = val
                        exp.p_activity_value = val
                    else:
                        print("Warning: missing Value field in: ", row)
                        exp.standard_activity_value = None
                        exp.p_activity_value = None
                elif pd.notna(act_type) and pd.notna(row['Value']):
                    raw = float(row['Value'])
                    exp.standard_activity_value = round(raw, 2) if act_type != 'Emax' else round(float(row['Emax (%)']), 2)
                    if pd.notna(row['Unit']) and raw != 0 and act_type!='Specific binding':
                        factor = unit_to_factor.get(str(row['Unit']).strip())
                        exp.p_activity_value = round(-math.log10(raw * factor), 2) if factor is not None else None
                    else:
                        exp.p_activity_value = None
                else:
                    exp.standard_activity_value = None
                    exp.p_activity_value = None

                exp.p_activity_ranges = None
                exp.standard_relation = row['Sign'] if pd.notna(row['Sign']) else None
                exp.value_type = row['Activity Type'] if pd.notna(row['Activity Type']) else None
                exp.assay_classification = Command.fetch_assay_classification(exp.value_type, row['Assay Type'])
                exp.source = 'Evolvus'
                exp.document_chembl_id = None
                exp.reference_ligand = row['Reference ligand name'] if pd.notna(row['Reference ligand name']) else None
                exp.qualitative_activity = row['Qualitative activity'] if pd.notna(row['Qualitative activity']) else None
                bioacts.append(exp)

                # stash publication to attach later (after PKs exist)
                if pd.notna(row.get('Reference')):
                    pub = Command.fetch_publication(row['Reference'])
                    # stash the publication for the current exp index within this batch
                    pub_links.append((len(bioacts) - 1, pub.id))

            else:
                if row['UniProt'] not in self.errors:
                    self.errors[row['UniProt']] = [row['Ligand Name']]
                else:
                    self.errors[row['UniProt']].append(row['Ligand Name'])

            # if (len(bioacts) == Command.bulk_size) or (index == bio_entries - 1):
            if (len(bioacts) == Command.bulk_size) or (index == bio_entries - 1):
                created = AssayExperiment.objects.bulk_create(bioacts)
                # Build through rows
                Through = AssayExperiment.publication.through
                through_rows = [
                    Through(assayexperiment_id=created[idx].id, publication_id=pub_id)
                    for idx, pub_id in pub_links
                    if created[idx].id is not None
                ]

                # Bulk create M2M links (ignore duplicates if needed)
                Through.objects.bulk_create(through_rows, ignore_conflicts=True)
                print("Inserted", index, "out of", bio_entries, "bioactivities")
                bioacts = []
                pub_links = []

            Command.assign_ligand_target_pairing(self.ligand_cache[ligand_label], receptor, None, row['Activity Type'])

    @staticmethod
    def calculate_potency_and_affinity():
        ligand_target_couples = AssayExperiment.objects.exclude(p_activity_value='None').exclude(value_type=None).exclude(p_activity_value=None).values_list('ligand_id',
                                                                                                     'protein_id',
                                                                                                     'value_type',
                                                                                                     'p_activity_value').distinct()
        connections = {}
        SI_dict = {}
        norm_to_raws = {}

        for pair in ligand_target_couples:
            raw_vt = pair[2]
            norm_vt = raw_vt if raw_vt.startswith(('p', 'P')) else 'p' + raw_vt
            norm_to_raws.setdefault(norm_vt, set()).add(raw_vt)
            if pair[0] not in connections:
                connections[pair[0]] = {}
            if pair[1] not in connections[pair[0]]:
                connections[pair[0]][pair[1]] = {}
            if norm_vt not in connections[pair[0]][pair[1]]:
                connections[pair[0]][pair[1]][norm_vt] = []
            connections[pair[0]][pair[1]][norm_vt].append(float(pair[3]))

        for ligand in connections:
            for target in connections[ligand]:
                for value in connections[ligand][target]:
                    connections[ligand][target][value] = round(statistics.mean(connections[ligand][target][value]),2)

        for ligand in connections:
            SI_dict[ligand] = {}
            all_vts = set(vt for tgt in connections[ligand] for vt in connections[ligand][tgt])
            for norm_vt in all_vts:
                targets_with_vt = [(tg, connections[ligand][tg][norm_vt])
                                   for tg in connections[ligand] if norm_vt in connections[ligand][tg]]
                vt_sort = sorted(targets_with_vt, key=lambda x: x[1], reverse=True)
                SI_dict[ligand][norm_vt] = {'Count': len(vt_sort)}
                if len(vt_sort) > 1:
                    for tg, val in targets_with_vt:
                        if tg == vt_sort[0][0]:
                            va = val - vt_sort[1][1]
                            SI_dict[ligand][norm_vt][tg] = int(10 ** abs(va))
                        else:
                            va = val - vt_sort[0][1]
                            SI_dict[ligand][norm_vt][tg] = -int(10 ** abs(va))

        vt_to_field = {}
        for cls in AssayClassification.objects.values('unit', 'assay_type'):
            unit = cls['unit']
            norm_unit = unit if unit.startswith(('p', 'P')) else 'p' + unit
            if 'Binding' in cls['assay_type']:
                vt_to_field[norm_unit] = 'affinity'
            elif 'Functional' in cls['assay_type']:
                vt_to_field[norm_unit] = 'potency'

        #This step is kinda slow and may be sped up using bulk_update
        #but I don't know if you can apply bulk_update with filters
        for lig in SI_dict:
            for norm_vt in SI_dict[lig]:
                field = vt_to_field.get(norm_vt)
                if field is None:
                    continue
                count = SI_dict[lig][norm_vt]['Count']
                raw_vts = norm_to_raws.get(norm_vt, [norm_vt])
                for tg, si_val in ((k, v) for k, v in SI_dict[lig][norm_vt].items() if k != 'Count'):
                    qs = AssayExperiment.objects.filter(
                        ligand_id=lig, protein_id=tg,
                        value_type__in=raw_vts,
                        p_activity_value__isnull=False,
                    )
                    if field == 'affinity':
                        qs.update(affinity=si_val, count_affinity_test=count)
                    else:
                        qs.update(potency=si_val, count_potency_test=count)

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
        effect_slug = None
        if unit in ['EC50', 'pEC50', 'Emax']:
            effect_slug = 'stimulatory'
        elif unit in ['pA2', 'pKB', 'IC50', 'pIC50', 'pKb', 'K(b)', 'pK(b)', 'KB']:
            effect_slug = 'inhibitory'
        elif unit in ['pKi', 'pKd', 'Ki', 'Kd', 'Potency', 'Specific binding', 'K(i)', 'K(d)', 'K(act)', 'pK(d)', 'pK(i)']:
            effect_slug = 'binding'
        elif unit in ['AC50','ED50']:
            desc = assay_description.lower()
            # stimulatory if “agonist” but not part of “antagonist”
            if 'agonist' in desc and 'antagonist' not in desc:
                effect_slug = 'stimulatory'
            # inhibitory if either “antagonist” or “inhibitor”
            elif 'antagonist' in desc or 'inhibitor' in desc:
                effect_slug = 'inhibitory'
            else:
                effect_slug = 'unknown'
        elif unit=='XC50':
            effect_slug = 'unknown'
        else:
            print(f"This value was not in the selection: {unit}")
        # 2. Lookup or fall back
        effect = None
        if effect_slug:
            try:
                effect = LigandEffect.objects.get(slug=effect_slug)
            except LigandEffect.DoesNotExist:
                print(f"No LigandEffect with slug {effect_slug} found; leaving effect NULL")
        # 3. Create + save
        with transaction.atomic():
            pairing = LigandTargetPairing(
                ligand=ligand,
                target=receptor,
                effect=effect,
            )
            pairing.save()

        return pairing

    @staticmethod
    def read_data(datadir, filename, sheet):
        source_file_path = os.sep.join([datadir, filename]).replace('//', '/')
        xls = pd.ExcelFile(source_file_path)
        df = pd.read_excel(xls, sheet, dtype=str)
        return df
