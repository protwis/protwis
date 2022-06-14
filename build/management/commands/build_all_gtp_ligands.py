from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from build.management.commands.build_ligand_functions import *

from protein.models import Protein, Species
from ligand.models import Ligand, LigandType, LigandRole
from common.models import WebLink, WebResource
from common.tools import get_or_create_url_cache, fetch_from_web_api

import logging
import time
import statistics
import os
import pandas as pd

#This command will build ALL the Guide to Pharmacology data starting
#from several data sheet. First implementation will feature all the
#downloaded file in the gtp_data folder, while second interaction
#will be structured to read data directly from the web given the
#correct URLs of the filenames

class Command(BaseBuild):
    help = 'Update the Ligand Model with the data from ALL Guide to Pharmacology ligands'

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
            #let's skip some horrible mistakes for now
        if options['purge']:
            try:
                self.purge_old_data()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
        self.analyse_rows()

# pylint: disable=R0201
    @staticmethod
    def purge_old_data():
        print("# Purging data")
        Ligand.objects.all().delete()
        print("# Old data removed")

    def analyse_rows(self):
        """
        Fetch data to models
        Saves to DB
        """
        print('\n---Starting---')

        print('\n#1 Reading info sheets from Guide to Pharmacology')
        gtp_uniprot_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/GtP_to_UniProt_mapping.csv", 7 * 24 * 3600)
        gtp_uniprot = pd.read_csv(gtp_uniprot_link, dtype=str, header = 1)

        gtp_complete_ligands_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/ligands.csv", 7 * 24 * 3600)
        gtp_complete_ligands = pd.read_csv(gtp_complete_ligands_link, dtype=str, header = 1)

        gtp_ligand_mapping_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/ligand_id_mapping.csv", 7 * 24 * 3600)
        gtp_ligand_mapping = pd.read_csv(gtp_ligand_mapping_link, dtype=str, header = 1)

        gtp_detailed_endogenous_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/detailed_endogenous_ligands.csv", 7 * 24 * 3600)
        gtp_detailed_endogenous = pd.read_csv(gtp_detailed_endogenous_link, dtype=str, header = 1)

        gtp_interactions_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/interactions.csv", 7 * 24 * 3600)
        gtp_interactions = pd.read_csv(gtp_interactions_link, dtype=str, header = 1)
        gtp_interactions.columns = gtp_interactions.columns.str.strip() # Fixing GtP whitespace issues in the headers

        gtp_peptides_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/peptides.csv", 7 * 24 * 3600)
        gtp_peptides = pd.read_csv(gtp_peptides_link, dtype=str, header = 1)

        print('\n#2 Retrieving IUPHAR ids from UniProt ids')
        iuphar_ids = self.compare_proteins(gtp_uniprot)

        print('\n#3 Retrieving ALL ligands from GTP associated to GPCRs')
        bioactivity_ligands_ids = self.obtain_ligands(gtp_interactions, iuphar_ids, ['Target ID','Ligand ID'])
        endogenous_ligands_ids = self.obtain_ligands(gtp_detailed_endogenous, iuphar_ids, ['Target ID','Ligand ID'])
        ligand_ids = list(set(bioactivity_ligands_ids + endogenous_ligands_ids))

        print('\n#4 Collating all info from GPCR related ligands in the GTP')
        ligand_data = self.get_ligands_data(ligand_ids, gtp_complete_ligands, gtp_ligand_mapping)

        print('\n#5 Saving the ligands in the models')
        self.save_the_ligands_save_the_world(ligand_data, gtp_peptides)

        print('\n---Finished---')

    @staticmethod
    def compare_proteins(gtp_data):
        gpcrdb_proteins = Protein.objects.filter(family__slug__startswith="00", sequence_type__slug="wt").values_list('entry_name','accession')
        entries = gtp_data.loc[gtp_data['UniProtKB ID'].isin([protein[1].split("-")[0] for protein in gpcrdb_proteins]), ['UniProtKB ID', 'GtoPdb IUPHAR ID']]
        return list(entries['GtoPdb IUPHAR ID'].unique())

    @staticmethod
    def obtain_ligands(data, compare_set, labels):
        interactions_targets = list(data[labels[0]].unique())
        targets_with_ligands = set(interactions_targets).intersection(set(compare_set))
        ligands = list(set(data.loc[data[labels[0]].isin(targets_with_ligands), labels[1]]))
        ligands = [x for x in ligands if x == x] #remove nan
        return ligands

    #gtp_complete_ligands, gtp_ligand_mapping
    @staticmethod
    def get_ligands_data(ligands, complete_ligands, ligand_mapping):
        full_info = ['Ligand id', 'Name','Species','Type','Approved','Withdrawn','Labelled','Radioactive', 'PubChem SID', 'PubChem CID',
                     'UniProt id','IUPAC name', 'INN', 'Synonyms','SMILES','InChIKey','InChI','GtoImmuPdb','GtoMPdb']
        ligand_data = complete_ligands.loc[complete_ligands['Ligand id'].isin(ligands), full_info]

        weblinks = ['Ligand id', 'ChEMBl ID','Chebi ID','CAS','DrugBank ID','Drug Central ID']
        ligand_weblinks = ligand_mapping.loc[ligand_mapping['Ligand id'].isin(ligands), weblinks]

        ligand_complete = ligand_data.merge(ligand_weblinks, on="Ligand id")
        ligand_complete = ligand_complete.rename(columns={"Ligand id": "Ligand ID"})

        return ligand_complete

    @staticmethod
    def save_the_ligands_save_the_world(lig_df, pep_df):
        types_dict = {'Inorganic': 'small-molecule',
                      'Metabolite': 'small-molecule',
                      'Natural product': 'small-molecule',
                      'Peptide': 'peptide',
                      'Synthetic organic': 'small-molecule',
                      'Antibody': 'protein',
                      None: 'na'}
        weblink_keys = {'smiles': 'SMILES',
                        'inchikey': 'InChIKey',
                        'pubchem': 'PubChem CID',
                        'gtoplig': 'Ligand ID',
                        'chembl_ligand': 'ChEMBl ID',
                        'drugbank': 'DrugBank ID',
                        'drug_central': 'Drug Central ID'}

        #remove radioactive ligands
        lig_df = lig_df.loc[lig_df['Radioactive'] != 'yes']

        #Process dataframes and replace nan with None
        lig_df = lig_df.astype(str)
        for column in lig_df:
            lig_df[column] = lig_df[column].replace({'nan':None})
        pep_df = pep_df.astype(str)
        for column in pep_df:
            pep_df[column] = pep_df[column].replace({'nan':None})

        #start parsing the GtP ligands
        issues = []
        for index, row in lig_df.iterrows():
            if row['Name']:
                ids = {}
                for key, value in weblink_keys.items():
                    if row[value] != None:
                        ids[key] = str(row[value])
                        if is_float(ids[key]):
                            ids[key] = str(int(float(ids[key])))

                # When protein or peptide => get UNIPROT AND FASTA Here
                uniprot_ids = []
                if types_dict[row['Type']] != "small-molecule":
                    if pep_df["Ligand id"].eq(row["Ligand ID"]).any():
                        peptide_entry = pep_df.loc[pep_df["Ligand id"] == row["Ligand ID"]]

                        sequence = peptide_entry["Single letter amino acid sequence"].item()
                        if sequence != None and sequence != "" and len(sequence) < 1000:
                            ids["sequence"] = sequence

                        uniprot = peptide_entry["UniProt id"].item()
                        if uniprot != None and uniprot != "":
                            # TODO - when multiple UniProt IDs => clone ligand add new UniProt IDs for other species
                            uniprot_ids = uniprot.split("|")
                            ids["uniprot"] = uniprot_ids[0].strip()

                # Create or get ligand
                ligand = get_or_create_ligand(row['Name'], ids, types_dict[row['Type']], True, False)
                if ligand == None:
                    print("Issue with", row['Name'])
                    exit()

                # Process multiple UniProt IDs
                if len(uniprot_ids) > 1:
                    ligand.uniprot = ",".join(uniprot_ids)
                    ligand.save()

        print("DONE - Ligands with issues:")
        print(issues)
