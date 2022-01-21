from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from build.management.commands.build_ligand_functions import *

from protein.models import Protein, Species
from ligand.models import Ligand, LigandProperties, LigandType, LigandRole, TestLigand
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
        TestLigand.objects.all().delete()
        LigandProperties.objects.all().delete()
        print("# Old data removed")

    def analyse_rows(self):
        """
        Fetch data to models
        Saves to DB
        """
        ligands_ids = []
        print('\n---Starting---')

        print('\n#1 Reading info sheets from Guide to Pharmacology')
        gtp_uniprot_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/GtP_to_UniProt_mapping.csv", 7 * 24 * 3600)
        gtp_uniprot = pd.read_csv(gtp_uniprot_link)

        gtp_complete_ligands_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/ligands.csv", 7 * 24 * 3600)
        gtp_complete_ligands = pd.read_csv(gtp_complete_ligands_link)

        gtp_ligand_mapping_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/ligand_id_mapping.csv", 7 * 24 * 3600)
        gtp_ligand_mapping = pd.read_csv(gtp_ligand_mapping_link)

        gtp_detailed_endogenous_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/detailed_endogenous_ligands.csv", 7 * 24 * 3600)
        gtp_detailed_endogenous = pd.read_csv(gtp_detailed_endogenous_link)

        gtp_interactions_link = get_or_create_url_cache("https://www.guidetopharmacology.org/DATA/interactions.csv", 7 * 24 * 3600)
        gtp_interactions = pd.read_csv(gtp_interactions_link)

        print('\n#2 Retrieving IUPHAR ids from UniProt ids')
        all_data_dataframe, iuphar_ids = self.compare_proteins(gtp_uniprot)

        print('\n#3 Retrieving ALL ligands from GTP associated to GPCRs')
        ligands_ids = self.obtain_ligands(ligands_ids, gtp_interactions, iuphar_ids, ['target_id','ligand_id']) #4181
        ligands_ids = self.obtain_ligands(ligands_ids, gtp_detailed_endogenous, iuphar_ids, ['Target ID','Ligand ID']) #4325
        ligands_ids = [x for x in ligands_ids if x == x] #remove nan

        print('\n#4 Mapping ligands to the tested receptors')
        all_data_dataframe = self.add_ligands_to_dataframe(all_data_dataframe, iuphar_ids, ligands_ids, gtp_interactions, gtp_detailed_endogenous)

        print('\n#5 Collating all info from GPCR related ligands in the GTP')
        complete_data, ligand_data = self.get_ligands_data(all_data_dataframe, ligands_ids, gtp_complete_ligands, gtp_ligand_mapping)

        print('\n#6 Saving the ligands in the models')
        self.save_the_ligands_save_the_world(ligand_data)

        print('\n---Finished---')

    @staticmethod
    def compare_proteins(gtp_data):
        data = pd.DataFrame(columns=['Name','UniProt','IUPHAR'])
        gpcrdb_proteins = Protein.objects.filter(family__slug__startswith="00", sequence_type__slug="wt").values_list('entry_name','accession')
        # FOR DEBUGGING
        #gpcrdb_proteins = Protein.objects.filter(family__slug__startswith="001_002", sequence_type__slug="wt").values_list('entry_name','accession')
        for protein in gpcrdb_proteins:
            try:
                new_row = [protein[0], protein[1]]
                new_row.append(gtp_data.loc[gtp_data['uniprot_id'] == protein[1], 'iuphar_id'].values[0])
                data.loc[len(data)] = new_row
            except:
                continue
        iuphar_ids = list(data['IUPHAR'].unique())
        return data, iuphar_ids

    @staticmethod
    def obtain_ligands(ligands, data, compare_set, labels):
        interactions_targets = list(data[labels[0]].unique())
        targets_with_ligands = set(interactions_targets).intersection(set(compare_set))
        for target_id in targets_with_ligands:
            ligands = ligands + list(data.loc[data[labels[0]] == target_id, labels[1]])
        ligands = list(set(ligands))
        return ligands

    #all_data_dataframe, iuphar_ids, ligands_ids, gtp_interactions, gtp_detailed_endogenous
    @staticmethod
    def add_ligands_to_dataframe(df, ids, ligands, interactions, endogenous):
        uniprots = df['UniProt'].tolist()
        new_df = pd.DataFrame(columns=['Receptor', 'IUPHAR', 'UniProt', 'Ligand ID'])
        for target in uniprots:
            inter_ligs = interactions.loc[interactions['target_uniprot'] == target, 'ligand_id'].tolist()
            endo_ligs = endogenous.loc[endogenous['Target UniProt ID'] == target, 'Ligand ID'].tolist()
            tot_ligs = list(set(inter_ligs + endo_ligs))
            if len(tot_ligs) > 0:
                for lig in tot_ligs:
                    row = list(df.loc[df['UniProt'] == target, ['Name','IUPHAR']].values[0])
                    row = row + [target, lig]
                    new_df.loc[len(new_df)] = row
        return new_df

    #gtp_complete_ligands, gtp_ligand_mapping
    @staticmethod
    def get_ligands_data(data, ligands, complete_ligands, ligand_mapping):
        cols = data.columns.tolist()
        full_info = ['Name','Species','Type','Approved','Withdrawn','Labelled','Radioactive', 'PubChem SID', 'PubChem CID',
                     'UniProt id','IUPAC name', 'INN', 'Synonyms','SMILES','InChIKey','InChI','GtoImmuPdb','GtoMPdb']
        weblinks = ['ChEMBl ID','Chebi ID','CAS','DrugBank ID','Drug Central ID']
        cols = cols + full_info + weblinks
        data = data.reindex(columns=cols)
        for lig in ligands:
            for info in full_info:
                data.loc[data['Ligand ID'] == lig, info] = complete_ligands.loc[complete_ligands['Ligand id'] ==  lig, info].values[0]
            for link in weblinks:
                data.loc[data['Ligand ID'] == lig, link] = ligand_mapping.loc[ligand_mapping['Ligand id'] ==  lig, link].values[0]
        ligs_info = data.drop(columns=['Receptor','IUPHAR','UniProt']).drop_duplicates()
        return data, ligs_info

    @staticmethod
    def save_the_ligands_save_the_world(lig_df):
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

        #Process dataframe and replace nan with None
        lig_df = lig_df.astype(str)
        for column in lig_df:
            lig_df[column] = lig_df[column].replace({'nan':None})

        #start parsing the GtP ligands
        issues = []
        for index, row in lig_df.iterrows():
            if row['Name']:
                ids = {}
                for key, value in weblink_keys.items():
                    if row[value] != None:
                        ids[key] = str(row[value])

                if types_dict[row['Type']] == "small-molecule":
                    ligand = get_or_create_smallmolecule(row['Name'], ids)
                    if not ligand:
                        issues.append([row['Name'], row["Ligand ID"]])
                else:
                    # Handling peptide or protein => skipping InChIKey and SMILES
                    lt = LigandType.objects.get(slug=types_dict[row['Type']])
                    lig_props = LigandProperties()
                    lig_props.ligand_type = lt
                    lig_props.sequence = None # TODO => update sequence for peptides
                    lig_props.save()

                    # Weblinks for each resource (PubChem, GtP, ChEMBL, ..)
                    for key in ids:
                        if key not in ["smiles", "inchikey"] and ids[key]:
                            value = ids[key]
                            if value.isnumeric():
                                value = int(value)

                            wr = WebResource.objects.get(slug=key)
                            wl, created = WebLink.objects.get_or_create(index=ids[key], web_resource=wr)
                            lig_props.web_links.add(wl)

                    ligand = TestLigand()
                    ligand.name = row['Name']
                    ligand.type = lt
                    ligand.properties = lig_props
                    ligand.canonical = True
                    ligand.ambiguous_alias = False
                    ligand.save()

        print("DONE - Ligands with issues:")
        print(issues)
