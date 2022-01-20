from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db import IntegrityError
from django.utils.text import slugify
from django.http import HttpResponse, JsonResponse
from build.management.commands.base_build import Command as BaseBuild
from protein.models import Protein, Species
from ligand.models import Endogenous_GTP, Ligand, LigandProperties, LigandType, LigandRole, TestLigand
from common.models import WebLink, WebResource, Publication
import logging
import time
import requests
import statistics
import os
import pandas as pd
import numpy as np

#defining globals and URLs
missing_info = []
not_commented = []

#This command will build ALL the Guide to Pharmacology data starting
#from several data sheet. First implementation will feature all the
#downloaded file in the gtp_data folder, while second interaction
#will be structured to read data directly from the web given the
#correct URLs of the filenames

class Command(BaseBuild):
    mylog = logging.getLogger(__name__)
    mylog.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
    file_handler = logging.FileHandler('biasDataTest.log')
    file_handler.setLevel(logging.ERROR)
    file_handler.setFormatter(formatter)
    mylog.addHandler(file_handler)
    endogenous_data_dir = os.sep.join([settings.DATA_DIR, 'ligand_data', 'gtp_data'])

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
    def purge_old_data(self):
        print("\n# Purging data")
        TestLigand.objects.all().delete()
        LigandProperties.objects.all().delete()
        print("\n# Old data removed")

    def analyse_rows(self):
        """
        Fetch data to models
        Saves to DB
        """
        ligands_ids = []
        gtp_uniprot_link = "https://www.guidetopharmacology.org/DATA/GtP_to_UniProt_mapping.csv"
        gtp_complete_ligands_link = "https://www.guidetopharmacology.org/DATA/ligands.csv"
        gtp_ligand_mapping_link = "https://www.guidetopharmacology.org/DATA/ligand_id_mapping.csv"
        gtp_detailed_endogenous_link = "https://www.guidetopharmacology.org/DATA/detailed_endogenous_ligands.csv"
        gtp_interactions_link = "https://www.guidetopharmacology.org/DATA/interactions.csv"
        print('---Starting---')
        print('\n#1 Reading info sheets from Guide to Pharmacology')
        gtp_uniprot = pd.read_csv(gtp_uniprot_link)
        gtp_complete_ligands = pd.read_csv(gtp_complete_ligands_link)
        gtp_ligand_mapping = pd.read_csv(gtp_ligand_mapping_link)
        gtp_detailed_endogenous = pd.read_csv(gtp_detailed_endogenous_link)
        gtp_interactions = pd.read_csv(gtp_interactions_link)
        print('\n#2 Retrieving IUPHAR ids from UniProt ids')
        all_data_dataframe, iuphar_ids = self.compare_proteins(gtp_uniprot)
        print('\n#3 Retrieving ALL ligands from GTP associated to GPCRs')
        ligands_ids = self.obtain_ligands(ligands_ids, gtp_interactions, iuphar_ids, ['target_id','ligand_id']) #4181
        ligands_ids = self.obtain_ligands(ligands_ids, gtp_detailed_endogenous, iuphar_ids, ['Target ID','Ligand ID']) #4325
        #remove nan
        ligands_ids = [x for x in ligands_ids if x == x]
        print('\n#4 Mapping ligands to the tested receptors')
        all_data_dataframe = self.add_ligands_to_dataframe(all_data_dataframe, iuphar_ids, ligands_ids, gtp_interactions, gtp_detailed_endogenous)
        print('\n#5 Collating all info from GPCR related ligands in the GTP')
        complete_data, ligand_data = self.get_ligands_data(all_data_dataframe, ligands_ids, gtp_complete_ligands, gtp_ligand_mapping)
        print('\n#6 Saving the ligands in the models')
        self.save_the_ligands_save_the_world(ligand_data)
        print('\n---Finished---')

    # @staticmethod
    def compare_proteins(self, gtp_data):
        data = pd.DataFrame(columns=['Name','UniProt','IUPHAR'])
        gpcrdb_proteins = Protein.objects.filter(family__slug__startswith="00", sequence_type__slug="wt").values_list('entry_name','accession')
        for protein in gpcrdb_proteins:
            try:
                new_row = [protein[0], protein[1]]
                new_row.append(gtp_data.loc[gtp_data['uniprot_id'] == protein[1], 'iuphar_id'].values[0])
                data.loc[len(data)] = new_row
            except:
                continue
        iuphar_ids = list(data['IUPHAR'].unique())
        return data, iuphar_ids

    # @staticmethod
    def obtain_ligands(self, ligands, data, compare_set, labels):
        interactions_targets = list(data[labels[0]].unique())
        targets_with_ligands = set(interactions_targets).intersection(set(compare_set))
        for target_id in targets_with_ligands:
            ligands = ligands + list(data.loc[data[labels[0]] == target_id, labels[1]])
        ligands = list(set(ligands))
        return ligands

    #all_data_dataframe, iuphar_ids, ligands_ids, gtp_interactions, gtp_detailed_endogenous
    # @staticmethod
    def add_ligands_to_dataframe(self, df, ids, ligands, interactions, endogenous):
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
    # @staticmethod
    def get_ligands_data(self, data, ligands, complete_ligands, ligand_mapping):
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
        ligs_info = data.drop(['Receptor','IUPHAR','UniProt'], 1).drop_duplicates()
        return data, ligs_info

    # @staticmethod
    def save_the_ligands_save_the_world(self, lig_df):
        types_dict = {'Inorganic': 'small-molecule',
                      'Metabolite': 'small-molecule',
                      'Natural product': 'small-molecule',
                      'Peptide': 'peptide',
                      'Synthetic organic': 'small-molecule',
                      'Antibody': 'protein',
                      None: 'na'}
        weblink_keys = {'pubchem': 'PubChem CID',
                        'gtoplig': 'Ligand ID',
                        'chembl_ligand': 'ChEMBl ID',
                        'drugbank': 'DrugBank ID',
                        'drug_central': 'Drug Central ID'}

        #remove radioactive ligands
        lig_df = lig_df.loc[lig_df['Radioactive'] != 'yes']

        #convert dataframe into list of dictionaries
        lig_df = lig_df.astype(str)
        for column in lig_df:
            lig_df[column] = lig_df[column].replace({'nan':None})
        lig_dicts = lig_df.to_dict("records")

        #start parsing the list of dictionaties
        for row in lig_dicts:
            if row['Name']:
                lt = LigandType.objects.get(slug=types_dict[row['Type']])
                #check if inchikey exists in LigandProperties
                check_inchi = None
                if row['InChIKey']:
                    check_inchi = LigandProperties.objects.filter(inchikey=row['InChIKey']).first()

                if check_inchi:
                    ligand = TestLigand()
                    ligand.name = row['Name']
                    ligand.type = lt
                    ligand.properties = check_inchi
                    ligand.canonical = True
                    ligand.ambiguous_alias = False
                else:
                    lig_props = LigandProperties()
                    lig_props.smiles = row['SMILES']
                    lig_props.inchikey = row['InChIKey']
                    lig_props.mw = None
                    lig_props.rotatable_bonds = None
                    lig_props.hacc = None
                    lig_props.hdon = None
                    lig_props.logp = None
                    lig_props.ligand_type = lt
                    lig_props.sequence = None
                    lig_props.save()
                    for key, value in weblink_keys.items():
                        if row[value]:
                            try:
                                index = int(row[value])
                            except (ValueError,TypeError):
                                index = row[value]
                    #weblinks for each resource (PubChem, GtP, ChEMBL, )
                            resource = WebResource.objects.get(slug=key)
                            wl, created = WebLink.objects.get_or_create(index=index, web_resource=resource)
                            lig_props.web_links.add(wl)
                    ligand = TestLigand()
                    ligand.name = row['Name']
                    ligand.type = lt
                    ligand.properties = lig_props
                    ligand.canonical = True
                    ligand.ambiguous_alias = False
                ligand.save()
