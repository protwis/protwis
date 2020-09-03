from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.utils.text import slugify
from django.conf import settings
from django.db import IntegrityError
from ligand.functions import get_or_make_ligand
from common.models import WebLink, WebResource, Publication
from protein.models import Protein
from ligand.models import Ligand, LigandProperities, LigandRole, LigandType, ChemblAssay, AssayExperiment
from ligand.models import LigandVendorLink, LigandVendors
from ligand.functions import get_or_make_ligand
from common.tools import fetch_from_web_api
from collections import defaultdict
import requests
from optparse import make_option
import logging
import shlex
import csv
import os
import json
from collections import OrderedDict
import datetime
"""
building ligads assays without vendor info
"""
class Command(BaseBuild):
    help = 'Reads source data and creates links to other databases'

    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
            type=int,
            action='store',
            dest='proc',
            default=1,
            help='Number of processes to run')
        parser.add_argument('--filename', action='append', dest='filename',
            help='Filename to import. Can be used multiple times')
        parser.add_argument('--test_run', action='store_true', help='Skip this during a test run', default=False)

    logger = logging.getLogger(__name__)
    ligand_cache = {}
    publication_cache = {}
    # source file directory
    links_data_dir = os.sep.join([settings.DATA_DIR, 'ligand_assay'])

    # wr = WebResource.objects.get(slug='chembl_ligand')
    # wr_pubchem = WebResource.objects.get(slug='pubchem')

    def handle(self, *args, **options):
        if options['test_run']:
            print('Skipping in test run')
            return
        if options['filename']:
            filenames = options['filename']
        else:
            filenames = os.listdir(self.links_data_dir)
        print("size of current dict",filenames)
        self.data = self.load_data(filenames) #reading the entry file from the data filder.
        print("Total rows in dataset",len(self.data))
        # Load all ligands ( possible to skip if believed to be included or already imported )
        #self.prepare_input(options['proc'], self.chembl_mol_ids, 0)
        # Insert the actual data points
        self.prepare_input(options['proc'], self.data ,1)

    def fetch_protein(self, protein_json):
        """
        fetch receptor with Protein model
        requires: protein id, source
        """
        protein = None
        if Protein.objects.get(entry_name = protein_json):
            protein = Protein.objects.get(entry_name=protein_json)

        return protein

    def fetch_publication(self, publication_doi):
        """
        fetch publication with Publication model
        requires: publication doi or pmid
        """
        pub = None
        if Publication.objects.get(reference=publication_doi):
            pub = Publication.objects.get(reference=publication_doi)
        return pub

    def fetch_ligand(self, ligand_id, smiles):
        """
        fetch ligands with Ligand model
        requires: ligand id, ligand id type, ligand name
        requires: source_file name
        """
        l = None
        try:
            if ligand_id in self.ligand_cache:
                l = self.ligand_cache[ligand_id]
            else:
                l = Ligand.objects.get(name=ligand_id)
                if l:
                    return l
                else:
                    l = get_or_make_ligand(smiles, 'SMILES', ligand_id,  )
        except Exception as msg:
            l = None
            # print('ligand_id---',l,'\n end')
        return l

    def fetch_assay(self, assay_id):
        assay = None
        try:
            if ChemblAssay.objects.get(assay_id = assay_id):
                assay = ChemblAssay.objects.get(assay_id = assay_id)
        except:
            assay = None
        return assay

    def main_func(self,positions, iteration,count,lock):
        #####Create chembl compound link and connect it to the corresponding ligand/cid#####
        print('---main func--')

            # First load makes sure ligands are there
        for i in self.data:
            ligand	= self.fetch_ligand(i['ligand'],i['smiles'])
            if not ligand:
                print(i['ligand'])
            protein	= self.fetch_protein(i['protein'])
            assay	= self.fetch_assay(i['assay'])
            publication = None

            if i['publication'] != None:
                publication	= self.fetch_publication(i['publication'])

            assay_experiment = AssayExperiment(
                ligand	= ligand,
                protein	= protein,
                assay	= assay,
                publication = publication,
                assay_type	= i['assay_type'],
                assay_description= i['assay_description'],
                pchembl_value	= i['pchembl_value'],
                published_value	= i['published_value'],
                published_relation	= i['published_relation'],
                published_type	= i['published_type'],
                published_units	= i['published_units'],
                standard_value	= i['standard_value'],
                standard_relation	= i['standard_relation'],
                standard_type	= i['standard_type'],
                standard_units	= i['standard_units'],
                chembl	= i['chembl'],
                smiles	= i['smiles'],
                activity	= i['activity'],
                document_chembl_id	= i['document_chembl_id'],
                cell_line	= 'delete_me'
            )
            assay_experiment.save()
            print('saved')
            count.value +=1


    ##read pre-generated file and extract the chembl_ids
    def load_data(self, filenames=False):
        for filename in filenames:
            full_file = os.sep.join([self.links_data_dir,filename])
            data = []
            if filename[-2:]=='gz':
                import gzip, io
                with gzip.open(full_file, "r") as fii:
                    data_test = json.load(fii)
                    for x,row in enumerate(data_test):
                        data.append(row)

        return data
