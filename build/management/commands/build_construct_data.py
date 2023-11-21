from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.utils.text import slugify
from django.db import IntegrityError
from common.tools import test_model_updates
from protein.models import Protein, ProteinConformation
from residue.models import Residue
from structure.models import Structure
from construct.models import (Construct,Crystallization,CrystallizationLigandConc,ChemicalType,Chemical,ChemicalConc,ChemicalList,
CrystallizationMethods,CrystallizationTypes,ChemicalListName,ContributorInfo,ConstructMutation,ConstructInsertion,ConstructInsertionType,
ConstructDeletion,ConstructModification,CrystalInfo,ExpressionSystem,Solubilization,PurificationStep,Purification)
from construct.functions import add_construct, fetch_pdb_info

from ligand.models import Ligand, LigandType, LigandRole

from optparse import make_option
import logging
import csv
import os
import json
import datetime
import django.apps

class Command(BaseCommand):
    help = 'Build construct data'

    def add_arguments(self, parser):
        parser.add_argument('--filename', action='append', dest='filename',
            help='Filename to import. Can be used multiple times')
        parser.add_argument('--local', action='store_true', dest='local', default=False,
            help='Read local construct files')
        parser.add_argument('--purge', help='Purge all existing records', default=False, action='store_true')

    logger = logging.getLogger(__name__)

        # source file directory
    construct_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data','construct_data'])
    construct_data_local_dir = "../files/construct_data"

    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    def handle(self, *args, **options):
        if options['filename']:
            filenames = options['filename']
        else:
            filenames = False

        if options['local']:
            local_fill = True
        else:
            local_fill = False

        if options['purge']:
            self.purge_construct_data()
            self.tracker = {}
            test_model_updates(self.all_models, self.tracker, initialize=True)

        if not local_fill:
            self.create_construct_data(filenames)
            test_model_updates(self.all_models, self.tracker, check=True)
        else:
            self.create_construct_local_data()
            test_model_updates(self.all_models, self.tracker, check=True)
        # except Exception as msg:
        #     print("ERROR: "+str(msg))
        #     self.logger.error(msg)


    def purge_construct_data(self):
        Construct.objects.all().delete()
        Crystallization.objects.all().delete()
        ChemicalConc.objects.all().delete()
        Chemical.objects.all().delete()
        ChemicalType.objects.all().delete()
        ChemicalList.objects.all().delete()
        CrystallizationLigandConc.objects.all().delete()
        CrystallizationMethods.objects.all().delete()
        CrystallizationTypes.objects.all().delete()
        ChemicalListName.objects.all().delete()
        ContributorInfo.objects.all().delete()
        ConstructMutation.objects.all().delete()
        ConstructDeletion.objects.all().delete()
        ConstructInsertion.objects.all().delete()
        ConstructInsertionType.objects.all().delete()
        ConstructModification.objects.all().delete()
        CrystalInfo.objects.all().delete()
        ExpressionSystem.objects.all().delete()
        Solubilization.objects.all().delete()
        Purification.objects.all().delete()
        PurificationStep.objects.all().delete()

    # def create_construct_local_data(self, filenames=False):
    #     self.logger.info('ADDING EXPERIMENTAL CONSTRUCT DATA')

    #     # read source files
    #     if not filenames:
    #         #delete existing if nothing specific is defined
    #         self.purge_construct_data()
    #         filenames = os.listdir(self.construct_data_dir)

    #     for filename in sorted(filenames):
    #         print('dealing with',filename)
    #         if filename[-4:]!='json':
    #             continue
    #         filepath = os.sep.join([self.construct_data_dir, filename])
    #         with open(filepath) as json_file:
    #             d = json.load(json_file)
    #             add_construct(d)

    #     filenames = os.listdir(self.construct_data_local_dir)

    #     for filename in sorted(filenames):
    #         print('dealing with',filename)
    #         if filename[-4:]!='json':
    #             continue
    #         filepath = os.sep.join([self.construct_data_local_dir, filename])
    #         with open(filepath) as json_file:
    #             d = json.load(json_file)
    #             add_construct(d)

    #     if not filenames:
    #         structures = Structure.objects.all()
    #         for s in structures:
    #             pdbname = str(s)
    #             print(pdbname)
    #             try:
    #                 exists = Construct.objects.filter(structure__pdb_code__index=pdbname).exists()
    #                 if not exists:
    #                     print(pdbname)
    #                     protein = Protein.objects.filter(entry_name=pdbname.lower()).get()
    #                     d = fetch_pdb_info(pdbname,protein)
    #                     add_construct(d)
    #                 else:
    #                     print("Entry for",pdbname,"already there")
    #             except:
    #                 print(pdbname,'failed')

    def create_construct_data(self, filenames=False):
        self.logger.info('ADDING EXPERIMENTAL CONSTRUCT DATA')

        # read source files
        do_all = False
        if not filenames:
            do_all = True
            # self.purge_construct_data()
            # filenames = os.listdir(self.construct_data_dir)

        if filenames:
            for filename in filenames:
                if filename[-4:]!='json':
                    continue
                filepath = os.sep.join([self.construct_data_dir, filename])
                print('Adding '+filepath)
                with open(filepath) as json_file:
                    d = json.load(json_file)
                    add_construct(d)

        if do_all:
            structures = Structure.objects.all().exclude(structure_type__slug__startswith='af-')
            for s in structures:
                pdbname = str(s)
                try:
                    exists = Construct.objects.filter(structure__pdb_code__index=pdbname).exists()
                    if not exists:
                        # print(pdbname)
                        protein = Protein.objects.filter(entry_name=pdbname.lower()).get()
                        d = fetch_pdb_info(pdbname,protein)
                        add_construct(d)
                    else:
                        # pass
                        print("Entry for",pdbname,"already there")
                except:
                    print(pdbname,'failed')

        self.logger.info('COMPLETED CREATING EXPERIMENTAL CONSTRUCT DATA')
