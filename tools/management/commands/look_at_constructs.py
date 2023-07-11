from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from django.db.models import Q
from django.template.loader import render_to_string
from structure.models import Structure
from common import definitions
from common.selection import SelectionItem
from common.alignment_gpcr import Alignment
import xlsxwriter, xlrd
from construct.models import *
from construct.functions import *
from construct.tool import *


import logging, json, os
from collections import OrderedDict



class Command(BaseCommand):

    help = "Output all uniprot mappings"

    logger = logging.getLogger(__name__)



    def handle(self, *args, **options):

        for s in Structure.objects.all().exclude(structure_type__slug__startswith='af-'):
            slug = str(s)
            # if slug != '5N2R':
            #     continue
            print(s)
            protein = Protein.objects.filter(entry_name=slug.lower()).get()

            d = fetch_pdb_info(slug,protein)

            #delete before adding new
            #print(d['construct_crystal']['pdb_name'])
            #Construct.objects.filter(name__iexact=d['construct_crystal']['pdb_name']).delete()
            #add_construct(d)

            #cache.delete(d['construct_crystal']['pdb_name']+'_schematics')
            #cache.delete(d['construct_crystal']['pdb_name']+'_snake')
            #print(d['construct_sequences'])
            for seg,v in d['construct_sequences'].items():
                if 'b562' in seg:
                    print(slug,seg,v)
            # print(d)
