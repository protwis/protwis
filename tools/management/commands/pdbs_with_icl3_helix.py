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

from residue.models import *
from protein.models import *


import logging, json, os
from collections import OrderedDict



class Command(BaseCommand):

    help = "Output all uniprot mappings"

    logger = logging.getLogger(__name__)



    def handle(self, *args, **options):

        for s in Structure.objects.all().exclude(structure_type__slug__startswith='af-').order_by('protein_conformation__protein__parent__entry_name'):
            slug = str(s)
            pc = ProteinConformation.objects.filter(protein__entry_name=slug.lower()).get()
            rs = pc.residue_set.filter(generic_number__label='34x50')
            if not len(rs):
                print(slug, pc.protein.parent.entry_name, s.state, s.representative)
