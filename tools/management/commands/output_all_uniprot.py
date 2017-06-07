from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from django.db.models import Q
from django.template.loader import render_to_string
from protein.models import Protein
from residue.models import ResidueGenericNumber, ResidueGenericNumberEquivalent
from common import definitions
from common.selection import SelectionItem
from common.alignment_gpcr import Alignment
import xlsxwriter, xlrd


import logging, json, os



class Command(BaseCommand):

    help = "Output all uniprot mappings"

    logger = logging.getLogger(__name__)


    def handle(self, *args, **options):

        #Get the proteins
        f = open('uniprot.json', 'w')
        ps = Protein.objects.filter(Q(source__name='SWISSPROT') | Q(source__name='TREMBL'),web_links__web_resource__slug='uniprot').all().prefetch_related('web_links__web_resource')
        print('total:',len(ps))
        mapping = {}
        for p in ps:
            uniprot = p.web_links.filter(web_resource__slug='uniprot').values_list('index', flat = True)
            mapping[p.entry_name] = list(uniprot)


        json.dump(mapping,f, indent=4, separators=(',', ': '))

        # print("Seqs: {}\tNot matching: {}".format(num_of_sequences, num_of_non_matching_sequences))
        # open("uniprot.txt", "w").write()
