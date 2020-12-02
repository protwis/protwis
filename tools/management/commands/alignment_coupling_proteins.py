from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from django.db.models import Q
from django.template.loader import render_to_string
from django.shortcuts import render
from django.template import loader
from django.template import Template, Context
from protein.models import Protein, ProteinConformation, ProteinAlias, ProteinSegment, ProteinFamily, Gene, ProteinGProtein, ProteinGProteinPair
from residue.models import ResidueGenericNumber, ResidueGenericNumberEquivalent
from common import definitions
from common.selection import SelectionItem
from common.alignment_gpcr import Alignment
import xlsxwriter, xlrd


import logging, json, os
from collections import OrderedDict



class Command(BaseCommand):

    help = "Output all uniprot mappings"

    logger = logging.getLogger(__name__)


    def handle(self, *args, **options):


        ## For GPCRs
        # couplings = ProteinGProteinPair.objects.filter(source='Aska').distinct('protein').prefetch_related('protein')
        # print(len(couplings))

        # segments = ProteinSegment.objects.filter(proteinfamily="GPCR")
        # for c in couplings:
        #     p = c.protein
        #     protein_orthologues = Protein.objects.filter(family=p.family, parent=None, source__name='SWISSPROT').all()
        #     print(p,len(protein_orthologues))

        #     # create an alignment object
        #     a = Alignment()
        #     a.show_padding = False

        #     # load data from selection into the alignment
        #     a.load_proteins(protein_orthologues)
        #     a.load_segments(segments)

        #     # build the alignment data matrix
        #     a.build_alignment()

        #     string = render_to_string('alignment/alignment_fasta.html',{'a': a})
        #     f = open('manbir/{}.fasta'.format(c.protein.entry_name), 'w')
        #     f.write(string)
        #     f.close()
        #     # break


        # For Gproteins
        couplings = ProteinGProteinPair.objects.filter(source='Aska').distinct('g_protein_subunit').prefetch_related('g_protein_subunit')
        print(len(couplings))

        segments = ProteinSegment.objects.filter(proteinfamily="Alpha")
        for c in couplings:
            p = c.g_protein_subunit
            protein_orthologues = Protein.objects.filter(family=p.family, parent=None, source__name='SWISSPROT').all()
            print(p,len(protein_orthologues))

            # create an alignment object
            a = Alignment()
            a.show_padding = False

            # load data from selection into the alignment
            a.load_proteins(protein_orthologues)
            a.load_segments(segments)

            # build the alignment data matrix
            a.build_alignment()

            string = render_to_string('alignment/alignment_fasta.html',{'a': a})
            f = open('manbir/{}.fasta'.format(p.entry_name), 'w')
            f.write(string)
            f.close()
            # break

            


