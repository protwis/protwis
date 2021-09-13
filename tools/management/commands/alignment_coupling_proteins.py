from django.core.management.base import BaseCommand
from django.template.loader import render_to_string
from protein.models import Protein, ProteinAlias, ProteinSegment, ProteinFamily, Gene, ProteinCouplings
from residue.models import ResidueGenericNumberEquivalent
from common.alignment_gpcr import Alignment

import logging

class Command(BaseCommand):

    help = "Output all uniprot mappings"

    logger = logging.getLogger(__name__)


    def handle(self, *args, **options):
        # For Gproteins
        couplings = ProteinCouplings.objects.filter(source='Inoue').distinct('g_protein_subunit').prefetch_related('g_protein_subunit')
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
