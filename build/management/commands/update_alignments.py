from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from protein.models import (Protein, ProteinConformation, ProteinSequenceType, ProteinSegment,
    ProteinConformationTemplateStructure)
from structure.models import Structure
from common.alignment import Alignment

import logging


class Command(BaseCommand):
    help = 'Updates protein alignments based on structure data'

    logger = logging.getLogger(__name__)

    def handle(self, *args, **options):
        # find protein templates
        try:
            self.update_alignments()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def update_alignments(self):
        self.logger.info('UPDATING PROTEIN ALIGNMENTS') 

        # segments, move to settings
        segment_slugs = settings.COMPARISON_SEGMENTS
        segments = ProteinSegment.objects.filter(slug__in=segment_slugs)

        # fetch all conformations of wild-type proteins
        protein_sequence_type = ProteinSequenceType.objects.get(slug='wt')
        pconfs = ProteinConformation.objects.filter(protein__sequence_type=protein_sequence_type).select_related(
            'protein')

        # fetch wild-type sequences of receptors with available structures
        structures = Structure.objects.order_by('protein_conformation__protein__parent', 'resolution').distinct(
            'protein_conformation__protein__parent').select_related('protein_conformation__protein__parent')
        sps = []
        for structure in structures:
            sps.append(structure.protein_conformation.protein.parent) # use the wild-type sequence

        # find templates
        # for pconf in pconfs:
        #     # overall
        #     print("Overall template for {}".format(pconf))
        #     template = self.find_segment_template(pconf, sps, segments)
        #     pconf.template_structure = self.fetch_template_structure(structures, template.protein.entry_name)
        #     pconf.save()

        #     # for each segment
        #     for segment in segments:
        #         print("{}".format(segment.slug))
        #         template = self.find_segment_template(pconf, sps, [segment])
        #         template_structure = self.fetch_template_structure(structures, template.protein.entry_name)
        #         pcts, created = ProteinConformationTemplateStructure.objects.get_or_create(protein_conformation=pconf,
        #             protein_segment=segment, defaults={'structure': template_structure})
        #         if pcts.structure != template_structure:
        #             pcts.structure = template_structure
        #             pcts.save()

        self.logger.info('COMPLETED UPDATING PROTEIN ALIGNMENTS')