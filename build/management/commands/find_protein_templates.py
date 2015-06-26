from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from protein.models import (Protein, ProteinConformation, ProteinSequenceType, ProteinSegment,
    ProteinConformationTemplateStructure)
from structure.models import Structure
from common.alignment import Alignment

import logging


class Command(BaseCommand):
    help = 'Goes though all protein records and finds the best structural templates for the whole 7TM domain, and ' \
        + 'each helix individually'

    logger = logging.getLogger(__name__)

    def handle(self, *args, **options):
        # find protein templates
        try:
            self.find_protein_templates()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def find_protein_templates(self):
        self.logger.info('ASSIGNING STRUCTURE TEMPLATES FOR PROTEINS')

        # segments
        segments = ProteinSegment.objects.filter(slug__in=settings.REFERENCE_POSITIONS)

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
        for pconf in pconfs:
            # overall
            template = self.find_segment_template(pconf, sps, segments)
            template_structure = self.fetch_template_structure(structures, template.protein.entry_name)
            pconf.template_structure = template_structure
            pconf.save()
            self.logger.info("Assigned {} as overall template for {}".format(template_structure, pconf))

            # for each segment
            for segment in segments:
                template = self.find_segment_template(pconf, sps, [segment])
                template_structure = self.fetch_template_structure(structures, template.protein.entry_name)
                pcts, created = ProteinConformationTemplateStructure.objects.get_or_create(protein_conformation=pconf,
                    protein_segment=segment, defaults={'structure': template_structure})
                if pcts.structure != template_structure:
                    pcts.structure = template_structure
                    pcts.save()
                self.logger.info("Assigned {} as {} template for {}".format(template_structure, segment, pconf))

        self.logger.info('COMPLETED ASSIGNING STRUCTURE TEMPLATES FOR PROTEINS')

    def find_segment_template(self, pconf, sconfs, segments):
            a = Alignment()
            a.load_reference_protein(pconf.protein)
            a.load_proteins(sconfs)
            a.load_segments(segments)
            a.build_alignment()
            a.calculate_similarity()

            return a.proteins[1]

    def fetch_template_structure(self, structures, template_protein):
        for structure in structures:
            structure_protein = structure.protein_conformation.protein.parent.entry_name
            if structure_protein == template_protein:
                return structure