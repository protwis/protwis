from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from protein.models import ProteinConformation, ProteinSegment
from structure.models import StructureSegment
from residue.models import Residue
from residue.functions import *

import logging


class Command(BaseCommand):
    help = 'Updates protein alignments based on structure data'

    logger = logging.getLogger(__name__)

    generic_numbers_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'generic_numbers'])
    ref_position_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'reference_positions'])

    def handle(self, *args, **options):
        # find protein templates
        try:
            self.update_alignments()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def update_alignments(self):
        self.logger.info('UPDATING PROTEIN ALIGNMENTS')

        schemes = parse_scheme_tables(self.generic_numbers_source_dir)

        # fetch protein conformations
        segments = ProteinSegment.objects.filter(partial=False)
        pconfs = ProteinConformation.objects.filter(protein__sequence_type__slug='wt').select_related(
            'protein__residue_numbering_scheme__parent', 'template_structure')
        for pconf in pconfs:
            sequence_number_counter = 0
            
            # read reference positions for this protein
            ref_position_file_path = os.sep.join([self.ref_position_source_dir, pconf.protein.entry_name + '.yaml'])
            ref_positions = load_reference_positions(ref_position_file_path)

            # determine segment ranges, and update residues residues
            nseg = len(segments)
            for i, segment in enumerate(segments):
                self.logger.info("Updating segment borders for {} of {}, using template {}".format(segment.slug, pconf,
                    pconf.template_structure))

                # find template segment (for segments borders)
                try:
                    tss = StructureSegment.objects.get(structure=pconf.template_structure, protein_segment=segment)
                except StructureSegment.DoesNotExist:
                    continue

                if segment.slug in settings.REFERENCE_POSITIONS:
                    # get reference positions of this segment (e.g. 1x50)
                    segment_ref_position = settings.REFERENCE_POSITIONS[segment.slug]

                    # template segment reference residue number
                    tsrrn = Residue.objects.get(protein_conformation=pconf.template_structure.protein_conformation,
                        generic_number__label=segment_ref_position)

                    # number of residues before and after the reference position
                    tpl_res_before_ref = tsrrn.sequence_number - tss.start
                    tpl_res_after_ref = tss.end - tsrrn.sequence_number

                    # FIXME check whether this segments actually needs an update

                    segment_start = (ref_positions[segment_ref_position]
                        - tpl_res_before_ref)
                    segment_end = (ref_positions[segment_ref_position]
                        + tpl_res_after_ref)
                else:
                    segment_start = sequence_number_counter + 1
                    segment_end = tss.end - tss.start

                # update residues for this segment
                create_or_update_residues_in_segment(pconf, segment, segment_start, segment_end, schemes,
                    ref_positions)
                
                sequence_number_counter = segment_end

        self.logger.info('COMPLETED UPDATING PROTEIN ALIGNMENTS')