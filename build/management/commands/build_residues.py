from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from protein.models import Protein, ProteinConformation, ProteinSegment
from residue.functions import *

import logging
import os
import yaml


class Command(BaseCommand):
    help = 'Creates residue records'

    logger = logging.getLogger(__name__)

    generic_numbers_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'generic_numbers'])
    ref_position_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'reference_positions'])
    default_segment_length_file_path = os.sep.join([settings.DATA_DIR, 'residue_data', 'default_segment_length.yaml'])

    def handle(self, *args, **options):
        # create residue records for all proteins
        try:
            self.create_residues(args)
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def create_residues(self, args):
        self.logger.info('CREATING RESIDUES')

        schemes = parse_scheme_tables(self.generic_numbers_source_dir)

        # default segment length
        with open(self.default_segment_length_file_path, 'r') as default_segment_length_file:
            segment_length = yaml.load(default_segment_length_file)

        # fetch protein conformations
        segments = ProteinSegment.objects.filter(partial=False)
        pcs = ProteinConformation.objects.filter(protein__sequence_type__slug='wt').select_related(
            'protein__residue_numbering_scheme__parent')
        for pc in pcs:
            sequence_number_counter = 0
            # read reference positions for this protein
            ref_position_file_path = os.sep.join([self.ref_position_source_dir, pc.protein.entry_name + '.yaml'])
            ref_positions = load_reference_positions(ref_position_file_path)

            # check whether all segments have annotated reference positions
            if len(ref_positions) == len(settings.REFERENCE_POSITIONS):
                assign_generic_numbers = True
            else:
                self.logger.warning('Missing reference positions for {}'.format(pc.protein))
                assign_generic_numbers = False

            # determine segment ranges, and create residues
            nseg = len(segments)
            for i, segment in enumerate(segments):
                if (segment.slug in settings.REFERENCE_POSITIONS and
                    settings.REFERENCE_POSITIONS[segment.slug] in ref_positions):
                    segment_start = (ref_positions[settings.REFERENCE_POSITIONS[segment.slug]]
                        - segment_length[segment.slug]['before'])
                    segment_end = (ref_positions[settings.REFERENCE_POSITIONS[segment.slug]]
                        + segment_length[segment.slug]['after'])
                else:
                    segment_start = sequence_number_counter + 1
                    
                    if (i+1) < nseg and assign_generic_numbers:
                        next_segment = segments[i+1]
                        if next_segment.slug in settings.REFERENCE_POSITIONS:
                            segment_end = (ref_positions[settings.REFERENCE_POSITIONS[next_segment.slug]]
                            - segment_length[next_segment.slug]['before'] - 1)
                        else: 
                            raise Exception('Not enough reference data to build alignment, aborting!')
                    else:
                        segment_end = len(pc.protein.sequence)

                # create residues for this segment
                create_or_update_residues_in_segment(pc, segment, segment_start, segment_end, schemes,
                    ref_positions, [])

                sequence_number_counter = segment_end

        self.logger.info('COMPLETED CREATING RESIDUES')