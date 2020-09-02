from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db import IntegrityError

from build.management.commands.base_build import Command as BaseBuild
from protein.models import Protein, ProteinConformation, ProteinSegment, ProteinFamily
from residue.functions import *

import os
import yaml
from collections import OrderedDict
import copy

class Command(BaseBuild):
    help = 'Creates residue records for human receptors'

    generic_numbers_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'generic_numbers'])
    ref_position_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'reference_positions'])
    auto_ref_position_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'auto_reference_positions'])
    default_segment_length_file_path = os.sep.join([settings.DATA_DIR, 'residue_data', 'default_segment_length.yaml'])

    segments = ProteinSegment.objects.filter(partial=False)
    all_segments = {ps.slug: ps for ps in ProteinSegment.objects.all()}  # all segments dict for faster lookups

    pconfs = ProteinConformation.objects.filter(protein__species__id=1).prefetch_related(
        'protein__residue_numbering_scheme__parent')

    schemes = parse_scheme_tables(generic_numbers_source_dir)

    # default segment length
    with open(default_segment_length_file_path, 'r') as default_segment_length_file:
        segment_length = yaml.load(default_segment_length_file, Loader=yaml.FullLoader)

    def handle(self, *args, **options):
        try:
            self.logger.info('CREATING RESIDUES')

            # run the function twice (second run for proteins without reference positions)
            iterations = 2
            for i in range(1,iterations+1):
                self.prepare_input(options['proc'], self.pconfs, i)

            self.logger.info('COMPLETED CREATING RESIDUES')
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def main_func(self, positions, iteration):
        # pconfs
        if not positions[1]:
            pconfs = self.pconfs[positions[0]:]
        else:
            pconfs = self.pconfs[positions[0]:positions[1]]

        for pconf in pconfs:
            # read reference positions for this protein
            ref_position_file_path = os.sep.join([self.ref_position_source_dir, pconf.protein.entry_name + '.yaml'])
            ref_positions = load_reference_positions(ref_position_file_path)

            # look for automatically generated ref positions if annotations are not found
            if not ref_positions:
                auto_ref_position_file_path = os.sep.join([self.auto_ref_position_source_dir,
                    pconf.protein.entry_name + '.yaml'])
                ref_positions = load_reference_positions(auto_ref_position_file_path)

            # if auto refs are not found, generate them
            if not ref_positions:
                # is this the second iteration of this function? We want all proteins with reference positions to be
                # processed before looking for positions for those who lack them
                if iteration == 2:
                    self.logger.info("Reference positions for {} not annotated, looking for a template".format(
                        pconf.protein))

                    # required information about this protein
                    up = {}
                    up['entry_name'] = pconf.protein.entry_name
                    up['sequence'] = pconf.protein.sequence

                    # find closest protein (by family) to get ref positions
                    # - level 3 parent family
                    # - - level2 parent family
                    # - - - level1 parent family
                    # - - - - current proteins family
                    # - - - - - current protein
                    template_found = False

                    # try level1 families first, then level2, then level3
                    parent_family_levels = [pconf.protein.family.parent, pconf.protein.family.parent.parent,
                        pconf.protein.family.parent.parent.parent]
                    for parent_family in parent_family_levels:
                        if template_found:
                            break

                        # find sub families
                        related_families = ProteinFamily.objects.filter(parent=parent_family)

                        # loop through families and search for proteins to use as template
                        for family in related_families:
                            if template_found:
                                break
                            proteins = Protein.objects.filter(family=family)
                            if not proteins:
                                proteins = Protein.objects.filter(family__parent=family)
                                if not proteins:
                                    proteins = Protein.objects.filter(family__parent__parent=family)
                            for p in proteins:
                                tpl_ref_position_file_path = os.sep.join([self.ref_position_source_dir,
                                    p.entry_name + '.yaml'])
                                tpl_ref_positions = load_reference_positions(tpl_ref_position_file_path)
                                if tpl_ref_positions:
                                    self.logger.info("Found template {}".format(p))
                                    ref_positions = align_protein_to_reference(up, tpl_ref_position_file_path, p)
                                    # write reference positions to a file
                                    with open(auto_ref_position_file_path, "w") as auto_ref_position_file:
                                        yaml.dump(ref_positions, auto_ref_position_file, default_flow_style=False)
                                    template_found = True
                                    break
                    else:
                        if not template_found:
                            self.logger.error('No template reference positions found for {}'.format(pconf.protein))
                else:
                    continue
            elif iteration == 2:
                # proteins with ref positions have already been processed in the first iteration
                continue

            # remote empty ref positions
            ref_positions_copy = copy.deepcopy(ref_positions)
            for position, position_value in ref_positions_copy.items():
                if position_value == '-':
                    del ref_positions[position]

            # determine segment ranges, and create residues
            nseg = self.segments.count()
            sequence_number_counter = 0
            for i, segment in enumerate(self.segments):
                # should this segment be aligned? This value is updated below
                unaligned_segment = True

                # next segment (for checking start positions)
                if (i+1) < nseg:
                    next_segment = self.segments[i+1]
                else:
                    next_segment = False

                # is this an alignable segment?
                if segment.slug in settings.REFERENCE_POSITIONS:
                    # is there a reference position available?
                    if ref_positions and settings.REFERENCE_POSITIONS[segment.slug] in ref_positions:
                        # mark segment as aligned
                        unaligned_segment = False

                        segment_start = (ref_positions[settings.REFERENCE_POSITIONS[segment.slug]]
                            - self.segment_length[segment.slug]['before'])
                        aligned_segment_start = segment_start
                        segment_end = (ref_positions[settings.REFERENCE_POSITIONS[segment.slug]]
                            + self.segment_length[segment.slug]['after'])
                        aligned_segment_end = segment_end

                        # is this segment is not fully aligned, find the start and stop (not just the aligned start
                        # and stop)
                        if not segment.fully_aligned:
                            segment_start = sequence_number_counter + 1
                            if next_segment:
                                next_segment_start = (ref_positions[settings.REFERENCE_POSITIONS[next_segment.slug]]
                                - self.segment_length[next_segment.slug]['before'])
                                segment_end = next_segment_start - 1

                                if (next_segment.slug in settings.REFERENCE_POSITIONS and ref_positions and
                                    settings.REFERENCE_POSITIONS[next_segment.slug] in ref_positions):
                                    next_segment_start = (
                                        ref_positions[settings.REFERENCE_POSITIONS[next_segment.slug]]
                                        - self.segment_length[next_segment.slug]['before'])
                                    segment_end = next_segment_start - 1
                                else:
                                    self.logger.warning('A non-fully aligned segment {} is followed a unaligned' \
                                        + 'segment {}. Skipping.'.format(segment.slug, next_segment.slug))
                                    continue
                            else:
                                segment_end = len(pconf.protein.sequence)

                            if next_segment:
                                if (next_segment.slug in settings.REFERENCE_POSITIONS and ref_positions and
                                    settings.REFERENCE_POSITIONS[next_segment.slug] in ref_positions):
                                    segment_end = (ref_positions[settings.REFERENCE_POSITIONS[next_segment.slug]]
                                    - self.segment_length[next_segment.slug]['before'] - 1)
                                    next_ref_found = True
                                else:
                                    continue
                            else:
                                # for the last segment, the end is the last residue of the sequence
                                segment_end = len(pconf.protein.sequence)

                    else:
                        if segment.fully_aligned:
                            # stop processing this segment
                            self.logger.error('Reference position missing for fully aligned segment {} in {},' \
                                + ' skipping'.format(segment, pconf))
                            continue
                        else:
                            # log the missing reference position
                            self.logger.warning('Reference position missing for segment {} in {}'.format(segment,
                                pconf))

                if unaligned_segment:
                    segment_start = sequence_number_counter + 1

                    # if this is not the last segment, find next segments reference position
                    if next_segment:
                        if (next_segment.slug in settings.REFERENCE_POSITIONS and ref_positions and
                            settings.REFERENCE_POSITIONS[next_segment.slug] in ref_positions):
                            segment_end = (ref_positions[settings.REFERENCE_POSITIONS[next_segment.slug]]
                            - self.segment_length[next_segment.slug]['before'] - 1)
                        else:
                            continue
                    else:
                        # for the last segment, the end is the last residue of the sequence
                        segment_end = len(pconf.protein.sequence)

                    aligned_segment_start = None
                    aligned_segment_end = None

                # skip if the segment ends before it starts (can happen if the next segment is long)
                if segment_start > segment_end:
                    self.logger.warning('Start of segment {} is larger than its end'.format(segment))
                    continue

                # create residues for this segment
                create_or_update_residues_in_segment(pconf, segment, segment_start, aligned_segment_start,
                    segment_end, aligned_segment_end, self.schemes, ref_positions, [], True)

                sequence_number_counter = segment_end
