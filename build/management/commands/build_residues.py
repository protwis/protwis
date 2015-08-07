from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from protein.models import Protein, ProteinConformation, ProteinSegment, ProteinFamily
from residue.functions import *

import logging
import os
import yaml


class Command(BaseCommand):
    help = 'Creates residue records'

    logger = logging.getLogger(__name__)

    generic_numbers_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'generic_numbers'])
    ref_position_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'reference_positions'])
    auto_ref_position_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'auto_reference_positions'])
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
        pcs = ProteinConformation.objects.all().select_related(
            'protein__residue_numbering_scheme__parent')
        
        # run the loop twice, once for annotated proteins, and again for those without annotations (using the closest
        # annotated protein as a reference)
        pclists = [pcs, []]
        for pclist in pclists:
            for pc in pclist:
                sequence_number_counter = 0
                # read reference positions for this protein
                ref_position_file_path = os.sep.join([self.ref_position_source_dir, pc.protein.entry_name + '.yaml'])
                ref_positions = load_reference_positions(ref_position_file_path)

                # look for automatically generated ref positions if annotations are not found
                if not ref_positions:
                    auto_ref_position_file_path = os.sep.join([self.auto_ref_position_source_dir,
                        pc.protein.entry_name + '.yaml'])
                    ref_positions = load_reference_positions(auto_ref_position_file_path)

                # if auto refs are not found, generate them
                if not ref_positions:
                    # is this protein in the "not-annotated list"?
                    if pc in pclists[1]:
                        self.logger.info("Reference positions for {} not annotated, looking for a template".format(
                            pc.protein))
                        
                        # required information about this protein
                        up = {}
                        up['entry_name'] = pc.protein.entry_name
                        up['sequence'] = pc.protein.sequence

                        # find closest protein (by family) to get ref positions
                        # - level 3 parent family
                        # - - level2 parent family
                        # - - - level1 parent family
                        # - - - - current proteins family
                        # - - - - - current protein
                        template_found = False

                        # try level1 families first, then level2, then level3
                        parent_family_levels = [pc.protein.family.parent, pc.protein.family.parent.parent,
                            pc.protein.family.parent.parent.parent]
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
                                self.logger.error('No template reference positions found for {}'.format(pc.protein))
                    else:
                        pclists[1].append(pc)
                        continue

                # determine segment ranges, and create residues
                nseg = len(segments)
                for i, segment in enumerate(segments):
                    # is this an alignable segment?
                    if segment.slug in settings.REFERENCE_POSITIONS:
                        # is there a reference position available?
                        if settings.REFERENCE_POSITIONS[segment.slug] in ref_positions:
                            segment_start = (ref_positions[settings.REFERENCE_POSITIONS[segment.slug]]
                                - segment_length[segment.slug]['before'])
                            segment_end = (ref_positions[settings.REFERENCE_POSITIONS[segment.slug]]
                                + segment_length[segment.slug]['after'])
                        else:
                            # skip this segment if the reference position is missing
                            self.logger.warning('Reference position missing for segment {} in {}'.format(segment, pc))
                            continue
                    else:
                        segment_start = sequence_number_counter + 1
                        
                        # if this is not the last segment, find next segments reference position
                        if (i+1) < nseg:
                            next_segment = segments[i+1]
                            if (next_segment.slug in settings.REFERENCE_POSITIONS and 
                                settings.REFERENCE_POSITIONS[next_segment.slug] in ref_positions):
                                segment_end = (ref_positions[settings.REFERENCE_POSITIONS[next_segment.slug]]
                                - segment_length[next_segment.slug]['before'] - 1)
                                next_ref_found = True
                            else: 
                                continue
                        else:
                            # for the last segment, the end is the last residue of the sequence
                            segment_end = len(pc.protein.sequence)

                        # skip if the segment ends before it starts (can happen if the next segment is long)
                        if segment_start > segment_end:
                            self.logger.warning('Start of segment {} is larger than its end'.format(segment))
                            continue

                    # create residues for this segment
                    create_or_update_residues_in_segment(pc, segment, segment_start, segment_end, schemes,
                        ref_positions, [], True)

                    sequence_number_counter = segment_end

        self.logger.info('COMPLETED CREATING RESIDUES')