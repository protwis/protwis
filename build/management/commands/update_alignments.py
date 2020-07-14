from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from build.management.commands.base_build import Command as BaseBuild
from protein.models import ProteinConformation, ProteinSegment, ProteinAnomaly, ProteinConformationTemplateStructure
from structure.models import StructureSegment
from residue.models import Residue
from residue.functions import *
from common.alignment import Alignment

import os
from collections import OrderedDict
import copy


class Command(BaseBuild):
    help = 'Updates protein alignments based on structure data'

    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
            type=int,
            action='store',
            dest='proc',
            default=1,
            help='Number of processes to run')

    generic_numbers_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'generic_numbers'])
    ref_position_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'reference_positions'])
    auto_ref_position_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'auto_reference_positions'])
    default_segment_length_file_path = os.sep.join([settings.DATA_DIR, 'residue_data', 'default_segment_length.yaml'])

    # default segment length
    with open(default_segment_length_file_path, 'r') as default_segment_length_file:
        segment_length = yaml.load(default_segment_length_file, Loader=yaml.FullLoader)

    pconfs = ProteinConformation.objects.order_by('protein__parent', 'id').prefetch_related(
            'protein__residue_numbering_scheme__parent', 'protein__genes', 'template_structure')

    def handle(self, *args, **options):
        try:
            self.logger.info('UPDATING PROTEIN ALIGNMENTS')
            self.prepare_input(options['proc'], self.pconfs)
            self.logger.info('COMPLETED UPDATING PROTEIN ALIGNMENTS')
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def main_func(self, positions, iteration):
        # pconfs
        if not positions[1]:
            pconfs = self.pconfs[positions[0]:]
        else:
            pconfs = self.pconfs[positions[0]:positions[1]]

        schemes = parse_scheme_tables(self.generic_numbers_source_dir)

        # pre-fetch protein anomalies and associated rules
        anomaly_rule_sets = {}
        anomalies = {}
        pas = ProteinAnomaly.objects.all().prefetch_related(
            'rulesets__protein_anomaly__generic_number__protein_segment', 'rulesets__rules')
        for pa in pas:
            segment = pa.generic_number.protein_segment
            if segment.slug not in anomaly_rule_sets:
                anomaly_rule_sets[segment.slug] = {}
            anomaly_label = pa.generic_number.label
            anomalies[anomaly_label] = pa
            if anomaly_label not in anomaly_rule_sets[segment.slug]:
                anomaly_rule_sets[segment.slug][anomaly_label] = []
            for pars in pa.rulesets.all():
                anomaly_rule_sets[segment.slug][anomaly_label].append(pars)

        # pre-fetch protein conformations
        segments = ProteinSegment.objects.filter(partial=False)

        for pconf in pconfs:
            # skip protein conformations without a template (consensus sequences)
            if not pconf.template_structure:
                continue
            else:
                template_structure = pconf.template_structure

            # get the sequence number of the first residue for this protein conformation
            pconf_residues = Residue.objects.filter(protein_conformation=pconf)
            if pconf_residues.exists():
                sequence_number_counter = pconf_residues[0].sequence_number - 1
            else:
                sequence_number_counter = 0

            # read reference positions for this protein
            ref_position_file_paths = [
                # canonical ref positions
                os.sep.join([self.ref_position_source_dir, pconf.protein.entry_name + '.yaml']),
                # auto-generated ref positions
                os.sep.join([self.auto_ref_position_source_dir, pconf.protein.entry_name + '.yaml']),
            ]
            if pconf.protein.parent:
                parent_ref_position_file_paths = [
                    # parent ref positions
                    os.sep.join([self.ref_position_source_dir, pconf.protein.parent.entry_name + '.yaml']),
                    # parent auto-generated ref positions
                    os.sep.join([self.auto_ref_position_source_dir, pconf.protein.parent.entry_name + '.yaml']),
                ]
                ref_position_file_paths += parent_ref_position_file_paths

            for file_path in ref_position_file_paths:
                ref_positions = load_reference_positions(file_path)
                if ref_positions:
                    self.logger.info("Reference positions for {} found in {}".format(pconf.protein, file_path))
                    break
            else:
                self.logger.error("No reference positions found for {}, skipping".format(pconf.protein))
                continue

            # remove empty values from reference positions
            ref_positions_copy = copy.deepcopy(ref_positions)
            for position, position_value in ref_positions_copy.items():
                if position_value == '-':
                    del ref_positions[position]

            # protein anomalies in main template
            main_tpl_pas = template_structure.protein_anomalies.all()
            main_tpl_pa_labels = []
            for main_tpl_pa in main_tpl_pas:
                main_tpl_pa_labels.append(main_tpl_pa.generic_number.label)

            # dictionary of updated segment values
            update_segments = []

            # determine segment ranges, and update residues residues
            nseg = len(segments)
            for i, segment in enumerate(segments):
                self.logger.info("Updating segment borders for {} of {}, using template {}".format(segment.slug, pconf,
                    template_structure))

                # segment template structure (only differs from the main template if segment is missing in main
                # structure, and segment is not fully aligned)
                segment_template_structure = template_structure

                # should this segment be aligned? This value is updated below
                unaligned_segment = True

                # add segment to updated values
                update_segments.append({'segment': segment})

                # protein anomalies to include
                protein_anomalies = []

                if segment.slug in settings.REFERENCE_POSITIONS:
                    # are all template requirements satisfied? Updated below
                    templates_found = True

                    # find template segment (for segments borders)
                    try:
                        main_tpl_ss = StructureSegment.objects.get(structure=segment_template_structure,
                            protein_segment=segment)
                    except StructureSegment.DoesNotExist:
                        self.logger.info('Segment records not found for {} in template structure {}, looking for alternatives'.format(
                            segment, segment_template_structure))
                        if not segment.fully_aligned:
                            segment_tpl = ProteinConformationTemplateStructure.objects.get(protein_conformation=pconf,
                                protein_segment=segment)
                            try:
                                main_tpl_ss = StructureSegment.objects.get(structure=segment_tpl.structure,
                                    protein_segment=segment)
                                segment_template_structure = segment_tpl.structure
                                self.logger.info('Using structure {} for {} in {}'.format(segment_tpl.structure,
                                    segment, pconf))
                            except:
                                templates_found = False
                                self.logger.warning('No template found for {} in {}, skipping'.format(segment, pconf))
                        else:
                            templates_found = False
                            self.logger.warning('No template found for {} in {}, skipping'.format(segment, pconf))

                    # get reference positions of this segment (e.g. 1x50)
                    segment_ref_position = settings.REFERENCE_POSITIONS[segment.slug]

                    # is there a defined reference position for this protein and segment
                    if segment_ref_position not in ref_positions:
                        self.logger.warning("{} missing definition for {}".format(pconf, segment_ref_position))
                        templates_found = False

                    # template segment reference residue number
                    try:
                        tsrrn = Residue.objects.get(
                            protein_conformation=segment_template_structure.protein_conformation,
                            generic_number__label=segment_ref_position)
                    except Residue.DoesNotExist:
                        self.logger.info("Template residues for {} in {} not found, looking for alternatives!".format(
                            segment, pconf))
                        if not segment.fully_aligned:
                            segment_tpl = ProteinConformationTemplateStructure.objects.get(protein_conformation=pconf,
                                protein_segment=segment)
                            try:
                                tsrrn = Residue.objects.get(
                                    protein_conformation=segment_tpl.structure.protein_conformation,
                                    generic_number__label=segment_ref_position)
                                main_tpl_ss = StructureSegment.objects.get(structure=segment_tpl.structure,
                                    protein_segment=segment)
                                self.logger.info('Using residues from {} for {} in {}'.format(segment_tpl.structure,
                                    segment, pconf))
                            except:
                                self.logger.warning("No template residues for {} in {} not found, skipping!".format(
                                    segment, pconf))
                                templates_found = False
                        else:
                            self.logger.warning("No template residues for {} in {} not found, skipping!".format(
                                segment, pconf))
                            templates_found = False

                    # if template requirements are fulfilled, treat as an aligned segment
                    if templates_found:
                        unaligned_segment = False
                    # if not, stop here for fully aligned segments, otherwise treat as unaligned
                    else:
                        if segment.fully_aligned:
                            continue

                if not unaligned_segment:
                    # number of residues before and after the reference position
                    tpl_res_before_ref = tsrrn.sequence_number - main_tpl_ss.start
                    tpl_res_after_ref = main_tpl_ss.end - tsrrn.sequence_number
                    aligned_segment_start = ref_positions[segment_ref_position] - tpl_res_before_ref
                    aligned_segment_end = ref_positions[segment_ref_position] + tpl_res_after_ref

                    # protein anomaly rules
                    if segment.slug in anomaly_rule_sets:
                        # if there exists a structure for this particular protein, don't use the rules
                        ignore_rules = False

                        if pconf.protein.parent:
                            # use parent protein for constructs and other non wild-type sequences
                            current_protein = pconf.protein.parent
                        else:
                            current_protein = pconf.protein
                        current_protein_genes = current_protein.genes.order_by('position')
                        if current_protein_genes.count():
                            current_gene = current_protein_genes[0]
                            template_structure_protein = segment_template_structure.protein_conformation.protein.parent
                            template_structure_gene = template_structure_protein.genes.order_by('position')[0]
                            if current_gene.name.lower() == template_structure_gene.name.lower():
                                ignore_rules = True
                                self.logger.info('Ignoring anomaly rules because of {} structure'
                                    .format(template_structure_protein))

                        # get a list of generic numbers in main template
                        main_tpl_gn_labels = Residue.objects.filter(
                            protein_conformation=segment_template_structure.protein_conformation, generic_number__isnull=False,
                            protein_segment=segment).values_list('generic_number__label', flat=True)

                        for pa, parss in anomaly_rule_sets[segment.slug].items():
                            # check whether this anomaly is inside the segment borders
                            numbers_within_segment = generic_number_within_segment_borders(pa, main_tpl_gn_labels)
                            if not numbers_within_segment:
                                self.logger.info("Anomaly {} excluded for {} (outside segment borders)".format(pa,
                                    pconf))
                                continue

                            # use similarity to decide on anomaly, rules can override this
                            use_similarity = True

                            # go through rule sets
                            if not ignore_rules:
                                for pars in parss:
                                    # fetch rules in this rule set
                                    rules = pars.rules.all()
                                    for rule in rules:
                                        # fetch the residue in question
                                        try:
                                            r = Residue.objects.get(protein_conformation=pconf,
                                                generic_number=rule.generic_number)
                                        except Residue.DoesNotExist:
                                            self.logger.warning('Residue {} in {} not found, skipping'.format(
                                                rule.generic_number.label, pconf.protein.entry_name))
                                            continue

                                        # does the rule break the set? Then go to next set..
                                        if ((r.amino_acid == rule.amino_acid and rule.negative) or
                                            (r.amino_acid != rule.amino_acid and not rule.negative)):
                                            break
                                    # if the loop was not broken, all rules passed, and the set controls the anomaly
                                    else:
                                        # do not use similarity, since a rule matched
                                        use_similarity = False

                                        # add the anomaly to the list for this segment (if the rule set is not
                                        # exclusive)
                                        if not pars.exclusive:
                                            protein_anomalies.append(anomalies[pa])

                                        # break the set loop, because one set match is enough for a decision
                                        break

                            # use similarity?
                            if use_similarity:
                                # does the template have the anomaly in question?
                                if pa in segment_template_structure.protein_anomalies.all().values_list(
                                    'generic_number__label', flat=True):

                                    # add it to the list of anomalies for this segment
                                    protein_anomalies.append(anomalies[pa])
                                    self.logger.info("Anomaly {} included for {} (similarity to {})".format(pa,
                                        pconf, segment_template_structure))
                                else:
                                    self.logger.info("Anomaly {} excluded for {} (similarity to {})".format(pa,
                                        pconf, segment_template_structure))
                            else:
                                if anomalies[pa] in protein_anomalies:
                                    self.logger.info("Anomaly {} included for {} (rule)".format(pa, pconf))
                                else:
                                    self.logger.info("Anomaly {} excluded for {} (rule)".format(pa, pconf))

                    # update start and end positions based on anomalies in this protein
                    pa_labels = []
                    ref_generic_index = int(segment_ref_position.split("x")[1])
                    for pa in protein_anomalies:
                        # does the anomaly belong to this segment?
                        if pa.generic_number.protein_segment != segment:
                            continue

                        # Add bulge to protein_protein_anomalies
                        pconf.protein_anomalies.add(pa)

                        # add to list of anomaly labels (to compare with template below)
                        pa_labels.append(pa.generic_number.label)

                        # do change segment borders if this anomaly is in the main template
                        if pa.generic_number.label in main_tpl_pa_labels:
                            continue

                        # generic number without the prime for bulges
                        pa_generic_index = int(pa.generic_number.label.split("x")[1][:2])
                        if (pa_generic_index > ref_generic_index and pa.anomaly_type.slug == 'bulge'):
                            aligned_segment_end += 1
                        elif (pa_generic_index > ref_generic_index and pa.anomaly_type.slug == 'constriction'):
                            aligned_segment_end -= 1
                        elif (pa_generic_index < ref_generic_index and pa.anomaly_type.slug == 'bulge'):
                            aligned_segment_start -= 1
                        elif (pa_generic_index < ref_generic_index and pa.anomaly_type.slug == 'constriction'):
                            aligned_segment_start += 1
                    # update start and end positions based on anomalies in the template
                    for pa in main_tpl_pas:
                        # does the anomaly belong to this segment?
                        if pa.generic_number.protein_segment != segment:
                            continue

                        # do change segment borders if this anomaly is in the current protein
                        if pa.generic_number.label in pa_labels:
                            continue

                        # generic number without the prime for bulges
                        pa_generic_index = int(pa.generic_number.label.split("x")[1][:2])
                        if (pa_generic_index > ref_generic_index and pa.anomaly_type.slug == 'bulge'):
                            aligned_segment_end -= 1
                        elif (pa_generic_index > ref_generic_index and pa.anomaly_type.slug == 'constriction'):
                            aligned_segment_end += 1
                        elif (pa_generic_index < ref_generic_index and pa.anomaly_type.slug == 'bulge'):
                            aligned_segment_start += 1
                        elif (pa_generic_index < ref_generic_index and pa.anomaly_type.slug == 'constriction'):
                            aligned_segment_start -= 1

                    # set start and end positions (not just the aligned)
                    if segment.fully_aligned:
                        segment_start = aligned_segment_start
                        segment_end = aligned_segment_end
                    else:
                        segment_start = sequence_number_counter + 1
                        segment_end = 0 # will be updated in the next iteration
                else:
                    segment_start = sequence_number_counter + 1
                    if i == (nseg-1):
                        segment_end = len(pconf.protein.sequence)
                    else:
                        segment_end = 0 # will be updated in the next iteration

                    aligned_segment_start = None
                    aligned_segment_end = None

                update_segments[i]['start'] = segment_start
                update_segments[i]['aligned_start'] = aligned_segment_start
                update_segments[i]['end'] = segment_end
                update_segments[i]['aligned_end'] = aligned_segment_end
                update_segments[i]['protein_anomalies'] = protein_anomalies
                if segment_end:
                    sequence_number_counter = segment_end

                # update previous segment end if needed
                if (i and 'end' in update_segments[i-1] and (not update_segments[i-1]['end']
                    or update_segments[i-1]['end'] != update_segments[i-1]['aligned_end'])):
                    update_segments[i-1]['end'] = segment_start - 1

                # check whether minimum length is fulfilled for last segment
                if (i and update_segments[i-1]['segment'].slug in self.segment_length
                    and 'min' in self.segment_length[update_segments[i-1]['segment'].slug]
                    and 'end' in update_segments[i-1]):
                    # if it is not, find out how many residues are missing
                    last_min_segment_length = self.segment_length[update_segments[i-1]['segment'].slug]['min']
                    # +1 because a segment starting at 6 and ending at 10 is 5 positions, but 10-6 is 4
                    last_segment_length = update_segments[i-1]['end'] - update_segments[i-1]['start'] + 1
                    if last_segment_length < 0:
                        last_segment_length = 0
                    if last_segment_length < last_min_segment_length:
                        missing_residues = last_min_segment_length - last_segment_length

                        # take the missing residues from the segments before and after
                        add_residues_before = round(missing_residues / 2) + missing_residues % 2
                        add_residues_after = round(missing_residues / 2)
                        update_segments[i-1]['start'] -= add_residues_before
                        update_segments[i-1]['end'] = update_segments[i-1]['start'] + last_min_segment_length - 1

                        # update aligned start and end if they exceed the updated start and stop values
                        if (update_segments[i-1]['aligned_start']
                            and update_segments[i-1]['aligned_start'] < update_segments[i-1]['start']):
                            update_segments[i-1]['aligned_start'] = update_segments[i-1]['start']
                        if (update_segments[i-1]['aligned_end']
                            and update_segments[i-1]['aligned_end'] > update_segments[i-1]['end']):
                            update_segments[i-1]['aligned_end'] = update_segments[i-1]['end']

                        # update this segment's start
                        update_segments[i]['start'] = update_segments[i-1]['end'] + 1
                        if update_segments[i]['aligned_start'] < update_segments[i]['start']:
                            update_segments[i]['aligned_start'] = update_segments[i]['start']

            for us in update_segments:
                if 'start' in us and 'end' in us and us['end']:
                    create_or_update_residues_in_segment(pconf, us['segment'], us['start'], us['aligned_start'],
                        us['end'], us['aligned_end'], schemes, ref_positions, us['protein_anomalies'], False)
