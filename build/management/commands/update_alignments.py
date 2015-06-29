from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from protein.models import ProteinConformation, ProteinSegment, ProteinAnomaly, ProteinConformationTemplateStructure
from structure.models import StructureSegment
from residue.models import Residue
from residue.functions import *

import logging
import os
from collections import OrderedDict


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
        pconfs = ProteinConformation.objects.filter(protein__sequence_type__slug='wt').select_related(
            'protein__residue_numbering_scheme__parent', 'template_structure')
        
        for pconf in pconfs:
            sequence_number_counter = 0
            
            # read reference positions for this protein
            ref_position_file_path = os.sep.join([self.ref_position_source_dir, pconf.protein.entry_name + '.yaml'])
            if not os.path.isfile(ref_position_file_path):
                self.logger.error("File {} not found, skipping!".format(ref_position_file_path))
                continue
            
            ref_positions = load_reference_positions(ref_position_file_path)

            # check whether all segments have annotated reference positions
            if len(ref_positions) != len(settings.REFERENCE_POSITIONS):
                self.logger.error('Missing reference positions for {}'.format(pconf))

            # protein anomalies in main template
            main_tpl_pas = pconf.template_structure.protein_anomalies.all().values_list(
                    'generic_number__label', flat=True)

            # dicitonary of updated segment values
            update_segments = []

            # determine segment ranges, and update residues residues
            nseg = len(segments)
            for i, segment in enumerate(segments):
                self.logger.info("Updating segment borders for {} of {}, using template {}".format(segment.slug, pconf,
                    pconf.template_structure))

                # add segment to updated values
                update_segments.append({'segment': segment})

                # protein anomalies to include
                protein_anomalies = []

                # find template segment (for segments borders)
                try:
                    main_tpl_ss = StructureSegment.objects.get(structure=pconf.template_structure,
                        protein_segment=segment)
                except StructureSegment.DoesNotExist:
                    self.logger.warning('Segment records not found for {} in template structure {}, skipping!'.format(
                        segment, pconf.template_structure))
                    continue

                if segment.slug in settings.REFERENCE_POSITIONS:
                    # get reference positions of this segment (e.g. 1x50)
                    segment_ref_position = settings.REFERENCE_POSITIONS[segment.slug]

                    # is there a defined reference position for this protein and segment
                    if segment_ref_position not in ref_positions:
                        self.logger.warning("{} missing definition for {}".format(pconf, segment_ref_position))
                        continue

                    # protein anomaly rules
                    if segment.slug in anomaly_rule_sets:
                        for pa, parss in anomaly_rule_sets[segment.slug].items():
                            # use similarity to decide on anomaly, rules can override this
                            use_similarity = True

                            # go through rule sets
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
                                    
                                    # add the anomaly to the list for this segment (if the rule set is not exclusive)
                                    if not pars.exclusive:
                                        protein_anomalies.append(anomalies[pa])
                                    
                                    # break the set loop, because one set match is enough for a decision
                                    break

                            # use similarity?
                            if use_similarity:
                                # template anomalies
                                tplpas = ProteinConformationTemplateStructure.objects.get(
                                    protein_conformation=pconf, protein_segment=segment)
                                
                                # does the template have the anomaly in question?
                                if pa in tplpas.structure.protein_anomalies.all().values_list(
                                    'generic_number__label', flat=True):
                                    
                                    # add it to the list of anomalies for this segment
                                    protein_anomalies.append(anomalies[pa])
                                    self.logger.info("Anomaly {} included by similarity to {} in {}".format(pa,
                                        tplpas.structure, pconf))
                                else:
                                    self.logger.info("Anomaly {} excluded by similarity to {} in {}".format(pa,
                                        tplpas.structure, pconf))
                            else:
                                if anomalies[pa] in protein_anomalies:
                                    self.logger.info("Anomaly {} included by rule in {}".format(pa, pconf))
                                else:
                                    self.logger.info("Anomaly {} excluded by rule in {}".format(pa, pconf))

                    # template segment reference residue number
                    try:
                        tsrrn = Residue.objects.get(
                            protein_conformation=pconf.template_structure.protein_conformation,
                            generic_number__label=segment_ref_position)
                    except Residue.DoesNotExist:
                        self.logger.warning("Template residues for {} in {} not found, skipping!".format(segment,
                            pconf))
                        continue

                    # number of residues before and after the reference position
                    tpl_res_before_ref = tsrrn.sequence_number - main_tpl_ss.start
                    tpl_res_after_ref = main_tpl_ss.end - tsrrn.sequence_number

                    # FIXME check whether this segments actually needs an update

                    segment_start = ref_positions[segment_ref_position] - tpl_res_before_ref
                    segment_end = ref_positions[segment_ref_position] + tpl_res_after_ref
                    
                    # update start and end positions based on anomalies
                    ref_generic_index = int(segment_ref_position.split("x")[1])
                    for pa in protein_anomalies:
                        # Add bulge to protein_protein_anomalies
                        pconf.protein_anomalies.add(pa)

                        # do change segment borders if  this anomaly is in the main template
                        if pa.generic_number.label in main_tpl_pas:
                            continue
                        
                        # generic number without the prime for bulges
                        pa_generic_index = int(pa.generic_number.label.split("x")[1][:2])
                        if (pa_generic_index > ref_generic_index and pa.anomaly_type.slug == 'bulge'):
                            segment_end += 1
                        elif (pa_generic_index > ref_generic_index and pa.anomaly_type.slug == 'constriction'):
                            segment_end -= 1
                        elif (pa_generic_index < ref_generic_index and pa.anomaly_type.slug == 'bulge'):
                            segment_start -= 1
                        elif (pa_generic_index < ref_generic_index and pa.anomaly_type.slug == 'constriction'):
                            segment_start += 1

                    # update previous segment end
                    update_segments[i-1]['end'] = segment_start - 1
                else:
                    segment_start = sequence_number_counter + 1
                    if i == (nseg-1):
                        segment_end = len(pconf.protein.sequence)
                    else:
                        segment_end = 0 # will be updated in the next iteration
                
                update_segments[i]['start'] = segment_start
                update_segments[i]['end'] = segment_end
                update_segments[i]['protein_anomalies'] = protein_anomalies
                sequence_number_counter = segment_end

            for us in update_segments:
                # update residues for this segment
                if 'end' in us and us['end']: # FIXME take split segments into account
                    create_or_update_residues_in_segment(pconf, us['segment'], us['start'], us['end'], schemes,
                        ref_positions, us['protein_anomalies'])
                

        self.logger.info('COMPLETED UPDATING PROTEIN ALIGNMENTS')