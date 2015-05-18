from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from protein.models import ProteinConformation, ProteinSegment, ProteinAnomaly
from structure.models import StructureSegment
from residue.models import Residue
from residue.functions import *

import logging
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
        anomalies = {}
        pas = ProteinAnomaly.objects.all().prefetch_related(
            'rulesets__protein_anomaly__generic_number__protein_segment', 'rulesets__rules')
        for pa in pas:
            segment = pa.generic_number.protein_segment
            if segment.slug not in anomalies:
                anomalies[segment.slug] = {}
            anomaly_label = pa.generic_number.label
            if anomaly_label not in anomalies[segment.slug]:
                anomalies[segment.slug][anomaly_label] = []
            for pars in pa.rulesets.all():
                anomalies[segment.slug][anomaly_label].append(pars)

        # pref-fetch protein conformations
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

                # protein anomalies to include
                protein_anomalies = []

                # find template segment (for segments borders)
                try:
                    main_tpl_ss = StructureSegment.objects.get(structure=pconf.template_structure, protein_segment=segment)
                except StructureSegment.DoesNotExist:
                    self.logger.error('Segment records not found for template structure {}, skipping!'.format(
                        pconf.template_structure))
                    continue

                # find template segments (for anomalies)
                # tpl_ss = .. FIXME write

                if segment.slug in settings.REFERENCE_POSITIONS:
                    # protein anomaly rules
                    if segment.slug in anomalies:
                        for pa, parss in anomalies[segment.slug].items():
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
                                        self.logger.error('Residue {} in {} not found, skipping'.format(
                                            rule.generic_number.label, pconf.protein.entry_name))
                                        rule_set_met
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
                                        protein_anomalies.append(pa)
                                    
                                    # break the set loop, because one set match is enough for a decision
                                    break

                            # use similarity?
                            if use_similarity:
                                # template anomalies
                                # tplpas = .. # FIXME write
                                print(pconf, pa, "use similarity")
                            else:
                                if pa in protein_anomalies:
                                    print(pconf, pa, "Included by rule")
                                else:
                                    print(pconf, pa, "Excluded by rule")
                    continue # FIXME temp

                    # get reference positions of this segment (e.g. 1x50)
                    segment_ref_position = settings.REFERENCE_POSITIONS[segment.slug]

                    # template segment reference residue number
                    tsrrn = Residue.objects.get(protein_conformation=pconf.template_structure.protein_conformation,
                        generic_number__label=segment_ref_position)

                    # number of residues before and after the reference position
                    tpl_res_before_ref = tsrrn.sequence_number - main_tpl_ss.start
                    tpl_res_after_ref = main_tpl_ss.end - tsrrn.sequence_number

                    # FIXME check whether this segments actually needs an update

                    segment_start = (ref_positions[segment_ref_position]
                        - tpl_res_before_ref)
                    segment_end = (ref_positions[segment_ref_position]
                        + tpl_res_after_ref)
                else:
                    segment_start = sequence_number_counter + 1
                    segment_end = main_tpl_ss.end - main_tpl_ss.start

                # update residues for this segment
                create_or_update_residues_in_segment(pconf, segment, segment_start, segment_end, schemes,
                    ref_positions)
                
                sequence_number_counter = segment_end

        self.logger.info('COMPLETED UPDATING PROTEIN ALIGNMENTS')