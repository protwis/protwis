from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.utils.text import slugify

from build.management.commands.build_human_proteins import Command as BuildHumanProteins
from residue.functions import *
from protein.models import Protein, ProteinConformation, ProteinFamily, ProteinSegment, ProteinSequenceType
from common.alignment import Alignment
from common.tools import test_model_updates
from alignment.models import AlignmentConsensus

import os
import yaml
import pickle
import django.apps
from collections import OrderedDict

class Command(BuildHumanProteins):
    help = 'Builds consensus sequences for human proteins in all families'

    ref_position_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'reference_positions'])
    auto_ref_position_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'auto_reference_positions'])
    generic_numbers_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'generic_numbers'])

    # numbering scheme tables
    schemes = parse_scheme_tables(generic_numbers_source_dir)

    # get all segments
    segments = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR')

    # fetch families
    families = ProteinFamily.objects.filter(slug__startswith='00').all()

    #Setting the variables for the test tracking of the model upadates
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc', type=int, action='store', dest='proc', default=1, help='Number of processes to run')
        parser.add_argument('--signprot', type=str, action='store', dest='signprot', default=False, help='Only run for either G proteins or arrestins')
        parser.add_argument('--input-slug', type=str, action='store', dest='input-slug', default=False, help='Run only on a family slug from ProteinFamily table')
        parser.add_argument('--purge', action='store_true', dest='purge', default=False, help='Purge all consensus data')

    def handle(self, *args, **options):
        try:
            self.signprot = options['signprot']
            self.input_slug = options['input-slug']
            if options['purge']:
                self.purge_consensus_sequences()
                self.tracker = {}
                test_model_updates(self.all_models, self.tracker, initialize=True)
            self.logger.info('CREATING CONSENSUS SEQUENCES')
            self.prepare_input(options['proc'], self.families)
            test_model_updates(self.all_models, self.tracker, check=True)
            self.logger.info('COMPLETED CREATING CONSENSUS SEQUENCES')
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def purge_consensus_sequences(self):
        Protein.objects.filter(sequence_type__slug='consensus').delete()
        AlignmentConsensus.objects.all().delete()

    def get_segment_residue_information(self, consensus_sequence):
        ref_positions = dict()
        segment_starts = dict()
        segment_aligned_starts = dict()
        segment_ends = dict()
        segment_aligned_ends = dict()
        sequence_num = 1
        unaligned_prefixes = ['00', '01', 'zz']
        for segment_slug, s in consensus_sequence.items():
            if self.signprot and self.signprot in ['Alpha','Arrestin']:
                segment = ProteinSegment.objects.get(slug=segment_slug, proteinfamily=self.signprot)
            elif self.signprot and self.signprot not in ['Alpha','Arrestin']:
                raise AssertionError('{} not supported as signprot type'.format(self.signprot))
            else:
                segment = ProteinSegment.objects.get(slug=segment_slug)
            i = 1
            for gn, aa in s.items():
                if segment_slug in settings.REFERENCE_POSITIONS and gn[-2:] == '50':
                    ref_positions[gn] = sequence_num
                if i == 1:
                    segment_starts[segment_slug] = sequence_num
                if i == len(s):
                    segment_ends[segment_slug] = sequence_num

                # aligned start and end
                if gn[:2] not in unaligned_prefixes:
                    if segment_slug not in segment_aligned_starts:
                        segment_aligned_starts[segment_slug] = sequence_num
                    segment_aligned_ends[segment_slug] = sequence_num

                sequence_num += 1
                i += 1
        for segment_slug in consensus_sequence:
            if segment_slug not in segment_aligned_starts:
                segment_aligned_starts[segment_slug] = None
            if segment_slug not in segment_aligned_ends:
                segment_aligned_ends[segment_slug] = None

        return ref_positions, segment_starts, segment_aligned_starts, segment_ends, segment_aligned_ends

    def main_func(self, positions, iteration,count,lock):
        # families
        # if not positions[1]:
        #     families = self.families[positions[0]:]
        # else:
        #     families = self.families[positions[0]:positions[1]]
        if self.signprot:
            signprot_fam = ProteinFamily.objects.get(name=self.signprot)
            families = ProteinFamily.objects.filter(slug__startswith=signprot_fam.slug+'_').all() # The '_' at the end is needed to skip the Alpha and Arrestin consensus sequences
            self.segments = ProteinSegment.objects.filter(partial=False, proteinfamily=self.signprot)
        else:
            families = self.families

        if self.input_slug:
            families = ProteinFamily.objects.filter(slug__startswith=self.input_slug)


        while count.value<len(families):
            with lock:
                # Check if the count changed before having obtained the lock
                if count.value < len(families):
                    family = families[count.value]
                    count.value +=1
                else:
                    continue

        # for family in families:
            # get proteins in this family
            proteins = Protein.objects.filter(family__slug__startswith=family.slug, sequence_type__slug='wt',
                species__common_name="Human").prefetch_related('species', 'residue_numbering_scheme')

            # if family does not have human equivalents, like Class D1
            if len(proteins)==0:
                proteins = Protein.objects.filter(family__slug__startswith=family.slug, sequence_type__slug='wt',).prefetch_related('species', 'residue_numbering_scheme')

            if proteins.count() <= 1:
                continue
            self.logger.info('Building alignment for {}'.format(family))
            # create alignment
            a = Alignment()
            a.load_proteins(proteins)
            a.load_segments(self.segments)
            a.build_alignment()
            a.calculate_statistics()

            try:
                # Save alignment
                AlignmentConsensus.objects.create(slug=family.slug, alignment=pickle.dumps(a))

                # Load alignment to ensure it works
                a = pickle.loads(AlignmentConsensus.objects.get(slug=family.slug).alignment)
                self.logger.info('Succesfully pickled {}'.format(family))
            except:
                self.logger.error('Failed pickle for {}'.format(family))

            self.logger.info('Completed building alignment for {}'.format(family))

            # get (forced) consensus sequence from alignment object
            family_consensus = str()
            for segment, s in a.forced_consensus.items():
                for gn, aa in s.items():
                    family_consensus += aa

            # create sequence type 'consensus'
            sequence_type, created = ProteinSequenceType.objects.get_or_create(slug='consensus',
                defaults={'name': 'Consensus',})
            if created:
                self.logger.info('Created protein sequence type {}'.format(sequence_type.name))

            # create a protein record
            consensus_name = family.name + " consensus"
            residue_numbering_scheme = proteins[0].residue_numbering_scheme
            up = dict()
            up['entry_name'] = slugify(consensus_name)
            if Protein.objects.filter(entry_name=up['entry_name']).exists():
                up['entry_name'] += "-" + family.slug.split('_')[0]
            up['source'] = "OTHER"
            up['species_latin_name'] = proteins[0].species.latin_name
            up['species_common_name'] = proteins[0].species.common_name
            up['sequence'] = family_consensus

            up['names'] = up['genes'] = []
            self.create_protein(consensus_name, family, sequence_type, residue_numbering_scheme, False, up)

            # get protein anomalies in family
            all_constrictions = []
            constriction_freq = dict()
            consensus_pas = dict() # a constriction has to be in all sequences to be included in the consensus
            pcs = ProteinConformation.objects.filter(protein__in=proteins,
                state__slug=settings.DEFAULT_PROTEIN_STATE).prefetch_related('protein_anomalies')
            for pc in pcs:
                pas = pc.protein_anomalies.all().prefetch_related('generic_number__protein_segment', 'anomaly_type')
                for pa in pas:
                    pa_label = pa.generic_number.label
                    pa_type = pa.anomaly_type.slug
                    pa_segment_slug = pa.generic_number.protein_segment.slug

                    # bulges are directly added to the consensus list
                    if pa_type == 'bulge':
                        if pa_segment_slug not in consensus_pas:
                            consensus_pas[pa_segment_slug] = []
                        if pa not in consensus_pas[pa_segment_slug]:
                            consensus_pas[pa_segment_slug].append(pa)

                    # a constriction's frequency is counted
                    else:
                        if pa not in all_constrictions:
                            all_constrictions.append(pa)
                        if pa_label in constriction_freq:
                            constriction_freq[pa_label] += 1
                        else:
                            constriction_freq[pa_label] = 1

            # go through constrictions to see which ones should be included in the consensus
            for pa in all_constrictions:
                pa_label = pa.generic_number.label
                pa_segment_slug = pa.generic_number.protein_segment.slug
                freq = constriction_freq[pa_label]

                # is the constriction in all sequences?
                if freq == len(all_constrictions):
                    if pa_segment_slug not in consensus_pas:
                        consensus_pas[pa_segment_slug] = []
                    consensus_pas[pa_segment_slug].append(pa)

            # create residues
            pc = ProteinConformation.objects.get(protein__entry_name=up['entry_name'],
                state__slug=settings.DEFAULT_PROTEIN_STATE)
            segment_info = self.get_segment_residue_information(a.forced_consensus)
            ref_positions, segment_starts, segment_aligned_starts, segment_ends, segment_aligned_ends = segment_info
            for segment_slug, s in a.forced_consensus.items():
                if self.signprot:
                    segment = ProteinSegment.objects.get(slug=segment_slug, proteinfamily=self.signprot)
                else:
                    segment = ProteinSegment.objects.get(slug=segment_slug)
                if segment_slug in consensus_pas:
                    protein_anomalies = consensus_pas[segment_slug]
                else:
                    protein_anomalies = []
                if segment_slug in segment_starts:
                    if self.signprot:
                        create_or_update_residues_in_segment(pc, segment, segment_starts[segment_slug],
                            segment_aligned_starts[segment_slug], segment_ends[segment_slug],
                            segment_aligned_ends[segment_slug], self.schemes, ref_positions, protein_anomalies, True, self.signprot)
                    else:
                        create_or_update_residues_in_segment(pc, segment, segment_starts[segment_slug],
                            segment_aligned_starts[segment_slug], segment_ends[segment_slug],
                            segment_aligned_ends[segment_slug], self.schemes, ref_positions, protein_anomalies, True)
