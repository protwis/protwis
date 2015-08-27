from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.utils.text import slugify

from build.management.commands.build_human_proteins import Command as BuildHumanProteins
from residue.functions import *
from protein.models import Protein, ProteinConformation, ProteinFamily, ProteinSegment, ProteinSequenceType
from common.alignment import Alignment

import logging
import os
import yaml


class Command(BuildHumanProteins):
    help = 'Builds consensus sequences for human proteins in all families'

    ref_position_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'reference_positions'])
    auto_ref_position_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'auto_reference_positions'])
    generic_numbers_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'generic_numbers'])

    def handle(self, *args, **options):
        # create proteins
        try:
            self.purge_consensus_sequences()
            self.create_consensus_sequences()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def purge_consensus_sequences(self):
        Protein.objects.filter(sequence_type__slug='consensus').delete()

    def get_segment_residue_information(self, consensus_sequence):
        ref_positions = dict()
        segment_starts = dict()
        segment_ends = dict()
        sequence_num = 1
        for segment_slug, s in consensus_sequence.items():
            segment = ProteinSegment.objects.get(slug=segment_slug)
            i = 1
            for gn, aa in s.items():
                if segment_slug in settings.REFERENCE_POSITIONS and gn[-2:] == '50':
                    ref_positions[gn] = sequence_num
                if i == 1:
                    segment_starts[segment_slug] = sequence_num
                if i == len(s):
                    segment_ends[segment_slug] = sequence_num
                sequence_num += 1
                i += 1
        return ref_positions, segment_starts, segment_ends

    def create_consensus_sequences(self):
        self.logger.info('CREATING CONSENSUS SEQUENCES')

        # numbering scheme tables
        schemes = parse_scheme_tables(self.generic_numbers_source_dir)

        # fetch families
        families = ProteinFamily.objects.all()
        for family in families:
            # get proteins in this family
            proteins = Protein.objects.filter(family__slug__startswith=family.slug, sequence_type__slug='wt',
                species__id=1)

            if len(proteins) <= 1:
                continue

            # get all segments
            segments = ProteinSegment.objects.filter(partial=False)

            # create alignment
            a = Alignment()
            a.load_proteins(proteins)
            a.load_segments(segments)
            a.build_alignment()
            a.calculate_statistics()

            # get (forced) consensus sequence from alignment object
            family_consensus = str()
            for segment, s in a.forced_consensus.items():
                for gn, aa in s.items():
                    family_consensus += aa

            # create sequence type 'consensus'
            sequence_type, created = ProteinSequenceType.objects.get_or_create(slug='consensus',
                defaults={
                'name': 'Consensus',
                })
            if created:
                self.logger.info('Created protein sequence type {}'.format(sequence_type.name))

            # create a protein record
            consensus_name = family.name + " consensus"
            if Protein.objects.filter(entry_name=consensus_name).exists():
                consensus_name += " " + family.slug.split('_')[0]
            residue_numbering_scheme = proteins[0].residue_numbering_scheme
            up = dict()
            up['entry_name'] = slugify(consensus_name)
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
            pcs = ProteinConformation.objects.filter(protein__in=proteins, state__slug=settings.DEFAULT_PROTEIN_STATE)
            for pc in pcs:
                pas = pc.protein_anomalies.all()
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
            ref_positions, segment_starts, segment_ends = self.get_segment_residue_information(a.forced_consensus)
            for segment_slug, s in a.forced_consensus.items():
                segment = ProteinSegment.objects.get(slug=segment_slug)
                if segment_slug in consensus_pas:
                    protein_anomalies = consensus_pas[segment_slug]
                else:
                    protein_anomalies = []
                if segment_slug in segment_starts:
                    create_or_update_residues_in_segment(pc, segment, segment_starts[segment_slug],
                        segment_ends[segment_slug], schemes, ref_positions, protein_anomalies, True)

        self.logger.info('COMPLETED CREATING CONSENSUS SEQUENCES')