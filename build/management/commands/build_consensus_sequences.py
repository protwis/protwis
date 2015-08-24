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

    def handle(self, *args, **options):
        # create proteins
        try:
            self.create_consensus_sequences()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def create_consensus_sequences(self):
        self.logger.info('CREATING CONSENSUS SEQUENCES')

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
            a.calculate_statistics

            # FIXME use this sequence
            # get (forced) consensus sequence from alignment object
            # family_consensus = str()
            # for segment in a.forced_consensus:
            #     for aa in segment:
            #         family_consensus += aa[1]

            # create sequence type 'consensus'
            sequence_type, created = ProteinSequenceType.objects.get_or_create(slug='consensus',
                defaults={
                'name': 'Consensus',
                })
            if created:
                self.logger.info('Created protein sequence type {}'.format(sequence_type.name))

            # create a protein record
            consensus_name = family.name + " consensus"
            residue_numbering_scheme = proteins[0].residue_numbering_scheme
            up = dict()
            up['entry_name'] = slugify(consensus_name)
            up['source'] = "OTHER"
            up['species_latin_name'] = proteins[0].species.latin_name
            up['species_common_name'] = proteins[0].species.common_name
            up['sequence'] = proteins[0].sequence

            up['names'] = up['genes'] = []
            self.create_protein(consensus_name, family, sequence_type, residue_numbering_scheme, False, up)

            # create residues FIXME use actual consensus
            pc = ProteinConformation.objects.get(protein__entry_name=up['entry_name'],
                state__slug=settings.DEFAULT_PROTEIN_STATE)
            ppc = ProteinConformation.objects.get(protein=proteins[0], state__slug=settings.DEFAULT_PROTEIN_STATE)
            prs = Residue.objects.filter(protein_conformation=ppc).prefetch_related(
                        'protein_conformation__protein', 'protein_segment', 'generic_number',
                        'display_generic_number__scheme', 'alternative_generic_numbers__scheme')
            for pr in prs:
                r = Residue()
                r.protein_conformation = pc
                r.generic_number = pr.generic_number
                r.display_generic_number = pr.display_generic_number
                r.sequence_number = pr.sequence_number
                r.protein_segment = pr.protein_segment
                r.amino_acid = pr.amino_acid

                # save residue before populating M2M relations
                r.save()

                # alternative generic numbers
                agns = pr.alternative_generic_numbers.all()
                for agn in agns:
                    r.alternative_generic_numbers.add(agn)

        self.logger.info('COMPLETED CREATING CONSENSUS SEQUENCES')