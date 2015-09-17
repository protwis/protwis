from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.utils.text import slugify
from django.utils.html import strip_tags
from django.db import IntegrityError

from build.management.commands.base_build import Command as BaseBuild
from protein.models import (Protein, ProteinConformation, ProteinState, ProteinSequenceType, ProteinSegment,
ProteinFusion, ProteinFusionProtein, ProteinSource)
from residue.models import Residue

from optparse import make_option
from datetime import datetime
import logging, os
import yaml


class Command(BaseBuild):
    help = 'Reads source data and creates protein records for constructs'
    
    def add_arguments(self, parser):
        parser.add_argument('--njobs', action='store', dest='njobs', help='Number of jobs to run')
        parser.add_argument('--filename', action='append', dest='filename',
            help='Filename to import. Can be used multiple times')
        parser.add_argument('--purge', action='store_true', dest='purge', default=False,
            help='Purge existing construct records')

    # source file directory
    construct_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'constructs'])

    # source files
    filenames = os.listdir(construct_data_dir)

    def handle(self, *args, **options):
        # delete any existing construct data
        if options['purge']:
            try:
                self.purge_constructs()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
        
        # how many jobs to run?
        if 'njobs' in options and options['njobs']:
            njobs = int(options['njobs'])
        else:
            njobs = 1

        # where filenames specified?
        if options['filename']:
            self.filenames = options['filename']

        try:
            self.logger.info('CREATING CONSTRUCTS')
            self.prepare_input(njobs, self.filenames)
            self.logger.info('COMPLETED CREATING CONSTRUCTS')
        except Exception as msg:
            print(msg)
            self.logger.error(msg)
    
    def purge_constructs(self):
        try:
            pst = ProteinSequenceType.objects.get(slug='mod')
            Protein.objects.filter(sequence_type=pst).delete()
        except ProteinSequenceType.DoesNotExist:
            self.logger.warning('ProteinSequenceType mod not found: nothing to delete.')

    def main_func(self, positions):
        # filenames
        if not positions[1]:
            filenames = self.filenames[positions[0]:]
        else:
            filenames = self.filenames[positions[0]:positions[1]]

        # parse files
        for source_file in filenames:
            source_file_path = os.sep.join([self.construct_data_dir, source_file])
            if os.path.isfile(source_file_path) and source_file[0] != '.':
                self.logger.info('Reading file {}'.format(source_file_path))
                # read the yaml file
                with open(source_file_path, 'r') as f:
                    sd = yaml.load(f)

                    # is a protein specified?
                    if 'protein' not in sd:
                        self.logger.error('Protein not specified for construct, skipping')
                        continue

                    # fetch the parent protein
                    try:
                        ppc = ProteinConformation.objects.prefetch_related('protein__family', 'protein__species',
                            'protein__residue_numbering_scheme').get(protein__entry_name=sd['protein'],
                            state__slug=settings.DEFAULT_PROTEIN_STATE)
                    except ProteinConformation.DoesNotExist:
                        # abort if parent protein is not found
                        self.logger.error('Parent protein {} for construct {} not found, aborting!'.format(
                            sd['protein'], sd['name']))
                        continue

                    # sequence type
                    try:
                        sequence_type, created = ProteinSequenceType.objects.get_or_create(slug='mod',
                            defaults={'name': 'Modified'})
                        if created:
                            self.logger.info('Created sequence type {}'.format(sequence_type))
                    except IntegrityError:
                        sequence_type = ProteinSequenceType.objects.get(slug='mod')

                    # protein source
                    try:
                        protein_source, created = ProteinSource.objects.get_or_create(name='OTHER')
                        if created:
                            self.logger.info('Created protein source {}'.format(protein_source))
                    except IntegrityError:
                        protein_source = ProteinSource.objects.get(name='OTHER')

                    # create a protein record
                    p = Protein()
                    p.parent = ppc.protein
                    p.family = ppc.protein.family
                    p.species = ppc.protein.species
                    p.residue_numbering_scheme = ppc.protein.residue_numbering_scheme
                    p.sequence_type= sequence_type
                    p.source = protein_source
                    p.entry_name = slugify(strip_tags(sd['name']))
                    p.name = sd['name']
                    p.sequence = ppc.protein.sequence
                    
                    # save protein (construct)
                    try:
                        p.save()
                        self.logger.info('Created construct {} with parent protein {}'.format(p.name,
                            ppc.protein.entry_name))
                    except:
                        self.logger.error('Failed creating construct {} with parent protein {}'.format(p.name,
                            ppc.protein.entry_name))
                        continue

                    # create protein conformation record
                    pc = ProteinConformation()
                    pc.protein = p
                    pc.state = ProteinState.objects.get(slug=settings.DEFAULT_PROTEIN_STATE)
                    try:
                        pc.save()
                        self.logger.info('Created conformation {} of protein {}'.format(pc.state.name, p.name))
                    except:
                        self.logger.error('Failed creating conformation {} of protein {}'.format(pc.state.name,
                            p.entry_name))

                    # create residue records
                    truncations = []
                    if 'truncations' in sd and sd['truncations']:
                        for t in sd['truncations']:
                            truncations += list(range(t[0],t[1]+1))

                    mutations = {}
                    if 'mutations' in sd and sd['mutations']:
                        for m in sd['mutations']:
                            res_num = m[1:-1]
                            mutations[res_num] = {
                                'wt_res': m[0],
                                'mut_res': m[-1],
                                'full': m,
                            }

                    # fusion proteins
                    split_segments = {}
                    if 'fusion_proteins' in sd and sd['fusion_proteins']:
                        for fp in sd['fusion_proteins']:
                            fp_start = Residue.objects.get(protein_conformation=ppc,
                                sequence_number=fp['positions'][0])
                            fp_end = Residue.objects.get(protein_conformation=ppc, sequence_number=fp['positions'][1])
                            # if the fusion protein is inserted within only one segment (the usual case), split that
                            # segment into two segments
                            if fp_start and fp_start.protein_segment == fp_end.protein_segment:
                                # get/create split protein segments
                                slug_1 = fp_start.protein_segment.slug + "_1"
                                try:
                                    segment_before, created = ProteinSegment.objects.get_or_create(slug=slug_1,
                                        defaults={'name': fp_start.protein_segment.name,
                                        'category': fp_start.protein_segment.category, 'partial': True})
                                    if created:
                                        self.logger.info('Created protein segment {}'.format(segment_before))
                                except IntegrityError:
                                    segment_before = ProteinSegment.objects.get(slug=slug_1)

                                slug_2 = fp_start.protein_segment.slug + "_2"
                                try:
                                    segment_after, created = ProteinSegment.objects.get_or_create(slug=slug_2,
                                        defaults={'name': fp_start.protein_segment.name,
                                        'category': fp_start.protein_segment.category, 'partial': True})
                                    if created:
                                        self.logger.info('Created protein segment {}'.format(segment_after))
                                except IntegrityError:
                                    segment_after = ProteinSegment.objects.get(slug=slug_2)

                                # keep track of  information about split segments
                                split_segments[fp_start.protein_segment.slug] = {
                                    'start': {
                                        'sequence_number': fp['positions'][0],
                                        'segment': segment_before,
                                    },
                                    'end': {
                                        'sequence_number': fp['positions'][1],
                                        'segment': segment_after,
                                    },
                                }

                            # get/insert fusion protein
                            fusion, create = ProteinFusion.objects.get_or_create(name=fp['name'], defaults={
                                'sequence': fp['sequence']})

                            # create relationship with protein
                            ProteinFusionProtein.objects.create(protein=p, protein_fusion=fusion,
                                segment_before=segment_before, segment_after=segment_after)

                    prs = Residue.objects.filter(protein_conformation=ppc).prefetch_related(
                        'protein_conformation__protein', 'protein_segment', 'generic_number',
                        'display_generic_number__scheme', 'alternative_generic_numbers__scheme')
                    updated_sequence = ''
                    for pr in prs:
                        if pr.sequence_number not in truncations:
                            r = Residue()
                            r.protein_conformation = pc
                            r.generic_number = pr.generic_number
                            r.display_generic_number = pr.display_generic_number
                            r.sequence_number = pr.sequence_number
                            
                            # check for split segments
                            if pr.protein_segment.slug in split_segments:
                                rsns = split_segments[pr.protein_segment.slug]['start']['sequence_number']
                                rsne = split_segments[pr.protein_segment.slug]['end']['sequence_number']
                                if r.sequence_number <= rsns:
                                    r.protein_segment = split_segments[pr.protein_segment.slug]['start']['segment']
                                elif r.sequence_number >= rsne:
                                    r.protein_segment = split_segments[pr.protein_segment.slug]['end']['segment']
                            else:
                                r.protein_segment = pr.protein_segment

                            # amino acid, check for mutations
                            if r.sequence_number in mutations:
                                if mutations[r.sequence_number]['wt_res'] == pr.amino_acid:
                                    r.amino_acid = mutations[r.sequence_number]['mut_res']
                                else:
                                    self.logger.error('Mutation {} in construct {} does not match wild-type sequence' \
                                        + ' of {}'.format(mutations[r.sequence_number]['full'], pc.protein.name,
                                        ppc.protein.entry_name))
                            else:
                                r.amino_acid = pr.amino_acid

                            # save amino acid to updated sequence
                            updated_sequence += r.amino_acid

                            # save residue before populating M2M relations
                            r.save()

                            # alternative generic numbers
                            agns = pr.alternative_generic_numbers.all()
                            for agn in agns:
                                r.alternative_generic_numbers.add(agn)
                    
                    # update sequence
                    p.sequence = updated_sequence
                    p.save()
