from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from common.models import WebResource
from protein.models import ProteinSegment
from residue.models import ResidueNumberingScheme

import logging
import shlex
import os

class Command(BaseCommand):
    help = 'Reads source data and creates common database tables'

    logger = logging.getLogger(__name__)

    resource_source_file = os.sep.join([settings.DATA_DIR, 'common_data', 'resources.txt'])
    segment_source_file = os.sep.join([settings.DATA_DIR, 'protein_data', 'segments.txt'])
    residue_number_scheme_source_file = os.sep.join([settings.DATA_DIR, 'residue_data', 'generic_numbers',
        'schemes.txt'])

    def handle(self, *args, **options):
        functions = [
            'create_resources',
            'create_protein_segments',
            'create_residue_numbering_schemes',
        ]
        
        # execute functions
        for f in functions:
            try:
                getattr(self, f)()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)

    def create_resources(self):
        self.logger.info('CREATING RESOURCES')
        self.logger.info('Parsing file ' + self.resource_source_file)

        with open(self.resource_source_file, "r", encoding='UTF-8') as resource_source_file:
            for row in resource_source_file:
                split_row = shlex.split(row)

                # create resource
                try:
                    defaults = {
                        'name': split_row[1],
                        'url': split_row[2]
                    }

                    wr, created = WebResource.objects.get_or_create(slug=split_row[0], defaults=defaults)

                    if created:
                        self.logger.info('Created resource ' + wr.slug)
                    
                except:
                    self.logger.error('Failed creating resource ' + split_row[0])
                    continue

        self.logger.info('COMPLETED CREATING RESOURCES')

    def create_protein_segments(self):
        self.logger.info('CREATING PROTEIN SEGMENTS')
        self.logger.info('Parsing file ' + self.segment_source_file)

        with open(self.segment_source_file, "r", encoding='UTF-8') as segment_file:
            for row in segment_file:
                split_row = shlex.split(row)

                # create segment
                try:
                    defaults={
                        'category': split_row[1],
                        'name': split_row[2]
                    }

                    s, created = ProteinSegment.objects.get_or_create(slug=split_row[0], defaults=defaults)

                    if created:
                        self.logger.info('Created protein segment ' + s.name)
                except:
                    self.logger.error('Failed creating protein segment {}'.format(split(row[0])))
                    continue

        self.logger.info('COMPLETED CREATING PROTEIN SEGMENTS')

    def create_residue_numbering_schemes(self):
        self.logger.info('CREATING RESIDUE NUMBERING SCHEMES')
        self.logger.info('Parsing file ' + self.residue_number_scheme_source_file)

        with open(self.residue_number_scheme_source_file, "r", encoding='UTF-8') as residue_number_scheme_source_file:
            for row in residue_number_scheme_source_file:
                split_row = shlex.split(row)

                # create scheme
                try:
                    defaults={
                        'short_name': split_row[1],
                        'name': split_row[2]
                    }
                    if len(split_row) == 4:
                        try:
                            prns = split_row[3]
                            parent = ResidueNumberingScheme.objects.get(slug=prns)
                        except ResidueNumberingScheme.DoesNotExists:
                            raise Exception('Parent scheme {} does not exists, aborting!'.format(prns))
                        defaults['parent'] = parent

                    s, created = ResidueNumberingScheme.objects.get_or_create(slug=split_row[0], defaults=defaults)
                    
                    if created:
                        self.logger.info('Created residue numbering scheme ' + s.short_name)
                except:
                    self.logger.error('Failed creating residue numbering scheme {}'.format(split_row[0]))
                    continue

        self.logger.info('COMPLETED CREATING RESIDUE NUMBERING SCHEMES')