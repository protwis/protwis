from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from common.models import WebResource

import logging
import shlex
import os

class Command(BaseCommand):
    help = 'Reads source data and creates common database tables'

    logger = logging.getLogger(__name__)

    resource_source_file = os.sep.join([settings.DATA_DIR, 'common_data', 'resources.txt'])

    def handle(self, *args, **options):
        functions = [
            'create_resources',
        ]
        
        # execute functions
        for f in functions:
            try:
                getattr(self, f)()
            except Exception as msg:
                # print(msg)
                self.logger.error(msg)

    def create_resources(self):
        self.logger.info('Parsing file ' + self.resource_source_file)
        self.logger.info('CREATING RESOURCES')

        with open(self.resource_source_file, "r", encoding='UTF-8') as resource_source_file:
            for row in resource_source_file:
                split_row = shlex.split(row)

                wr = WebResource()
                wr.slug = split_row[0]
                wr.name = split_row[1]
                wr.url = split_row[2]

                try:
                    wr.save()
                    self.logger.info('Created resource ' + wr.name)
                except:
                    self.logger.error('Failed creating resource ' + wr.name)
                    continue

        self.logger.info('COMPLETED CREATING RESOURCES')