from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from common.models import Publication,PublicationJournal
import logging
import os
import yaml

class Command(BaseCommand):
    help = 'Exports reference positions for all proteins'

    logger = logging.getLogger(__name__)

    export_dir_path = os.sep.join([settings.DATA_DIR, 'publications_data'])
    if not os.path.exists(export_dir_path):
        os.makedirs(export_dir_path)
    def handle(self, *args, **options):
        functions = [
            'export_publications',
        ]
        
        # execute functions
        for f in functions:
            try:
                getattr(self, f)()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)

    def export_publications(self):
        self.logger.info('EXPORTING PUBLICATIONS')

        # fetch all ligands
        ps = Publication.objects.all().prefetch_related('journal','web_link')
        export_data = []

        for p in ps:
            export = {}
            export['title'] = p.title
            export['authors'] = p.authors
            export['year'] = p.year
            export['reference'] = p.reference
            export['journal_slug'] = p.journal.slug
            export['journal_name'] = p.journal.name
            export['weblink_index'] = p.web_link.index
            export['weblink_resource'] = p.web_link.web_resource.slug
            export_data.append(export)

            # open export file for writing
        export_file_path = os.sep.join([self.export_dir_path,'publications.yaml'])
        with open(export_file_path, 'w') as export_file:
            yaml.dump(export_data, export_file, default_flow_style=False)

        self.logger.info('COMPLETED EXPORTING PUBLICATIONS - Entries:'+str(len(export_data)))