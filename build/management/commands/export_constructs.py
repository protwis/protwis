from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from construct.models import *
import logging
import os
import json

class Command(BaseCommand):
    help = 'Exports all construct json files'

    logger = logging.getLogger(__name__)

    export_dir_path = os.sep.join([settings.DATA_DIR, 'construct_export'])
    if not os.path.exists(export_dir_path):
        os.makedirs(export_dir_path)
    def handle(self, *args, **options):
        functions = [
            'export_constructs',
        ]
        
        # execute functions
        for f in functions:
            try:
                getattr(self, f)()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)

    def export_constructs(self):
        self.logger.info('EXPORTING CONSTRUCTS')

        # fetch all ligands
        cs = Construct.objects.all()
            # open export file for writing
        for c in cs:
            export_file_path = os.sep.join([self.export_dir_path,c.name + '.json'])
            f = open(export_file_path,"w") #opens file with name of "test.txt"
            f.write(c.json)
            f.close()
            # with open(export_file_path, 'w') as export_file:
            #     json.dump(c.json, export_file, indent=4, separators=(',', ': '))

        self.logger.info('COMPLETED EXPORTING CONSTRUCT - Entries:'+str(len(cs)))