from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection

import datetime
import logging
from optparse import make_option
import os
import shutil


class Command(BaseCommand):
    help = 'Basic functions for build scrips'

    logger = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument('-f', '--filename',
            action='store',
            dest='filename',
            help='Path to Uniprot text file')

    local_uniprot_dir = os.sep.join([settings.DATA_DIR, 'protein_data', 'uniprot'])

    def handle(self, *args, **options):
        if 'filename' in options and options['filename']:
            filename = options['filename']
            if os.path.isfile(filename):
                self.separate(filename)
        else:
            self.separate(os.sep.join([settings.DATA_DIR, 'protein_data', 'uniprot_all_proteins.txt']))
            # self.logger.error('No filename specified, aborting')

    def separate(self, input_filename):
        print(input_filename)
        with open(input_filename) as input_file:
            line_buffer = []
            acs = []
            for line in input_file:
                if line.startswith('//'):
                    # Only pick the first one
                    ac = acs[1]
                    output_filename = ac + '.txt'
                    output_file_path = os.sep.join([self.local_uniprot_dir, output_filename])
                    with open(output_file_path, "w") as output_file:
                        for bline in line_buffer:
                            output_file.write(bline)
                        output_file.write(line)
                    line_buffer = []
                    acs = []
                else:
                    line_buffer.append(line)
                    if line.startswith('AC'):
                        sline = line.split()
                        for ac in sline:
                            acs.append(ac.strip(';'))