from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from build.management.commands.base_build import Command as BaseBuild

import os


class Command(BaseBuild):
    help = 'Reads custom Structure and related objects based on input'

    def add_arguments(self, parser):
        parser.add_argument('--input_folder',
            dest='input_folder',
            help='Add custom structure based on input',
            nargs='+')
        
    def handle(self, *args, **options):
        input_folder = options['input_folder'][0]
        pdb_id = input_folder.split('/')[-1].split('_')[0]

        ### Build structure object and ligand interactions
        call_command('build_structures', structure=pdb_id.upper(), custom=[os.sep.join([input_folder, f'{pdb_id}.json'])])
        files = os.listdir(input_folder)
        signprot = [i for i in files if 'signprot' in i][0]
        if len(signprot)>0:
            ### Build signprot_complex object - connecting receptor with signprot
            call_command('build_signprot_complex', custom=os.sep.join([input_folder, signprot]))
            ### Build G protein structure related objects
            call_command('build_g_protein_structures', s=[pdb_id.upper()], debug=True)
            ### ToDo signprot
            ### build_complex_interactions - within receptor and with signprot
        else:
            pass
            ### ToDo non-signprot

        ### ToDo either
        ### build_structure_angles - for structure comparison tool