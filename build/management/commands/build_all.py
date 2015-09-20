from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command

import datetime


class Command(BaseCommand):
    help = 'Runs all build functions'

    def add_arguments(self, parser):
        parser.add_argument('--njobs', action='store', dest='njobs', help='Number of jobs to run')

    def handle(self, *args, **options):
        # how many jobs to run?
        if 'njobs' in options and options['njobs']:
            njobs = int(options['njobs'])
        else:
            njobs = 1

        commands = [
            ['build_common'],
            ['build_human_proteins'],
            ['build_human_residues', {'njobs': njobs}],
            ['build_blast_database'],
            ['build_other_proteins'],
            ['build_other_residues', {'njobs': njobs}],
            ['build_blast_database'],
            ['build_links'],
            ['build_construct_proteins', {'njobs': njobs}],
            ['build_structures', {'njobs': njobs}],
            ['build_mutant_data'],
            ['find_protein_templates', {'njobs': njobs}],
            ['update_alignments', {'njobs': njobs}],
            ['build_protein_sets'],
            ['build_consensus_sequences', {'njobs': njobs}],
        ]

        for c in commands:
            print('{} Running {}'.format(
                datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'), c[0]))
            if len(c) > 1:
                call_command(c[0], **c[1])
            else:
                call_command(c[0])

        print('{} Build completed'.format(datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')))