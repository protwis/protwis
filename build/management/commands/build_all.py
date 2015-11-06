from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command

import datetime


class Command(BaseCommand):
    help = 'Runs all build functions'

    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
                            type=int,
                            action='store',
                            dest='proc',
                            default=1,
                            help='Number of processes to run')
        parser.add_argument('-t', '--test',
                            action='store_true',
                            dest='test',
                            default=False,
                            help='Include only a subset of data for testing')

    def handle(self, *args, **options):
        if options['test']:
            print('Running in test mode')

        commands = [
            ['build_common'],
            ['build_human_proteins', {'test': options['test']}],
            ['build_human_residues', {'proc': options['proc']}],
            ['build_blast_database'],
            ['build_other_proteins', {'constructs_only': options['test']}], # build only constructs in test mode
            ['build_other_residues', {'proc': options['proc']}],
            ['build_blast_database'],
            ['build_links'],
            ['build_construct_proteins', {'proc': options['proc']}],
            ['build_structures', {'proc': options['proc']}],
            ['build_mutant_data'],
            ['find_protein_templates', {'proc': options['proc']}],
            ['update_alignments', {'proc': options['proc']}],
            ['build_protein_sets'],
            ['build_consensus_sequences', {'proc': options['proc']}],
            ['build_text'],
            ['build_release_notes'],
        ]

        for c in commands:
            print('{} Running {}'.format(
                datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'), c[0]))
            if len(c) > 1:
                call_command(c[0], **c[1])
            else:
                call_command(c[0])

        print('{} Build completed'.format(datetime.datetime.strftime(
            datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')))
