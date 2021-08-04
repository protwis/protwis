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
        # parser.add_argument('--hommod',
        #                     action='store_true',
        #                     dest='hommod',
        #                     default=False,
        #                     help='Include build of homology models')
        parser.add_argument('--phase',
                            type=int,
                            action='store',
                            dest='phase',
                            default=None,
                            help='Specify build phase to run (1 or 2, default: None)')

    def handle(self, *args, **options):
        if options['test']:
            print('Running in test mode')

        phase1 = [
            ['upload_excel_bias'],
            ['build_bias_data-inferred'],
            ['build_bias_data-inferred-subtypes'],
            ['build_bias_data-predicted']
        ]

        if options['phase']:
            if options['phase']==1:
                commands = phase1
        else:
            commands = phase1

        for c in commands:
            print('\n{} Running {}\n'.format(
                datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'), c[0]))
            call_command(c[0])

        print('{} Build completed'.format(datetime.datetime.strftime(
            datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')))
