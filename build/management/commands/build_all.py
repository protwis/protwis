from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command

import datetime


class Command(BaseCommand):
    help = 'Runs all build functions'

    def handle(self, *args, **options):
        commands = [
            ['build_common'],
            ['build_proteins'],
            ['build_residues', 'all'],
            ['build_constructs'],
            ['build_structures'],
        ]

        for c in commands:
            print('{} Running {}'.format(
                datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'), c[0]))
            if len(c) > 1:
                call_command(c[0],c[1])
            else:
                call_command(c[0])