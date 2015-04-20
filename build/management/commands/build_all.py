from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command


class Command(BaseCommand):
    help = 'Runs all build functions'

    def handle(self, *args, **options):
        commands = [
            'build_common',
            'build_proteins',
            # 'build_residues',
            # 'build_constructs',
            # 'build_structures',
        ]

        for c in commands:
            print('Running ' + c)
            call_command(c)