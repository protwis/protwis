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

    def handle(self, *args, **options):
        if options['test']:
            print('Running in test mode')

        commands = [
            ['build_common'],
            ['build_human_proteins'],
            ['build_blast_database'],
            ['build_other_proteins', {'constructs_only': options['test'] ,'proc': options['proc']}], # build only constructs in test mode
            ['build_annotation', {'proc': options['proc']}],
            ['build_blast_database'],
            ['build_links'],
            ['build_construct_proteins'],
            ['build_structures', {'proc': options['proc']}],
            ['build_endogenous_ligands'],
            ['build_structure_angles'],
            ['build_distance_representative'],
            ['build_contact_representative'],
            ['build_construct_data'],
            ['update_construct_mutations'],
            ['build_ligands_from_cache', {'proc': options['proc'], 'test_run': options['test']}],
            ['build_ligand_assays', {'proc': options['proc'], 'test_run': options['test']}],
            ['build_mutant_data', {'proc': options['proc'], 'test_run': options['test']}],
            ['build_protein_sets'],
            ['build_consensus_sequences', {'proc': options['proc']}],
            ['build_g_proteins'],
            ['build_arrestins'],
            ['build_signprot_complex'],
            ['build_g_protein_structures'],
            ['build_structure_extra_proteins'],
            ['build_drugs'],
            ['build_nhs'],
            ['build_mutational_landscape'],
            ['build_residue_sets'],
            ['build_dynamine_annotation', {'proc': options['proc']}],
            ['build_blast_database'],
            ['build_complex_interactions'],
            # ['build_homology_models', ['--update', '-z'], {'proc': options['proc'], 'test_run': options['test']}],
            ['build_text'],
            ['build_release_notes'],
        ]

        for c in commands:
            print('{} Running {}'.format(
                datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'), c[0]))
            if len(c) == 2:
                call_command(c[0], **c[1])
            elif len(c) == 3:
                call_command(c[0], *c[1], **c[2])
            else:
                call_command(c[0])

        print('{} Build completed'.format(datetime.datetime.strftime(
            datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')))
