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

        if options['proc']>4:
            safe_proc_num = 4
        else:
            safe_proc_num = options['proc']

        phase1 = [
            ['clear_cache'],
            ['build_common'],
            ['build_citations'],
            ['build_human_proteins'],
            ['build_blast_database'],
            ['build_other_proteins', {'constructs_only': options['test'] ,'proc': options['proc']}], # build only constructs in test mode
            ['build_annotation', {'proc': options['proc']}],
            ['build_blast_database'],
            ['build_links'],
            ['build_construct_proteins'],
            ['build_experimental_data', {'test_run': options['test']}],
            # ['build_all_gtp_ligands', {'test_run': options['test']}],
            # ['build_endogenous_data_from_gtp_source', {'test_run': options['test']}],
            ['build_bias_preprocess_data', {'test_run': options['test']}],
            #['build_balanced_ligands', {'test_run': options['test']}],
            # ['build_chembl_data', {'test_run': options['test']}],
            ['build_mutant_data', {'test_run': options['test']}],
            ['build_structures', {'proc': safe_proc_num, 'skip_cn': options['test']}],
            ['build_consensus_sequences', {'proc': options['proc']}],
            ['build_g_proteins'],
            ['build_consensus_sequences', {'proc': options['proc'], 'signprot': 'Alpha'}],
            ['build_arrestins'],
            ['build_coupling_data'],
            ['build_consensus_sequences', {'proc': options['proc'], 'signprot': 'Arrestin'}],
            ['build_signprot_complex'],
            ['build_g_protein_structures', {'proc': options['proc']}],
            ['build_arrestin_structures'],
            ['build_structure_extra_proteins'],
            ['build_structure_model_rmsd'],
            ['build_af_models', {'proc': safe_proc_num}],
            ['build_blast_database']
        ]
        phase2 = [
            ['build_structure_angles', {'proc': options['proc']}],
            ['build_construct_data'],
            ['update_construct_mutations'],
            ['build_protein_sets'],
            ['build_drugs'],
            ['build_mutational_landscape'],
            ['build_residue_sets'],
            ['build_dynamine_annotation', {'proc': options['proc']}],
            ['build_complex_interactions'],
            ['assign_structure_states'],
            ['build_contact_representative'],
            ['build_mammalian_representative'],
            ['upload_excel_bias_pathways'],
            ['build_text'],
            ['build_release_notes'],
        ]

        if options['phase']:
            if options['phase']==1:
                commands = phase1
            elif options['phase']==2:
                commands = phase2
        else:
            commands = phase1+phase2

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
