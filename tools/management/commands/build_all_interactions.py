from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from structure.models import Structure

from contactnetwork.cube import *
from contactnetwork.cube import compute_interactions
from build import ligand_imports

import logging, json, os

class Command(BaseBuild):

    # help = "Function to calculate interaction for all GPCR structures."
    help = ("Recompute the interactions of all experimental GPCR structures: the "
            "intra-receptor contact network, and the ligand interactions imported "
            "from the Schrodinger deliveries as build_all imports them. Do not change "
            "the deliveries while it runs.")

    logger = logging.getLogger(__name__)
    pdbs = Structure.objects.filter(structure_type__origin='experiment').values_list('pdb_code__index', flat=True)


    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
            type=int,
            action='store',
            dest='proc',
            default=1,
            help='Number of processes to run')
        ligand_imports.add_arguments(parser)

    def handle(self, *args, **options):
        # The ligand imports of build_all. Maps and dry runs come before the long
        # contact-network pass, so a problem stops the command early. The contact
        # network writes only InteractingResiduePair and Interaction, which the
        # imports neither read nor write, so the dry runs still hold after it.
        before, imports = ligand_imports.split(
            ligand_imports.with_tests(ligand_imports.steps(options))
        )
        ligand_imports.check_deliveries(options, [c for c, _o in before + imports])
        ligand_imports.run(before)
        try:
            self.logger.info('CREATING ALL INTERACTIONS')
            self.prepare_input(options['proc'], self.pdbs)
        except Exception as msg:
            print(msg)
            self.logger.error(msg)
        ligand_imports.run(imports)
        self.logger.info('COMPLETED ALL INTERACTIONS')

    def main_func(self, positions, iteration,count,lock):
        pdbs = self.pdbs
        while count.value<len(pdbs):
            with lock:
                pdb = pdbs[count.value]
                count.value +=1
            try:
                # compute_interactions(pdb, True)
                # The contact network only: the receptor x peptide pairs are imported
                # in handle() (import_schrodinger_peptides, via ligand_imports.run).
                # compute_interactions(pdb, do_interactions=True, do_peptide_ligand=True, save_to_db=True)
                compute_interactions(pdb, do_interactions=True, do_peptide_ligand=False, save_to_db=True)
            except:
                print('Issue making interactions for',pdb)
