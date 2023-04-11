from django.core.management.base import BaseCommand
from build.management.commands.base_build import Command as BaseBuild

from django.db import transaction

from contactnetwork.models import *
from residue.models import Residue

import contactnetwork.interaction as ci

import logging
import datetime

from contactnetwork.cube import compute_interactions

from django.contrib.contenttypes.models import ContentType

class Command(BaseBuild):
    help = 'Compute interactions for all available crystals.'

    logger = logging.getLogger(__name__)
    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
            type=int,
            action='store',
            dest='proc',
            default=1,
            help='Number of processes to run')

    def handle(self, *args, **options):
        self.delete_all()
        self.structures = Structure.objects.all().exclude(structure_type__slug__startswith='af-')
        self.prepare_input(options['proc'], self.structures)
        self.logger.info('Finished building crystal interaction data for all PDBs!')

    def delete_all(self):
        VanDerWaalsInteraction.objects.all().delete()
        HydrophobicInteraction.objects.all().delete()
        PolarBackboneSidechainInteraction.objects.all().delete()
        PolarSidechainSidechainInteraction.objects.all().delete()
        FaceToFaceInteraction.objects.all().delete()
        FaceToEdgeInteraction.objects.all().delete()
        PiCationInteraction.objects.all().delete()
        InteractingResiduePair.objects.all().delete()
        self.logger.info('Deleted crystal interactions data all PDBs...')

    # @transaction.atomic
    def main_func(self, positions, iteration,count,lock):
        while count.value<len(self.structures):
            with lock:
                s = self.structures[count.value]
                pdb_code = s.protein_conformation.protein.entry_name
                count.value +=1
                self.logger.info('Generating crystal interactions data for PDB \'{}\'... ({} out of {})'.format(pdb_code, count.value, len(self.structures)))

            try:
                # interacting_pairs = compute_interactions(pdb_code)
                interacting_pairs = compute_interactions(pdb_code, do_interactions=True)
            except:
                self.logger.error('Error with computing interactions (%s)' % (pdb_code))
                continue

            for p in interacting_pairs:
                # Create the pair
                res1_seq_num = p.get_residue_1().id[1]
                res2_seq_num = p.get_residue_2().id[1]
                conformation = s.protein_conformation

                # Get the residues
                try:
                    res1 = Residue.objects.get(sequence_number=res1_seq_num, protein_conformation=conformation)
                    res2 = Residue.objects.get(sequence_number=res2_seq_num, protein_conformation=conformation)
                except Residue.DoesNotExist:
                    self.logger.warning('Error with pair between %s and %s (%s)' % (res1_seq_num,res2_seq_num,conformation))
                    continue

                # Save the pair
                pair = InteractingResiduePair()
                pair.res1 = res1
                pair.res2 = res2
                pair.referenced_structure = s
                pair.save()

                # Add the interactions to the pair
                for i in p.get_interactions():
                    if type(i) is ci.VanDerWaalsInteraction:
                        ni = VanDerWaalsInteraction()
                        ni.interacting_pair = pair
                        ni.save()
                    elif type(i) is ci.HydrophobicInteraction:
                        ni = HydrophobicInteraction()
                        ni.interacting_pair = pair
                        ni.save()
                    elif type(i) is ci.PolarSidechainSidechainInteraction:
                        ni = PolarSidechainSidechainInteraction()
                        ni.interacting_pair = pair
                        ni.is_charged_res1 = i.is_charged_res1
                        ni.is_charged_res2 = i.is_charged_res2
                        ni.save()
                    elif type(i) is ci.PolarBackboneSidechainInteraction:
                        ni = PolarBackboneSidechainInteraction()
                        ni.interacting_pair = pair
                        ni.is_charged_res1 = i.is_charged_res1
                        ni.is_charged_res2 = i.is_charged_res2
                        ni.res1_is_sidechain = False
                        ni.save()
                    elif type(i) is ci.PolarSideChainBackboneInteraction:
                        ni = PolarBackboneSidechainInteraction()
                        ni.interacting_pair = pair
                        ni.is_charged_res1 = i.is_charged_res1
                        ni.is_charged_res2 = i.is_charged_res2
                        ni.res1_is_sidechain = True
                        ni.save()
                    elif type(i) is ci.FaceToFaceInteraction:
                        ni = FaceToFaceInteraction()
                        ni.interacting_pair = pair
                        ni.save()
                    elif type(i) is ci.FaceToEdgeInteraction:
                        ni = FaceToEdgeInteraction()
                        ni.interacting_pair = pair
                        ni.res1_has_face = True
                        ni.save()
                    elif type(i) is ci.EdgeToFaceInteraction:
                        ni = FaceToEdgeInteraction()
                        ni.interacting_pair = pair
                        ni.res1_has_face = False
                        ni.save()
                    elif type(i) is ci.PiCationInteraction:
                        ni = PiCationInteraction()
                        ni.interacting_pair = pair
                        ni.res1_has_pi = True
                        ni.save()
                    elif type(i) is ci.CationPiInteraction:
                        ni = PiCationInteraction()
                        ni.interacting_pair = pair
                        ni.res1_has_pi = False
                        ni.save()

            self.logger.info('Generated crystal interactions data for PDB \'{}\'...'.format(pdb_code))
