from django.core.management.base import BaseCommand

from django.db import transaction

from contactnetwork.models import *
from residue.models import Residue

import contactnetwork.interaction as ci

import logging

from contactnetwork.cube import compute_interactions

class Command(BaseCommand):
    help = 'Compute interactions for all available crystals.'

    logger = logging.getLogger(__name__)

    def handle(self, *args, **options):
        self.delete_all()
        self.compute_all()

    def compute_all(self):
        structures = Structure.objects.all()
        for s in structures:
            self.build_interactions(s)
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

    @transaction.atomic
    def build_interactions(self, s):
        pdb_code = s.protein_conformation.protein.entry_name

        interacting_pairs = compute_interactions(pdb_code)

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
