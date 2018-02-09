from structure.models import Structure

from django.db import models

from polymorphic.models import PolymorphicModel

class InteractingResiduePair(PolymorphicModel):
    referenced_structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    res1 = models.ForeignKey('residue.Residue', related_name='residue1', on_delete=models.CASCADE)
    res2 = models.ForeignKey('residue.Residue', related_name='residue2', on_delete=models.CASCADE)

    class Meta():
        db_table = 'interacting_residue_pair'


class Interaction(PolymorphicModel):
    interacting_pair = models.ForeignKey('contactnetwork.InteractingResiduePair', on_delete=models.CASCADE)

    class Meta():
        db_table = 'interaction'


class VanDerWaalsInteraction(Interaction):
    def get_name(self):
        return 'van-der-waals'

    class Meta():
        db_table = 'interaction_van_der_waals'


class HydrophobicInteraction(Interaction):
    def get_name(self):
        return 'hydrophobic'

    class Meta():
        db_table = 'interaction_hydrophobic'


class PolarInteraction(Interaction):
    is_charged_res1 = models.BooleanField()
    is_charged_res2 = models.BooleanField()

    class Meta():
        db_table = 'interaction_polar'


class PolarSidechainSidechainInteraction(PolarInteraction):
    def get_name(self):
        return 'polar-sidechain-sidechain'

    class Meta():
        db_table = 'interaction_polar_sidechain_sidechain'


class PolarBackboneSidechainInteraction(PolarInteraction):
    res1_is_sidechain = models.BooleanField()

    def get_name(self):
        return 'polar-backbone-sidechain'

    class Meta():
        db_table = 'interaction_polar_backbone_sidechain'


class AromaticInteraction(Interaction):
    class Meta():
        db_table = 'interaction_aromatic'


class FaceToFaceInteraction(AromaticInteraction):
    def get_name(self):
        return 'aromatic-face-to-face'

    class Meta():
        db_table = 'interaction_aromatic_face_face'


class FaceToEdgeInteraction(AromaticInteraction):
    res1_has_face = models.BooleanField()

    def get_name(self):
        return 'aromatic-face-to-edge'

    class Meta():
        db_table = 'interaction_aromatic_face_edge'


class PiCationInteraction(AromaticInteraction):
    res1_has_pi = models.BooleanField()

    def get_name(self):
        return 'aromatic-pi-cation'

    class Meta():
        db_table = 'interaction_aromatic_pi_cation'
