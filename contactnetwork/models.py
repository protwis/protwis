from structure.models import Structure

from django.db import models


class InteractingResiduePair(models.Model):
    referenced_structure = models.ForeignKey('structure.Structure')
    res1 = models.ForeignKey('residue.Residue', related_name='residue1')
    res2 = models.ForeignKey('residue.Residue', related_name='residue2')

    class Meta():
        db_table = 'interacting_residue_pair'


class Interaction(models.Model):
    interacting_pair = models.ForeignKey('contactnetwork.InteractingResiduePair')

    class Meta():
        abstract = True


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
        abstract = True


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
        abstract = True


class FaceToFaceInteraction(AromaticInteraction):
    def get_name(self):
        return 'face-to-face'

    class Meta():
        db_table = 'interaction_aromatic_face_face'


class FaceToEdgeInteraction(AromaticInteraction):
    res1_has_face = models.BooleanField()

    def get_name(self):
        return 'face-to-edge'

    class Meta():
        db_table = 'interaction_aromatic_face_edge'


class PiCationInteraction(AromaticInteraction):
    res1_has_pi = models.BooleanField()

    def get_name(self):
        return 'pi-cation'

    class Meta():
        db_table = 'interaction_aromatic_pi_cation'
