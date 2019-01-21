from structure.models import Structure

from django.db import models

class InteractingResiduePair(models.Model):
    referenced_structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    res1 = models.ForeignKey('residue.Residue', related_name='residue1', on_delete=models.CASCADE)
    res2 = models.ForeignKey('residue.Residue', related_name='residue2', on_delete=models.CASCADE)

    class Meta():
        db_table = 'interacting_residue_pair'


class Interaction(models.Model): 
    interacting_pair = models.ForeignKey('contactnetwork.InteractingResiduePair', on_delete=models.CASCADE)
    interaction_type = models.CharField(max_length=100)
    specific_type = models.CharField(max_length=100, null=True)

    class Meta():
        db_table = 'interaction'


class Distance(models.Model):
    interacting_pair = models.ForeignKey('contactnetwork.InteractingResiduePair', on_delete=models.CASCADE)
    distance = models.FloatField()

    class Meta():
        db_table = 'distance'