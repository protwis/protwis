from django.db import models

class Fragment(models.Model):
    residue = models.ForeignKey('residue.Residue')
    protein = models.ForeignKey('protein.Protein')
    crystalstructure = models.ForeignKey('crystalstructure.Crystalstructure')
    ligand = models.ForeignKey('ligand.Ligand')
    interaction = models.ForeignKey('Interaction')
    filename = models.CharField(max_length=50)
    
    def __str__(self):
        return self.filename
    
    class Meta():
        db_table = 'fragment'

 
class Interaction(models.Model):
    interaction_type = models.CharField(max_length=10)
    description = models.CharField(max_length = 200)
    
    class Meta():
        db_table = 'interaction'

