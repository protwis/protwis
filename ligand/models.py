from django.db import models

class Ligand(models.Model):
    name = models.CharField(max_length=50)
    role = models.ForeignKey('LigandRole')
    xray_name = models.CharField(max_length=3, null=True) #residue name from pdb file
    smiles = models.TextField()
    
    class Meta():
        db_table = 'ligand'
        
        
class LigandRole(models.Model):
    role = models.CharField(max_length=20)
    
    def __str__(self):
        return self.role
    
    class Meta():
        db_table = 'ligand_role'
    