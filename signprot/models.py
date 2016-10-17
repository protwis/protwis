from django.db import models

# Create your models here.

class SignprotStructure(models.Model):
    origin = models.ManyToManyField('protein.Protein')
    
    PDB_code = models.CharField(max_length=4)
    resolution = models.DecimalField(max_digits=5, decimal_places=3)

    def __str__(self):
        return self.PDB_code

    class Meta():
        db_table = 'signprot_structure'
