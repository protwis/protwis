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

class SignprotBarcode(models.Model):

    protein = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    residue = models.ForeignKey('residue.Residue', on_delete=models.CASCADE)
    seq_similarity = models.DecimalField(max_digits=5, decimal_places=4)
    seq_identity = models.DecimalField(max_digits=5, decimal_places=4)
    paralog_score = models.DecimalField(max_digits=5, decimal_places=4)

    def __str__(self):
        return self.protein.name

    class Meta():
        db_table = 'gprotein_barcode'
