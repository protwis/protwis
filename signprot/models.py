from django.db import models

from protein.models import Protein
from structure.models import Structure

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

    protein = models.ForeignKey('protein.Protein')
    residue = models.ForeignKey('residue.Residue')
    seq_similarity = models.DecimalField(max_digits=5, decimal_places=4)
    seq_identity = models.DecimalField(max_digits=5, decimal_places=4)
    paralog_score = models.DecimalField(max_digits=5, decimal_places=4)

    def __str__(self):
        return self.protein.name

    class Meta():
        db_table = 'gprotein_barcode'

class SignprotComplex(models.Model):
    protein = models.ForeignKey('protein.Protein')
    structure = models.ForeignKey('structure.Structure')
    chain = models.CharField(max_length=1)

    def __str__(self):
        return '<SignprotComplex: {} {}>'.format(self.protein.entry_name, self.structure.pdb_code.index)

    class Meta():
        db_table = 'signprot_complex'