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

    protein = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    residue = models.ForeignKey('residue.Residue', on_delete=models.CASCADE)
    seq_similarity = models.DecimalField(max_digits=5, decimal_places=4)
    seq_identity = models.DecimalField(max_digits=5, decimal_places=4)
    paralog_score = models.DecimalField(max_digits=5, decimal_places=4)

    def __str__(self):
        return self.protein.name

    class Meta():
        db_table = 'gprotein_barcode'

class SignprotComplex(models.Model):
    protein = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    alpha = models.CharField(max_length=1)
    beta_chain = models.CharField(max_length=1, null=True)
    beta_protein = models.ForeignKey('protein.Protein', related_name='beta_protein', on_delete=models.CASCADE, null=True)
    gamma_chain = models.CharField(max_length=1, null=True)
    gamma_protein = models.ForeignKey('protein.Protein', related_name='gamma_protein', on_delete=models.CASCADE, null=True)

    def __str__(self):
        return '<SignprotComplex: {} {}>'.format(self.protein.entry_name, self.structure.pdb_code.index)

    class Meta():
        db_table = 'signprot_complex'

class SignprotInteractions(models.Model):

    gpcr_residue = models.ForeignKey('residue.Residue', on_delete=models.CASCADE, related_name='gpcr_residue')
    signprot_residue = models.ForeignKey('residue.Residue', on_delete=models.CASCADE, related_name='signprot_residue')
    interaction_type = models.CharField(max_length=200)
    structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    
    def __str__(self):
        return '{} between {} and {}'.format(self.interaction_type, self.gpcr_residue, self.signprot_residue)
