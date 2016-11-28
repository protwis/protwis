from django.db import models

# Create your models here.

class NaturalMutations(models.Model):

    protein = models.ForeignKey('protein.Protein')
    residue = models.ForeignKey('residue.Residue')
    amino_acid = models.CharField(max_length=1)
    allele_frequency = models.FloatField()
    allele_count = models.FloatField()
    allele_number = models.FloatField()
    number_homozygotes = models.FloatField()

    def __str__(self):
        return self.protein.name

    class Meta():
        db_table = 'mutation_natural'
        unique_together = ('protein','residue','amino_acid')

class CancerMutations(models.Model):

    protein = models.ForeignKey('protein.Protein')
    residue = models.ForeignKey('residue.Residue')
    amino_acid = models.CharField(max_length=1)
    cancer_type = models.CharField(max_length=100)

    def __str__(self):
        return self.protein.name

    class Meta():
        db_table = 'mutation_cancer'
        unique_together = ('protein','residue','amino_acid')

class DiseaseMutations(models.Model):

    protein = models.ForeignKey('protein.Protein')
    residue = models.ForeignKey('residue.Residue')
    amino_acid = models.CharField(max_length=1)
    disease = models.CharField(max_length=100)

    def __str__(self):
        return self.protein.name

    class Meta():
        db_table = 'mutation_disease'
        unique_together = ('protein','residue','amino_acid')
