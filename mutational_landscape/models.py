from django.db import models

# Create your models here.

class NaturalMutations(models.Model):

    protein = models.ForeignKey('protein.Protein')
    residue = models.ForeignKey('residue.Residue')
    amino_acid = models.CharField(max_length=3)
    type = models.CharField(max_length=30)
    allele_frequency = models.FloatField()
    allele_count = models.IntegerField()
    allele_number = models.IntegerField()
    number_homozygotes = models.IntegerField()
    sift_score = models.FloatField()
    polyphen_score = models.FloatField()

    def __str__(self):
        return self.protein.name + '_' + str(self.residue.sequence_number) + '_' + self.residue.amino_acid

    class Meta():
        db_table = 'mutation_natural'
        # unique_together = ('protein','residue','amino_acid','allele_frequency')

class CancerMutations(models.Model):

    protein = models.ForeignKey('protein.Protein')
    residue = models.ForeignKey('residue.Residue')
    amino_acid = models.CharField(max_length=1)
    cancer_type = models.CharField(max_length=100)
    # allele_frequency = models.FloatField()
    # allele_count = models.IntegerField()

    def __str__(self):
        return self.protein.name

    class Meta():
        db_table = 'mutation_cancer'
        # unique_together = ('protein','residue','amino_acid')

class DiseaseMutations(models.Model):

    protein = models.ForeignKey('protein.Protein')
    residue = models.ForeignKey('residue.Residue')
    amino_acid = models.CharField(max_length=1)
    disease = models.CharField(max_length=100)

    def __str__(self):
        return self.protein.name

    class Meta():
        db_table = 'mutation_disease'
        # unique_together = ('protein','residue','amino_acid')
