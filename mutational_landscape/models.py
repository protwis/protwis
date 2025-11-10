from django.db import models

# Create your models here.

class NaturalMutations(models.Model):

    protein = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    residue = models.ForeignKey('residue.Residue', on_delete=models.CASCADE)
    amino_acid = models.CharField(max_length=20)
    type = models.CharField(max_length=30)
    allele_frequency = models.FloatField()
    allele_count = models.IntegerField()
    allele_number = models.IntegerField()
    number_homozygotes = models.IntegerField()
    sift_score = models.FloatField(null=True)
    polyphen_score = models.FloatField(null=True)

    def __str__(self):
        return self.protein.name + '_' + str(self.residue.sequence_number) + '_' + self.residue.amino_acid

    class Meta():
        db_table = 'mutation_natural'
        # unique_together = ('protein','residue','amino_acid','allele_frequency')



class PTMs(models.Model):

    protein = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    residue = models.ForeignKey('residue.Residue', on_delete=models.CASCADE)
    modification = models.CharField(max_length=40)
    # source = models.CharField(max_length=30)

    def __str__(self):
        return self.protein.name + '_' + str(self.residue.sequence_number) + '_' + self.residue.amino_acid

    class Meta():
        db_table = 'residue_ptm'

# class PTMsType(models.Model):
#     modification = models.CharField(max_length=100, unique=True)
#
#     def __str__(self):
#         return self.modification
#
#     class Meta():
#         db_table = 'residue_ptm_types'
