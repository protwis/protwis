from django.db import models

# Create your models here.

# ==========

# class DrugsIndication(models.Model):
#     indications = models.ManyToManyField('indication', through='ProteinDrug')

#     def __str__(self):
#         return self.drug_name

#     class Meta():
#         db_table = 'protein_drug_indication'

class Drugs(models.Model):
    target = models.ManyToManyField('protein.Protein')
    name = models.CharField(max_length=100)
    
    synonym = models.CharField(max_length=100, null=True) # addName as list, make as ManyField

    drugtype = models.CharField(max_length=50, null=True)

    status = models.CharField(max_length=30)
    indication = models.CharField(max_length=200, null=True)
    approval = models.IntegerField()

    externallink = models.CharField(max_length=200, null=True) # Link Framework

    novelty = models.CharField(max_length=50) #Boolean

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'drugs'
