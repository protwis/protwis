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
    name = models.CharField(max_length=75)
    
    synonym = models.CharField(max_length=300, null=True) # addName as list, make as ManyField

    drugtype = models.CharField(max_length=40, null=True)

    status = models.CharField(max_length=15)
    indication = models.CharField(max_length=150, null=True)
    approval = models.CharField(max_length=5)

    externallink = models.CharField(max_length=150, null=True) # Link Framework

    novelty = models.CharField(max_length=15) #Boolean

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'drugs'
