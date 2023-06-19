from django.db import models


class Drugs(models.Model):
    target = models.ManyToManyField('protein.Protein')
    name = models.CharField(max_length=100)
    synonym = models.CharField(max_length=300, null=True) # addName as list, make as ManyField
    drugtype = models.CharField(max_length=40, null=True)
    status = models.CharField(max_length=15)
    indication = models.CharField(max_length=150, null=True)
    approval = models.CharField(max_length=5)
    clinicalstatus = models.CharField(max_length=30, null=True, default=False)
    phasedate = models.CharField(max_length=15, null=True, default=True)
    phase = models.CharField(max_length=1, null=True, default=True)
    moa = models.CharField(max_length=30, null=True, default=True)
    targetlevel = models.CharField(max_length=10, null=True, default='NA')
    externallink = models.CharField(max_length=150, null=True) # Link Framework
    novelty = models.CharField(max_length=15) #Boolean
    references = models.CharField(max_length=180) #Boolean
    publication = models.ManyToManyField('common.Publication')


    def __str__(self):
        return self.name

    class Meta():
        db_table = 'drugs'
