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

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'drugs'


class Drugs2024(models.Model):
    target = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    indication = models.ForeignKey('Indication', on_delete=models.CASCADE)
    ligand = models.ForeignKey('ligand.Ligand', on_delete=models.CASCADE)
    charge = models.CharField(max_length=2)
    complexity = models.CharField(max_length=10)
    tpsa = models.CharField(max_length=10)
    drug_status = models.CharField(max_length=15)
    approval_year = models.CharField(max_length=5)
    indication_max_phase = models.CharField(max_length=1)
    moa = models.CharField(max_length=30, null=True)
    affected_pathway = models.CharField(max_length=30)
    somatic_mutation = models.CharField(max_length=30)
    similarity_to_model = models.CharField(max_length=30)
    novelty_score = models.CharField(max_length=5)
    publication = models.ManyToManyField('common.Publication')

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'drugs'

class Indication(models.Model):
    name =  models.CharField(max_length=55)
    code =  models.CharField(max_length=10)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'indication'
