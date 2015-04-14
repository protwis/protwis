from django.db import models


class Ligand(models.Model):
    web_link = models.ManyToManyField('common.WebLink')
    ligand_type = models.ForeignKey('LigandType')
    name = models.TextField()
    smiles = models.TextField(null=True)
    inchikey = models.CharField(max_length=50, null=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'ligand'


class LigandType(models.Model):
    slug = models.SlugField(max_length=20)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name
    
    class Meta():
        db_table = 'ligand_type'


class LigandAlias(models.Model):
    ligand = models.ForeignKey('Ligand')
    name = models.TextField()

    def __str__(self):
        return self.name
    
    class Meta():
        db_table = 'ligand_alias'


class LigandRole(models.Model):
    slug = models.SlugField(max_length=10)
    name = models.CharField(max_length=100)
    
    def __str__(self):
        return self.name
    
    class Meta():
        db_table = 'ligand_role'