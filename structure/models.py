from django.db import models

class Structure(models.Model):
    protein = models.ForeignKey('protein.Protein')
    structure_type = models.ForeignKey('StructureType')
    pdb_code = models.ForeignKey('common.WebLink')
    preferred_chain = models.CharField(max_length=20)
    resolution = models.DecimalField(max_digits=5, decimal_places=3)
    endogenous_ligand = models.ForeignKey('ligand.Ligand', related_name='endogenous_ligand', null=True)
    xray_ligand = models.ForeignKey('ligand.Ligand', related_name='xray_ligand', null=True)
    publication = models.CharField(max_length=20)
    pdb_publication_date = models.DateField()
    stabilizing_agents = models.ManyToManyField('StructureStabilizingAgent', null=True)

    class Meta():
        db_table = 'structure'


class StructureType(models.Model):
    slug = models.CharField(max_length=20)
    description = models.TextField()

    class Meta():
        db_table = "structure_type"

class StructureStabilizingAgnet(models.Model):
    slug = models.CharField(max_length=20)

    class Meta():
        db_table = "structure_stabilizing_agent"