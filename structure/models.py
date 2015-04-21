from django.db import models


class Structure(models.Model):
    protein = models.ForeignKey('protein.Protein')
    structure_type = models.ForeignKey('StructureType')
    pdb_code = models.ForeignKey('common.WebLink')
    state = models.ForeignKey('StructureState')
    ligands = models.ManyToManyField('ligand.Ligand', through='interaction.StructureLigandInteraction')
    protein_anomalies = models.ManyToManyField('protein.ProteinAnomaly')
    stabilizing_agents = models.ManyToManyField('StructureStabilizingAgent')
    preferred_chain = models.CharField(max_length=20)
    resolution = models.DecimalField(max_digits=5, decimal_places=3)
    publication = models.ForeignKey('common.Publication', null=True)
    publication_date = models.DateField()

    def __str__(self):
        return self.pdb_code.index

    class Meta():
        db_table = 'structure'


class StructureModel(models.Model):
    protein = models.ForeignKey('protein.Protein')

    def __str__(self):
        return self.protein.entry_name

    class Meta():
        db_table = 'structure_model'


class StructureType(models.Model):
    slug = models.SlugField(max_length=20)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = "structure_type"


class StructureState(models.Model):
    slug = models.SlugField(max_length=20)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = "structure_state"


class StructureStabilizingAgent(models.Model):
    slug = models.SlugField(max_length=20)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = "structure_stabilizing_agent"
