from django.db import models

class Drugs(models.Model):
    target = models.ForeignKey('protein.Protein', on_delete=models.CASCADE, null=True)
    indication = models.ForeignKey('Indication', on_delete=models.CASCADE, null=True)
    ligand = models.ForeignKey('ligand.Ligand', on_delete=models.CASCADE, null=True)
    charge = models.CharField(max_length=5, null=True)
    complexity = models.FloatField(max_length=4, null=True)
    tpsa = models.CharField(max_length=10, null=True)
    drug_status = models.CharField(max_length=15, null=True)
    approval_year = models.IntegerField(null=True)
    indication_max_phase = models.IntegerField(null=True)
    indication_status = models.CharField(max_length=15, null=True)
    moa = models.ForeignKey('ligand.LigandRole', on_delete=models.CASCADE, null=True)
    genetic_association = models.CharField(max_length=30, null=True)
    affected_pathway = models.CharField(max_length=30, null=True)
    somatic_mutation = models.CharField(max_length=30, null=True)
    similarity_to_model = models.FloatField(max_length=4, null=True)
    novelty_score = models.FloatField(max_length=4, null=True)
    publication_count = models.IntegerField(null=True)
    target_level = models.CharField(max_length=30, null=True)
    disease_association = models.ForeignKey('IndicationAssociation', on_delete=models.CASCADE, null=True)
    reference = models.ManyToManyField('common.Publication')


    def __str__(self):
        return self.ligand.name

    class Meta():
        db_table = 'drugs'

#Do we need to fix this model structure?
class Indication(models.Model):
    title =  models.CharField(max_length=255)
    code =  models.CharField(max_length=8, null=True)
    slug = models.CharField(max_length=255, unique=True, null=True)
    uri =  models.ForeignKey('common.WebLink', on_delete=models.CASCADE, null=True)
    parent = models.ForeignKey('self', on_delete=models.CASCADE, null=True)
    level = models.IntegerField(null=True)

    def get_level_0(self):
        current = self
        while current.parent is not None:
            current = current.parent
        return current

    def get_level_1(self):
        current = self
        while current.parent.parent is not None:
            current = current.parent
        return current

    def __str__(self):
        return self.title

    class Meta():
        db_table = 'indication'

class IndicationAssociation(models.Model):
    target = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    indication = models.ForeignKey('Indication', on_delete=models.CASCADE)
    association_score = models.FloatField(max_length=7, null=True)
    ot_genetics_portal = models.FloatField(max_length=7, null=True)
    gene_burden = models.FloatField(max_length=7, null=True)
    clingen = models.FloatField(max_length=7, null=True)
    gene2phenotype = models.FloatField(max_length=7, null=True)
    orphanet = models.FloatField(max_length=7, null=True)
    genomics_england = models.FloatField(max_length=7, null=True)
    uniprot_literature = models.FloatField(max_length=7, null=True)
    uniprot_variants = models.FloatField(max_length=7, null=True)
    sysbio = models.FloatField(max_length=7, null=True)
    cancer_gene_census = models.FloatField(max_length=7, null=True)
    cancer_biomarkers = models.FloatField(max_length=7, null=True)
    intogen = models.FloatField(max_length=7, null=True)
    eva = models.FloatField(max_length=7, null=True)
    eva_somatic = models.FloatField(max_length=7, null=True)
    chembl = models.FloatField(max_length=7, null=True)
    slapenrich = models.FloatField(max_length=7, null=True)
    crispr = models.FloatField(max_length=7, null=True)
    crispr_screen = models.FloatField(max_length=7, null=True)
    reactome = models.FloatField(max_length=7, null=True)
    europepmc = models.FloatField(max_length=7, null=True)
    expression_atlas = models.FloatField(max_length=7, null=True)
    impc = models.FloatField(max_length=7, null=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'indication_association'

class ATCCodes(models.Model):
    ligand =  models.ForeignKey('ligand.Ligand', on_delete=models.CASCADE, null=True)
    code =  models.ForeignKey('common.WebLink', on_delete=models.CASCADE, null=True)
    level_0 = models.CharField(max_length=1, null=True)
    level_1 = models.CharField(max_length=3, null=True)
    name_0 = models.CharField(max_length=100, null=True)
    name_1 = models.CharField(max_length=100, null=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'atc_codes'
