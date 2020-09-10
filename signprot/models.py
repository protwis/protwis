from django.db import models
from protein.models import Protein, ProteinConformation
from structure.models import Structure, StructureType, StructureExtraProteins, StructureStabilizingAgent
from common.models import WebLink, Publication


class SignprotStructure(models.Model):
    protein = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    pdb_code = models.ForeignKey('common.WebLink', on_delete=models.CASCADE)
    extra_proteins = models.ManyToManyField('SignprotStructureExtraProteins', related_name='extra_proteins')
    structure_type = models.ForeignKey('structure.StructureType', on_delete=models.CASCADE)
    publication_date = models.DateField()
    publication = models.ForeignKey('common.Publication', null=True, on_delete=models.CASCADE)
    stabilizing_agents = models.ManyToManyField('structure.StructureStabilizingAgent')
    resolution = models.DecimalField(max_digits=5, decimal_places=3)

    def __str__(self):
        return self.pdb_code.index

    class Meta():
        db_table = 'signprot_structure'


class SignprotStructureExtraProteins(models.Model):
    structure = models.ForeignKey('SignprotStructure', on_delete=models.CASCADE, null=True)
    wt_protein = models.ForeignKey('protein.Protein', on_delete=models.CASCADE, null=True)
    protein_conformation = models.ForeignKey('protein.ProteinConformation', on_delete=models.CASCADE, null=True)
    display_name = models.CharField(max_length=20)
    note = models.CharField(max_length=50, null=True)
    chain = models.CharField(max_length=1)
    category = models.CharField(max_length=20)
    wt_coverage = models.IntegerField(null=True)

    def __str__(self):
        return self.display_name

    class Meta():
        db_table = "signprot_extra_proteins"


class SignprotBarcode(models.Model):
    protein = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    residue = models.ForeignKey('residue.Residue', on_delete=models.CASCADE)
    seq_similarity = models.DecimalField(max_digits=5, decimal_places=4)
    seq_identity = models.DecimalField(max_digits=5, decimal_places=4)
    paralog_score = models.DecimalField(max_digits=5, decimal_places=4)

    def __str__(self):
        return self.protein.name

    class Meta():
        db_table = 'gprotein_barcode'


class SignprotComplex(models.Model):
    protein = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    alpha = models.CharField(max_length=1)
    beta_chain = models.CharField(max_length=1, null=True)
    beta_protein = models.ForeignKey('protein.Protein', related_name='beta_protein', on_delete=models.CASCADE, null=True)
    gamma_chain = models.CharField(max_length=1, null=True)
    gamma_protein = models.ForeignKey('protein.Protein', related_name='gamma_protein', on_delete=models.CASCADE, null=True)

    def __str__(self):
        return '<SignprotComplex: {} {}>'.format(self.protein.entry_name, self.structure.pdb_code.index)

    class Meta():
        db_table = 'signprot_complex'


class SignprotInteractions(models.Model):
    gpcr_residue = models.ForeignKey('residue.Residue', on_delete=models.CASCADE, related_name='gpcr_residue')
    signprot_residue = models.ForeignKey('residue.Residue', on_delete=models.CASCADE, related_name='signprot_residue')
    interaction_type = models.CharField(max_length=200)
    structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)

    def __str__(self):
        return '{} between {} and {}'.format(self.interaction_type, self.gpcr_residue, self.signprot_residue)
