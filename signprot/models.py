from django.db import models
from protein.models import Protein, ProteinConformation
from structure.models import StructureStabilizingAgent, PdbData
from common.models import WebLink, Publication


class SignprotStructure(models.Model):
    protein = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    pdb_code = models.ForeignKey('common.WebLink', on_delete=models.CASCADE)
    structure_type = models.ForeignKey('structure.StructureType', on_delete=models.CASCADE)
    publication_date = models.DateField()
    publication = models.ForeignKey('common.Publication', null=True, on_delete=models.CASCADE)
    stabilizing_agents = models.ManyToManyField('structure.StructureStabilizingAgent')
    resolution = models.DecimalField(max_digits=5, decimal_places=3)
    pdb_data = models.ForeignKey('structure.PdbData', null=True, on_delete=models.CASCADE)

    def __str__(self):
        return self.pdb_code.index

    def get_cleaned_pdb(self, remove_waters=True, ligands_to_keep=None, remove_aux=False, aux_range=5.0):

        tmp = []
        for line in self.pdb_data.pdb.split('\n'):
            save_line = True
            if remove_waters and line.startswith('HET') and line[17:20] == 'HOH':
                save_line = False
            if ligands_to_keep and line.startswith('HET'):
                if line[17:20] != 'HOH' and line[17:20] in ligands_to_keep:
                    save_line = True
                elif line[17:20] != 'HOH':
                    save_line=False
            if save_line:
                tmp.append(line)

        return '\n'.join(tmp)

    class Meta():
        db_table = 'signprot_structure'

class SignprotStructureExtraProteins(models.Model):
    structure = models.ForeignKey('SignprotStructure', on_delete=models.CASCADE, null=True, related_name='extra_proteins')
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
