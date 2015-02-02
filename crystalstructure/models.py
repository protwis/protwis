from django.db import models

class CrystalStructure(models.Model):
    protein_class = models.ForeignKey('protein.ProteinFamily')
    receptor = models.ForeignKey('protein.Protein')
    pdb_code = models.CharField(max_length=4)
    preferred_chain = models.CharField(max_length=20)
    resolution = models.DecimalField(max_digits=5, decimal_places=3)
    endogenous_ligand = models.ForeignKey('ligand.Ligand', related_name='endogenous_ligand', null=True)
    xray_ligand = models.ForeignKey('ligand.Ligand', related_name='xray_ligand', null=True)
    pmid = models.CharField(max_length=20)
    publication_date = models.DateField()
    n_terminus = models.TextField()
    icl1 = models.TextField()
    ecl1 = models.TextField()
    icl2 = models.TextField()
    ecl2_1 = models.TextField()
    ecl2_2 = models.TextField()
    icl3 = models.TextField()
    ecl3 = models.TextField()
    c_terminus = models.TextField()
    
    class Meta():
        db_table = 'crystalstructure'
    