from django.db import models


class ResidueFragmentInteraction(models.Model):
    structure = models.ForeignKey('structure.Structure')
    residue = models.ForeignKey('residue.Residue')
    ligand = models.ForeignKey('ligand.Ligand')
    interaction_type = models.ForeignKey('ResidueFragmentInteractionType')
    fragment = models.TextField()

    def __str__(self):
        return "{} {} {}".format(self.structure.pdb_code, self.residue.generic_number.label, self.ligand.name)

    class Meta():
        db_table = 'interaction_residue_fragment'


class ResidueFragmentInteractionType(models.Model):
    slug = models.SlugField(max_length=10)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'interaction_type_residue_fragment'


class StructureLigandInteraction(models.Model):
    structure = models.ForeignKey('structure.Structure')
    ligand = models.ForeignKey('ligand.Ligand')
    ligand_role = models.ForeignKey('ligand.LigandRole')
    pdb_reference = models.CharField(max_length=3, null=True)

    def __str__(self):
        return "{} {}".format(self.structure.pdb_code, self.ligand.name)

    class Meta():
        db_table = 'interaction_structure_ligand'


class ProteinLigandInteraction(models.Model):
    protein = models.ForeignKey('protein.Protein')
    ligand = models.ForeignKey('ligand.Ligand')

    def __str__(self):
        return "{} {}".format(self.protein.entry_name, self.ligand.name)

    class Meta():
        db_table = 'interaction_protein_ligand'