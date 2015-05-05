from django.db import models


class ResidueFragmentInteraction(models.Model):

    structure = models.ForeignKey('structure.Structure') #redundant via rotamer
    residue = models.ForeignKey('residue.Residue') #redundant via rotamer
    rotamer = models.ForeignKey('structure.Rotamer')
    fragment = models.ForeignKey('structure.Fragment')
    ligand = models.ForeignKey('ligand.Ligand') #redundant via fragment
    interaction_type = models.ForeignKey('ResidueFragmentInteractionType')

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
    pdb_file = models.ForeignKey('structure.PdbData')

    def __str__(self):
        return "{} {}".format(self.structure.pdb_code, self.ligand.name)

    class Meta():
        db_table = 'interaction_structure_ligand'


class ProteinLigandInteraction(models.Model):
    protein = models.ForeignKey('protein.ProteinConformation')
    ligand = models.ForeignKey('ligand.Ligand')

    def __str__(self):
        return "{} {}".format(self.protein.entry_name, self.ligand.name)

    class Meta():
        db_table = 'interaction_protein_ligand'

class ResidueFragmentAtom(models.Model):
    structureligandpair = models.ForeignKey('StructureLigandInteraction')
    interaction = models.ForeignKey('ResidueFragmentInteraction', null=True)
    atomtype = models.CharField(max_length = 20)
    atomnr = models.SmallIntegerField()
    atomclass = models.CharField(max_length = 20)
    residuename = models.CharField(max_length = 20)
    chain = models.CharField(max_length = 20)
    residuenr = models.SmallIntegerField()
    x = models.DecimalField(max_digits=6, decimal_places=3)
    y = models.DecimalField(max_digits=6, decimal_places=3)
    z = models.DecimalField(max_digits=6, decimal_places=3)
    occupancy = models.DecimalField(max_digits=6, decimal_places=2)
    temperature = models.DecimalField(max_digits=6, decimal_places=2)
    element_name = models.CharField(max_length = 20)


    class Meta():
        db_table = 'interaction_residue_fragment_atoms'

    def __str__(self):
        return self.atomtype
