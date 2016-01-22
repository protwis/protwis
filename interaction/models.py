from django.db import models


class ResidueFragmentInteraction(models.Model):

    structure_ligand_pair = models.ForeignKey('StructureLigandInteraction')
    rotamer = models.ForeignKey('structure.Rotamer')
    fragment = models.ForeignKey('structure.Fragment')
    interaction_type = models.ForeignKey('ResidueFragmentInteractionType')

    def __str__(self):
        if self.rotamer.residue.display_generic_number is not None:
            return "{!s} {!s} {!s}".format(self.structure_ligand_pair.structure.pdb_code.index, self.rotamer.residue.display_generic_number.label, self.structure_ligand_pair.ligand.name)
        else:
            return "{!s} {!s} {!s}".format(self.structure_ligand_pair.structure.pdb_code.index, self.rotamer.residue, self.structure_ligand_pair.ligand.name)
    class Meta():
        db_table = 'interaction_residue_fragment'

    def get_pdbdata(self):
        return "{!s}\n{!s}".format(self.rotamer.pdbdata, self.fragment.pdbdata)


    def generate_filename(self):

        if self.rotamer.residue.display_generic_number is not None:
            generic_num = self.rotamer.residue.display_generic_number.label
        else:
            generic_num = self.rotamer.residue.sequence_number
        res_name = self.rotamer.residue.amino_acid
        prot_entry_name = str(self.structure_ligand_pair.structure.protein_conformation.protein.parent.entry_name)
        pdb_code = self.structure_ligand_pair.structure.pdb_code.index
        interaction = self.interaction_type.slug

        return "{}_{}_{}_{}_{}.pdb".format(generic_num.replace('.','_'), res_name, prot_entry_name, pdb_code, interaction)


class ResidueFragmentInteractionType(models.Model):
    slug = models.SlugField(max_length=40)
    name = models.CharField(max_length=100)
    type = models.CharField(max_length=50, null=True)
    direction = models.CharField(max_length=30, null=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'interaction_type_residue_fragment'


class StructureLigandInteraction(models.Model):
    structure = models.ForeignKey('structure.Structure')
    ligand = models.ForeignKey('ligand.Ligand')
    ligand_role = models.ForeignKey('ligand.LigandRole')
    pdb_reference = models.CharField(max_length=3, null=True)
    pdb_file = models.ForeignKey('structure.PdbData', null=True)
    annotated = models.BooleanField(default=False)

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
