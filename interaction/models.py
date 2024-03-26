from django.db import models


class ResidueFragmentInteraction(models.Model):

    structure_ligand_pair = models.ForeignKey('StructureLigandInteraction', on_delete=models.CASCADE)
    rotamer = models.ForeignKey('structure.Rotamer', on_delete=models.CASCADE)
    fragment = models.ForeignKey('structure.Fragment', on_delete=models.CASCADE)
    interaction_type = models.ForeignKey('ResidueFragmentInteractionType', on_delete=models.CASCADE)

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
    slug = models.SlugField(max_length=40, unique=True)
    name = models.CharField(max_length=100)
    type = models.CharField(max_length=50, null=True)
    direction = models.CharField(max_length=30, null=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'interaction_type_residue_fragment'


class StructureLigandInteraction(models.Model):
    structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE, null=True)
    ligand = models.ForeignKey('ligand.Ligand', on_delete=models.CASCADE)
    ligand_role = models.ForeignKey('ligand.LigandRole', on_delete=models.CASCADE)
    pdb_reference = models.CharField(max_length=3, null=True)
    pdb_file = models.ForeignKey('structure.PdbData', null=True, on_delete=models.CASCADE)
    annotated = models.BooleanField(default=False)

    def __str__(self):
        return "{} {}".format(self.structure.pdb_code, self.ligand.name)

    class Meta():
        db_table = 'interaction_structure_ligand'
