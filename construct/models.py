from django.db import models

from protein.models import Protein
from residue.models import Residue
from mutation.models import Mutation


class ConstructMutation(models.Model):
    construct = models.ForeignKey('protein.ProteinConformation')
    sequence_number = models.SmallIntegerField()
    wild_type_amino_acid = models.CharField(max_length=1)
    mutated_amino_acid = models.CharField(max_length=1)

    def __str__(self):
        return '{} {}{}{}'.format(self.construct.protein.entry_name, self.wild_type_amino_acid, self.sequence_number,
            self.mutated_amino_acid)

    class Meta():
        db_table = 'construct_mutation'


class ConstructDeletion(models.Model):
    construct = models.ForeignKey('protein.ProteinConformation')
    start = models.SmallIntegerField()
    end = models.SmallIntegerField()

    def __str__(self):
        return '{} {}-{}'.format(self.construct.protein.entry_name, self.start, self.end)

    class Meta():
        db_table = 'construct_deletion'


class ConstructInsertion(models.Model):
    construct = models.ForeignKey('protein.ProteinConformation')
    protein_type = models.ForeignKey('ConstructInsertionType')
    name = models.CharField(max_length=100, null=True)
    uniprot_id = models.CharField(max_length=20)
    sequence = models.TextField(null=True)
    deletions = models.TextField(max_length=100, null=True)
    position = models.TextField(max_length=20, null=True)
    remarks = models.TextField(null=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'construct_aux_protein'


class ConstructInsertionType(models.Model):
    # values like FUsion, tag, linkers,signal_peptide,prtn_Cleavage_site, aa_modfcn (comma separated data)
    name = models.TextField(max_length=50, null=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'construct_aux_protein_type'


class Chemical(models.Model):
    chemical_type = models.ForeignKey('ChemicalType')
    name = models.CharField(max_length=200, null=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'construct_chemical'


class ChemicalType(models.Model):
    # type like detergent, lipid
    name = models.CharField(max_length=100, null=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'construct_chemical_type'


class ChemicalConc(models.Model):
    chemical = models.ForeignKey('Chemical')
    concentration = models.TextField(null=True)

    def __str__(self):
        return self.chemical.name

    class Meta():
        db_table = 'construct_chemical_conc'


# includes all chemicals, type & concentration. Chemicals can be from LCPlipid, detergent, etc.
class ChemicalList(models.Model):
    chemicals = models.ManyToManyField('ChemicalConc')

    def __str__(self):
        return self.chemicals.concentration

    class Meta():
        db_table = 'construct_chemical_list'


class Expression(models.Model):
    # since Expression construct can vary from  crystallization construct so the experiments refer to construct
    # instead
    construct = models.ForeignKey('protein.ProteinConformation')
    expression_system = models.ForeignKey('ExpressionSystem')
    remarks = models.TextField(null=True)

    def __str__(self):

        return self.construct + ' (' + expression_sys + ')'

    class Meta():
        db_table = 'construct_expression'


class ExpressionSystem (models.Model):
    expression_method = models.CharField(max_length=100)
    host_cell_type = models.CharField(max_length=100)
    host_cell = models.CharField(max_length=100)

    def __str__(self):

        return self.host_cell

    class Meta():
        db_table = 'construct_expression_system'



class Purification(models.Model):
    construct = models.ForeignKey('protein.ProteinConformation')
    remarks = models.TextField(null=True)

    def __str__(self):
        return self.construct.protein_conformation.protein.entry_name

    class Meta():
        db_table = 'construct_purification'


class PurificationStep(models.Model):
    purification = models.ForeignKey('Purification')
    purification_type = models.ForeignKey('PurificationStepType')
    # IMAC/Flag tag/1D4/Alprenolol/Strep tag/Ligand affinity
    description = models.TextField(null=True)

    def __str__(self):

        return self.description

    class Meta():
        db_table = 'construct_purification_step'


class PurificationStepType(models.Model):
    # Chromatography/Enzyme modification/Chem modification IMAC/Flag tag/1D4/Alprenolol/Strep tag/Ligand affinity
    name = models.TextField(null=True)

    def __str__(self):

        return self.name

    class Meta():
        db_table = 'construct_purification_step_type'


class Solubilization(models.Model):
    construct = models.ForeignKey('protein.ProteinConformation')
    # includes chem name, type[detergent & solubilisation_lipid] and concentration
    chemical_list = models.ForeignKey('ChemicalList')
    remarks = models.TextField(null=True)

    def __str__(self):
        return self.chemical_list.chemicals.concentration

    class Meta():
        db_table = 'construct_solubilization'


class ChemicalModification(models.Model):
    construct_solubilization = models.ForeignKey('Solubilization')
    description =  models.TextField()

    def __str__(self):
        return self.description

    class Meta():
        db_table = 'construct_chemical_modification'


class Crystallization(models.Model):
    construct = models.ForeignKey('protein.ProteinConformation')
    crystal_type = models.ForeignKey('CrystallizationMethodTypes', null=True)
    # chemical type= LCPlipid for LCP exp, else type= lipid for in surfo exp includes info on lcp_lipidic_condition,
    # protein_component etc.
    chemical_list =  models.ForeignKey('ChemicalList')
    method = models.TextField(max_length=100,null=True)
    settings = models.TextField(max_length=100,null=True)
    remarks = models.TextField(max_length=1000, null=True)
    protein_conc =  models.SlugField(max_length=20,blank=True)
    ligands = models.ManyToManyField('ligand.Ligand', through='CrystallizationLigandConc')
    aqueous_solution_lipid_ratio=  models.SlugField(max_length=20, null=True)
    lcp_bolus_volume =  models.SlugField(max_length=20,null=True)
    precipitant_solution_volume =  models.SlugField(max_length=20, null=True)
    temp = models.CharField(max_length=5, null=True)
    ph = models.TextField(max_length=10, null=True)

    def __str__(self):
        return self.method

    class Meta():
        db_table = 'construct_crystallization'


class CrystallizationLigandConc(models.Model):
    construct_crystallization = models.ForeignKey('Crystallization')
    ligand = models.ForeignKey('ligand.Ligand')
    ligand_conc = models.TextField(null=True)

    def __str__(self):
        return self.ligand_conc

    class Meta():
        db_table = 'construct_crystallization_ligand_conc'


class CrystallizationMethodTypes(models.Model):
    # LCP/ in surfo/bicelles
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'construct_crystallization_method_types'
