from django.db import models

from protein.models import Protein
from residue.models import Residue
from mutation.models import Mutation 


class Construct(models.Model):

    parent = models.ForeignKey('self', null=True)
    protein_conformation = models.ForeignKey('protein.ProteinConformation') 

    deletions = models.TextField(max_length=100, null=True)

    def __str__(self):
        return self.protein_conformation.protein.entry_name
    
    class Meta():
        db_table = 'construct'


class AuxProtein(models.Model):
    
    construct = models.ForeignKey('Construct')
    protein_type = models.ForeignKey('AuxProteinType')
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


class AuxProteinType(models.Model):

    name = models.TextField(max_length=50, null=True) #values like FUsion, tag, linkers,signal_peptide,prtn_Cleavage_site, aa_modfcn  #comma separated data

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

    name = models.CharField(max_length=100, null=True) # type like detergent, lipid 

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


class ChemicalList(models.Model):   ## includes all chemicals, type & concentration. Chemicals can be from LCPlipid, detergent, etc.

    chemicals = models.ManyToManyField('ChemicalConc')        

    def __str__(self):
        return self.chemicals.concentration

    class Meta():
        db_table = 'construct_chemical_list'


class ChemicalModification(models.Model):
    
    construct_solubilization = models.ForeignKey('ConstructSolubilization')
    description =  models.TextField()

    def __str__(self):
        return self.description

    class Meta():
        db_table = 'construct_chemical_modification'


class PurificationStep(models.Model):
   
    purification = models.ForeignKey('ConstructPurification')
    purification_type = models.ForeignKey('PurificationStepType')
    description = models.TextField(null=True) # IMAC/Flag tag/1D4/Alprenolol/Strep tag/Ligand affinity


    def __str__(self):            

        return self.description

    class Meta():
        db_table = 'construct_purification_step'



class PurificationStepType(models.Model):

    name = models.TextField(null=True) # Chromatography/Enzyme modification/Chem modification    IMAC/Flag tag/1D4/Alprenolol/Strep tag/Ligand affinity

    def __str__(self):            

        return self.name

    class Meta():
        db_table = 'construct_purification_step_type'




class ConstructExpression(models.Model):

    construct = models.ForeignKey('Construct') #since Expression construct can vary from  crystallization construct so the experiments refer to construct instead
    expression_system = models.ForeignKey('ConstructExpressionSystem')
    remarks = models.TextField(null=True)

    def __str__(self):            

        return self.construct + ' (' + expression_sys + ')'

    class Meta():
        db_table = 'construct_expression'


class ConstructExpressionSystem (models.Model):
     
    expression_method = models.CharField(max_length=100)
    host_cell_type = models.CharField(max_length=100)
    host_cell = models.CharField(max_length=100)

    def __str__(self):            

        return self.host_cell

    class Meta():
        db_table = 'construct_expression_system'


class ConstructSolubilization(models.Model):

    construct = models.ForeignKey('Construct')
    chemical_list = models.ForeignKey('ChemicalList') #includes chem name, type[detergent & solubilisation_lipid]  and concentration
    remarks = models.TextField(null=True)

    def __str__(self):
        return self.chemical_list.chemicals.concentration

    class Meta():
        db_table = 'construct_solubilization'


class ConstructPurification(models.Model):

    construct = models.ForeignKey('Construct')
    remarks = models.TextField(null=True)

    def __str__(self):
        return self.construct.protein_conformation.protein.entry_name

    class Meta():
        db_table = 'construct_purification'


class ConstructCrystallization(models.Model):

    construct = models.ForeignKey('Construct')
    crystal_type = models.ForeignKey('CrystallizationMethodTypes') 
    chemical_list =  models.ForeignKey('ChemicalList') #chemical type= LCPlipid for LCP exp, else type= lipid for in surfo exp --- includes info on lcp_lipidic_condition, protein_component etc.
    method = models.TextField(max_length=100,null=True) 
    settings = models.TextField(max_length=100,null=True)
    remarks = models.TextField(max_length=1000, null=True)
    protein_conc =  models.SlugField(max_length=20,blank=True)
    ligands = models.ManyToManyField('ligand.Ligand', through='ConstructCrystallizationLigandConc')
    aqueous_solution_lipid_ratio=  models.SlugField(max_length=20,null=True)
    lcp_bolus_volume =  models.SlugField(max_length=20,null=True)
    precipitant_solution_volume =  models.SlugField(max_length=20,null=True)
    temp = models.CharField(max_length=5, null=True)
    ph = models.TextField(max_length=10, null=True)

    def __str__(self):
        return self.method

    class Meta():
        db_table = 'construct_crystallization'


class ConstructCrystallizationLigandConc(models.Model):

    construct_crystallization = models.ForeignKey('ConstructCrystallization')
    ligand = models.ForeignKey('ligand.Ligand')
    ligand_conc = models.TextField(null=True)
    
    def __str__(self):
        return self.ligand_conc

    class Meta():
        db_table = 'construct_ligand_conc_of_crystallization'

class CrystallizationMethodTypes(models.Model):

    name = models.CharField(max_length=100,null=True) #LCP/ in surfo/bicelles ... can be null if that construct did not crystallize

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'construct_crystallization_method_types'  


