from django.db import models

from protein.models import Protein
from residue.models import Residue
from mutation.models import Mutation ## check ??

# include models for Construct, prtn_order, exprssn, solubilisation, purifcn, and crystallisation 

class Construct(models.Model):

    parent = models.ForeignKey('self', null=True)
    protein = models.ForeignKey('protein.Protein') 
#    name = models.CharField('protein.Protein',max_length=200) # can get details like name,seq,pdbcode, ligand etc. from protein model
    mutation = models.TextField(max_length=100) # how to show [position & reasons]?  ###comma separated data for listing many mutations
####    mutation_to = models.ForeignKey('mutation.amino_acid') ## linked from Protein model
####    mutation_reason = models.ForeignKey('mutation.type')

    deletion = models.TextField(max_length=100)  #comma separated data
    aux_protein = models.ForeignKey('AuxProtein') 
#   expression = models.ForeignKey('ConstructExpression') 
 #  solubilization = models.ForeignKey('ConstructSolubilization')
  # purification = models.ForeignKey('ConstructPurification')
   #crystallization = models.ForeignKey('ConstructCrystallization')## NOTE: If you need to create a relationship on a model that has not yet been defined, you can use the name of the model, rather than the model object itself

     def __str__(self):
        return self.protein.slug
    
    class Meta():
        db_table = 'construct'


class AuxProtein(models.Model):

    ptn_type = models.ForeignKey('AuxProteinType')
    name = models.CharField(max_length=100, blank=True)
    sequence = models.TextField(blank=True)
    mutation = models.TextField(max_length=100, blank=True) #comma separated data
####    mutation_to = models.ForeignKey('mutation.mutation_to')
####    mutation_reason = models.ForeignKey('mutation.MutationType')
    deletion = models.TextField(max_length=100, blank=True) #comma separated data
    position = models.TextField(max_length=20, blank=True) # aa_modfcn like position of disulfide bond is imp to show (is it shown already in GPCRdb Sites page??) For 3 linkers, values will be 1,2,3
    remarks = models.TextField()

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'aux_protein'


class AuxProteinType(models.Model):

    name = models.TextField(max_length=50) #values like FUsion, tag, linkers,signal_peptide,prtn_Cleavage_site, aa_modfcn  #comma separated data

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'aux_protein_type'


# how to define protein order field?? Input as string of types separated by ":" ? and checked via ptn_type field  ::: NOt needed to define in models.py as it will used only by "build.py" script while making the order of construct

#==============================================================================


class Chemical(models.Model):

    chemical_type = models.ForeignKey('ChemicalType')
    name = models.CharField(max_length=200)

    def __str__(self):
        return self.name
    
    class Meta():
        db_table = 'chemical'


class ChemicalType(models.Model):

    name = models.CharField(max_length=100) # type like detergent, lipid  ## If there is no type present in database then Django creates another type

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'chemical_type'


class ChemicalConc(models.Model):
  
    chemical = models.ForeignKey('Chemical')
    concentration = models.TextField(null=True)

    def __str__(self):
        return self.chemical.name

    class Meta():
        db_table = 'chemical_conc'


class ChemicalList(models.Model):   ## includes all chemicals, type & concentration. Chemicals can be from LCPlipid, detergent, etc.

    chemicals = models.ManyToManyField('ChemicalConc')        

    def __str__(self):
        return self.chemicals.concentration

    class Meta():
        db_table = 'chemical_list'


class ChemicalModification(models.Model):  # can have more than one chem_modfcn as comma separated fields

    description =  models.TextField()

    def __str__(self):
        return self.description

    class Meta():
        db_table = 'chemical_modification'


class EnzymeModification(models.Model): # can have more than one enzymatic_modfcn as comma separated fields
    description =  models.TextField()

    def __str__(self):
        return self.description

    class Meta():
        db_table = 'enzyme_modification'


class ChromatographyType(models.Model):

    description = models.TextField() # can have more than one chromatography method as comma separated fields

    def __str__(self):
        return self.description

    class Meta():
        db_table = 'chromatography_method'


class ConstructExpression(models.Model):

    construct_sequence = models.ForeignKey('Construct') #since Expression construct can vary from  crystallization construct so the experiments refer to construct instead
    expression_method = models.CharField(max_length=100)
    host_cell_type = models.CharField(max_length=100)
    host_cell = models.CharField(max_length=100)
    remarks = models.TextField()

    def __str__(self):
        return self.host_cell

    class Meta():
        db_table = 'construct_expression'


class ConstructSolubilization(models.Model):

    construct_sequence = models.ForeignKey('Construct')
    chemical_list = models.ForeignKey('ChemicalList') #includes chem name, type[detergent & solubilisation_lipid]  and concentration
    chemical_modification = models.ForeignKey('ChemicalModification')#more than one chem_modfcn ## remove field, not easily searchable?
    remarks = models.TextField()

    def __str__(self):
        return self.chemical_list.chemicals.concentration

    class Meta():
        db_table = 'construct_solubilization'


class ConstructPurification(models.Model):

    construct_sequence = models.ForeignKey('Construct')
    chromatography_type = models.ForeignKey('ChromatographyType') 
    enzyme_modification = models.ForeignKey('EnzymeModification')
    chemical_modification = models.ForeignKey('ChemicalModification')
    remarks = models.TextField()

    def __str__(self):
        return self.chromatography_type.description

    class Meta():
        db_table = 'construct_purification'


class ConstructCrystallization(models.Model):

    construct_sequence = models.ForeignKey('Construct')
    crystal_type = models.ForeignKey('CrystallizationMethodTypes') # if Many crystallography types for one construct :named as separate Xtal Exp
    chemical_list =  models.ForeignKey('ChemicalList') #chemical type= LCPlipid for LCP exp, else type= lipid for in surfo exp --- includes info on lcp_lipidic_condition, protein_component etc.
    method = models.TextField(max_length=100) # more than one solved by comma separated
    settings = models.TextField(max_length=100)# can have more than one
    remarks = models.TextField()
    protein_conc =  models.SlugField(max_length=20,blank=True)
    ligands = models.ManyToManyField('ligand.Ligand', through='ConstructCrystallizationLigandConc')
    aqueous_solution_lipid_ratio=  models.SlugField(max_length=20,blank=True)
    lcp_bolus_volume =  models.SlugField(max_length=20,blank=True)
    precipitant_solution_volume =  models.SlugField(max_length=20,blank=True)
#   condition = models.ForeignKey('XtalCondition')

    def __str__(self):
        return self.method

    class Meta():
        db_table = 'construct_crystallization'


class ConstructCrystallizationLigandConc(models.Model):

    ligand_conc = models.TextField()
    
    def __str__(self):
        return self.ligand_conc

    class Meta():
        db_table = 'ligand_conc_of_crystallization'

class CrystallizationMethodTypes(models.Model):

    name = models.CharField(max_length=100) #LCP/ in surfo/bicelles

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'crystallization_method_types'  


## Construct refers to Expression Construct (longer seq) and is parent class. Crystallization construct is often shorter (removal of linkers, tags) and maybe "child" of the Construct class thus inheriting all features of it but present only if crystallizationConstruct sequence different than that of ExpressionConstruct sequence

class CrystallizationCondition(models.Model):   #called one or many times when explicitly fetched within ConstructCrystallization

    experiment = models.ForeignKey('ConstructCrystallisation')
    temp = models.CharField(max_length=5)
    ph = models.DecimalField(max_digits=2,decimal_places=2)
    chemical_list =  models.ForeignKey('ChemicalList') # can have many chemicals and their concentration listed

    def __str__(self):
        return self.temp

    class Meta():
        db_table = 'crystallization_condition'
