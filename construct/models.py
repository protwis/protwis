from django.db import models

from protein.models import Protein
from residue.models import Residue
from mutation.models import Mutation ## check ??

# include models for Construct, prtn_order, exprssn, solubilisation, purifcn, and crystallisation 

class Construct(models.Model):

    parent = models.ForeignKey('self', null=True)
    protein_conformation = models.ForeignKey('protein.ProteinConformation') 
## Can get details like name,seq,pdbcode, ligand etc. from protein model
## Mutation data can also be extracted from protein entry [since here 2 protein entries will be made during db buildup : Expression construct protein n Cyrstallzn construct protein] both position & reasons of mutation can be accessed from "mutation.Mutation()" class/model in VIEWS.py! 
####    mutation_to = models.ForeignKey('mutation.amino_acid') ## linked from Protein model
####    mutation_reason = models.ForeignKey('mutation.type')

    deletions = models.TextField(max_length=100, null=True)  #comma separated data
#    aux_protein = models.ForeignKey('AuxProtein') 
#    expression = models.ForeignKey('ConstructExpression') 
 #  solubilization = models.ForeignKey('ConstructSolubilization')
  # purification = models.ForeignKey('ConstructPurification')
   #crystallization = models.ForeignKey('ConstructCrystallization')
## NOTE: If you need to create a relationship on a model that has not yet been defined, you can use the name of the model, rather than the model object itself

    def __str__(self):
        return self.protein.slug
    
    class Meta():
        db_table = 'construct'


class AuxProtein(models.Model):
    
    construct = models.ForeignKey('Construct')
    protein_type = models.ForeignKey('AuxProteinType')
    name = models.CharField(max_length=100, null=True)
    uniprot_id = models.CharField(max_length=20)
    sequence = models.TextField(null=True)
## AuxProtein mutation info also accessed from Protein model  ??? YES
#    mutation = models.TextField(max_length=100, blank=True) #comma separated data
####    mutation_to = models.ForeignKey('mutation.mutation_to')
####    mutation_reason = models.ForeignKey('mutation.MutationType')
    deletions = models.TextField(max_length=100, null=True) #comma separated data
    position = models.TextField(max_length=20, null=True) # aa_modfcn like position of disulfide bond is imp to show (is it shown already in GPCRdb Sites page??) For 3 linkers, values will be 1,2,3
    remarks = models.TextField(null=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'aux_protein'


class AuxProteinType(models.Model):

    name = models.TextField(max_length=50, null=True) #values like FUsion, tag, linkers,signal_peptide,prtn_Cleavage_site, aa_modfcn  #comma separated data

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'aux_protein_type'

            

# how to define protein order field?? Input as string of types separated by ":" ? and checked via ptn_type field  ::: NOt needed to define in models.py as it will used only by "build.py" script while making the order of construct

#==============================================================================


class Chemical(models.Model):

    chemical_type = models.ForeignKey('ChemicalType')
    name = models.CharField(max_length=200, null=True)

    def __str__(self):
        return self.name
    
    class Meta():
        db_table = 'chemical'


class ChemicalType(models.Model):

    name = models.CharField(max_length=100, null=True) # type like detergent, lipid  ## If there is no type present in database then Django creates another type

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
    
    construct_solubilization = models.ForeignKey('ConstructSolubilization')
    description =  models.TextField()

    def __str__(self):
        return self.description

    class Meta():
        db_table = 'chemical_modification'


#class ChromatographyType(models.Model):

#   description = models.TextField() # can have more than one chromatography method as comma separated fields

#    def __str__(self):            

#        return self.description

#    class Meta():
#        db_table = 'chromatography_method'


class PurificationStep(models.Model):
   
    purification = models.ForeignKey('ConstructPurification')
    purification_type = models.ForeignKey('PurificationStepType')
    description = models.TextField(null=True) # IMAC/Flag tag/1D4/Alprenolol/Strep tag/Ligand affinity


    def __str__(self):            

        return self.description

    class Meta():
        db_table = 'purification_step'



class PurificationStepType(models.Model):

    name = models.TextField(null=True) # Chromatography/Enzyme modification/Chem modification    IMAC/Flag tag/1D4/Alprenolol/Strep tag/Ligand affinity

    def __str__(self):            

        return self.name

    class Meta():
        db_table = 'purification_step_type'




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
#    chemical_modification = models.ForeignKey('ChemicalModification')#more than one chem_modfcn ## remove field, not easily searchable?
    remarks = models.TextField(null=True)

    def __str__(self):
        return self.chemical_list.chemicals.concentration

    class Meta():
        db_table = 'construct_solubilization'


class ConstructPurification(models.Model):

    construct = models.ForeignKey('Construct')
#    purification_type = models.ForeignKey('PurificationStep')
#   chromatography_type = models.ForeignKey('ChromatographyType') 
#   enzyme_modification = models.ForeignKey('EnzymeModification')
#   chemical_modification = models.ForeignKey('ChemicalModification')
    remarks = models.TextField(null=True)

    def __str__(self):
        return self.construct.protein_conformation.protein.entry_name

    class Meta():
        db_table = 'construct_purification'


class ConstructCrystallization(models.Model):

    construct = models.ForeignKey('Construct')
    crystal_type = models.ForeignKey('CrystallizationMethodTypes') # if Many crystallography types for one construct :named as separate Xtal Exp
    chemical_list =  models.ForeignKey('ChemicalList') #chemical type= LCPlipid for LCP exp, else type= lipid for in surfo exp --- includes info on lcp_lipidic_condition, protein_component etc.
    method = models.TextField(max_length=100,null=True) # more than one solved by comma separated
    settings = models.TextField(max_length=100,null=True)# can have more than one
    remarks = models.TextField(max_length=1000, null=True)
    protein_conc =  models.SlugField(max_length=20,blank=True)
    ligands = models.ManyToManyField('ligand.Ligand', through='ConstructCrystallizationLigandConc')
    aqueous_solution_lipid_ratio=  models.SlugField(max_length=20,null=True)
    lcp_bolus_volume =  models.SlugField(max_length=20,null=True)
    precipitant_solution_volume =  models.SlugField(max_length=20,null=True)
#   condition = models.ForeignKey('XtalCondition')
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
        db_table = 'ligand_conc_of_crystallization'

class CrystallizationMethodTypes(models.Model):

    name = models.CharField(max_length=100,null=True) #LCP/ in surfo/bicelles ... can be null if that construct did not crystallize

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'crystallization_method_types'  


## Construct refers to Expression Construct (longer seq) and is parent class. Crystallization construct is often shorter (removal of linkers, tags) and maybe "child" of the Construct class thus inheriting all features of it but present only if crystallizationConstruct sequence different than that of ExpressionConstruct sequence

#class CrystallizationCondition(models.Model):   #called one or many times when explicitly fetched within ConstructCrystallization

 #   experiment = models.ForeignKey('ConstructCrystallization')
  #  temp = models.CharField(max_length=5)
   # ph = models.DecimalField(max_digits=2,decimal_places=2)
    #chemical_list =  models.ForeignKey('ChemicalList') # can have many chemicals and their concentration listed

    #def __str__(self):
     #   return self.temp

    #class Meta():
     #   db_table = 'crystallization_condition'
