from django.db import models

# include models for Construct, prtn_order, exprssn, solubilisation, purifcn, and crystallisation 

class Construct(models.Model):

    protein = models.ForeignKey('protein.Protein')
    name = models.CharField('protein.Protein',max_length=200)
    sequence = models.TextField('protein.Protein')
    source = models.ForeignKey('protein.ProteinSource')
    uniprot_id = models.CharField(max_length=10)
    pdbcode = models.CharField(max_length=4)
    resolution = models.DecimalField(max_digits=3, decimal_places=2)
    ligand = models.ForeignKey('mutation.MutationLigand')
    mutation = models.TextField(max_length=100) # how to show [position & reasons]?  ###comma separated data for listing many mutations
    deletion = models.TextField(max_length=100)  #comma separated data
    aux_ptn = models.ForeignKey('AuxProtein') 
    expression = models.ForeignKey('ConstructExpression') 
    solubilisation = models.ForeignKey('ConstructSolubilisation')
    purification = models.ForeignKey('ConstructPurification')
    crystallization = models.ForeignKey('ConstructCrystallisation')## NOTE: If you need to create a relationship on a model that has not yet been defined, you can use the name of the model, rather than the model object itself

    class Meta():
        db_table = 'construct'


class AuxProtein(models.Model):

    ptn_type = models.ForeignKey('AuxProteinType')
    name = models.CharField(max_length=100, blank=True)
    sequence = models.TextField(blank=True)
    fusion_mutation = models.TextField(max_length=100, blank=True) #comma separated data
    fusion_deletion = models.TextField(max_length=100, blank=True) #comma separated data
    position = models.TextField(max_length=20, blank=True) # aa_modfcn like position of disulfide bond is imp to show (is it shown already in GPCRdb Sites page??) For 3 linkers, values will be 1,2,3
    remarks = models.TextField()

    class Meta():
        db_table = 'aux_protein'


class AuxProteinType(models.Model):

    name = models.TextField(max_length=50) #values like FUsion, tag, linkers,signal_peptide,prtn_Cleavage_site, aa_modfcn  #comma separated data

    class Meta():
        db_table = 'aux_protein_type'


# how to define protein order field?? Input as string of types separated by ":" ? and checked via ptn_type field  ::: NOt needed to define in models.py as it will used only by "build.py" script while making the order of construct

#==============================================================================


class Chemical(models.Model):

    chemical_type = models.ForeignKey('ChemicalType')
    name = models.CharField(max_length=100)
    
    class Meta():
        db_table = 'chemical'


class ChemicalType(models.Model):

    name = models.CharField(max_length=100) # type like detergent, lipid  ## If there is no type present in database then Django creates another type

    class Meta():
        db_table = 'chemical_type'


class ChemicalConc(models.Model):
  
    chemical = models.ForeignKey('Chemical')
    concentration = models.TextField(null=True)

    class Meta():
        db_table = 'chemical_conc'


class ChemicalList(models.Model):   ## includes all chemicals, type & concentration. Chemicals can be from LCPlipid, detergent, etc.

    chemicals = models.ManyToManyField('ChemicalConc')        

    class Meta():
        db_table = 'chemical_list'


class ChemicalModification(models.Model):  # can have more than one chem_modfcn as comma separated fields

    description =  models.TextField()

    class Meta():
        db_table = 'chemical_modification'


class EnzymaticModification(models.Model): # can have more than one enzymatic_modfcn as comma separated fields
    description =  models.TextField()

    class Meta():
        db_table = 'enzymatic_modification'


class ChromatographyType(models.Model):

    description = models.TextField() # can have more than one enzymatic_modfcn as comma separated fields

    class Meta():
        db_table = 'chromatography_method'


class ConstructExpression(models.Model):

    expression_method = models.CharField(max_length=100)
    host_cell_type = models.CharField(max_length=100)
    host_cell = models.CharField(max_length=100)
    remarks = models.TextField()

    class Meta():
        db_table = 'construct_expression'


class ConstructSolubilisation(models.Model):

    chemical_list = models.ForeignKey('ChemicalList') #includes chem name, type[detergent & solubilisation_lipid]  and concentration
    chemical_modification = models.ForeignKey('ChemicalModification')#more than one chem_modfcn ## remove field, not easily searchable?
    remarks = models.TextField()

    class Meta():
        db_table = 'construct_solubilisation'


class ConstructPurification(models.Model):

    chromatography_type = models.ForeignKey('ChromatographyType') 
    enzyme_modification = models.ForeignKey('EnzymaticModification')
    chemical_modification = models.ForeignKey('ChemicalModification')
    remarks = models.TextField()

    class Meta():
        db_table = 'construct_purification'


class ConstructCrystallisation(models.Model):

    crystal_type = models.ForeignKey('CrystalMethodTypes') # if Many crystallography types for one construct :named as separate Xtal Exp
    chemical_list =  models.ForeignKey('ChemicalList') #chemical type= LCPlipid for LCP exp, else type= lipid for in surfo exp --- includes info on lcp_lipidic_condition, protein_component etc.
    method = models.TextField(max_length=100) # more than one solved by comma separated
    settings = models.TextField(max_length=100)# can have more than one
    remarks = models.TextField()
    ptn_conc =  models.SlugField(max_length=20,blank=True)
    ligands = models.ManyToManyField('ligand.Ligand', through='ConstructCrystallizationLigandConc')
    aqueous_soln_lipid_ratio=  models.SlugField(max_length=20,blank=True)
    lcp_bolus_vol =  models.SlugField(max_length=20,blank=True)
    precipitant_soln_vol =  models.SlugField(max_length=20,blank=True)
#   condition = models.ForeignKey('XtalCondition')


    class Meta():
        db_table = 'construct_crystallisation'


class ConstructCrystallizationLigandConc(models.Model):

    ligand_conc = models.TextField()

    class Meta():
        db_table = 'ligand_conc_of_crystallisation'

class CrystalMethodTypes(models.Model):

    name = models.CharField(max_length=100) #LCP/ in surfo/bicelles

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'crystal_method_types'  


class CrystalCondition(models.Model):   #called one or many times when explicitly fetched within ConstructCrystallisation

    experiment = models.ForeignKey('ConstructCrystallisation')
    temp = models.CharField(max_length=5)
    pH = models.DecimalField(max_digits=2,decimal_places=2)
    chemical_list =  models.ForeignKey('ChemicalList') # can have many chemicals and their concentration listed

    class Meta():
        db_table = 'crystal_condition'
