from django.db import models

# include models for Construct, prtn_order, exprssn, solubilisation, purifcn, and crystallisation 

class Construct(models.Model):

       
    protein = models.ForeignKey('protein.Protein')
    name = models.CharField('protein.Protein',max_length=200)
    sequence = models.TextField('protein.Protein')
    source = models.ForeignKey('protein.ProteinSource')
    uniprot_id = models.CharField(max_length=10)
    pdbcode = models.CharField(max_length=4)
    resolutn = models.DecimalField(max_digits=3, decimal_places=2)
    ligand = models.ForeignKey('mutation.MutationLigand')
    ligand_class = models.ForeignKey('mutation.MutationLigandClass')
    mutation = models.TextField(max_length=100) # how to show [position & reasons]?  ###comma separated data for listing many mutations
    deletion = models.TextField(max_length=100)  #comma separated data
    aux_ptn = models.ForeignKey('AuxProtein') 
    expression = models.ForeignKey('ConstructExprsn') 
    solubilisation = models.ForeignKey('ConstructSolblsn')
    purification = models.ForeignKey('ConstructPurifcn')
    crystallization = models.ForeignKey('ConstructXtalzn')## NOTE: If you need to create a relationship on a model that has not yet been defined, you can use the name of the model, rather than the model object itself

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

        name = models.CharField(max_length=100)
        chemical_type = models.ForeignKey('ChemicalType')

        class Meta():
                db_table = 'chemical'


class ChemicalType(models.Model):

        name = models.CharField(max_length=100) # type like detergent, lipid  ## If there is no type present in database then Django creates another type

        class Meta():
                db_table = 'chemical_type'


class ChemicalConc(models.Model):

        concentration = models.TextField(blank=True)
        chemical = models.ForeignKey('Chemical')

        class Meta():
                db_table = 'chemical_conc'


class ConstructExprsn(models.Model):

        expression_method = models.CharField(max_length=100)
        host_cell_type = models.CharField(max_length=100)
        host_cell = models.CharField(max_length=100)
        remarks = models.TextField()

        class Meta():
                db_table = 'construct_exprsn'

class ConstructSolblsn(models.Model):

#        detergent = models.ForeignKey('SolblsnDetergent') # many detergent across solubilisation model ? or all constructs? since same detergent for many constructs (detergent = models.ManyToManyField('')) ?? Include detergent conc by "intermediate model"??
 #       lipid = models.CharField(max_length=100) # ManyToManyField('')) ?? Include lipid conc by "intermediate model"??
        
        chemical_name = models.ForeignKey('Chemical') #includes chem name, type[detergent & solubilisation_lipid]  and concentration
        chem_modification = models.TextField(max_length=100)# can have more than one chem_modfcn 
        remarks = models.TextField()

        class Meta():
                db_table = 'construct_solblsn'


class ConstructPurifcn(models.Model):

        chromat_type = models.TextField(max_length=100) # Many chromatography types for one construct and one type used for many construct? (ManyToManyField(chromatography_method)) ???
        enzyme_modfcn = models.TextField(max_length=100)
        chem_modfcn = models.TextField(max_length=100)# can have more than one chem_modfcn as comma separated fields
        remarks = models.TextField()

        class Meta():
                db_table = 'construct_purifcn'

class ConstructXtalzn(models.Model):

        crystal_type = models.ManyToManyField('XtalMthdTypes') # Many crystallography types for one construct and one type used for many constructs
        method = models.TextField(max_length=100) # more than one solved by comma separated
        settings = models.TextField(max_length=100)# can have more than one
        lcp_lipidic_comp = models.ForeignKey('LCPlipidicComponent',blank=True)
        ptn_component = models.ForeignKey('ProteinComponent')
        condition1 = models.ForeignKey('XtalCond1')
        condition2 = models.ForeignKey('XtalCond2',blank=True)
        remarks = models.TextField()

        class Meta():
                db_table = 'construct_crystallisation'


class XtalMthdTypes(models.Model)

        name = models.CharField(max_length=100) #LCP/ in surfo/bicelles

        def __str__(self):
            return self.name

        class Meta():
                db_table = 'crystal_method_types'  


class LCPlipidicComponent(models.Model):

       
       lcp_lipid = models.CharField(max_length=50) # OR lcp_lipid = models.ForeignKey('Chemical') => to include chemical name/type to search/create in db
       lcp_lipid_conc = models.ForeignKey('ChemicalConc')
       many_lipids_one_exp = models.ForeignKey('ConstructXtalzn') #doesn't it create a loop while calling ??
       
       class Meta():
                db_table = 'LCP_lipidic_component'


class ProteinComponent(models.Model):

       detergent = models.CharField(max_length=50)
       detergent_conc =  models.ForeignKey('ChemicalConc')
       lipid = models.CharField(max_length=50,blank=True)
       lipid_conc =  models.ForeignKey('ChemicalConc')
       ptn_conc =  models.SlugField(max_length=20,blank=True)
       ligand = models.CharField(max_length=50,blank=True)
       ligand_conc =  models.ForeignKey('ChemicalConc')
       aqueous_soln_lipid_ratio=  models.SlugField(max_length=20,blank=True)
       lcp_bolus_vol =  models.SlugField(max_length=20,blank=True)
       precipitant_soln_vol =  models.SlugField(max_length=20,blank=True)


       class Meta():
                db_table = 'protein_component'

class XtalCond1(models.Model):

       temp = models.CharField(max_length=5)
       pH = models.DecimalField(max_digits=2,decimal_places=2)
       chem_component = models.CharField(max_length=50) 
       chem_component_conc =  models.ForeignKey('ChemicalConc') # can have many chemicals listed

       class Meta():
                db_table = 'xtal_condition1'

class XtalCond2(models.Model):

       temp = models.CharField(max_length=5)
       pH = models.DecimalField(max_digits=2,decimal_places=2)
       chem_component = models.CharField(max_length=50)
       chem_component_conc =  models.ForeignKey('ChemicalConc') # can have many chemicals listed

       class Meta():
                db_table = 'xtal_condition2'

