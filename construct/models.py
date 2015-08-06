from django.db import models

# include models like structure/ or mutation/ folders 

class Construct(models.Model):

       
    protein = models.ForeignKey('protein.Protein') ## need to show prtn name,deletions, mutations [position & reasons]
    name = models.CharField('protein.Protein',max_length=200)
    sequence = models.TextField('protein.Protein')
    source = models.ForeignKey('protein.ProteinSource')
    unipro_id = models.CharField(max_length=10)
    pdbcode = models.CharField(max_length=4)
    resolutn = models.DecimalField(max_digits=3, decimal_places=2)
    ligand = models.ForeignKey('mutation.MutationLigand')
    ligand_class = models.ForeignKey('mutation.MutationLigandClass')
    mutation = models.CharField(max_length=100) # How to show many mutations within a field? new table? 
    deletion = models.CharField(max_length=100)
    aux_ptn = models.ForeignKey('AuxProtein') ## 
    exprsn = models.ForeignKey('ConstructExprsn') 
    solblsn = models.ForeignKey('ConstructSolblsn')
    purifcn = models.ForeignKey('ConstructPurifcn')
    xtalzn = models.ForeignKey('ConstructXtalzn')## NOTE: If you need to create a relationship on a model that has not yet been defined, you can use the name of the model, rather than the model object itself

    class Meta():
        db_table = 'construct'


class ConstructExprsn(models.Model):

        exprsn_mtd = models.CharField(max_length=100)
        host_cell_type = models.CharField(max_length=100)
        host_cell = models.CharField(max_length=100)
        remarks = models.TextField()

        class Meta():
                db_table = 'construct_exprsn'

class ConstructSolblsn(models.Model):

        detergent = models.ForeignKey('SolblsnDetergent') # many detergent across solubilisation model ? or all constructs? since same detergent for many constructs (detergent = models.ManyToManyField('')) ?? Include detergent conc by "intermediate model"??
        lipid = models.CharField(max_length=100) # ManyToManyField('')) ?? Include lipid conc by "intermediate model"??
        chem_modfcn = models.TextField(max_length=100)# can have more than one chem_modfcn 
        remarks = models.TextField()

        class Meta():
                db_table = 'construct_solblsn'


class ConstructPurifcn(models.Model):

        chromat_type = models.ForeignKey('ChromatTypes') # Many chromatography types for one construct and one type used for many construct? (ManyToManyField(chromatography_method)) ???
        enzyme_modfcn = models.CharField(max_length=100) # more than one?
        chem_modfcn = models.TextField(max_length=100)# can have more than one chem_modfcn 
        remarks = models.TextField()

        class Meta():
                db_table = 'construct_purifcn'

class AuxProtein(models.Model):

        ptn_type = models.ForeignKey('AuxProteinType')
        name = models.CharField(max_length=100, blank=True)
        sequence = models.TextField(blank=True)
        fusion_mutation = models.CharField(max_length=100, blank=True) # How to show many mutations within a field? new table? 
        fusion_deletion = models.CharField(max_length=100, blank=True)
        position = models.SlugField(max_length=20, blank=True) # aa_modfcn like position of disulfide bond is imp to show (is it shown already in GPCRdb Sites page??) For 3 linkers, values will be 1,2,3
        remarks = models.TextField()

        class Meta():
                db_table = 'aux_protein'

class AuxProtein(models.Model):
        fusion = 'FP'
        linker = 'LN'
        tag = 'TG'
        sigpep = 'SP'
        pcsite = 'PC'
        aamod = 'AM'
        Aux_ptn_choice = (
           (fusion, 'Fusion_protein'),
           (linker, 'linkers'),
           (tag, 'tags'),
           (sigpep, 'signal_peptide'),
           (pcsite, 'protein_cleavage_site'),
           (aamod, 'Amino_Acid_modification'),
        )
        name = models.CharField(max_length=2, choices=Aux_ptn_choice) #values like FUsion, tag, linkers,signal_peptide,prtn_Cleavage_site, aa_modfcn

        class Meta():
                db_table = 'aux_protein_type'


# how to define protein order field?? Input as string of types separated by ":" ? and checked via ptn_type field
