from django.db import models
from django.core.cache import cache

from protein.models import Protein
from residue.models import Residue
from mutation.models import Mutation

from construct.schematics import generate_schematic
from common.diagrams_gpcr import DrawSnakePlot

import pickle

class Construct(models.Model): 
    #overall class 
    name = models.TextField(max_length=100, unique=True)
    protein = models.ForeignKey('protein.Protein')
    contributor = models.ForeignKey('ContributorInfo')
    #Modifications
    mutations = models.ManyToManyField('ConstructMutation')
    deletions = models.ManyToManyField('ConstructDeletion')
    modifications = models.ManyToManyField('ConstructModification')
    insertions = models.ManyToManyField('ConstructInsertion')

    expression = models.ForeignKey('ExpressionSystem', null=True)  #method description if present
    solubilization = models.ForeignKey('Solubilization', null=True)  #method description if present
    purification = models.ForeignKey('Purification', null=True)  #method description if present
    crystallization = models.ForeignKey('Crystallization', null=True)  #method description if present
    crystal = models.ForeignKey('CrystalInfo', null=True) #might not exist, if failed
    structure = models.ForeignKey('structure.Structure', null=True) #might not exist, if failed
    schematics = models.BinaryField(null=True)
    snakecache = models.BinaryField(null=True)

    #Back up of original entry
    json = models.TextField(null=True)

    def fusion(self):
        list_of_none_fusion = ['Expression tag','Linker','Not_Observed','Engineered mutation','','Conflict','Insertion','S-arrestin']
        list_of_comfirmed_fusion = ['C8TP59','Q0SXH8','Q9V2J8','Soluble cytochrome b562','Endolysin','Rubredoxin','Lysozyme']
        #ODD Rubredoxin
        #Q9V2J8 GlgA glycogen synthase  auto_4ZJ8
        #C8TP59  Cytochrome b562
        # Q0SXH8 Cytochrome b(562)
        result = []
        position = None
        for insert in self.insertions.all():
            if insert.insert_type.subtype in list_of_none_fusion:
                continue
            if insert.insert_type.name!='fusion' and insert.insert_type.name!='auto':
                continue
            #print(insert.insert_type.name, insert.insert_type.subtype)
            #print(insert.position, self.name)
            confirmed = False
            if insert.insert_type.name=='fusion' or insert.insert_type.subtype in list_of_comfirmed_fusion:
                confirmed = True
                if insert.position.startswith('N-term'):
                    position = 'nterm'
                else:
                    position = 'icl3'
            result.append([confirmed,insert.insert_type.name, insert.insert_type.subtype,insert.position])
        return position,result

    def schematic(self):
        ## Use cache if possible
        temp = self.schematics
        if temp==None:
            # print(self.name+'_schematics no cache')
            temp = generate_schematic(self)
            self.schematics = pickle.dumps(temp)
            self.save()
        else:
            # print(self.name+'_schematics used cache')
            temp = pickle.loads(temp)
        return temp

    def snake(self):
        ## Use cache if possible
        temp = self.snakecache
        if temp==None:
            # print(self.name+'_snake no cache')
            residues = Residue.objects.filter(protein_conformation__protein=self.protein).order_by('sequence_number').prefetch_related(
                'protein_segment', 'generic_number', 'display_generic_number')
            temp = DrawSnakePlot(residues,self.protein.get_protein_class(),str(self.protein),nobuttons = True)
            self.snakecache = pickle.dumps(temp)
            self.save()
        else:
            temp = pickle.loads(temp)
        return temp

class CrystalInfo(models.Model):
    resolution = models.DecimalField(max_digits=5, decimal_places=3) #probably want more values
    pdb_data = models.ForeignKey('structure.PdbData', null=True) #if exists
    pdb_code = models.TextField(max_length=10, null=True) #if exists
    #No not include ligands here, as they should be part of crystalization 


class ContributorInfo(models.Model):
    name = models.TextField(max_length=50)
    pi_email = models.TextField(max_length=50)
    pi_name = models.TextField(max_length=50)
    urls = models.TextField() #can be comma seperated if many
    date = models.DateField() 
    address = models.TextField()


class ConstructMutation(models.Model):
    sequence_number = models.SmallIntegerField()
    wild_type_amino_acid = models.CharField(max_length=1)
    mutated_amino_acid = models.CharField(max_length=1)
    mutation_type = models.CharField(max_length=30, null=True)

    def __str__(self):
        return '{}{}{}'.format(self.wild_type_amino_acid, self.sequence_number,
            self.mutated_amino_acid)

    class Meta():
        db_table = 'construct_mutation'


class ConstructDeletion(models.Model):
    start = models.SmallIntegerField()
    end = models.SmallIntegerField()

    def __str__(self):
        return '{}-{}'.format(self.start, self.end)

    class Meta():
        db_table = 'construct_deletion'


class ConstructModification(models.Model):
    modification = models.TextField(max_length=50)
    position_type = models.TextField(max_length=20)
    pos_start = models.SmallIntegerField()
    pos_end = models.SmallIntegerField()
    remark = models.TextField()

    def __str__(self):
        return '{}'.format(self.modification)

    class Meta():
        db_table = 'construct_modification'


class ConstructInsertion(models.Model):
    insert_type = models.ForeignKey('ConstructInsertionType')
    position = models.TextField(max_length=20, null=True) #N-term1, N-term2 etc -- to track order
    presence = models.TextField(max_length=20, null=True) #YES or NO (presence in crystal)
    start = models.SmallIntegerField(null=True) #pos if in recepter
    end = models.SmallIntegerField(null=True) #pos if in recepter
    remarks = models.TextField(null=True)

    def __str__(self):
        return 'Protein: {}<br> Name: {}<br>Presence in Crystal: {}'.format(self.insert_type.name,self.insert_type.subtype,self.presence)

    class Meta():
        db_table = 'construct_insertion'


class ConstructInsertionType(models.Model):
    # values like FUsion, tag, linkers,signal_peptide,prtn_Cleavage_site, aa_modfcn (comma separated data)
    name = models.TextField(max_length=50, null=True)
    subtype = models.TextField(max_length=50, null=True)
    sequence = models.TextField(null=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'construct_insertion_type'


class Chemical(models.Model):
    chemical_type = models.ForeignKey('ChemicalType')
    name = models.CharField(max_length=200, null=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'construct_chemical'
        unique_together = ('name', 'chemical_type')


class ChemicalType(models.Model):
    # type like detergent, lipid
    name = models.CharField(max_length=100, null=True, unique=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'construct_chemical_type'


class ChemicalConc(models.Model):
    chemical = models.ForeignKey('Chemical')
    # concentration = models.FloatField()
    concentration = models.CharField(max_length=200)
    concentration_unit = models.TextField(null=True)

    def __str__(self):
        return self.chemical.name

    class Meta():
        db_table = 'construct_chemical_conc'    


# includes all chemicals, type & concentration. Chemicals can be from LCPlipid, detergent, etc.
class ChemicalList(models.Model):
    chemicals = models.ManyToManyField('ChemicalConc')
    name = models.ForeignKey('ChemicalListName',null=True)

    def __str__(self):
        return self.name.name

    class Meta():
        db_table = 'construct_chemical_list'

class ChemicalListName(models.Model):
    name = models.CharField(max_length=100, null=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'construct_chemical_list_name'


class ExpressionSystem (models.Model):
    expression_method = models.CharField(max_length=100)
    host_cell_type = models.CharField(max_length=100)
    host_cell = models.CharField(max_length=100)
    remarks = models.TextField(null=True)

    def __str__(self):

        return self.host_cell

    class Meta():
        db_table = 'construct_expression_system'



class Purification(models.Model):
    steps = models.ManyToManyField('PurificationStep')
    remarks = models.TextField(null=True)

    def __str__(self):
        return self.steps.name

    class Meta():
        db_table = 'construct_purification'


class PurificationStep(models.Model):
    name = models.TextField()
    # IMAC/Flag tag/1D4/Alprenolol/Strep tag/Ligand affinity
    description = models.TextField(null=True)

    def __str__(self):

        return self.name

    class Meta():
        db_table = 'construct_purification_step'



class Solubilization(models.Model):
    # includes chem name, type[detergent & solubilisation_lipid] and concentration
    chemical_list = models.ForeignKey('ChemicalList')
    remarks = models.TextField(null=True)

    def __str__(self):
        return self.chemical_list.chemicals.concentration

    class Meta():
        db_table = 'construct_solubilization'


class Crystallization(models.Model):
    crystal_type = models.ForeignKey('CrystallizationTypes', null=True)
    # chemical type= LCPlipid for LCP exp, else type= lipid for in surfo exp includes info on lcp_lipidic_condition,
    # protein_component etc.
    crystal_method = models.ForeignKey('CrystallizationMethods', null=True)
    remarks = models.TextField(max_length=1000, null=True)

    chemical_lists =  models.ManyToManyField('ChemicalList') #LCP list, detergent list, lip conc and other chem components
    #protein_conc =  models.FloatField()
    protein_conc = models.CharField(max_length=200,null=True)
    protein_conc_unit =  models.TextField(max_length=10,null=True)
    ligands = models.ManyToManyField('CrystallizationLigandConc')

    #temp = models.FloatField()
    temp = models.CharField(max_length=200)
    ph_start = models.DecimalField(max_digits=3, decimal_places=1) #probably want more values
    ph_end = models.DecimalField(max_digits=3, decimal_places=1) #probably want more values

    def __str__(self):
        return self.crystal_method.name

    def chemicals_total(self):
        c = 0
        for cl in self.chemical_lists.all():
            c += cl.chemicals.count()
        return c

    class Meta():
        db_table = 'construct_crystallization'


class CrystallizationLigandConc(models.Model):
    ligand = models.ForeignKey('ligand.Ligand')
    ligand_role = models.ForeignKey('ligand.LigandRole')
    #ligand_conc = models.FloatField(null=True)
    ligand_conc = models.CharField(max_length=200,null=True)
    ligand_conc_unit = models.TextField(null=True)

    def __str__(self):
        return self.ligand_conc

    class Meta():
        db_table = 'construct_crystallization_ligand_conc'


class CrystallizationTypes(models.Model):
    # LCP/ in surfo/bicelles
    name = models.CharField(max_length=100)
    sub_name = models.CharField(max_length=100,null=True) #for type of LCP method -- maybe other uses too

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'construct_crystallization_types'

class CrystallizationMethods(models.Model):
    # microbatch
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'construct_crystallization_methods'
