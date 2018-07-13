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
    # linked onto the WT protein
    protein = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    contributor = models.ForeignKey('ContributorInfo', on_delete=models.CASCADE)
    #Modifications
    # mutations = models.ManyToManyField('ConstructMutation') ## MOVE TO ONETOMANY FIELDS in children
    # deletions = models.ManyToManyField('ConstructDeletion')
    # modifications = models.ManyToManyField('ConstructModification')
    # insertions = models.ManyToManyField('ConstructInsertion')

    expression = models.ForeignKey('ExpressionSystem', null=True, on_delete=models.CASCADE)  #method description if present
    solubilization = models.ForeignKey('Solubilization', null=True, on_delete=models.CASCADE)  #method description if present
    purification = models.ForeignKey('Purification', null=True, on_delete=models.CASCADE)  #method description if present
    crystallization = models.ForeignKey('Crystallization', null=True, on_delete=models.CASCADE)  #method description if present
    crystal = models.ForeignKey('CrystalInfo', null=True, on_delete=models.CASCADE) #might not exist, if failed
    structure = models.ForeignKey('structure.Structure', null=True, on_delete=models.CASCADE, related_name='construct') #might not exist, if failed
    schematics = models.BinaryField(null=True)
    snakecache = models.BinaryField(null=True)

    #Back up of original entry
    json = models.TextField(null=True)

    def fusion(self):
        list_of_none_fusion = ['Expression tag','Linker','Not_Observed','Engineered mutation','','Conflict','Insertion','S-arrestin','3A Arrestin']
        list_of_comfirmed_fusion = ['C8TP59','Q0SXH8','Q9V2J8','Soluble cytochrome b562','Endolysin','Rubredoxin','Lysozyme','Flavodoxin']
        #ODD Rubredoxin
        #Q9V2J8 GlgA glycogen synthase  auto_4ZJ8
        #C8TP59  Cytochrome b562
        # Q0SXH8 Cytochrome b(562)
        result = []
        position = None
        linker = {'before':'','after':''}
        for insert in self.insertions.all():
            if not insert.presence=='YES':
                continue
            if insert.insert_type.subtype in list_of_none_fusion:
                continue
            if insert.insert_type.name!='fusion' and insert.insert_type.name!='auto':
                continue
            #print(insert.insert_type.name, insert.insert_type.subtype)
            confirmed = False
            if insert.insert_type.name=='fusion' or insert.insert_type.subtype in list_of_comfirmed_fusion:
                confirmed = True
                # if position != None:
                #     print("new fusion??",position,insert.position,self.name)
                if insert.position.startswith('N-term'):
                    if position:
                        position += '_nterm'
                    else:
                        position = 'nterm'
                else:
                    if position:
                        position += '_icl3'
                    else:
                        position = 'icl3'
            result.append([confirmed,insert.insert_type.name, insert.insert_type.subtype,insert.position,insert.start,insert.end,'',''])
        
        if position:
            for insert in self.insertions.all():
                if insert.presence=='YES' and insert.insert_type.name=='linker':
                    if result[0][3].split("_")[0] == insert.position.split("_")[0]:
                        if result[0][4] is None or insert.start is None or abs(result[0][4]-insert.start)<len(insert.insert_type.subtype)+5:
                            # pass
                            i_relative = 'after'
                            if int(insert.position.split("_")[-1])<int(result[0][3].split("_")[-1]):
                                i_relative = 'before'
                            linker[i_relative] = insert.insert_type.subtype
                            print("LINKER around fusion",self.structure, self.protein.entry_name,insert.position,insert.insert_type.subtype,result)

        if len(result)>1:
            print(position,self.structure, self.protein.entry_name,result)
        return position,result,linker

    def cons_schematic(self):
        cache_key = self.name + "_cons_schematic"
        schematic = cache.get(cache_key)
        if schematic==None:
            temp = self.schematic()
            cache_key = self.name + "_wt_schematic"
            schematic = temp['schematic_2_wt']
            cache.set(cache_key,schematic,60*60*24*7)
            cache_key = self.name + "_chem_summary"
            summary = temp['summary']
            cache.set(cache_key,summary,60*60*24*7)
            schematic = temp['schematic_2_c']
            cache.set(cache_key,schematic,60*60*24*7)

        return schematic

    def wt_schematic(self):
        cache_key = self.name + "_wt_schematic"
        schematic = cache.get(cache_key)
        if schematic==None:
            temp = self.schematic()
            cache_key = self.name + "_cons_schematic"
            schematic = temp['schematic_2_c']
            cache.set(cache_key,schematic,60*60*24*7)
            cache_key = self.name + "_chem_summary"
            summary = temp['summary']
            cache.set(cache_key,summary,60*60*24*7)
            schematic = temp['schematic_2_wt']
            cache.set(cache_key,schematic,60*60*24*7)


        return schematic

    def chem_summary(self):
        cache_key = self.name + "_chem_summary"
        summary = cache.get(cache_key)
        if summary==None:
            temp = self.schematic()
            summary = temp['summary']
            cache.set(cache_key,summary,60*60*24*7)
            cache_key = self.name + "_cons_schematic"
            schematic = temp['schematic_2_c']
            cache.set(cache_key,schematic,60*60*24*7)
            cache_key = self.name + "_wt_schematic"
            schematic = temp['schematic_2_wt']
            cache.set(cache_key,schematic,60*60*24*7)

        return summary

    def invalidate_schematics(self):
        cache.delete_many([self.name + "_cons_schematic",self.name + "_wt_schematic",self.name + "_chem_summary"])
        self.schematics = None
        self.save()


    def schematic(self):

        cache_key = self.name + "_all_schematics"
        try:
            temp = cache.get(cache_key)
        except:
            temp = None

        if temp==None:
            temp = generate_schematic(self)
            cache.set(cache_key,temp,60*60*24*7)
            cache_key = self.name + "_cons_schematic"
            schematic = temp['schematic_2_c']
            cache.set(cache_key,schematic,60*60*24*7)
            cache_key = self.name + "_wt_schematic"
            schematic = temp['schematic_2_wt']
            cache.set(cache_key,schematic,60*60*24*7)
            cache_key = self.name + "_chem_summary"
            summary = temp['summary']
            cache.set(cache_key,summary,60*60*24*7)
        return temp

    def snake(self):
        ## Use cache if possible
        temp = self.snakecache
        if temp==None:
            print(self.name+'_snake no cache')
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
    pdb_data = models.ForeignKey('structure.PdbData', null=True, on_delete=models.CASCADE) #if exists
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
    # construct = models.ManyToManyField('ConstructMutation')
    construct = models.ForeignKey('Construct', related_name = 'mutations', on_delete=models.CASCADE)
    sequence_number = models.SmallIntegerField()
    wild_type_amino_acid = models.CharField(max_length=1)
    mutated_amino_acid = models.CharField(max_length=1)
    # mutation_type = models.CharField(max_length=30, null=True)
    remark = models.TextField(null=True)
    residue = models.ForeignKey('residue.Residue', null=True, on_delete=models.CASCADE)
    effects = models.ManyToManyField('ConstructMutationType', related_name="bar")

    def get_res(self):
        '''Retrieve the residue connected to this mutation, and save it as a FK field.'''
        # try:
        #     construct = self.construct_set.get().structure.protein_conformation.protein
        # except Construct.DoesNotExist:
        #     print('no construct for this mutation')
        #     return None
        seq_no = self.sequence_number
        try:
            # res_cons = Residue.objects.get(protein_conformation__protein=construct, sequence_number=seq_no)
            res_wt = Residue.objects.get(protein_conformation__protein=self.construct.structure.protein_conformation.protein.parent, sequence_number=seq_no)
            if res_wt.amino_acid != self.wild_type_amino_acid:
                pass
                # print('aa dont match',construct.name,seq_no,"annotated wt:", self.wild_type_amino_acid, "DB wt:",res_wt.amino_acid, "Annotated Mut",self.mutated_amino_acid)
            #     print('records wt',res_wt.amino_acid,'construct res',res_cons.amino_acid)
            return res_wt
        except Residue.DoesNotExist:
            print('no residue for',construct,seq_no)
            return None



    # def save(self, *args, **kwargs):
    #     '''Modify save function to automatically get the associated residue, should it exist'''
    #     self.residue_id = self.get_res()
    #     super(ConstructMutation, self).save(*args, **kwargs)


    def __str__(self):
        return '{}{}{}'.format(self.wild_type_amino_acid, self.sequence_number,
                               self.mutated_amino_acid)

    class Meta():
        db_table = 'construct_mutation'

class ConstructMutationType(models.Model):
    slug = models.SlugField(max_length=100)
    name = models.CharField(max_length=100)
    effect = models.CharField(max_length=100, null=True)

    def __str__(self):
        return self.slug

    class Meta():
        db_table = 'construct_mutation_type'
        unique_together = ('slug', 'name','effect')


class ConstructDeletion(models.Model):
    construct = models.ForeignKey('Construct', related_name = 'deletions', on_delete=models.CASCADE)
    start = models.SmallIntegerField()
    end = models.SmallIntegerField()

    def __str__(self):
        return '{}-{}'.format(self.start, self.end)

    class Meta():
        db_table = 'construct_deletion'


class ConstructModification(models.Model):
    construct = models.ForeignKey('Construct', related_name = 'modifications', on_delete=models.CASCADE)
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
    construct = models.ForeignKey('Construct', related_name = 'insertions', on_delete=models.CASCADE)
    insert_type = models.ForeignKey('ConstructInsertionType', on_delete=models.CASCADE)
    position = models.TextField(max_length=20, null=True) #N-term1, N-term2 etc -- to track order
    presence = models.TextField(max_length=20, null=True) #YES or NO (presence in crystal)
    start = models.SmallIntegerField(null=True) #pos if in recepter
    end = models.SmallIntegerField(null=True) #pos if in recepter
    remarks = models.TextField(null=True)

    def __str__(self):
        return 'Protein: {}<br> Name: {}<br>Presence in Crystal: {}'.format(self.insert_type.name,self.insert_type.subtype,self.presence)

    def autotype(self):
        t = self.insert_type.subtype
        if t=='Expression tag':
            return 'tag'
        elif t=='Linker':
            return 'linker'
        else:
            return self.insert_type.name

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
    chemical_type = models.ForeignKey('ChemicalType', on_delete=models.CASCADE)
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
    chemical = models.ForeignKey('Chemical', on_delete=models.CASCADE)
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
    name = models.ForeignKey('ChemicalListName',null=True, on_delete=models.CASCADE)

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
    expression_time = models.CharField(max_length=100,null=True)
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
    chemical_list = models.ForeignKey('ChemicalList', on_delete=models.CASCADE)
    remarks = models.TextField(null=True)

    def __str__(self):
        return self.chemical_list.chemicals.concentration

    class Meta():
        db_table = 'construct_solubilization'


class Crystallization(models.Model):
    crystal_type = models.ForeignKey('CrystallizationTypes', null=True, on_delete=models.CASCADE)
    # chemical type= LCPlipid for LCP exp, else type= lipid for in surfo exp includes info on lcp_lipidic_condition,
    # protein_component etc.
    crystal_method = models.ForeignKey('CrystallizationMethods', null=True, on_delete=models.CASCADE)
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
    ligand = models.ForeignKey('ligand.Ligand', on_delete=models.CASCADE)
    ligand_role = models.ForeignKey('ligand.LigandRole', on_delete=models.CASCADE)
    #ligand_conc = models.FloatField(null=True)
    ligand_conc = models.CharField(max_length=200,null=True)
    ligand_conc_unit = models.TextField(null=True)

    def __str__(self):
        return self.ligand_conc

    class Meta():
        db_table = 'construct_crystallization_ligand_conc'


class CrystallizationTypes(models.Model):
    # LCP/ in surfo/bicelles
    name = models.CharField(max_length=300)
    sub_name = models.CharField(max_length=300,null=True) #for type of LCP method -- maybe other uses too

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
