from django.db import models
from common.models import Publication
from ligand.models import Ligand, LigandRole
import ast

# Create your models here.
class Mutation(models.Model):
    protein = models.ForeignKey('protein.Protein') 
    residue = models.ForeignKey('residue.Residue', null=True) #If auxilliary it will be null
    mutation_type = models.ForeignKey('MutationType', null=True)

    amino_acid = models.CharField(max_length=1) #amino acid one-letter

    class Meta():
        db_table = 'mutation'


class MutationType(models.Model):

    type = models.CharField(max_length=100) #type of mutation, eg thermostabilization, expression-increasing..

    class Meta():
        db_table = 'mutation_type'

class MutationExperiment(models.Model):

	#links
    refs = models.ForeignKey('common.Publication', null=True) #Change to a common model?
    protein = models.ForeignKey('protein.Protein')
    residue = models.ForeignKey('residue.Residue')
    mutation = models.ForeignKey('Mutation')
    ligand = models.ForeignKey('ligand.Ligand', null=True, related_name='ligand') #Change to a ligand model?
    ligand_role = models.ForeignKey('ligand.LigandRole', null=True) #Change to a ligand model?
    ligand_ref = models.ForeignKey('ligand.Ligand', null=True, related_name='reference_ligand') #Change to a ligand model?
    raw = models.ForeignKey('MutationRaw')
    optional = models.ForeignKey('MutationOptional', null=True)
    exp_type = models.ForeignKey('MutationExperimentalType', null=True)
    exp_func= models.ForeignKey('MutationFunc', null=True)
    exp_measure = models.ForeignKey('MutationMeasure', null=True)
    exp_qual = models.ForeignKey('MutationQual', null=True)

    #Values
    wt_value = models.DecimalField(max_digits=10, decimal_places=2)
    wt_unit = models.CharField(max_length=10)
    mu_value = models.DecimalField(max_digits=10, decimal_places=2)
    mu_sign = models.CharField(max_length=2)
    foldchange = models.FloatField()

    def citation(self):

        mainauthor = ast.literal_eval(self.refs.authors)[0]
        return mainauthor + " et al ("+str(self.refs.year)+")"

    def getCalculation(self):

        temp = ("Type: "+ self.exp_measure.measure + " <br> Measure: "+self.exp_type.type+" <br> Unit: " + str(self.wt_unit) +  " <br> WT: " + str(self.wt_value) + " <br> Mu: "+ str(self.mu_value) +" <br> Foldchange: "+str(self.foldchange))
        
        # if ($this->mut_effect_qual_id!=0) {
        #     $temp .= "\n".$this->mut_effect_qual->effect_qual. " ". $this->mut_effect_qual->effect_prop;    

        # }
        return temp

    def getFoldorQual(self):
        temp = self.exp_measure.measure
        if self.wt_value>0:
            if self.exp_measure.measure=='Activity/affinity':
                temp = round(self.mu_value/self.wt_value,2)
            elif self.exp_measure.measure=='Fold effect (mut/wt)':
                temp = round(self.mu_value,2)
            sign = ''
            if self.mu_sign!="=": sign = self.mu_sign
            
            if self.exp_measure.measure!='Qualitative effect': temp = round(self.foldchange,2) #use saved foldchange instaed
            
            if self.exp_measure.measure!='Qualitative effect' and temp!=0:
                
                if temp>1: 
                    temp =  "<font color='green'>"+sign + str(temp) + "↑</font>"
                elif temp<1:
                    temp =  "<font color='red'>"+sign + str(-temp) + "↓</font>"
                
            if self.exp_qual.qual:
                temp = self.exp_qual.qual +  " " +  self.exp_qual.prop  

        elif self.exp_qual.qual: #only display those with qual_id
            temp = self.exp_qual.qual +  " " + self.exp_qual.prop  
        
        return temp
    
    
    
    class Meta():
        db_table = 'mutation_experiment'



class MutationRefs(models.Model):

	year = models.SmallIntegerField()
	journal = models.CharField(max_length=100)
	title = models.TextField()
	citation = models.TextField()
	link = models.URLField()
	ref_type = models.CharField(max_length=100)
	reference  = models.CharField(max_length=100)

	def __str__(self):
		return self.link

	class Meta():
		db_table = 'mutation_refs'
		
class MutationLigand(models.Model):

	idtype = models.CharField(max_length=100)
	name = models.CharField(max_length=100)
	idid = models.CharField(max_length=100)
	longseq = models.TextField()

	class Meta():
		db_table = 'mutation_ligands'


class MutationOptional(models.Model):

	type = models.CharField(max_length=100)

	wt = models.DecimalField(max_digits=10, decimal_places=2)
	mu = models.DecimalField(max_digits=10, decimal_places=2)
	sign = models.CharField(max_length=2)
	percentage = models.DecimalField(max_digits=10, decimal_places=2)
	qual = models.CharField(max_length=100)
	agonist = models.CharField(max_length=100)

	class Meta():
		db_table = 'mutation_opt'

class MutationRaw(models.Model):


	reference = models.CharField(max_length=100)
	protein = models.CharField(max_length=100)
	mutation_pos = models.SmallIntegerField()
	mutation_from = models.CharField(max_length=1)
	mutation_to = models.CharField(max_length=1)

	ligand_name = models.CharField(max_length=100)
	ligand_idtype = models.CharField(max_length=100)
	ligand_id = models.CharField(max_length=100)
	ligand_class = models.CharField(max_length=100)

	exp_type = models.CharField(max_length=100)
	exp_func = models.CharField(max_length=100)

	exp_wt_value = models.DecimalField(max_digits=10, decimal_places=2)
	exp_wt_unit = models.CharField(max_length=10)

	exp_mu_effect_type = models.CharField(max_length=100)
	exp_mu_effect_sign = models.CharField(max_length=2)
	exp_mu_effect_value = models.DecimalField(max_digits=10, decimal_places=2)
	exp_mu_effect_qual = models.CharField(max_length=100)
	exp_mu_effect_ligand_prop = models.CharField(max_length=100)
	exp_mu_ligand_ref = models.CharField(max_length=100)

	opt_type = models.CharField(max_length=100)
	opt_wt = models.DecimalField(max_digits=10, decimal_places=2)
	opt_mu = models.DecimalField(max_digits=10, decimal_places=2)
	opt_sign = models.CharField(max_length=5)
	opt_percentage  = models.DecimalField(max_digits=10, decimal_places=2)
	opt_qual = models.CharField(max_length=100)
	opt_agonist = models.CharField(max_length=100)

	added_by = models.CharField(max_length=100)
	added_date = models.DateField()



	class Meta():
		db_table = 'mutation_raw'


class MutationMeasure(models.Model):

	measure = models.CharField(max_length=100)

	class Meta():
		db_table = 'mutation_measure'


class MutationQual(models.Model):

	qual = models.CharField(max_length=100)
	prop = models.CharField(max_length=100)

	class Meta():
		db_table = 'mutation_qual'


class MutationFunc(models.Model):

	func = models.CharField(max_length=100)

	class Meta():
		db_table = 'mutation_func'


class MutationExperimentalType(models.Model):

	type = models.CharField(max_length=100)

	class Meta():
		db_table = 'mutation_experimental_type'


class MutationLigandClass(models.Model):

	classname = models.CharField(max_length=100)

	class Meta():
		db_table = 'mutation_ligand_class'

		
class MutationLigandRef(models.Model):

	reference = models.CharField(max_length=100)

	class Meta():
		db_table = 'mutation_ligand_reference'