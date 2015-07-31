from django.db import models

# Create your models here.
class Mutation(models.Model):

	#links
    refs = models.ForeignKey('MutationRefs', null=True) #Change to a common model?
    protein = models.ForeignKey('protein.Protein')
    residue = models.ForeignKey('residue.Residue')
    ligand = models.ForeignKey('MutationLigand', null=True) #Change to a ligand model?
    ligand_class = models.ForeignKey('MutationLigandClass') #Change to a ligand model?
    ligand_ref = models.ForeignKey('MutationLigandRef') #Change to a ligand model?
    raw = models.ForeignKey('MutationRaw')
    optional = models.ForeignKey('MutationOptional')
    exp_type = models.ForeignKey('MutationType')
    exp_func= models.ForeignKey('MutationFunc')
    exp_measure = models.ForeignKey('MutationMeasure')
    exp_qual = models.ForeignKey('MutationQual')


    #Values
    mutation_to = models.CharField(max_length=1)
    wt_value = models.DecimalField(max_digits=10, decimal_places=2)
    wt_unit = models.CharField(max_length=10)
    mu_value = models.DecimalField(max_digits=10, decimal_places=2)
    mu_sign = models.CharField(max_length=2)
    foldchange = models.FloatField()

    def citation(self):

        temp = self.refs.citation.split(',')
        
        return temp[0] + " et al"
    
    class Meta():
        db_table = 'mutation'



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


class MutationType(models.Model):

	type = models.CharField(max_length=100)

	class Meta():
		db_table = 'mutation_type'


class MutationLigandClass(models.Model):

	classname = models.CharField(max_length=100)

	class Meta():
		db_table = 'mutation_ligand_class'

		
class MutationLigandRef(models.Model):

	reference = models.CharField(max_length=100)

	class Meta():
		db_table = 'mutation_ligand_reference'