from django.db import models
from common.models import Publication
from ligand.models import Ligand, LigandRole
import ast

# Create your models here.
class Mutation(models.Model):
    protein = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    residue = models.ForeignKey('residue.Residue', null=True, on_delete=models.CASCADE) #If auxilliary it will be null
    amino_acid = models.CharField(max_length=1) #amino acid one-letter

    class Meta():
        db_table = 'mutation'
        unique_together = ('protein','residue','amino_acid')


class MutationExperiment(models.Model):

    #Data-submitting group
    submitting_group = models.CharField(max_length=200, null=True)

    #links
    refs = models.ForeignKey('common.Publication', null=True, on_delete=models.CASCADE) #Change to a common model?
    data_container = models.CharField(max_length=200, null=True) # is the data from a table, figure, Supl. info or inside the text in reference
    data_container_number = models.CharField(max_length=20, null=True) # if the data above is from a table/figure, what is the number

    review = models.ForeignKey('common.Publication', null=True, related_name='review', on_delete=models.CASCADE) #Change to a common model?
    protein = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    residue = models.ForeignKey('residue.Residue', on_delete=models.CASCADE)
    mutation = models.ForeignKey('Mutation', on_delete=models.CASCADE)
    ligand = models.ForeignKey('ligand.Ligand', null=True, related_name='ligand', on_delete=models.CASCADE) #Change to a ligand model?
    ligand_role = models.ForeignKey('ligand.LigandRole', null=True, on_delete=models.CASCADE) #Change to a ligand model?
    ligand_ref = models.ForeignKey('ligand.Ligand', null=True, related_name='reference_ligand', on_delete=models.CASCADE) #Change to a ligand model?
    raw = models.ForeignKey('MutationRaw', null=True, on_delete=models.CASCADE)
    exp_type = models.ForeignKey('MutationExperimentalType', null=True, on_delete=models.CASCADE)
    exp_func= models.ForeignKey('MutationFunc', null=True, on_delete=models.CASCADE)
    exp_qual = models.ForeignKey('MutationQual', null=True, on_delete=models.CASCADE)

    #Values
    wt_value = models.FloatField()
    wt_unit = models.CharField(max_length=10)
    mu_value = models.FloatField()
    mu_sign = models.CharField(max_length=2)
    foldchange = models.FloatField()

    #Opt
    opt_receptor_expression = models.FloatField(null=True)
    opt_basal_activity = models.FloatField(null=True)
    opt_gain_of_activity = models.CharField(max_length=100, null=True)
    opt_ligand_emax = models.FloatField(null=True)
    opt_agonist = models.CharField(max_length=100, null=True)

    def citation(self):

        if self.refs.authors:
            try:
                mainauthor = ast.literal_eval(self.refs.authors)[0]
                return mainauthor + " et al ("+str(self.refs.year)+")"
            except:
                # print(self.refs.authors)
                # print(self.refs.authors.split(','))
                mainauthor = self.refs.authors.split(',')[0]
                return mainauthor + " et al ("+str(self.refs.year)+")"
        else:
            return "N/A"

    def review_citation(self):

        if self.review:
            if self.review.year:
                try:
                    mainauthor = ast.literal_eval(self.review.authors)[0]
                    return mainauthor + " et al ("+str(self.review.year)+")"
                except:
                    #print(self.refs.authors.split(','))
                    mainauthor = self.review.authors.split(',')[0]
                    return mainauthor + " et al ("+str(self.review.year)+")"
            else:
                return self.review.web_link.index
        else:
            return ''
        #return  " et al ("+str(self.refs.year)+")"

    def getCalculation(self):

        if self.foldchange and self.exp_type and self.wt_value:
            # if self.wt_unit=='%':
            #     # Check calc is done right, since current error
            #     temp = round(self.mu_value/self.wt_value,3);
            #     if temp<1 and temp!=0:
            #         temp = -1/temp
            #     temp = -round(temp,3)
            #     if temp != self.foldchange:
            #         self.foldchange = temp
            #         self.save()
            temp = (" Measure: "+self.exp_type.type+" <br> Unit: " + str(self.wt_unit) +  " <br> WT: " + str(self.wt_value) + " <br> Mu: "+ str(self.mu_value) +" <br> Foldchange: "+str(self.foldchange))
        else:
            temp = "No information"
        # if ($this->mut_effect_qual_id!=0) {
        #     $temp .= "\n".$this->mut_effect_qual->effect_qual. " ". $this->mut_effect_qual->effect_prop;

        # }
        return temp

    def getFoldorQual(self):
        if self.foldchange!=0:
            temp = self.foldchange
            sign = ''
            if self.mu_sign!="=":
                sign = self.mu_sign

            # CHECK FOR % CALC
            # if self.wt_unit=='%' and self.wt_value:
            #     # Check calc is done right, since current error
            #     temp = round(self.mu_value/self.wt_value,3);
            #     if temp<1:
            #         temp = -1/temp
            #     temp = -round(temp,3)
            #     if temp != self.foldchange:
            #         self.foldchange = temp
            #         self.save()

            if temp>1:
                temp =  "<font color='green'>"+sign + str(temp) + "↑</font>"
            elif temp<1:
                temp =  "<font color='red'>"+sign + str(-temp) + "↓</font>"
            if self.exp_qual:
                temp = self.exp_qual.qual +  " " +  self.exp_qual.prop

        elif self.exp_qual: #only display those with qual_id
            temp = self.exp_qual.qual +  " " + self.exp_qual.prop
        else:
            temp = "N/A"
        return temp



    class Meta():
        db_table = 'mutation_experiment'


class MutationRaw(models.Model):


    submitting_group = models.CharField(max_length=200, null=True)
    reference = models.CharField(max_length=100)
    review = models.CharField(max_length=100, null=True)
    data_container = models.CharField(max_length=200, null=True) # is the data from a table, figure, Supl. info or inside the text in reference
    data_container_number = models.CharField(max_length=20, null=True) # if the data above is from a table/figure, what is the number
    protein = models.CharField(max_length=100)
    mutation_pos = models.SmallIntegerField()
    mutation_from = models.CharField(max_length=1)
    mutation_to = models.CharField(max_length=1)

    ligand_name = models.CharField(max_length=100)
    ligand_idtype = models.CharField(max_length=100)
    ligand_id = models.CharField(max_length=300)
    ligand_class = models.CharField(max_length=100)

    exp_type = models.CharField(max_length=100)
    exp_func = models.CharField(max_length=100)

    exp_wt_value = models.FloatField()
    exp_wt_unit = models.CharField(max_length=10)

    exp_mu_effect_type = models.CharField(max_length=100)
    exp_mu_effect_sign = models.CharField(max_length=2)
    exp_mu_effect_value = models.FloatField()
    exp_fold_change = models.FloatField()
    exp_mu_effect_qual = models.CharField(max_length=100)
    exp_mu_effect_ligand_prop = models.CharField(max_length=100)
    exp_mu_ligand_ref = models.CharField(max_length=100)

    opt_receptor_expression = models.CharField(max_length=100)
    opt_basal_activity = models.CharField(max_length=100)
    opt_gain_of_activity = models.CharField(max_length=100)
    opt_ligand_emax = models.CharField(max_length=100)
    opt_agonist = models.CharField(max_length=100)

    # opt_type = models.CharField(max_length=100)
    # opt_wt = models.FloatField()
    # opt_mu = models.FloatField()
    # opt_sign = models.CharField(max_length=5)
    # opt_percentage  = models.FloatField()
    # opt_qual = models.CharField(max_length=100)
    # opt_agonist = models.CharField(max_length=100)

    added_by = models.CharField(max_length=100)
    added_date = models.DateField()



    class Meta():
        db_table = 'mutation_raw'

    def __iter__(self):
        for field_name in self._meta.get_fields():
            try:
                value = getattr(self, field_name)
            except:
                value = None
            yield (field_name, value)


class MutationQual(models.Model):

    qual = models.CharField(max_length=100)
    prop = models.CharField(max_length=100)

    class Meta():
        db_table = 'mutation_qual'
        unique_together = ('qual','prop')



class MutationFunc(models.Model):

    func = models.CharField(max_length=100, unique=True)

    class Meta():
        db_table = 'mutation_func'


class MutationExperimentalType(models.Model):

    type = models.CharField(max_length=100, unique=True)

    class Meta():
        db_table = 'mutation_experimental_type'
