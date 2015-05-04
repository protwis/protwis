from django.db import models

class Fragment(models.Model):
    #Basic set 
    residue = models.ForeignKey('residue.Residue')
    protein = models.ForeignKey('protein.Protein')
    structure = models.ForeignKey('structure.Structure', null=True)
    #Extension to interacting moiety
    ligand = models.ForeignKey('ligand.Ligand', null=True)
    interaction = models.ForeignKey('Interaction', null=True) #should this be M2M just in case?
    
    def __str__(self):
        return self.filename
    
    class Meta():
        db_table = 'fragment'

    def generate_file_name(self, pref_generic_numbering_scheme):
        generic_num = self.residue.generic_number.filter(scheme__slug=pref_generic_numbering_scheme)
        res_name = self.residue.amino_acid
        prot_entry_name = str(self.protein)
        pdb_code = self.structure.pdb_code.index
        if self.interaction is not null:
            interaction = str(self.interaction)
            return "{}_{}_{}_{}_{}.pdb".format(generic_num.replace('.','_'), res_name, prot_entry_name, prb_code, interaction)
        return "{}_{}_{}_{}.pdb".format(generic_num.replace('.','_'), res_name, prot_entry_name, prb_code)
        
 
class Interaction(models.Model):
    slug = models.CharField(max_length=10)
    description = models.CharField(max_length = 200)
    
    class Meta():
        db_table = 'interaction'

    def __str__(self):
        return self.slug


class FragmentAtom(models.Model):
    atomtype = models.CharField(max_length = 20)
    atomnr = models.SmallIntegerField()
    atomtype = models.CharField(max_length = 20)
    residuename = models.CharField(max_length = 20)
    chain = models.CharField(max_length = 20)
    residuenr = models.SmallIntegerField()
    x = models.DecimalField(max_digits=6, decimal_places=3)
    y = models.DecimalField(max_digits=6, decimal_places=3)
    z = models.DecimalField(max_digits=6, decimal_places=3)
    occupancy = models.DecimalField(max_digits=6, decimal_places=3)
    temperature = models.DecimalField(max_digits=6, decimal_places=3)
    element_name = models.CharField(max_length = 20)


    class Meta():
        db_table = 'fragment_atoms'

    def __str__(self):
        return self.atomtype

