from django.db import models
import common.definitions as definitions


class Residue(models.Model):
    protein_conformation = models.ForeignKey('protein.ProteinConformation')
    protein_segment = models.ForeignKey('protein.ProteinSegment', null=True)
    generic_number = models.ForeignKey('ResidueGenericNumber', related_name='compare', null=True)
    display_generic_number = models.ForeignKey('ResidueGenericNumber', related_name='display', null=True)
    alternative_generic_number = models.ManyToManyField('ResidueGenericNumber', related_name='alternative')
    sequence_number = models.SmallIntegerField()
    amino_acid = models.CharField(max_length=1)

    def __str__(self):
        return str(self.sequence_number) + self.amino_acid

    class Meta():
        db_table = 'residue'
        ordering = ['sequence_number']

    def three_letter(self):
        return definitions.AMINO_ACIDS[self.amino_acid]


class ResidueSet(models.Model):
    residue = models.ManyToManyField('Residue')
    name = models.CharField(max_length=50)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'residue_set'
    
    
class ResidueGenericNumber(models.Model):
    label = models.CharField(db_index=True, max_length=10)
    scheme = models.ForeignKey('ResidueNumberingScheme')
    protein_segment = models.ForeignKey('protein.ProteinSegment', null=True)

    def __str__(self):
        return self.label
    
    class Meta():
        db_table = 'generic_number'


class ResidueNumberingScheme(models.Model):
    slug = models.SlugField(max_length=20)
    short_name = models.CharField(max_length=20)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'residue_numbering_scheme'