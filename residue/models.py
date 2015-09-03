from django.db import models
import common.definitions as definitions


class Residue(models.Model):
    protein_conformation = models.ForeignKey('protein.ProteinConformation')
    protein_segment = models.ForeignKey('protein.ProteinSegment', null=True)
    generic_number = models.ForeignKey('ResidueGenericNumber', related_name='compare', null=True)
    display_generic_number = models.ForeignKey('ResidueGenericNumber', related_name='display', null=True)
    alternative_generic_numbers = models.ManyToManyField('ResidueGenericNumber', related_name='alternative')
    sequence_number = models.SmallIntegerField()
    amino_acid = models.CharField(max_length=1)

    def __str__(self):
        return self.amino_acid + str(self.sequence_number)

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
    scheme = models.ForeignKey('ResidueNumberingScheme')
    protein_segment = models.ForeignKey('protein.ProteinSegment', related_name='generic_numbers', null=True)
    label = models.CharField(db_index=True, max_length=10)

    def __str__(self):
        return self.label
    
    class Meta():
        db_table = 'residue_generic_number'
        unique_together = ('scheme', 'label')


class ResidueNumberingScheme(models.Model):
    parent = models.ForeignKey('self', null=True)
    slug = models.SlugField(max_length=20)
    short_name = models.CharField(max_length=20)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'residue_generic_numbering_scheme'


class ResidueGenericNumberEquivalent(models.Model):
    default_generic_number = models.ForeignKey('ResidueGenericNumber')
    scheme = models.ForeignKey('ResidueNumberingScheme')
    label = models.CharField(db_index=True, max_length=10)

    def __str__(self):
        return self.label
    
    class Meta():
        db_table = 'residue_generic_number_equivalent'
        unique_together = ('scheme', 'label')