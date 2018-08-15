from django.db import models
import common.definitions as definitions


class Residue(models.Model):
    protein_conformation = models.ForeignKey('protein.ProteinConformation', on_delete=models.CASCADE)
    protein_segment = models.ForeignKey('protein.ProteinSegment', null=True, on_delete=models.CASCADE)
    generic_number = models.ForeignKey('ResidueGenericNumber', related_name='compare', null=True, on_delete=models.CASCADE)
    display_generic_number = models.ForeignKey('ResidueGenericNumber', related_name='display', null=True, on_delete=models.CASCADE)
    alternative_generic_numbers = models.ManyToManyField('ResidueGenericNumber', related_name='alternative')
    sequence_number = models.SmallIntegerField()
    amino_acid = models.CharField(max_length=1)

    def __str__(self):
        return self.amino_acid + str(self.sequence_number)

    def short_display_generic_number(self):
        if self.display_generic_number:
            return '%sx%s' % (self.display_generic_number.label.split(".")[0], self.display_generic_number.label.split("x")[1])
        else:
            return None

    class Meta():
        db_table = 'residue'
        ordering = ['sequence_number']

    def three_letter(self):
        if self.amino_acid in definitions.AMINO_ACIDS:
            return definitions.AMINO_ACIDS[self.amino_acid]
        else:
            return False

class ResidueDataType(models.Model):
    slug = models.SlugField(max_length=20)
    name = models.CharField(max_length=100)

class ResidueDataPoint(models.Model):
    data_type = models.ForeignKey('ResidueDataType', on_delete=models.CASCADE)
    residue = models.ForeignKey('Residue', related_name='annotations', on_delete=models.CASCADE)
    value = models.FloatField(null=True)
    value_text = models.CharField(max_length=50,null=True)

    def __str__(self):
        return '%s (%s)' % (self.data_type.name, str(self.value))

class ResidueSet(models.Model):
    residue = models.ManyToManyField('Residue')
    name = models.CharField(max_length=50)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'residue_set'


class ResidueGenericNumber(models.Model):
    scheme = models.ForeignKey('ResidueNumberingScheme', on_delete=models.CASCADE)
    protein_segment = models.ForeignKey('protein.ProteinSegment', related_name='generic_numbers', null=True, on_delete=models.CASCADE)
    label = models.CharField(db_index=True, max_length=12)

    def __str__(self):
        return self.label

    class Meta():
        db_table = 'residue_generic_number'
        unique_together = ('scheme', 'label')


class ResidueNumberingScheme(models.Model):
    parent = models.ForeignKey('self', null=True, on_delete=models.CASCADE)
    slug = models.SlugField(max_length=20)
    short_name = models.CharField(max_length=20)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'residue_generic_numbering_scheme'


class ResidueGenericNumberEquivalent(models.Model):
    default_generic_number = models.ForeignKey('ResidueGenericNumber', on_delete=models.CASCADE)
    scheme = models.ForeignKey('ResidueNumberingScheme', on_delete=models.CASCADE)
    label = models.CharField(db_index=True, max_length=12)

    def __str__(self):
        return self.label

    class Meta():
        db_table = 'residue_generic_number_equivalent'
        unique_together = ('scheme', 'label')


class ResiduePositionSet(models.Model):
    residue_position = models.ManyToManyField('ResidueGenericNumberEquivalent')
    name = models.CharField(max_length=50)
    set_type = models.CharField(max_length=100)
    protein_group = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'residue_position_set'
