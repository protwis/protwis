from django.db import models


class Residue(models.Model):
    protein = models.ForeignKey('protein.Protein')
    protein_segment = models.ForeignKey('protein.ProteinSegment', null=True)
    generic_number = models.ForeignKey('ResidueGenericNumber', related_name='compare', null=True)
    display_generic_number = models.ForeignKey('ResidueGenericNumber', related_name='display', null=True)
    alternative_generic_number = models.ManyToManyField('ResidueGenericNumber', related_name='alternative')
    sequence_number = models.SmallIntegerField()
    amino_acid = models.CharField(max_length=1)

    def __str__(self):
        return str(self.sequence_number) + self.amino_acid

    def three_letter(self):
        aa = AminoAcid.objects.get(one_letter_code=self.amino_acid)
        return aa.three_letter_code

    class Meta():
        db_table = 'residue'


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


class AminoAcid(models.Model):
    one_letter_code = models.CharField(max_length=1)
    three_letter_code = models.CharField(max_length=3)

    def __str__(self):
        return self.one_letter_code

    class Meta():
        db_table = 'amino_acid'


class AminoAcidSet(models.Model):
    amino_acids = models.ManyToManyField('AminoAcid')
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'amino_acid_set'
