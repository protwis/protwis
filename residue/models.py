from django.db import models


class Residue(models.Model):
    protein = models.ForeignKey('protein.Protein')
    protein_segment = models.ForeignKey('protein.ProteinSegment', null=True)
    generic_number = models.ForeignKey('ResidueGenericNumber', related_name='compare', null=True)
    display_generic_number = models.ForeignKey('ResidueGenericNumber', related_name='display', null=True)
    alternative_generic_number = models.ManyToManyField('ResidueGenericNumber', related_name='alternative', null=True)
    sequence_number = models.SmallIntegerField()
    amino_acid = models.CharField(max_length=1)

    def __str__(self):
        return str(self.sequence_number) + self.amino_acid

    def three_letter(self):
        amino_acids = {
            'A': 'Ala',
            'C': 'Cys',
            'D': 'Asp',
            'E': 'Glu',
            'F': 'Phe',
            'G': 'Gly',
            'H': 'His',
            'I': 'Ile',
            'K': 'Lys',
            'L': 'Leu',
            'M': 'Met',
            'N': 'Asn',
            'P': 'Pro',
            'Q': 'Gln',
            'R': 'Arg',
            'S': 'Ser',
            'T': 'Thr',
            'V': 'Val',
            'W': 'Trp',
            'Y': 'Tyr',
        }
        return amino_acids[self.amino_acid]

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
    label = models.CharField(max_length=10)
    scheme = models.ForeignKey('ResidueNumberingScheme')
    protein_segment = models.ForeignKey('protein.ProteinSegment', null=True)

    def __str__(self):
        return self.label
    
    class Meta():
        db_table = 'generic_number'


class ResidueNumberingScheme(models.Model):

    #oli->gpcrdb->b-w->b-s
    #NUMBERING_SCHEMES = (
    #    ('bw', 'Ballesteros-Weinstein'),
    #    ('gpcrdb', 'GPCRdb generic'),
    #    ('oli', "Oliveira"),
    #    ('baldwin', "Baldwin-Schwartz")
    #)
    slug = models.SlugField(max_length=20)
    short_name = models.CharField(max_length=20)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'residue_numbering_scheme'