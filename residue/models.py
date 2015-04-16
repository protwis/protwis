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

    class Meta():
        db_table = 'residue'

    def three_letter(self):
        return self.AMINO_ACIDS[self.amino_acid]

    AMINO_ACIDS = {
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

    AMINO_ACID_GROUPS = {
        'hp' :    ('A', 'C', 'F', 'I', 'L', 'M', 'P', 'V', 'W', 'Y'),
        'alhp':   ('A', 'C', 'I', 'L', 'M', 'P', 'V'),
        'arhp':   ('F', 'W', 'Y'),
        'ar':     ('F', 'H', 'W', 'Y'),
        'pol':    ('D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T'),
        'hbd':    ('H', 'K', 'N', 'Q', 'R', 'S', 'T', 'W', 'Y'),
        'hbd':    ('H', 'K', 'N', 'Q', 'R', 'S', 'T', 'W', 'Y'),
        'hba':    ('D', 'E', 'H', 'N', 'Q', 'S', 'T', 'Y'),
        'neg':    ('D', 'E'),
        'pos':    ('H', 'K', 'R'),
        'lar':    ('E', 'F', 'H', 'K', 'Q', 'R', 'W', 'Y'), 
        'sma':    ('A', 'C', 'D', 'G', 'N', 'P', 'S', 'T', 'V'), 
        'any':    (),
        'custom': (),
    }

    AMINO_ACID_GROUP_NAMES = {
        'hp' :    'Hydrophobic - All',
        'alhp':   'Hydrophobic - Aliphatic',
        'arhp':   'Hydrophobic - Aromatic',
        'ar' :    'Aromatic',
        'pol':    'Polar',
        'hbd':    'H-Bond Donor',
        'hba':    'H-Bond Acceptor',
        'neg':    'Negative charge',
        'pos':    'Positive charge',
        'lar':    'Large',
        'sma':    'Small',
        'any':    'Any feature',
        'custom': 'Custom',
    }


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