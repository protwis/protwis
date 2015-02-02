from django.db import models


class Residue(models.Model):
    protein = models.ForeignKey('protein.Protein')
    protein_segment = models.ForeignKey('protein.ProteinSegment', null=True)
    sequence_number = models.SmallIntegerField()
    generic_number = models.CharField(max_length=10, blank=True)
    generic_number_alt = models.CharField(max_length=10, blank=True)
    generic_number_alt2 = models.CharField(max_length=10, blank=True)
    generic_number_alt3 = models.CharField(max_length=10, blank=True)
    generic_number_alt4 = models.CharField(max_length=10, blank=True)
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
    
    
class ResidueNumber(models.Model):
    generic_number = models.CharField(max_length=10)

    def __str__(self):
        return self.generic_number
    
    class Meta():
        db_table = 'residue_number'