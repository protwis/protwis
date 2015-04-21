from django.db import models


class Protein(models.Model):
    parent = models.ForeignKey('self', null=True)
    family = models.ForeignKey('ProteinFamily')
    species = models.ForeignKey('Species')
    source = models.ForeignKey('ProteinSource', null=True)
    residue_numbering_scheme = models.ForeignKey('residue.ResidueNumberingScheme')
    sequence_type = models.ForeignKey('ProteinSequenceType')
    endogenous_ligand = models.ManyToManyField('ligand.Ligand')
    web_link = models.ManyToManyField('common.WebLink')
    entry_name = models.SlugField(max_length=100, db_index=True, null=True)
    accession = models.CharField(max_length=100, db_index=True, null=True)
    name = models.CharField(max_length=200)
    sequence = models.TextField()
    
    # non-database attributes
    identity = False # % identity to a reference sequence in an alignment
    similarity = False # % similarity to a reference sequence in an alignment (% BLOSUM62 score > 0)
    similarity_score = False # similarity score to a reference sequence in an alignment (sum of BLOSUM62 scores)
    alignment = False # residues formatted for use in an Alignment class

    def __str__(self):
        if not self.entry_name:
            return self.name
        else:
            return str(self.entry_name)
    
    class Meta():
        db_table = 'protein'


class Gene(models.Model):
    proteins = models.ManyToManyField('Protein')
    species = models.ForeignKey('Species')
    name = models.CharField(max_length=100)
    position = models.SmallIntegerField()

    def __str__(self):
        return self.name

    class Meta():
        ordering = ('position', )
        db_table = 'gene'


class Species(models.Model):
    latin_name = models.CharField(max_length=100)
    common_name = models.CharField(max_length=100, blank=True)

    def __str__(self):
        return self.latin_name

    class Meta():     
        db_table = 'species'


class ProteinAlias(models.Model):
    protein = models.ForeignKey('Protein')
    name = models.CharField(max_length=200)
    position = models.SmallIntegerField()

    def __str__(self):
        return self.name

    class Meta():
        ordering = ('position', )
        db_table = 'protein_alias'


class ProteinSet(models.Model):
    protein = models.ManyToManyField('Protein')
    name = models.CharField(max_length=50)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'protein_set'


class ProteinSegment(models.Model):
    slug = models.SlugField(max_length=100)
    name = models.CharField(max_length=50)
    category = models.CharField(max_length=50)
    position = models.SmallIntegerField()

    def __str__(self):
        return self.name

    class Meta():
        ordering = ('position', )
        db_table = 'protein_segment'


class ProteinSource(models.Model):
    name = models.CharField(max_length=20)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'protein_source'


class ProteinFamily(models.Model):
    parent = models.ForeignKey('self', null=True)
    slug = models.SlugField(max_length=100)
    name = models.CharField(max_length=200)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'protein_family'
        ordering = ['id']


class ProteinSequenceType(models.Model):
    slug = models.SlugField(max_length=20)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'protein_sequence_type'


class ProteinAnomaly(models.Model):
    anomaly_type = models.ForeignKey('ProteinAnomalyType')
    generic_number = models.ForeignKey('residue.ResidueGenericNumber')

    def __str__(self):
        return self.generic_number.label

    class Meta():
        db_table = 'protein_anomaly'

class ProteinAnomalyType(models.Model):
    slug = models.SlugField(max_length=20)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'protein_anomaly_type'


class ProteinAnomalyRuleSet(models.Model):
    protein_anomaly = models.ForeignKey('ProteinAnomaly')

    def __str__(self):
        return self.protein_anomaly.generic_number.label

    class Meta():
        db_table = 'protein_anomaly_rule_set'


class ProteinAnomalyRule(models.Model):
    rule_set = models.ForeignKey('ProteinAnomalyRuleSet')
    generic_number = models.ForeignKey('residue.ResidueGenericNumber')
    amino_acid = models.CharField(max_length=1)

    def __str__(self):
        return "{} {}".format(self.generic_number.label, self.amino.acid)

    class Meta():
        db_table = 'protein_anomaly_rule'
