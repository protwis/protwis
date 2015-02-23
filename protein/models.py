from django.db import models


class Protein(models.Model):
    family = models.ForeignKey('ProteinFamily')
    species = models.ForeignKey('Species')
    source = models.ForeignKey('ProteinSource')
    accession = models.CharField(max_length=100)
    entry_name = models.SlugField(max_length=100, unique=True)
    name = models.CharField(max_length=200)
    sequence = models.TextField()
    web_link = models.ManyToManyField('common.WebLink')

    def __str__(self):
        return self.entry_name
    
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
    slug = models.SlugField(max_length=100, unique=True)
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
    slug = models.SlugField(max_length=100, unique=True)
    name = models.CharField(max_length=200)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'protein_family'