from django.db import models
from common.diagrams_gpcr import DrawHelixBox, DrawSnakePlot
from residue.models import Residue

class Protein(models.Model):
    parent = models.ForeignKey('self', null=True)
    family = models.ForeignKey('ProteinFamily')
    species = models.ForeignKey('Species')
    source = models.ForeignKey('ProteinSource')
    residue_numbering_scheme = models.ForeignKey('residue.ResidueNumberingScheme')
    sequence_type = models.ForeignKey('ProteinSequenceType')
    states = models.ManyToManyField('ProteinState', through='ProteinConformation')
    endogenous_ligands = models.ManyToManyField('ligand.Ligand')
    web_links = models.ManyToManyField('common.WebLink')
    entry_name = models.SlugField(max_length=100, unique=True)
    accession = models.CharField(max_length=100, db_index=True, null=True)
    name = models.CharField(max_length=200)
    sequence = models.TextField()
    
    
    def __str__(self):
        if not self.entry_name:
            return self.name
        else:
            return self.entry_name
    
    class Meta():
        db_table = 'protein'

    def get_protein_class(self):
        tmp = self.family
        while tmp.parent.parent is not None:
            tmp = tmp.parent
        return tmp.name

    def get_helical_box(self):
        residuelist = Residue.objects.filter(protein_conformation__protein__entry_name=str(self)).prefetch_related('protein_segment','display_generic_number','generic_number')
        return DrawHelixBox(residuelist,self.get_protein_class(),str(self))

    def get_snake_plot(self):
        residuelist = Residue.objects.filter(protein_conformation__protein__entry_name=str(self)).prefetch_related('protein_segment','display_generic_number','generic_number')
        return DrawSnakePlot(residuelist,self.get_protein_class(),str(self))

    def get_protein_family(self):
        tmp = self.family
        while tmp.parent.parent.parent is not None:
            tmp = tmp.parent
        return tmp.name


class ProteinConformation(models.Model):
    protein = models.ForeignKey('Protein')
    state = models.ForeignKey('ProteinState')
    template_structure = models.ForeignKey('structure.Structure', null=True)
    protein_anomalies = models.ManyToManyField('protein.ProteinAnomaly')

    # non-database attributes
    identity = 0 # % identity to a reference sequence in an alignment
    similarity = 0 # % similarity to a reference sequence in an alignment (% BLOSUM62 score > 0)
    similarity_score = 0 # similarity score to a reference sequence in an alignment (sum of BLOSUM62 scores)
    alignment = 0 # residues formatted for use in an Alignment class
    alignment_list = 0 # FIXME redundant, remove when dependecies are removed

    def __str__(self):
        return self.protein.entry_name + " (" + self.state.slug + ")"

    class Meta():
        ordering = ('id', )
        db_table = "protein_conformation"


class ProteinState(models.Model):
    slug = models.SlugField(max_length=20, unique=True)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.slug

    class Meta():
        db_table = "protein_state"


class Gene(models.Model):
    proteins = models.ManyToManyField('Protein', related_name='genes')
    species = models.ForeignKey('Species')
    name = models.CharField(max_length=100)
    position = models.SmallIntegerField()

    def __str__(self):
        return self.name

    class Meta():
        ordering = ('position', )
        db_table = 'gene'


class Species(models.Model):
    latin_name = models.CharField(max_length=100, unique=True)
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
    proteins = models.ManyToManyField('Protein')
    name = models.CharField(max_length=50, unique=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'protein_set'


class ProteinSegment(models.Model):
    slug = models.SlugField(max_length=100, unique=True)
    name = models.CharField(max_length=50)
    category = models.CharField(max_length=50)
    fully_aligned = models.BooleanField(default=False)
    partial = models.BooleanField(default=False)

    def __str__(self):
        return self.slug

    class Meta():
        ordering = ('id', )
        db_table = 'protein_segment'


class ProteinSource(models.Model):
    name = models.CharField(max_length=20, unique=True)

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
        ordering = ('id', )


class ProteinSequenceType(models.Model):
    slug = models.SlugField(max_length=20, unique=True)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.slug

    class Meta():
        db_table = 'protein_sequence_type'


class ProteinAnomaly(models.Model):
    anomaly_type = models.ForeignKey('ProteinAnomalyType')
    generic_number = models.ForeignKey('residue.ResidueGenericNumber')

    def __str__(self):
        return self.generic_number.label

    class Meta():
        ordering = ('generic_number__label', )
        db_table = 'protein_anomaly'
        unique_together = ('anomaly_type', 'generic_number')

class ProteinAnomalyType(models.Model):
    slug = models.SlugField(max_length=20, unique=True)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.slug

    class Meta():
        db_table = 'protein_anomaly_type'


class ProteinAnomalyRuleSet(models.Model):
    protein_anomaly = models.ForeignKey('ProteinAnomaly', related_name='rulesets')
    exclusive = models.BooleanField(default=False)

    def __str__(self):
        return self.protein_anomaly.generic_number.label

    class Meta():
        db_table = 'protein_anomaly_rule_set'
        ordering = ('id', )


class ProteinAnomalyRule(models.Model):
    rule_set = models.ForeignKey('ProteinAnomalyRuleSet', related_name='rules')
    generic_number = models.ForeignKey('residue.ResidueGenericNumber')
    amino_acid = models.CharField(max_length=1)
    negative = models.BooleanField(default=False)

    def __str__(self):
        return "{} {}".format(self.generic_number.label, self.amino_acid)

    class Meta():
        db_table = 'protein_anomaly_rule'


class ProteinFusion(models.Model):
    proteins = models.ManyToManyField('Protein', through='ProteinFusionProtein')
    name = models.CharField(max_length=100, unique=True)
    sequence = models.TextField(null=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'protein_fusion'


class ProteinFusionProtein(models.Model):
    protein = models.ForeignKey('Protein')
    protein_fusion = models.ForeignKey('ProteinFusion')
    segment_before = models.ForeignKey('ProteinSegment', related_name='segment_before')
    segment_after = models.ForeignKey('ProteinSegment', related_name='segment_after')

    def __str__(self):
        return self.protein.name + " " + self.protein_fusion.name

    class Meta():
        db_table = 'protein_fusion_protein'


class ProteinConformationTemplateStructure(models.Model):
    protein_conformation = models.ForeignKey('ProteinConformation')
    protein_segment = models.ForeignKey('ProteinSegment')
    structure = models.ForeignKey('structure.Structure')

    def __str__(self):
        return self.protein_conformation.protein.name + " " + self.protein_segment.slug \
        + " " + self.structure.pdb_code.index

    class Meta():
        db_table = 'protein_conformation_template_structure'