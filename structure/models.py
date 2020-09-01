from django.db import models
from django.core.cache import cache

from io import StringIO
from Bio.PDB import PDBIO
import re
from protein.models import ProteinGProteinPair

class Structure(models.Model):
    # linked onto the Xtal ProteinConformation, which is linked to the Xtal protein
    protein_conformation = models.ForeignKey('protein.ProteinConformation', on_delete=models.CASCADE)
    structure_type = models.ForeignKey('StructureType', on_delete=models.CASCADE)
    pdb_code = models.ForeignKey('common.WebLink', on_delete=models.CASCADE)
    state = models.ForeignKey('protein.ProteinState', on_delete=models.CASCADE)
    publication = models.ForeignKey('common.Publication', null=True, on_delete=models.CASCADE)
    ligands = models.ManyToManyField('ligand.Ligand', through='interaction.StructureLigandInteraction')
    extra_proteins = models.ManyToManyField('StructureExtraProteins', related_name='extra_proteins')
    protein_anomalies = models.ManyToManyField('protein.ProteinAnomaly')
    stabilizing_agents = models.ManyToManyField('StructureStabilizingAgent')
    preferred_chain = models.CharField(max_length=20)
    resolution = models.DecimalField(max_digits=5, decimal_places=3)
    publication_date = models.DateField()
    pdb_data = models.ForeignKey('PdbData', null=True, on_delete=models.CASCADE) #allow null for now, since dump file does not contain.
    representative = models.BooleanField(default=False)
    distance_representative = models.BooleanField(default=True)
    contact_representative = models.BooleanField(default=False)
    contact_representative_score = models.DecimalField(max_digits=5, decimal_places=3, null=True)
    inactive_class_contacts_fraction = models.DecimalField(max_digits=5, decimal_places=3, null=True)
    active_class_contacts_fraction = models.DecimalField(max_digits=5, decimal_places=3, null=True)
    class_contact_representative = models.BooleanField(default=False)
    annotated = models.BooleanField(default=True)
    refined = models.BooleanField(default=False)
    distance = models.DecimalField(max_digits=5, decimal_places=2, null=True)
    tm6_angle = models.DecimalField(max_digits=5, decimal_places=2, null=True)
    gprot_bound_likeness = models.DecimalField(max_digits=5, decimal_places=2, null=True)
    sodium = models.BooleanField(default=False)
    signprot_complex = models.ForeignKey('signprot.SignprotComplex', null=True, on_delete=models.SET_NULL, related_name='signprot_complex')
    stats_text = models.ForeignKey('StatsText', null=True, on_delete=models.CASCADE)
    mammal = models.BooleanField(default=False) #whether the species of the structure is mammal
    closest_to_human = models.BooleanField(default=False) # A boolean to say if the receptor/state of this structure is the closest structure to human

    def __str__(self):
        return self.pdb_code.index

    def get_stab_agents_gproteins(self):
        objs = self.stabilizing_agents.all()
        elements = [element for obj in objs for element in obj.name.split(',') if re.match(".*G.*", element) and not re.match(".*thase.*|PGS", element)]
        if len(elements) > 0:
            return "\n".join(elements)
        else:
            return '-'

    def get_signprot_gprot_family(self):
        tmp = self.signprot_complex.protein.family
        while tmp.parent.parent.parent.parent is not None:
            tmp = tmp.parent
        return tmp.name

        return str(self.signprot_complex.protein)

    def get_cleaned_pdb(self, pref_chain=True, remove_waters=True, ligands_to_keep=None, remove_aux=False, aux_range=5.0):

        tmp = []
        for line in self.pdb_data.pdb.split('\n'):
            save_line = False
            if pref_chain:
                # or 'refined' bit needs rework, it fucks up the extraction
                if (line.startswith('ATOM') or line.startswith('HET')) and (line[21] == self.preferred_chain[0] or 'refined' in self.pdb_code.index):
                # if (line.startswith('ATOM') or line.startswith('HET')) and (line[21] == self.preferred_chain[0]):
                    save_line = True
            else:
                save_line = True
            if remove_waters and line.startswith('HET') and line[17:20] == 'HOH':
                save_line = False
            if ligands_to_keep and line.startswith('HET'):
                if pref_chain:
                    if line[17:20] != 'HOH' and line[17:20] in ligands_to_keep and line[21] == self.preferred_chain[0]:
                        save_line = True
                    elif line[17:20] != 'HOH':
                        save_line=False
                else:
                    if line[17:20] != 'HOH' and line[17:20] in ligands_to_keep:
                        save_line = True
                    elif line[17:20] != 'HOH':
                        save_line=False
            if save_line:
                tmp.append(line)

        return '\n'.join(tmp)

    def get_ligand_pdb(self, ligand):

        tmp = []
        for line in self.pdb_data.pdb.split('\n'):
            if line.startswith('HET') and line[21] == self.preferred_chain[0]:
                if line[17:20] != 'HOH' and line[17:20] == ligand:
                    tmp.append(line)
        return '\n'.join(tmp)

    def get_preferred_chain_pdb(self):

        tmp = []
        for line in self.pdb_data.pdb.split('\n'):
            # http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
            if (line.startswith('ATOM') or line.startswith('HET')) and line[21] == self.preferred_chain[0]:
                tmp.append(line)
        return '\n'.join(tmp)

    @property
    def is_refined(self):
        # Ugly way of speeding up -- preferable the DB should create a relationship between the entries.
        refined = cache.get(self.pdb_code.index+'_refined')
        if refined == None:
            s = Structure.objects.filter(refined=True, pdb_code__index=self.pdb_code.index+'_refined')
            if len(s)>0:
                refined = True
            else:
                refined = False
            cache.set(self.pdb_code.index+'_refined',refined, 24*60*60)
        return refined

    class Meta():
        db_table = 'structure'


class StructureComplexProtein(models.Model):
    structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    protein_conformation = models.ForeignKey('protein.ProteinConformation', on_delete=models.CASCADE)
    chain = models.CharField(max_length=1)

    def __repr__(self):
        return '<StructureComplexProtein: '+str(self.protein_conformation.protein)+'>'

    def __str__(self):
        return '<StructureComplexProtein: '+str(self.protein_conformation.protein)+'>'

    class Meta():
        db_table = 'structure_complex_protein'

class StructureVectors(models.Model):
    structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    translation = models.CharField(max_length=100, null=True)
    center_axis = models.CharField(max_length=100)

    class Meta():
        db_table = 'structure_vectors'


class StructureModel(models.Model):
    protein = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    state = models.ForeignKey('protein.ProteinState', on_delete=models.CASCADE)
    main_template = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    pdb_data = models.ForeignKey('PdbData', null=True, on_delete=models.CASCADE)
    version = models.DateField()
    stats_text = models.ForeignKey('StatsText', on_delete=models.CASCADE)

    def __repr__(self):
        return '<HomologyModel: '+str(self.protein.entry_name)+' '+str(self.state)+'>'

    def __str__(self):
        return '<HomologyModel: '+str(self.protein.entry_name)+' '+str(self.state)+'>'

    class Meta():
        db_table = 'structure_model'

    def get_cleaned_pdb(self):
        return self.pdb_data.pdb


class StructureComplexModel(models.Model):
    receptor_protein = models.ForeignKey('protein.Protein', related_name='+', on_delete=models.CASCADE)
    sign_protein = models.ForeignKey('protein.Protein', related_name='+', on_delete=models.CASCADE)
    main_template = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    pdb_data = models.ForeignKey('PdbData', null=True, on_delete=models.CASCADE)
    version = models.DateField()
    # prot_signprot_pair = models.ForeignKey('protein.ProteinGProteinPair', related_name='+', on_delete=models.CASCADE, null=True)
    stats_text = models.ForeignKey('StatsText', on_delete=models.CASCADE)

    def __repr__(self):
        return '<ComplexHomologyModel: '+str(self.receptor_protein.entry_name)+'-'+str(self.sign_protein.entry_name)+'>'

    def __str__(self):
        return '<ComplexHomologyModel: '+str(self.receptor_protein.entry_name)+'-'+str(self.sign_protein.entry_name)+'>'

    class Meta():
        db_table = 'structure_complex_model'

    def get_cleaned_pdb(self):
        return self.pdb_data.pdb

    def get_prot_gprot_pair(self):
        pgp = ProteinGProteinPair.objects.filter(protein=self.receptor_protein, g_protein__slug=self.sign_protein.family.parent.slug, source='GuideToPharma')
        if len(pgp)>0:
            return pgp[0].transduction
        else:
            return 'no evidence'


class StatsText(models.Model):
    stats_text = models.TextField()

    def __repr__(self):
        if self.stats_text and len(self.stats_text)>0:
            line = self.stats_text.split('\n')[0]
        else:
            line = 'empty object'
        return '<StatsText: >'.format(line)

    def __str__(self):
        if self.stats_text and len(self.stats_text)>0:
            line = self.stats_text.split('\n')[0]
        else:
            line = 'empty object'
        return '<StatsText: >'.format(line)

    class Meta():
        db_table = 'stats_text'


class StructureModelStatsRotamer(models.Model):
    homology_model = models.ForeignKey('structure.StructureModel', on_delete=models.CASCADE)
    residue = models.ForeignKey('residue.Residue', null=True, on_delete=models.CASCADE)
    rotamer_template = models.ForeignKey('structure.Structure', related_name='+', null=True, on_delete=models.CASCADE)
    backbone_template = models.ForeignKey('structure.Structure', related_name='+', null=True, on_delete=models.CASCADE)

    def __repr__(self):
        return '<StructureModelStatsRotamer: seqnum '+str(self.residue.sequence_number)+' hommod '+str(self.homology_model.protein)+'>'

    class Meta():
        db_table = 'structure_model_stats_rotamer'


class StructureComplexModelStatsRotamer(models.Model):
    homology_model = models.ForeignKey('structure.StructureComplexModel', on_delete=models.CASCADE)
    protein = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    residue = models.ForeignKey('residue.Residue', null=True, on_delete=models.CASCADE)
    rotamer_template = models.ForeignKey('structure.Structure', related_name='+', null=True, on_delete=models.CASCADE)
    backbone_template = models.ForeignKey('structure.Structure', related_name='+', null=True, on_delete=models.CASCADE)

    def __repr__(self):
        return '<StructureComplexModelStatsRotamer: seqnum '+str(self.residue.sequence_number)+' hommod '+str(self.homology_model.protein)+'>'

    class Meta():
        db_table = 'structure_complex_model_stats_rotamer'


class StructureRefinedStatsRotamer(models.Model):
    structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    residue = models.ForeignKey('residue.Residue', null=True, on_delete=models.CASCADE)
    rotamer_template = models.ForeignKey('structure.Structure', related_name='+', null=True, on_delete=models.CASCADE)
    backbone_template = models.ForeignKey('structure.Structure', related_name='+', null=True, on_delete=models.CASCADE)

    def __repr__(self):
        return '<StructureRefinedStatsRotamer: seqnum '+str(self.residue.sequence_number)+' '+str(self.structure.pdb_code.index)+'>'

    class Meta():
        db_table = 'structure_refined_stats_rotamer'


class StructureModelSeqSim(models.Model):
    homology_model = models.ForeignKey('structure.StructureModel', on_delete=models.CASCADE)
    template = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    similarity = models.IntegerField()

    def __repr__(self):
        return '<StructureModelSeqSim: {}>'.format(self.homology_model.protein.entry_name)

    class Meta():
        db_table = 'structure_model_seqsim'


class StructureComplexModelSeqSim(models.Model):
    homology_model = models.ForeignKey('structure.StructureComplexModel', on_delete=models.CASCADE)
    template = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    similarity = models.IntegerField()

    def __repr__(self):
        return '<StructureComplexModelSeqSim: {}>'.format(self.homology_model.protein.entry_name)

    class Meta():
        db_table = 'structure_complex_model_seqsim'


class StructureRefinedSeqSim(models.Model):
    structure = models.ForeignKey('structure.Structure', related_name='+', null=True, on_delete=models.CASCADE)
    template = models.ForeignKey('structure.Structure', related_name='+', null=True, on_delete=models.CASCADE)
    similarity = models.IntegerField()

    def __repr__(self):
        return '<StructureRefinedSeqSim: {}>'.format(self.structure.pdb_code.index)

    class Meta():
        db_table = 'structure_refined_seqsim'


class StructureModelRMSD(models.Model):
    homology_model = models.ForeignKey('structure.StructureModel', on_delete=models.CASCADE)
    pdb = models.CharField(max_length=4)
    service = models.CharField(max_length=10)
    version = models.DecimalField(max_digits=2, decimal_places=1)
    date = models.DateField(null=True)
    overall_all = models.DecimalField(null=True, max_digits=3, decimal_places=1)
    overall_backbone = models.DecimalField(null=True, max_digits=3, decimal_places=1)
    TM_all = models.DecimalField(null=True, max_digits=3, decimal_places=1)
    TM_backbone = models.DecimalField(null=True, max_digits=3, decimal_places=1)

    def __repr__(self):
        return '<StructureModelRMSD: '+self.service+' hommod '+str(self.homology_model.protein)+'>'

    class Meta():
        db_table = 'structure_model_rmsd'


class StructureType(models.Model):
    slug = models.SlugField(max_length=20, unique=True)
    name = models.CharField(max_length=100)

    def type_short(self):
        if self.name=="X-ray diffraction":
            return "X-ray"
        elif self.name=="Electron microscopy":
            return "cryo-EM"
        else:
            return self.name

    def __str__(self):
        return self.name

    class Meta():
        db_table = "structure_type"


class StructureExtraProteins(models.Model):
    structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE, null=True)
    wt_protein = models.ForeignKey('protein.Protein', on_delete=models.CASCADE, null=True)
    protein_conformation = models.ForeignKey('protein.ProteinConformation', on_delete=models.CASCADE, null=True)
    display_name = models.CharField(max_length=20)
    note = models.CharField(max_length=50, null=True)
    chain = models.CharField(max_length=1)
    category = models.CharField(max_length=20)
    wt_coverage = models.IntegerField(null=True)

    def __str__(self):
        return self.display_name

    class Meta():
        db_table = "extra_proteins"


class StructureStabilizingAgent(models.Model):
    slug = models.SlugField(max_length=75, unique=True)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = "structure_stabilizing_agent"


class PdbData(models.Model):
    pdb = models.TextField()

    def __str__(self):
        return self.pdb

    class Meta():
        db_table = "structure_pdb_data"


class Rotamer(models.Model):
    residue = models.ForeignKey('residue.Residue', on_delete=models.CASCADE)
    structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    pdbdata = models.ForeignKey('PdbData', on_delete=models.CASCADE)
    missing_atoms = models.BooleanField(default=False)
    # TODO
    # Values: Angles
    def __str__(self):
        return '{} {}{}'.format(self.structure.pdb_code.index, self.residue.amino_acid, self.residue.sequence_number)

    class Meta():
        db_table = "structure_rotamer"


class Fragment(models.Model):
    residue = models.ForeignKey('residue.Residue', on_delete=models.CASCADE)
    ligand = models.ForeignKey('ligand.Ligand', on_delete=models.CASCADE)
    structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    pdbdata = models.ForeignKey('PdbData', on_delete=models.CASCADE)

    def __str__(self):
        return '{} {}{} {}'.format(self.structure.pdb_code.index, self.residue.amino_acid,
            self.residue.sequence_number, self.ligand.name)

    class Meta():
        db_table = "structure_fragment"


class StructureSegment(models.Model):
    structure = models.ForeignKey('Structure', on_delete=models.CASCADE)
    protein_segment = models.ForeignKey('protein.ProteinSegment', on_delete=models.CASCADE)
    start = models.IntegerField()
    end = models.IntegerField()

    def __str__(self):
        return self.structure.pdb_code.index + " " + self.protein_segment.slug

    class Meta():
        db_table = "structure_segment"


class StructureSegmentModeling(models.Model):
    """Annotations of segment borders that are observed in exp. structures, and can be used for modeling.
    This class is indentical to StructureSegment, but is kept separate to avoid confusion."""
    structure = models.ForeignKey('Structure', on_delete=models.CASCADE)
    protein_segment = models.ForeignKey('protein.ProteinSegment', on_delete=models.CASCADE)
    start = models.IntegerField()
    end = models.IntegerField()

    def __str__(self):
        return self.structure.pdb_code.index + " " + self.protein_segment.slug

    class Meta():
        db_table = "structure_segment_modeling"


class StructureCoordinates(models.Model):
    structure = models.ForeignKey('Structure', on_delete=models.CASCADE)
    protein_segment = models.ForeignKey('protein.ProteinSegment', on_delete=models.CASCADE)
    description = models.ForeignKey('StructureCoordinatesDescription', on_delete=models.CASCADE)

    def __str__(self):
        return "{} {} {}".format(self.structure.pdb_code.index, self.protein_segment.slug, self.description.text)

    class Meta():
        db_table = "structure_coordinates"


class StructureCoordinatesDescription(models.Model):
    text = models.CharField(max_length=200, unique=True)

    def __str__(self):
        return self.text

    class Meta():
        db_table = "structure_coordinates_description"


class StructureEngineering(models.Model):
    structure = models.ForeignKey('Structure', on_delete=models.CASCADE)
    protein_segment = models.ForeignKey('protein.ProteinSegment', on_delete=models.CASCADE)
    description = models.ForeignKey('StructureEngineeringDescription', on_delete=models.CASCADE)

    def __str__(self):
        return "{} {} {}".format(self.structure.pdb_code.index, self.protein_segment.slug, self.description.text)

    class Meta():
        db_table = "structure_engineering"


class StructureEngineeringDescription(models.Model):
    text = models.CharField(max_length=200, unique=True)

    def __str__(self):
        return self.text

    class Meta():
        db_table = "structure_engineering_description"
