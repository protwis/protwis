from rest_framework import serializers

from protein.models import Protein, ProteinConformation, ProteinFamily, Species, ProteinSource, ProteinSegment
from residue.models import Residue, ResidueNumberingScheme, ResidueGenericNumber
from structure.models import Structure


class ProteinSerializer(serializers.ModelSerializer):
    family = serializers.SlugRelatedField(read_only=True, slug_field='slug')
    species = serializers.StringRelatedField(read_only=True)
    source = serializers.StringRelatedField(read_only=True)
    residue_numbering_scheme = serializers.SlugRelatedField(read_only=True, slug_field='short_name')
    class Meta:
        model = Protein
        fields = ('entry_name', 'name', 'accession', 'family', 'species', 'source', 'residue_numbering_scheme',
            'sequence')


class ProteinFromConformationSerializer(serializers.ModelSerializer):
    protein = serializers.SlugRelatedField(read_only=True, slug_field='entry_name')
    class Meta:
        model = ProteinConformation
        fields = ('protein', )


class ParentProteinFamilySerializer(serializers.ModelSerializer):
    class Meta:
        model = ProteinFamily
        fields = ('slug', 'name')


class ProteinFamilySerializer(serializers.ModelSerializer):
    parent = ParentProteinFamilySerializer(read_only=True)
    class Meta:
        model = ProteinFamily
        fields = ('slug', 'name', 'parent')


class SpeciesSerializer(serializers.ModelSerializer):
    class Meta:
        model = Species
        fields = ('latin_name', 'common_name')


class ResidueNumberingSchemeSerializer(serializers.ModelSerializer):
    class Meta:
        model = ResidueNumberingScheme
        fields = ('slug', 'short_name')


class ResidueGenericNumberSerializer(serializers.ModelSerializer):
    scheme = serializers.SlugRelatedField(read_only=True, slug_field='short_name')
    class Meta:
        model = ResidueGenericNumber
        fields = ('scheme', 'label')


class ResidueSerializer(serializers.ModelSerializer):
    protein_segment = serializers.StringRelatedField(read_only=True)
    display_generic_number = serializers.StringRelatedField(read_only=True)
    class Meta:
        model = Residue
        fields = ('sequence_number', 'amino_acid', 'protein_segment', 'display_generic_number')


class ResidueExtendedSerializer(serializers.ModelSerializer):
    protein_segment = serializers.StringRelatedField(read_only=True)
    display_generic_number = serializers.StringRelatedField(read_only=True)
    alternative_generic_numbers = ResidueGenericNumberSerializer(read_only=True, many=True)
    class Meta:
        model = Residue
        fields = ('sequence_number', 'amino_acid', 'protein_segment', 'display_generic_number',
            'alternative_generic_numbers')


class StructureSerializer(serializers.ModelSerializer):
    pdb_code = serializers.SlugRelatedField(read_only=True, slug_field='index')
    protein_conformation = ProteinFromConformationSerializer()
    class Meta:
        model = Structure
        fields = ('pdb_code', 'resolution', 'protein_conformation')