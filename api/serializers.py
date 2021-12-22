from rest_framework import serializers

from interaction.models import ResidueFragmentInteraction
from mutation.models import MutationRaw
from protein.models import Protein, ProteinConformation, ProteinFamily, Species, ProteinSource, ProteinSegment
from residue.models import Residue, ResidueNumberingScheme, ResidueGenericNumber
from structure.models import Structure


class ProteinSerializer(serializers.ModelSerializer):
    family = serializers.SlugRelatedField(read_only=True, slug_field='slug')
    species = serializers.StringRelatedField(read_only=True)
    source = serializers.StringRelatedField(read_only=True)
    residue_numbering_scheme = serializers.SlugRelatedField(read_only=True, slug_field='short_name')
    genes = serializers.StringRelatedField(many=True)
    class Meta:
        model = Protein
        fields = ('entry_name', 'name', 'accession', 'family', 'species', 'source', 'residue_numbering_scheme',
            'sequence','genes')


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


class ReceptorListSerializer(serializers.ModelSerializer):
    family = serializers.StringRelatedField(read_only=True)
    endogenous_ligands = serializers.SlugRelatedField(
        many=True,
        read_only=True,
        slug_field='name'
    )
    species = serializers.StringRelatedField(read_only=True)
    class Meta:
        model = Protein
        fields = ('entry_name', 'name', 'accession', 'family', 'endogenous_ligands', 'species')


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


class StructureLigandInteractionSerializer(serializers.ModelSerializer):
    pdb_code = serializers.ReadOnlyField(source='structure_ligand_pair.structure.pdb_code.index')
    ligand_name = serializers.ReadOnlyField(source='structure_ligand_pair.ligand.name')
    amino_acid = serializers.ReadOnlyField(source='fragment.residue.amino_acid')
    sequence_number = serializers.IntegerField(read_only=True, source='fragment.residue.sequence_number')
    display_generic_number = serializers.ReadOnlyField(source='fragment.residue.display_generic_number.label')
    interaction_type = serializers.ReadOnlyField(source='interaction_type.name')
    class Meta:
        model = ResidueFragmentInteraction
        fields = ('pdb_code', 'ligand_name', 'amino_acid', 'sequence_number',
                  'display_generic_number', 'interaction_type')


class MutationSerializer(serializers.ModelSerializer):
    class Meta:
        model = MutationRaw
        fields = ('reference', 'protein', 'mutation_pos', 'mutation_from', 'mutation_to',
            'ligand_name', 'ligand_idtype', 'ligand_id', 'ligand_class',
            'exp_type', 'exp_func',  'exp_wt_value',  'exp_wt_unit','exp_mu_effect_sign', 'exp_mu_effect_type', 'exp_mu_effect_value',
            'exp_fold_change',
            'exp_mu_effect_qual', 'exp_mu_effect_ligand_prop',  'exp_mu_ligand_ref', 'opt_receptor_expression', 'opt_basal_activity',
            'opt_gain_of_activity', 'opt_ligand_emax', 'opt_agonist')
