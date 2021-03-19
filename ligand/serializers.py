from django.db.models import IntegerField, FloatField
from django.db.models.functions import Cast
from django_filters import filters
from rest_framework_datatables.django_filters.filterset import DatatablesFilterSet
from rest_framework import serializers
from .models import AnalyzedExperiment

class AnalyzedExperimentFilter(DatatablesFilterSet):

    class_ = filters.CharFilter(
        field_name='receptor__family__parent__parent__parent_id',
        method='yadcf_multiple_choices_query', required=False
    )

    receptor_family = filters.CharFilter(
        field_name='receptor__family__parent_id',
        method='yadcf_multiple_choices_query', required=False
    )

    uniprot = filters.CharFilter(
        field_name='receptor_id',
        method='yadcf_multiple_choices_query', required=False
    )

    iuphar = filters.CharFilter(
        field_name='receptor_id',
        method='yadcf_multiple_choices_query', required=False
    )

    species = filters.CharFilter(
        field_name='receptor__species_id',
        method='yadcf_multiple_choices_query', required=False
    )

    endogenous_ligand = filters.CharFilter(
        field_name="endogenous_ligand_id",
        method='yadcf_multiple_choices_query', required=False
    )

    reference_ligand = filters.CharFilter(
        field_name="reference_ligand_id",
        method='yadcf_multiple_choices_query', required=False
    )

    ligand = filters.CharFilter(
        field_name="ligand_id",
        method='yadcf_multiple_choices_query', required=False
    )

    # We need to fix ordering for this field.
    vendor_quantity = filters.CharFilter(
        method='yadcf_range_filter_with_integer_cast', required=False
    )
    article_quantity = filters.CharFilter(
        method='yadcf_range_filter_with_integer_cast', required=False
    )
    labs_quantity = filters.CharFilter(
        method='yadcf_range_filter_with_integer_cast', required=False
    )

    primary = filters.CharFilter(
        # field_name='primary',
        method='transducers_multiple_choices_filter', required=False
    )
    secondary = filters.CharFilter(
        # field_name='primary',
        method='transducers_multiple_choices_filter', required=False
    )

    pathways_p1 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )
    pathways_p2 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )
    pathways_p3 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )
    pathways_p4 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )
    pathways_p5 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )

    opmodel_p2_p1 = filters.CharFilter(
        method='yadcf_range_filter_with_float_cast', required=False
    )
    opmodel_p3_p1 = filters.CharFilter(
        method='yadcf_range_filter_with_float_cast', required=False
    )
    opmodel_p4_p1 = filters.CharFilter(
        method='yadcf_range_filter_with_float_cast', required=False
    )
    opmodel_p5_p1 = filters.CharFilter(
        method='yadcf_range_filter_with_float_cast', required=False
    )

    lbf_p2_p1 = filters.CharFilter(
        method='yadcf_range_filter_with_float_cast', required=False
    )
    lbf_p3_p1 = filters.CharFilter(
        method='yadcf_range_filter_with_float_cast', required=False
    )
    lbf_p4_p1 = filters.CharFilter(
        method='yadcf_range_filter_with_float_cast', required=False
    )
    lbf_p5_p1 = filters.CharFilter(
        method='yadcf_range_filter_with_float_cast', required=False
    )

    potency_p2_p1 = filters.CharFilter(
        method='yadcf_range_filter', required=False
    )
    potency_p3_p1 = filters.CharFilter(
        method='yadcf_range_filter', required=False
    )
    potency_p4_p1 = filters.CharFilter(
        method='yadcf_range_filter', required=False
    )
    potency_p5_p1 = filters.CharFilter(
        method='yadcf_range_filter', required=False
    )

    activity_p1 = filters.CharFilter(
        method='yadcf_range_filter', required=False
    )
    activity_p2 = filters.CharFilter(
        method='yadcf_range_filter', required=False
    )
    activity_p3 = filters.CharFilter(
        method='yadcf_range_filter', required=False
    )
    activity_p4 = filters.CharFilter(
        method='yadcf_range_filter', required=False
    )
    activity_p5 = filters.CharFilter(
        method='yadcf_range_filter', required=False
    )

    emax_p1 = filters.CharFilter(
        method='yadcf_range_filter', required=False
    )
    emax_p2 = filters.CharFilter(
        method='yadcf_range_filter', required=False
    )
    emax_p3 = filters.CharFilter(
        method='yadcf_range_filter', required=False
    )
    emax_p4 = filters.CharFilter(
        method='yadcf_range_filter', required=False
    )
    emax_p5 = filters.CharFilter(
        method='yadcf_range_filter', required=False
    )

    tfactor_p1 = filters.CharFilter(
        method='yadcf_range_filter', required=False
    )
    tfactor_p2 = filters.CharFilter(
        method='yadcf_range_filter', required=False
    )
    tfactor_p3 = filters.CharFilter(
        method='yadcf_range_filter', required=False
    )
    tfactor_p4 = filters.CharFilter(
        method='yadcf_range_filter', required=False
    )
    tfactor_p5 = filters.CharFilter(
        method='yadcf_range_filter', required=False
    )

    assay_p1 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )
    assay_p2 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )
    assay_p3 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )
    assay_p4 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )
    assay_p5 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )

    cell_p1 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )
    cell_p2 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )
    cell_p3 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )
    cell_p4 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )
    cell_p5 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )

    time_p1 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )
    time_p2 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )
    time_p3 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )
    time_p4 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )
    time_p5 = filters.CharFilter(
        method='yadcf_multiple_choices_query', required=False
    )

    authors = filters.CharFilter(
        field_name='publication__authors',
        method='yadcf_multiple_choices_query', required=False
    )
    doi_reference = filters.CharFilter(
        field_name='publication__web_link_id',
        method='yadcf_multiple_choices_query', required=False
    )

    class Meta:
        queryset = AnalyzedExperiment.objects.all()
        fields = [
            'id','class_', 'receptor_family', 'uniprot', 'iuphar', 'species', 'endogenous_ligand',
            'reference_ligand', 'ligand_biased', 'vendor_quantity', 'article_quantity', 'ligand_labs',
            'primary', 'secondary',
            'pathways_p1', 'pathways_p2', 'pathways_p3', 'pathways_p4', 'pathways_p5',
            'opmodel_p2_p1', 'opmodel_p3_p1', 'opmodel_p4_p1', 'opmodel_p5_p1',
            'lbf_p2_p1', 'lbf_p3_p1', 'lbf_p4_p1', 'lbf_p5_p1',
            'potency_p2_p1', 'potency_p3_p1', 'potency_p4_p1', 'potency_p5_p1',
            'activity_p1', 'activity_p2', 'activity_p3', 'activity_p4', 'activity_p5',
            'emax_p1', 'emax_p2', 'emax_p3', 'emax_p4', 'emax_p5',
            'tfactor_p1', 'tfactor_p2', 'tfactor_p3', 'tfactor_p4', 'tfactor_p5',
            'assay_p1', 'assay_p2', 'assay_p3', 'assay_p4', 'assay_p5',
            'cell_p1', 'cell_p2', 'cell_p3', 'cell_p4', 'cell_p5',
            'time_p1', 'time_p2', 'time_p3', 'time_p4', 'time_p5',
            'authors', 'doi_reference','ligand_id','publication_link','quality_activity_p1',
            'quality_activity_p2','quality_activity_p3','quality_activity_p4','quality_activity_p5',
            'standard_type_p1','standard_type_p2','standard_type_p3','standard_type_p4','standard_type_p5',
            'ligand_source_id','ligand_source_type'
        ]
        read_only_fields = fields

    @staticmethod
    def yadcf_range_filter_with_integer_cast(queryset, field_name, value):
        min_value, max_value = value.split('-yadcf_delim-')
        try:
            min_value = float(min_value)
        except:# noqa
            min_value = None
        try:
            max_value = float(max_value)
        except:
            max_value = None

        queryset = queryset.annotate(**{f'{field_name}_as_integer': Cast(field_name, IntegerField())})
        if min_value is not None:
            queryset = queryset.filter(**{f'{field_name}_as_integer__gte': min_value})
        if max_value is not None:
            queryset = queryset.filter(**{f'{field_name}_as_integer__lte': max_value})
        return queryset

    @staticmethod
    def yadcf_range_filter_with_float_cast(queryset, field_name, value):
        min_value, max_value = value.split('-yadcf_delim-')

        try:
            min_value = float(min_value)
        except:
            min_value = None
        try:
            max_value = float(max_value)
        except:
            max_value = None

        queryset = queryset.annotate(**{f'{field_name}_as_float': Cast(field_name, FloatField())})
        if field_name in ['lbf_p2_p1', 'lbf_p3_p1', 'lbf_p4_p1', 'lbf_p5_p1']:
            # TODO: Come up with numeric values for:
            queryset = queryset.exclude(
                **{f'{field_name}__in': ['High Bias', 'Full Bias', 'Low Bias','']}
            )

        if min_value is not None:
            queryset = queryset.filter(**{f'{field_name}_as_float__gte': min_value})

        if max_value is not None:
            queryset = queryset.filter(**{f'{field_name}_as_float__lte': max_value})

        return queryset

    @staticmethod
    def yadcf_range_filter(queryset, field_name, value):
        min_value, max_value = value.split('-yadcf_delim-')

        try:
            min_value = float(min_value)
        except:
            min_value = None

        try:
            max_value = float(max_value)
        except:
            max_value = None

        if min_value:
            queryset = queryset.filter(**{f'{field_name}__gte': min_value})

        if max_value:
            queryset = queryset.filter(**{f'{field_name}__lte': max_value})

        return queryset

    @staticmethod
    def yadcf_multiple_choices_query(queryset, field_name, value):
        choices = value.replace('\\', '').split('|')

        return queryset.filter(**{f'{field_name}__in': choices})

    @staticmethod
    def transducers_multiple_choices_filter(queryset, field_name, value):
        choices = value.replace('\\', '').replace('_', ' ').split('|')

        return queryset.filter(**{f'{field_name}__in': choices})


class AnalyzedExperimentSerializer(serializers.ModelSerializer):
    # receptor
    class_ = serializers.SerializerMethodField()
    receptor_family = serializers.CharField(source='receptor.family.parent')
    uniprot = serializers.SerializerMethodField()
    iuphar = serializers.SerializerMethodField()
    species = serializers.CharField(source='receptor.species.common_name')
    endogenous_ligand = serializers.CharField()

    # Ligand
    reference_ligand = serializers.CharField()
    ligand = serializers.CharField(source='ligand.name')
    vendor_quantity = serializers.CharField()
    article_quantity = serializers.CharField()
    labs_quantity = serializers.CharField()

    # Receptor trunsducers
    primary = serializers.SerializerMethodField()
    secondary = serializers.SerializerMethodField()

    # Pathways
    pathways_p1 = serializers.CharField()
    pathways_p2 = serializers.CharField()
    pathways_p3 = serializers.CharField()
    pathways_p4 = serializers.CharField()
    pathways_p5 = serializers.CharField()

    # operational model
    opmodel_p2_p1 = serializers.CharField()
    opmodel_p3_p1 = serializers.CharField()
    opmodel_p4_p1 = serializers.CharField()
    opmodel_p5_p1 = serializers.CharField()

    # log bias factor
    lbf_p2_p1 = serializers.CharField()
    lbf_p3_p1 = serializers.CharField()
    lbf_p4_p1 = serializers.CharField()
    lbf_p5_p1 = serializers.CharField()

    # Potency ratio
    potency_p2_p1 = serializers.CharField()
    potency_p3_p1 = serializers.CharField()
    potency_p4_p1 = serializers.CharField()
    potency_p5_p1 = serializers.CharField()

    quality_activity_p1 = serializers.CharField()
    quality_activity_p2 = serializers.CharField()
    quality_activity_p3 = serializers.CharField()
    quality_activity_p4 = serializers.CharField()
    quality_activity_p5 = serializers.CharField()

    standard_type_p1 = serializers.CharField()
    standard_type_p2 = serializers.CharField()
    standard_type_p3 = serializers.CharField()
    standard_type_p4 = serializers.CharField()
    standard_type_p5 = serializers.CharField()

    # Potency
    activity_p1 = serializers.CharField()
    activity_p2 = serializers.CharField()
    activity_p3 = serializers.CharField()
    activity_p4 = serializers.CharField()
    activity_p5 = serializers.CharField()

    emax_p1 = serializers.CharField()
    emax_p2 = serializers.CharField()
    emax_p3 = serializers.CharField()
    emax_p4 = serializers.CharField()
    emax_p5 = serializers.CharField()

    tfactor_p1 = serializers.CharField()
    tfactor_p2 = serializers.CharField()
    tfactor_p3 = serializers.CharField()
    tfactor_p4 = serializers.CharField()
    tfactor_p5 = serializers.CharField()

    assay_p1 = serializers.CharField()
    assay_p2 = serializers.CharField()
    assay_p3 = serializers.CharField()
    assay_p4 = serializers.CharField()
    assay_p5 = serializers.CharField()

    cell_p1 = serializers.CharField()
    cell_p2 = serializers.CharField()
    cell_p3 = serializers.CharField()
    cell_p4 = serializers.CharField()
    cell_p5 = serializers.CharField()

    time_p1 = serializers.CharField()
    time_p2 = serializers.CharField()
    time_p3 = serializers.CharField()
    time_p4 = serializers.CharField()
    time_p5 = serializers.CharField()

    authors = serializers.CharField(source='publication.authors')
    doi_reference = serializers.CharField(source='publication.web_link.index')
    id = serializers.CharField()
    ligand_id =serializers.CharField(source='ligand.id')
    publication_link = serializers.CharField(source='publication.web_link')
    ligand_source_id= serializers.CharField()
    ligand_source_type= serializers.CharField()

    class Meta:
        model = AnalyzedExperiment
        fields = [
            'id','class_', 'receptor_family', 'uniprot', 'iuphar', 'species', 'endogenous_ligand',
            'reference_ligand', 'ligand', 'vendor_quantity', 'article_quantity', 'labs_quantity',
            'primary', 'secondary',
            'pathways_p1', 'pathways_p2', 'pathways_p3', 'pathways_p4', 'pathways_p5',
            'opmodel_p2_p1', 'opmodel_p3_p1', 'opmodel_p4_p1', 'opmodel_p5_p1',
            'lbf_p2_p1', 'lbf_p3_p1', 'lbf_p4_p1', 'lbf_p5_p1',
            'potency_p2_p1', 'potency_p3_p1', 'potency_p4_p1', 'potency_p5_p1',
            'activity_p1', 'activity_p2', 'activity_p3', 'activity_p4', 'activity_p5',
            'emax_p1', 'emax_p2', 'emax_p3', 'emax_p4', 'emax_p5',
            'tfactor_p1', 'tfactor_p2', 'tfactor_p3', 'tfactor_p4', 'tfactor_p5',
            'assay_p1', 'assay_p2', 'assay_p3', 'assay_p4', 'assay_p5',
            'cell_p1', 'cell_p2', 'cell_p3', 'cell_p4', 'cell_p5',
            'time_p1', 'time_p2', 'time_p3', 'time_p4', 'time_p5',
            'authors', 'doi_reference','ligand_id', 'publication_link','quality_activity_p1',
            'quality_activity_p2','quality_activity_p3','quality_activity_p4','quality_activity_p5',
            'standard_type_p1','standard_type_p2','standard_type_p3','standard_type_p4','standard_type_p5',
            'ligand_source_id','ligand_source_type'
        ]

    @staticmethod
    def get_class_(obj):
        class_receptor = obj.receptor.family.parent.parent.parent.name
        return class_receptor.replace('Class', '').strip()

    @staticmethod
    def get_iuphar(obj):
        iuphar_name = obj.receptor.name.split(' ', 1)[0].split('-adrenoceptor', 1)[0].strip()
        return iuphar_name

    @staticmethod
    def get_uniprot(obj):
        uniprot = obj.receptor.entry_short()

        return uniprot

    @staticmethod
    def get_primary(obj):
        return obj.primary.replace(' family,', '')

    @staticmethod
    def get_secondary(obj):
        return obj.secondary.replace(' family,', '')
