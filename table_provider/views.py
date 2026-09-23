#Django Framework imports
from django.http import JsonResponse
from django.core.exceptions import FieldError

#Rest Framework imports
from rest_framework import generics
from rest_framework.response import Response
from rest_framework.pagination import LimitOffsetPagination

#Table provider imports
from table_provider.configuration_provider.configuration_factory import ColumnConfigurationFactory, DataTablesConfigurationFactory
from table_provider.models import GpcrStructureStatisticsTable
from table_provider.serializers import GpcrStatisticsSummarySerializer
from table_provider.serverside.filters import FilterSet
from table_provider.serverside.querystringprocessing import DatatablesQueryStringProcessor
from structure.tables.structure_coverage_statistics_query import GpcrStructureCoverageStatisticsQuery

#Pagination handler for dataTables server side processing
class DataTablesLimitOffsetPagination(LimitOffsetPagination):
    default_limit = 30
    max_limit = 200
    limit_query_param = "length"
    offset_query_param = "start"

class GpcrStructureStatisticsSummaryTable(generics.ListCreateAPIView):
    """API endpoint for summary statistics for GPCR entries"""

    #############
    # API setup #
    #############

    serializer_class = GpcrStatisticsSummarySerializer
    pagination_class = DataTablesLimitOffsetPagination

    def get_queryset(self):
        filters = FilterSet(self.request, GpcrStatisticsSummarySerializer)
        queryset = self.statistics_table_fetch(filters)
        return queryset

    def list(self, request, *args, **kwargs):

        table_qs = GpcrStructureStatisticsTable.objects.all()
        
        if not table_qs.exists():
            self.build_statistics_summary_table()
        
        unfilteredCount = table_qs.count()

        queryset = self.get_queryset()
        url_params = DatatablesQueryStringProcessor(self.request).query_set

        # server side datatables processing response format
        if 'draw' in url_params:

            response = {}

            filteredCount = queryset.count()

            page = self.paginate_queryset(queryset)
            rows = page if page is not None else queryset

            serializer = self.get_serializer(rows, many=True)

            response["draw"] = url_params["draw"]
            response["recordsTotal"] = unfilteredCount
            response["recordsFiltered"] = filteredCount
            response["data"] = serializer.data

            return JsonResponse(response)

        else:
            # client side datatables processing response format
            serializer = self.get_serializer(queryset, many=True)
            # for consistency with server side processing response format,
            # wrap data in "data" key, even though client side processing doesn't require this
            response = {}
            response["data"] = serializer.data
            return Response(response)

    def statistics_table_fetch(self, filters):

        queryset = GpcrStructureStatisticsTable.objects.all()

        for query_filter in filters.get_filters():
            queryset = queryset.filter(query_filter.format_query())

        if filters.get_ordering():
            queryset = queryset.order_by(*filters.get_ordering())

        return queryset

    @staticmethod
    def get_select_options(request, column):
        #Map any complex table columns back to their primary data column for the query (e.g. weblinks with custom renderers)
        if column in GpcrStatisticsSummarySerializer.Meta.serializer_method_to_filter_field_map:
            column = GpcrStatisticsSummarySerializer.Meta.serializer_method_to_filter_field_map[column]

        try:
            queryset = GpcrStructureStatisticsTable.objects.values_list(column, flat=True).distinct().order_by(column)
        except FieldError:
            return JsonResponse({'error': f'Invalid column requested: {column}'}, status=400)

        return JsonResponse(list(queryset), safe=False)

    def build_statistics_summary_table(self):
        stat_data_model = [ GpcrStructureStatisticsTable(**data_item) for data_item in GpcrStructureCoverageStatisticsQuery() ]
        GpcrStructureStatisticsTable.objects.bulk_create(stat_data_model, batch_size=10000)


class ConfigurationFactoryView(generics.ListCreateAPIView):
    """Request handler for TableManager configuration factory.

    Loads configuration file from disk and returns a JSON response containing the requested configuration for a specific table and variant.
    """

    def get(self, request, *args, **kwargs):
        configuration_type = self.kwargs.get('configuration_type')
        table_name = self.kwargs.get('table_name')
        configuration_variant = self.kwargs.get('configuration_variant')

        if configuration_type == 'column':
            config_factory = ColumnConfigurationFactory(table_name, configuration_variant)
        elif configuration_type == 'datatable':
            config_factory = DataTablesConfigurationFactory(table_name, configuration_variant)
        else:
            return JsonResponse({'error': 'Invalid configuration type requested. Config: ' + configuration_type}, status=400)

        try:
            config = config_factory.fetch()
        except FileNotFoundError as e:
            return JsonResponse({'error': str(e)}, status=400)

        return JsonResponse(config, safe=False)