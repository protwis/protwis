from django.conf.urls import url
from table_provider.views import GpcrStructureStatisticsSummaryTable, ConfigurationFactoryView
from django.views.decorators.cache import cache_page

urlpatterns = [
    ###################################
    # GPCR Structure Statistics Table #
    ###################################
    
    #DataTables/TableManager data source - Non-cached URL for server-side processing
    url(r'^structure_statistics_table/gpcr/?$', (GpcrStructureStatisticsSummaryTable.as_view()), name='structure_statistics_table'),
    #DataTables/TableManager data source - Cached URL for client side processing
    url(r'^structure_statistics_table_cached/gpcr/?$', cache_page(60*60*24)(GpcrStructureStatisticsSummaryTable.as_view()), name='structure_statistics_table_cached'),  

    #AJAX source for dropdown select options for filtering various columns in the GPCR Structure Statistics Table
    url(r'^structure_statistics_table/gpcr/options/(?P<column>[^/]+)/?$', (GpcrStructureStatisticsSummaryTable.get_select_options), name='structure_statistics_table_options'),
    
    ######################################
    # TableManager configuration factory #
    ######################################

    url(r'^configuration/(?P<configuration_type>[^/]+)/(?P<table_name>[^/]+)/(?P<configuration_variant>[^/]+)/?$', ConfigurationFactoryView.as_view(), name='configuration_factory'),
]