from django.conf.urls import patterns, url
from structure.views import *
from django.conf import settings

urlpatterns = patterns('',
    url(r'^$', StructureBrowser.as_view(), name='structure_browser'),
    url(r'^statistics$', StructureStatistics.as_view(), name='structure_statistics'),
    url(r'^crystal_statistics_by_year$', StructureStatistics.get_crystalized_receptors_data, name='structure_statistics'),
    url(r'^generic_numbering_index', GenericNumberingIndex.as_view(), name='generic_numbering'),
    url(r'^generic_numbering_results', GenericNumberingResults.as_view(), name='generic_numbering'),
    url(r'superposition_workflow_index', SuperpositionWorkflowIndex.as_view(), name='superposition_workflow'),
    url(r'superposition_workflow_selection', SuperpositionWorkflowSelection.as_view(), name='superposition_workflow'),
    url(r'superposition_workflow_results', SuperpositionWorkflowResults.as_view(), name='superposition_workflow'),
    url(r'fragment_superposition_index', FragmentSuperpositionIndex.as_view(), name='fragment_superposition'),
    url(r'fragment_superposition_results', FragmentSuperpositionResults.as_view(), name='fragment_superposition'),
    url(r'^output/(?P<outfile>\w+.\w{3})/(?P<replacement_tag>\w+)$', ServePdbOutfile, name='structural_tools_result'),
    url(r'^zipoutput/(?P<outfile>\w+.\w{3})/(?P<replacement_tag>\w+)$', ServeZipOutfile, name='structural_tools_result'), 
)