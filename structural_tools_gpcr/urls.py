from django.conf.urls import patterns, url
from structural_tools_gpcr.views import *
#from structural_tools_gpcr.views import GenericNumberingIndex, GenericNumberingResults, SuperpositionWorkflowIndex, SuperpositionWorkflowSelection, SuperpositionWorkflowResults, ServePdbOutfile


urlpatterns = patterns('',
    url(r'^generic_numbering_index', GenericNumberingIndex.as_view(), name='generic_numbering'),
    url(r'^generic_numbering_results', GenericNumberingResults.as_view(), name='generic_numbering'),
    url(r'superposition_workflow_index', SuperpositionWorkflowIndex.as_view(), name='superposition_workflow'),
    url(r'superposition_workflow_selection', SuperpositionWorkflowSelection.as_view(), name='superposition_workflow'),
    url(r'superposition_workflow_results', SuperpositionWorkflowResults.as_view(), name='superposition_workflow'),
    url(r'^output/(?P<outfile>\w+.\w{3})/(?P<replacement_tag>\w+)$', ServePdbOutfile, name='structural_tools_result'),
)