from django.conf.urls import patterns, url
from structural_tools_gpcr.views import GenericNumberingIndex, GenericNumberingResults, SuperpositionWorkflowIndex, ServePdbOutfile


urlpatterns = patterns('',
    url(r'^generic_numbering_index', GenericNumberingIndex.as_view(), name='generic_numbering'),
    url(r'^generic_numbering_results', GenericNumberingResults.as_view(), name='generic_numbering'),
    url(r'superposition_workflow_index', SuperpositionWorkflowIndex.as_view(), name='superposition_workflow'),
    url(r'^output/(?P<outfile>\w+.\w{3})$', ServePdbOutfile, name='structural_tools_result'),
)