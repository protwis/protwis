from django.conf.urls import patterns, url
from structural_tools_gpcr.views import GenericNumberingIndex, GenericNumberingResults, ServeOutfile


urlpatterns = patterns('',
    url(r'^gn_uploadfile', GenericNumberingIndex.as_view(), name='gn_uploadfile'),
    url(r'^gn_results', GenericNumberingResults.as_view(), name='gn_results'),
    url(r'^/output/(?P<outfile>\w+)$', ServeOutfile, name='structural_tools_result'),
)