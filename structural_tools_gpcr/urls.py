from django.conf.urls import patterns, url

from structural_tools_gpcr.views import GenericNumberingStart


urlpatterns = patterns('',
    url(r'^gn_uploadfile', GenericNumberingStart.as_view(), name='gn_uploadfile'),
)