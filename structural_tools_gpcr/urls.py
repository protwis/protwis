from django.conf.urls import patterns, url

from structural_tools_gpcr import views


urlpatterns = patterns('',
    url(r'^gn_uploadfile', views.GenericNumberingStart.as_view(), name='gn_uploadfile'),
)