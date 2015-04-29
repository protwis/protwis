from django.conf.urls import patterns, url
from structure.views import StructureBrowser
from django.conf import settings

urlpatterns = patterns('',
    url(r'^$', StructureBrowser.as_view(), name='structure_browser'), 
)