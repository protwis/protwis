from django.conf.urls import url
from django.views.decorators.cache import cache_page

from protein import views


urlpatterns = [
    #url(r'^$', views.BrowseSelection.as_view(), name='index'),
    url(r'^isoforms[/]?$', views.isoforms, name='isoforms'),
    url(r'^isoform_lookup$', views.AlignIsoformWildtype, name='isoform_lookup'),
    url(r'^$', cache_page(60*60*24*7)(views.BrowseSelection.as_view()), name='index'),
    url(r'^autocomplete', views.SelectionAutocomplete, name='autocomplete'),
    url(r'^gproteins', views.g_proteins, name='g_proteins'),
    url(r'^(?P<slug>[-\w]+)/$', views.detail, name='detail'),
]
