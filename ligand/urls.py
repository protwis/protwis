from django.conf.urls import url
from django.views.decorators.cache import cache_page
from django.views.generic import TemplateView
from ligand.views import *

urlpatterns = [
    url(r'^$', LigandBrowser.as_view(), name='ligand_browser'),
    url(r'^target/(?P<slug>[-\w]+)/$',TargetDetails, name='ligand_target_detail'),
    url(r'^(?P<ligand_id>[-\w]+)/$',LigandDetails, name='ligand_detail'),
    url(r'^statistics', LigandStatistics.as_view(), name='ligand_statistics')
]
