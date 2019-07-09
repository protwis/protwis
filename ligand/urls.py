from django.conf.urls import url
from django.views.decorators.cache import cache_page
from django.views.generic import TemplateView
from ligand.views import *

urlpatterns = [
    url(r'^$', cache_page(3600*24*7)(LigandBrowser.as_view()), name='ligand_browser'),
    url(r'^target/all/(?P<slug>[-\w]+)/$',TargetDetails, name='ligand_target_detail'),
    url(r'^target/compact/(?P<slug>[-\w]+)/$',TargetDetailsCompact, name='ligand_target_detail_compact'),
    url(r'^targets$',TargetDetails, name='ligand_target_detail'),
    url(r'^targets_compact',TargetDetailsCompact, name='ligand_target_detail_compact'),
    url(r'^targets_purchasable',TargetPurchasabilityDetails, name='ligand_target_detail_purchasable'),
    url(r'^(?P<ligand_id>[-\w]+)/$',LigandDetails, name='ligand_detail'),
    url(r'^statistics', cache_page(3600*24*7)(LigandStatistics.as_view()), name='ligand_statistics'),
    #url(r'^bias', output_bias , name='biased_ligands'),
    ##url(r'^list', bias_list , name='biased_ligands_list'),
    url(r'^browser', bias_browser , name='bias_browser'),


]
