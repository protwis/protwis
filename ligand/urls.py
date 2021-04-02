from django.conf.urls import url
from django.urls import path
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
    # url(r'^(?P<ligand_id>[-\w]+)/$',LigandDetails, name='ligand_detail'),
    url(r'^statistics', cache_page(3600*24*7)(LigandStatistics.as_view()), name='ligand_statistics'),

    url(r'^biasedbrowser/$', BiasBrowser.as_view(), name='bias_browser-list'),
    url(r'^biasedsubtypesbrowser/$', TemplateView.as_view(template_name='bias_g_browser.html')),

    url(r'^biasedbrowser/experiment/(?P<pk>[-\w]+)/detail$', ExperimentEntryView.as_view()),
    url(r'^biasedsubtypesbrowser/experiment/(?P<pk>[-\w]+)/detail$', ExperimentEntryView.as_view()),

    url(r'^vendors$', test_link, name='test'),
    url(r'^browservendors$', BiasVendorBrowser.as_view(), name='browservendor'),
    url(r'^biasedpathways$', BiasPathways.as_view(), name='pathways'),
    url(r'^pathwaydata/(?P<pk>[-\w]+)/detail$', PathwayExperimentEntryView.as_view()),


]
