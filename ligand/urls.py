from django.conf.urls import url
from django.views.decorators.cache import cache_page
from ligand import views

urlpatterns = [
    url(r'^$', cache_page(3600*24*7)(views.LigandBrowser.as_view()), name='ligand_browser'),
    url(r'^target/all/(?P<slug>[-\w]+)/$', views.TargetDetails, name='ligand_target_detail'),
    url(r'^target/compact/(?P<slug>[-\w]+)/$', views.TargetDetailsCompact, name='ligand_target_detail_compact'),
    url(r'^targets$', views.TargetDetails, name='ligand_target_detail'),
    url(r'^targets_compact', views.TargetDetailsCompact, name='ligand_target_detail_compact'),
    url(r'^targets_purchasable', views.TargetPurchasabilityDetails, name='ligand_target_detail_purchasable'),
    url(r'^(?P<ligand_id>[-\w]+)/details$', views.LigandDetails, name='ligand_detail'),
    url(r'^statistics', cache_page(3600*24*7)(views.LigandStatistics.as_view()), name='ligand_statistics'),
    url(r'^bias_statistics', views.LigandBiasStatistics.as_view(), name='ligand_statistics'),
    url(r'^(?P<pk>[-\w]+)/info$', views.LigandInformationView.as_view()),

    url(r'^biased/$', views.CachedBiasBrowser, name='bias_browser-list'),
    #url(r'^biased/$', views.BiasBrowser.as_view(), name='bias_browser-list'),
    url(r'^biasedsubtypes/$',views.CachedBiasGBrowser, name='bias_browser-subtype'),
    #url(r'^biasedsubtypes/$',views.BiasGBrowser.as_view(), name='bias_browser-list'),
    url(r'^biasedbrowser',views.BiasTargetSelection.as_view(), name='bias_browser-list1'),
    url(r'^biasedsubtypesbrowser',views.BiasGTargetSelection.as_view(), name='bias_browser-list1'),

    url(r'^biased/experiment/(?P<pk>[-\w]+)/detail$', views.ExperimentEntryView.as_view()),
    url(r'^biasedsubtypes/experiment/(?P<pk>[-\w]+)/detail$', views.ExperimentEntryView.as_view()),

    url(r'^vendors$', views.test_link, name='test'),
    url(r'^browservendors$', views.BiasVendorBrowser.as_view(), name='browservendor'),
    url(r'^biasedpathways$', views.BiasPathways.as_view(), name='pathways'),
    url(r'^pathwaydata/(?P<pk>[-\w]+)/detail$', views.PathwayExperimentEntryView.as_view()),

]
