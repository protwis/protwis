from django.conf.urls import url
from django.views.decorators.cache import cache_page
from ligand import views
from django.views.generic import TemplateView

urlpatterns = [
    url(r'^$', cache_page(3600*24*7)(views.LigandBrowser.as_view()), name='ligand_browser'),
    url(r'^target/all/(?P<slug>[-\w]+)/$', views.TargetDetails, name='ligand_target_detail'),
    url(r'^target/compact/(?P<slug>[-\w]+)/$', views.TargetDetailsCompact, name='ligand_target_detail_compact'),
    url(r'^targets$', views.TargetDetails, name='ligand_target_detail'),
    url(r'^targets_compact', views.TargetDetailsCompact, name='ligand_target_detail_compact'),
    url(r'^targets_purchasable', views.TargetPurchasabilityDetails, name='ligand_target_detail_purchasable'),
    url(r'^(?P<ligand_id>[-\w]+)$', views.LigandDetails, name='ligand_detail'),
    url(r'^statistics', cache_page(3600*24*7)(views.LigandStatistics.as_view()), name='ligand_statistics'),
    url(r'^bias_statistics', cache_page(3600*24*7)(views.LigandBiasStatistics.as_view()), name='ligand_statistics'),

    url(r'^biasedbrowser/$', TemplateView.as_view(template_name='bias_browser.html')),
    url(r'^biasedsubtypesbrowser/$', TemplateView.as_view(template_name='bias_g_browser.html')),
    url(r'^bias/api/$', views.BiasAPI.as_view(), name='api-rest'),
    url(r'^biasg/api/$', views.GBiasAPI.as_view(), name='g-api-rest'),
    url(r'^biasedbrowser/experiment/(?P<pk>[-\w]+)/detail$', views.ExperimentEntryView.as_view()),
    url(r'^biasedsubtypesbrowser/experiment/(?P<pk>[-\w]+)/detail$', views.ExperimentEntryView.as_view()),
    url(r'^vendors$', views.test_link, name='test'),
    url(r'^browservendors$', views.BiasVendorBrowser.as_view(), name='browservendor'),

    url(r'^biasedpathways$', views.BiasPathways.as_view(), name='pathways'),
    url(r'^pathwaydata/(?P<pk>[-\w]+)/detail$', views.PathwayExperimentEntryView.as_view()),
]
