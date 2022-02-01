from django.conf.urls import url
from django.urls import path
from django.views.decorators.cache import cache_page
from ligand import views

urlpatterns = [
    #  url(r'^browser$', cache_page(3600*24*7)(views.LigandBrowser.as_view()), name='ligand_browser'),
    url(r'^$', views.LigandTargetSelection.as_view(), name='ligand_selection'),
    url(r'^browser$', views.LigandBrowser.as_view(), name='ligand_browser'),

    url(r'^target/all/(?P<slug>[-\w]+)/$', views.TargetDetails, name='ligand_target_detail'),
    url(r'^target/compact/(?P<slug>[-\w]+)/$', views.TargetDetailsCompact, name='ligand_target_detail_compact'),
    url(r'^targets$', views.TargetDetails, name='ligand_target_detail'),
    url(r'^targets_compact', views.TargetDetailsCompact, name='ligand_target_detail_compact'),
    url(r'^targets_purchasable', views.TargetPurchasabilityDetails, name='ligand_target_detail_purchasable'),
    url(r'^(?P<ligand_id>[-\w]+)/details$', views.LigandDetails, name='ligand_detail'),
    url(r'^statistics', cache_page(3600*24*7)(views.LigandStatistics.as_view()), name='ligand_statistics'),
    # url(r'^statistics', views.LigandStatistics.as_view(), name='ligand_statistics'),

    # BIASED LIGANDS
    url(r'^bias_statistics', cache_page(3600*24*7)(views.LigandStatistics.as_view(page='ligand_bias')), name='ligand_statistics'),
    # url(r'^bias_statistics', views.LigandStatistics.as_view(page='ligand_bias'), name='ligand_statistics'),
    path('emax_rank_order_selection', views.RankOrderSelection.as_view(), name='emax_ro_selection'),
    path('emax_rankorder', views.BiasedRankOrderOnTheFly.as_view(), name='biased_rank_order'),
    path('emax_path_profiles_selection', views.EmaxPathProfileSelection.as_view(), name='ema_pp_selection'),
    path('emax_path_profiles', views.BiasedRankOrderOnTheFly.as_view(page='pathwayprofiles'), name='biased_rank_order'),
    path('tau_rank_order_selection', views.TauRankOrderSelection.as_view(), name='tau_ro_selection'),
    path('tau_rankorder', views.BiasedRankOrderOnTheFly.as_view(label='tau'), name='biased_rank_order'),
    path('tau_path_profiles_selection', views.TauPathProfileSelection.as_view(), name='tau_pp_selection'),
    path('tau_path_profiles', views.BiasedRankOrderOnTheFly.as_view(page='pathwayprofiles', label='tau'), name='biased_rank_order'),

    # SUBTYPE BIASED
    url(r'^subtype_statistics', cache_page(3600*24*7)(views.LigandStatistics.as_view(page='subtype')), name='ligand_statistics'),
    # url(r'^subtype_statistics', views.LigandStatistics.as_view(page='subtype'), name='ligand_statistics'),
    path('subtype_emax_rankorder_selection', views.EmaxSubtypeRankOrderSelection.as_view(), name='ema_subtype_ro_selection'),
    path('subtype_emax_rankorder', views.BiasedRankOrderOnTheFly.as_view(subtype=True), name='biased_rank_order'),
    path('subtype_emax_path_profiles_selection', views.EmaxSubtypePathProfilesSelection.as_view(), name='ema_pathpref_pathprof_selection'),
    path('subtype_emax_path_profiles', views.BiasedRankOrderOnTheFly.as_view(page='pathwayprofiles', subtype=True), name='biased_rank_order'),
    path('subtype_tau_rank_order_selection', views.TauSubtypeRankOrderSelection.as_view(), name='tau_subtype_ro_selection'),
    path('subtype_tau_rankorder', views.BiasedRankOrderOnTheFly.as_view(label='tau', subtype=True), name='biased_rank_order'),
    path('subtype_tau_path_profiles_selection', views.TauSubtypePathProfileSelection.as_view(), name='tau_subtype_pp_selection'),
    path('subtype_tau_path_profiles', views.BiasedRankOrderOnTheFly.as_view(page='pathwayprofiles', subtype=True, label='tau'), name='biased_rank_order'),

    # PATHWAY PREFERRED
    url(r'^path_preference_statistics', cache_page(3600*24*7)(views.LigandStatistics.as_view(page='pathway_pref')), name='ligand_statistics'),
    # url(r'^path_preference_statistics', views.LigandStatistics.as_view(page='pathway_pref'), name='ligand_statistics'),
    path('path_preference_emax_rankorder_selection', views.EmaxPathPrefRankOrderSelection.as_view(), name='ema_pathpref_ro_selection'),
    path('path_preference_emax_rankorder', views.BiasedRankOrderOnTheFly.as_view(pathway=True), name='biased_rank_order'),
    path('path_preference_emax_path_profiles_selection', views.EmaxPathPrefPathProfilesSelection.as_view(), name='ema_pathpref_pathprof_selection'),
    path('path_preference_emax_path_profiles', views.BiasedRankOrderOnTheFly.as_view(page='pathwayprofiles', pathway=True), name='biased_rank_order'),

    url(r'^(?P<pk>[-\w]+)/info$', views.LigandInformationView.as_view()),
    #Browsers Cached
    url(r'^biased/$', views.CachedOTFBiasBrowser, name='bias_browser-list'),
    url(r'^userbiased/$', views.CachedOTFBiasBrowserUser, name='bias_browser-list'),
    url(r'^biasedsubtypes/$',views.CachedOTFBiasSubtypeBrowser, name='bias_browser-list'),
    url(r'^userbiasedsubtypes/$',views.CachedOTFBiasSubtypeBrowserUser, name='bias_browser-list'),
    url(r'^pathwaypreference/$',views.CachedOTFPathwayPrefBrowser, name='bias_browser-list'),
    #User selected calculations
    url(r'^userselectionbiased', views.UserBiased.as_view(), name='bias_browser-list'),
    url(r'^userselectionbiasedsubtype', views.UserBiasedSubtype.as_view(), name='bias_browser-list'),

    url(r'^biasedbrowser',views.BiasTargetSelection.as_view(), name='bias_browser-list1'),
    url(r'^pathwaypreferencebrowser',views.BiasPredictionTargetSelection.as_view(), name='bias_browser-list1'),
    url(r'^biasedsubtypesbrowser',views.BiasGTargetSelection.as_view(), name='bias_browser-list1'),

    url(r'^biased/experiment/(?P<pk>[-\w]+)/detail$', views.ExperimentEntryView.as_view()),
    url(r'^biasedsubtypes/experiment/(?P<pk>[-\w]+)/detail$', views.ExperimentEntryView.as_view()),
    url(r'^biasedpredicted/experiment/(?P<pk>[-\w]+)/detail$', views.ExperimentEntryView.as_view()),

    url(r'^vendors$', views.test_link, name='test'),
    url(r'^browservendors$', views.BiasVendorBrowser.as_view(), name='browservendor'),
    url(r'^biasedpathways$', views.BiasPathways.as_view(), name='pathways'),
    url(r'^pathwaydata/(?P<pk>[-\w]+)/detail$', views.PathwayExperimentEntryView.as_view()),

]
