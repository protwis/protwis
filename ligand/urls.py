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
    path('emax_rank_order_selection', views.BiasedSignallingSelection.as_view(way='EmaxRankOrder'), name='emax_ro_selection'),
    path('emax_rankorder', views.BiasedSignallingOnTheFlyCalculation.as_view(), name='biased_rank_order'),
    path('emax_rankorder_path_bias', views.BiasedSignallingOnTheFlyCalculation.as_view(balanced=True), name='biased_rank_order'),
    path('userbiased_emax_rank_order', views.BiasedSignallingOnTheFlyCalculation.as_view(user=True), name='biased_rank_order'),
    path('emax_path_profiles_selection', views.BiasedSignallingSelection.as_view(way='EmaxPathProfile'), name='ema_pp_selection'),
    path('emax_path_profiles', views.BiasedSignallingOnTheFlyCalculation.as_view(page='pathwayprofiles'), name='biased_rank_order'),
    path('emax_path_profiles_path_bias', views.BiasedSignallingOnTheFlyCalculation.as_view(page='pathwayprofiles', balanced=True), name='biased_rank_order'),
    path('userbiased_emax_path_profile', views.BiasedSignallingOnTheFlyCalculation.as_view(page='pathwayprofiles', user=True), name='biased_path_profile'),
    path('tau_rank_order_selection', views.BiasedSignallingSelection.as_view(way='TauRankOrder'), name='tau_ro_selection'),
    path('tau_rankorder', views.BiasedSignallingOnTheFlyCalculation.as_view(label='tau'), name='biased_rank_order'),
    path('tau_rankorder_path_bias', views.BiasedSignallingOnTheFlyCalculation.as_view(label='tau', balanced=True), name='biased_rank_order'),
    path('userbiased_tau_rank_order', views.BiasedSignallingOnTheFlyCalculation.as_view(label='tau', user=True), name='biased_rank_order'),
    path('tau_path_profiles_selection', views.BiasedSignallingSelection.as_view(way='TauPathProfile'), name='tau_pp_selection'),
    path('tau_path_profiles', views.BiasedSignallingOnTheFlyCalculation.as_view(page='pathwayprofiles', label='tau'), name='biased_rank_order'),
    path('tau_path_profiles_path_bias', views.BiasedSignallingOnTheFlyCalculation.as_view(page='pathwayprofiles', label='tau', balanced=True), name='biased_rank_order'),
    path('userbiased_tau_path_profile', views.BiasedSignallingOnTheFlyCalculation.as_view(page='pathwayprofiles', label='tau', user=True), name='biased_rank_order'),
    # SUBTYPE BIASED
    url(r'^subtype_statistics', cache_page(3600*24*7)(views.LigandStatistics.as_view(page='subtype')), name='ligand_statistics'),
    # url(r'^subtype_statistics', views.LigandStatistics.as_view(page='subtype'), name='ligand_statistics'),
    path('subtype_emax_rankorder_selection', views.BiasedSignallingSelection.as_view(subtype=True, way='EmaxRankOrderSubtype'), name='ema_subtype_ro_selection'),
    path('subtype_emax_rankorder', views.BiasedSignallingOnTheFlyCalculation.as_view(subtype=True), name='biased_rank_order'),
    path('subtype_emax_rankorder_path_bias', views.BiasedSignallingOnTheFlyCalculation.as_view(subtype=True, balanced=True), name='biased_rank_order'),
    path('userbiasedsubtypes_emax_rank_order', views.BiasedSignallingOnTheFlyCalculation.as_view(subtype=True, user=True), name='biased_rank_order'),
    path('subtype_emax_path_profiles_selection', views.BiasedSignallingSelection.as_view(subtype=True, way='EmaxPathProfileSubtype'), name='ema_pathpref_pathprof_selection'),
    path('subtype_emax_path_profiles', views.BiasedSignallingOnTheFlyCalculation.as_view(page='pathwayprofiles', subtype=True), name='biased_rank_order'),
    path('subtype_emax_path_profiles_path_bias', views.BiasedSignallingOnTheFlyCalculation.as_view(page='pathwayprofiles', subtype=True, balanced=True), name='biased_rank_order'),
    path('userbiasedsubtypes_emax_path_profile', views.BiasedSignallingOnTheFlyCalculation.as_view(page='pathwayprofiles', subtype=True, user=True), name='biased_path_profile'),
    path('subtype_tau_rank_order_selection', views.BiasedSignallingSelection.as_view(subtype=True, way='TauRankOrderSubtype'), name='tau_subtype_ro_selection'),
    path('subtype_tau_rankorder', views.BiasedSignallingOnTheFlyCalculation.as_view(label='tau', subtype=True), name='biased_rank_order'),
    path('subtype_tau_rankorder_path_bias', views.BiasedSignallingOnTheFlyCalculation.as_view(label='tau', subtype=True, balanced=True), name='biased_rank_order'),
    path('userbiasedsubtypes_tau_rank_order', views.BiasedSignallingOnTheFlyCalculation.as_view(label='tau', subtype=True, user=True), name='biased_rank_order'),
    path('subtype_tau_path_profiles_selection', views.BiasedSignallingSelection.as_view(subtype=True, way='TauPathProfileSubtype'), name='tau_subtype_pp_selection'),
    path('subtype_tau_path_profiles', views.BiasedSignallingOnTheFlyCalculation.as_view(page='pathwayprofiles', subtype=True, label='tau'), name='biased_rank_order'),
    path('subtype_tau_path_profiles_path_bias', views.BiasedSignallingOnTheFlyCalculation.as_view(page='pathwayprofiles', subtype=True, label='tau', balanced=True), name='biased_rank_order'),
    path('userbiasedsubtypes_tau_path_profile', views.BiasedSignallingOnTheFlyCalculation.as_view(page='pathwayprofiles', subtype=True, label='tau', user=True), name='biased_rank_order'),
    # PATHWAY PREFERRED
    url(r'^path_preference_statistics', cache_page(3600*24*7)(views.LigandStatistics.as_view(page='pathway_pref')), name='ligand_statistics'),
    # url(r'^path_preference_statistics', views.LigandStatistics.as_view(page='pathway_pref'), name='ligand_statistics'),
    path('path_preference_emax_rankorder_selection', views.BiasedSignallingSelection.as_view(pathway=True, way='EmaxRankOrderPathway'), name='ema_pathpref_ro_selection'),
    path('path_preference_emax_rankorder', views.BiasedSignallingOnTheFlyCalculation.as_view(pathway=True), name='biased_rank_order'),
    path('path_preference_emax_path_profiles_selection', views.BiasedSignallingSelection.as_view(pathway=True, way='EmaxPathProfilePathway'), name='ema_pathpref_pathprof_selection'),
    path('path_preference_emax_path_profiles', views.BiasedSignallingOnTheFlyCalculation.as_view(page='pathwayprofiles', pathway=True), name='biased_rank_order'),

    url(r'^(?P<pk>[-\w]+)/info$', views.LigandInformationView.as_view()),
    #Browsers Cached
    url(r'^userbiased/$', views.CachedOTFBiasBrowserUser, name='bias_browser-list'),
    url(r'^userbiasedsubtypes/$',views.CachedOTFBiasSubtypeBrowserUser, name='bias_browser-list'),
    url(r'^biased/$', views.CachedOTFBiasBrowser, name='bias_browser-list'),
    url(r'^biasedsubtypes/$',views.CachedOTFBiasSubtypeBrowser, name='bias_browser-list'),
    url(r'^pathwaypreference/$',views.CachedOTFPathwayPrefBrowser, name='bias_browser-list'),
    #User selected calculations
    #Biased Family
    path('userselectionbiased', views.UserBiased.as_view(way='Browser'), name='bias_browser-list'),
    path('userselectionbiased_emax_rank_order', views.UserBiased.as_view(way='EmaxRankOrder'), name='userbiased_emax_rank_order'),
    path('userselectionbiased_tau_rank_order', views.UserBiased.as_view(way='TauRankOrder'), name='userbiased_tau_rank_order'),
    path('userselectionbiased_emax_path_profile', views.UserBiased.as_view(way='EmaxPathProfile'), name='userbiased_emax_path_profile'),
    path('userselectionbiased_tau_path_profile', views.UserBiased.as_view(way='TauPathProfile'), name='userbiased_tau_path_profile'),
    #Biased subtype
    path('userselectionbiasedsubtype', views.UserBiased.as_view(subtype=True, way='BrowserSubtype'), name='bias_browser-list'),
    path('userselectionbiasedsubtype_emax_rank_order', views.UserBiased.as_view(subtype=True, way='EmaxRankOrderSubtype'), name='userbiasedsubtype_emax_rank_order'),
    path('userselectionbiasedsubtype_tau_rank_order', views.UserBiased.as_view(subtype=True, way='TauRankOrderSubtype'), name='userbiasedsubtype_tau_rank_order'),
    path('userselectionbiasedsubtype_emax_path_profile', views.UserBiased.as_view(subtype=True, way='EmaxPathProfileSubtype'), name='userbiasedsubtype_emax_path_profile'),
    path('userselectionbiasedsubtype_tau_path_profile', views.UserBiased.as_view(subtype=True, way='TauPathProfilesubtype'), name='userbiasedsubtype_tau_path_profile'),

    #Balanced reference calculations
    url(r'^pathwaybiased', views.CachedOTFBalancedBrowser, name='bias_browser-list'),
    url(r'^pathwaybiasedsubtype', views.CachedOTFBalancedSubtypeBrowser, name='bias_browser-list'),
    #Browsers selection pages
    url(r'^biasedbrowser',views.BiasedSignallingSelection.as_view(way='Browser'), name='bias_browser-list1'),
    url(r'^pathwaypreferencebrowser',views.BiasedSignallingSelection.as_view(pathway=True, way='BrowserPathway'), name='bias_browser-list1'),
    url(r'^biasedsubtypesbrowser',views.BiasedSignallingSelection.as_view(subtype=True, way='BrowserSubtype'), name='bias_browser-list1'),

    url(r'^vendors$', views.test_link, name='test'),
    url(r'^browservendors$', views.BiasVendorBrowser.as_view(), name='browservendor'),
    url(r'^biasedpathways$', views.BiasPathways.as_view(), name='pathways'),
    url(r'^pathwaydata/(?P<pk>[-\w]+)/detail$', views.PathwayExperimentEntryView.as_view()),

]
