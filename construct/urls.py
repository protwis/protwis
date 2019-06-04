from django.conf.urls import url
from django.views.decorators.cache import cache_page

from construct import views


urlpatterns = [
    # url(r'^$', views.ConstructBrowser.as_view(), name='browse'), #no cache, for dev
    url(r'^$', cache_page(3600*24)(views.ConstructBrowser.as_view()), name='browse'),
    url(r'^experiments[/]?$', views.ExperimentBrowser.as_view(), name='browse'), #no cache version
    # url(r'^experiments[/]?$', cache_page(3600*24)(views.ExperimentBrowser.as_view()), name='browse'), #cache
    # url(r'^analysis[/]?$', views.ConstructStatistics.as_view(), name='statistics'),
    url(r'^analysis[/]?$', cache_page(3600*24)(views.ConstructStatistics.as_view()), name='statistics'),
    # url(r'^mutations[/]?$', views.ConstructMutations.as_view(), name='mutations'),
    url(r'^mutations[/]?$', cache_page(3600*24)(views.ConstructMutations.as_view()), name='mutations'),
    url(r'^residuetable[/]?$', views.ConstructTable.as_view(), name='residuetable'),
    url(r'^auto_webform/(?P<slug>[-\w]+)/$', views.fetch_pdb_for_webform, name='fetch'),
    url(r'^auto/(?P<slug>[-\w]+)/$', views.fetch_pdb, name='fetch'),
    url(r'^auto_all$', views.fetch_all_pdb, name='fetch'),
    url(r'^align$', views.align, name='align'),
    url(r'^design[/]?$', views.design.as_view(), name='design'),
    # url(r'^tool/$', views.tool, name='tool'),
    url(r'^tool/$', views.new_tool, name='tool'),
    url(r'^tool/json/nterm/(?P<slug>[-\w]+)/$', views.json_nterm, name='nterm'),
    url(r'^tool/json/cterm/(?P<slug>[-\w]+)/$', views.json_cterm, name='cterm'),
    url(r'^tool/json/icl3/(?P<slug>[-\w]+)/$', views.json_icl3, name='icl3'),
    url(r'^tool/json/icl2/(?P<slug>[-\w]+)/$', views.json_icl2, name='icl2'),
    url(r'^tool/json/glyco/(?P<slug>[-\w]+)/$', views.json_glyco, name='glyco'),
    url(r'^tool/json/palmi/(?P<slug>[-\w]+)/$', views.json_palmi, name='palmi'),
    url(r'^tool/json/mutations/(?P<slug>[-\w]+)/$', views.mutations, name='mutations'),
    url(r'^tool/json/fusion/(?P<slug>[-\w]+)/$', views.json_fusion, name='fusion'),
    url(r'^tool/json/termo/(?P<slug>[-\w]+)/$', views.thermostabilising, name='termo'),
    url(r'^tool/json/struc_rules/(?P<slug>[-\w]+)/$', views.structure_rules, name='struc_rules'),
    url(r'^tool/json/cons_strucs/(?P<slug>[-\w]+)/$', views.cons_strucs, name='cons_strucs'),
    url(r'^tool/json/cons_rf/(?P<slug>[-\w]+)/$', views.cons_rf, name='cons_rf'),
    url(r'^tool/json/cons_rm_GP/(?P<slug>[-\w]+)/$', views.cons_rm_GP, name='cons_rm_GP'),
    url(r'^tool/json/cons_rf_and_class/(?P<slug>[-\w]+)/$', views.cons_rf_and_class, name='cons_rf_and_class'),
    url(r'^stabilisation[/]?$', views.stabilisation_browser, name='stabilisation'),
    url(r'^(?P<slug>[-\w]+)/$', views.detail, name='detail'),
]
