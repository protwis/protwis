from django.conf.urls import url
from django.views.decorators.cache import cache_page

from construct import views


urlpatterns = [
    url(r'^$', views.ConstructBrowser.as_view(), name='browse'), #no cache, for dev
    #url(r'^$', cache_page(60*60*24*7)(views.ConstructBrowser.as_view()), name='browse'),
    url(r'^statistics[/]?$', views.ConstructStatistics.as_view(), name='statistics'),
    url(r'^mutations[/]?$', views.ConstructMutations.as_view(), name='mutations'),
    url(r'^residuetable[/]?$', views.ConstructTable.as_view(), name='residuetable'),
    url(r'^auto_webform/(?P<slug>[-\w]+)/$', views.fetch_pdb_for_webform, name='fetch'),
    url(r'^auto/(?P<slug>[-\w]+)/$', views.fetch_pdb, name='fetch'),
    url(r'^auto_all$', views.fetch_all_pdb, name='fetch'),
    url(r'^align$', views.align, name='align'),
    url(r'^design[/]?$', views.design.as_view(), name='design'),
    url(r'^tool/$', views.tool, name='tool'),
    url(r'^tool/json/nterm/(?P<slug>[-\w]+)/$', views.json_nterm, name='nterm'),
    url(r'^tool/json/cterm/(?P<slug>[-\w]+)/$', views.json_cterm, name='cterm'),
    url(r'^tool/json/icl3/(?P<slug>[-\w]+)/$', views.json_icl3, name='icl3'),
    url(r'^tool/json/glyco/(?P<slug>[-\w]+)/$', views.json_glyco, name='glyco'),
    url(r'^tool/json/palmi/(?P<slug>[-\w]+)/$', views.json_palmi, name='palmi'),
    url(r'^tool/json/mutations/(?P<slug>[-\w]+)/$', views.mutations, name='mutations'),
    url(r'^tool/json/fusion/(?P<slug>[-\w]+)/$', views.json_fusion, name='fusion'),
    url(r'^tool/json/termo/(?P<slug>[-\w]+)/$', views.thermostabilising, name='termo'),
    url(r'^(?P<slug>[-\w]+)/$', views.detail, name='detail'),
]
