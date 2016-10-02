from django.conf.urls import url
from django.views.decorators.cache import cache_page

from mutation import views


urlpatterns = [
    url(r'^$', cache_page(60*60*24*7)(views.TargetSelection.as_view()), name='targetselection'),
    url(r'^import', views.importmutation, name='import'),
    url(r'^designpdb', cache_page(60*60*24*7)(views.designPDB.as_view()), name='design'),
    url(r'^design', cache_page(60*60*24*7)(views.design.as_view()), name='design'),
    url(r'^pocket', views.pocket, name='pocket'),
    url(r'^statistics', views.coverage, name='statistics'),
    url(r'^coverage', views.coverage, name='coverage'),
    url(r'^calculatepdb', views.showcalculationPDB, name='showcalculationPDB'),
    url(r'^calculate', views.showcalculation, name='showcalculation'),
    url(r'^targetselection', cache_page(60*60*24*7)(views.TargetSelection.as_view()), name='targetselection'),
    url(r'^segmentselection', views.SegmentSelection.as_view(), name='segmentselection'),
    url(r'^render', views.render_mutations, name='render'),
    url(r'^(?P<download>download)', views.render_mutations, name='render'),
    url(r'^protein/(?P<protein>[^/]*?)/$', views.render_mutations, name='render'),
    url(r'^protein/(?P<protein>[^/]*?)/(?P<download>download)$', views.render_mutations, name='render'),
    url(r'^list/(?P<receptor_class>[^/]*?)/(?P<gn>[^/]*?)/(?P<aa>[^/]*?)$', views.render_mutations, name='render'),
    url(r'^family/(?P<family>[^/]*?)/$', views.render_mutations, name='render'),
    url(r'^family/(?P<family>[^/]*?)/(?P<download>download)$', views.render_mutations, name='render'),
    url(r'^ajax/(?P<slug>[^/]*?)/$', views.ajax, name='ajax'),
    url(r'^ajax/(?P<slug>[^/]*?)/(?P<segments>.+)/$', views.ajaxSegments, name='ajaxSegments'),
]
