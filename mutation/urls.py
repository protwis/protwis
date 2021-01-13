from django.conf.urls import url
from django.urls import path
from django.views.decorators.cache import cache_page

from mutation import views


urlpatterns = [
    path('', cache_page(60*60*24*7)(views.TargetSelection.as_view()), name='targetselection'),
    path('import', views.importmutation, name='import'),
    path('designpdb', views.designPDB.as_view(), name='design'),
    path('design', cache_page(60*60*24*7)(views.design.as_view()), name='design'),


    path('state_stabilizing', views.designStateSelector.as_view(), name='design_state_selector'),
    path('state_stabilizing_<goal>', views.contactMutationDesign, name='design_state_mutations'),
    path('state_detail_gn', views.designStateDetailsGN, name='design_state_detail_gn'),

    path('gprot_coupling', views.designGprotSelector.as_view(), name='design_gprot_selector'),
    path('gprot_coupling_<goal>', views.gprotMutationDesign, name='design_gprot_mutations'),

    path('pocket', views.pocket, name='pocket'),
    path('statistics', views.coverage, name='statistics'),
    path('coverage', views.coverage, name='coverage'),
    path('calculatepdb', views.showcalculationPDB, name='showcalculationPDB'),
    path('calculate', views.showcalculation, name='showcalculation'),
    path('targetselection', cache_page(60*60*24*7)(views.TargetSelection.as_view()), name='targetselection'),
    path('segmentselection', views.SegmentSelection.as_view(), name='segmentselection'),
    path('render', views.render_mutations, name='render'),

    url(r'^(?P<download>download)', views.render_mutations, name='render'),
    url(r'^protein/(?P<protein>[^/]*?)/$', views.render_mutations, name='render'),
    url(r'^protein/(?P<protein>[^/]*?)/(?P<download>download)$', views.render_mutations, name='render'),
    url(r'^list/(?P<receptor_class>[^/]*?)/(?P<gn>[^/]*?)/(?P<aa>[^/]*?)$', views.render_mutations, name='render'),
    url(r'^family/(?P<family>[^/]*?)/$', views.render_mutations, name='render'),
    url(r'^family/(?P<family>[^/]*?)/(?P<download>download)$', views.render_mutations, name='render'),
    url(r'^ajax/(?P<slug>[^/]*?)/$', views.ajax, name='ajax'),
    url(r'^ajax/(?P<slug>[^/]*?)/(?P<segments>.+)/$', views.ajaxSegments, name='ajaxSegments'),
]
