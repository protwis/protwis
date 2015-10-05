from django.conf.urls import patterns, url

from mutation import views


urlpatterns = patterns('',
    url(r'^$', views.TargetSelection.as_view(), name='targetselection'),
    url(r'^import', views.importmutation, name='import'),
    url(r'^targetselection', views.TargetSelection.as_view(), name='targetselection'),
    url(r'^segmentselection', views.SegmentSelection.as_view(), name='segmentselection'),
    url(r'^render', views.render_mutations, name='render'),
    url(r'^(?P<download>download)', views.render_mutations, name='render'),
    url(r'^protein/(?P<protein>[^/]*?)/$', views.render_mutations, name='render'),
    url(r'^protein/(?P<protein>[^/]*?)/(?P<download>download)$', views.render_mutations, name='render'),
    url(r'^family/(?P<family>[^/]*?)/$', views.render_mutations, name='render'),
    url(r'^family/(?P<family>[^/]*?)/(?P<download>download)$', views.render_mutations, name='render'),
    url(r'^ajax/(?P<slug>[^/]*?)/$', views.ajax, name='ajax'),
    url(r'^ajax/(?P<slug>[^/]*?)/(?P<segments>.+)/$', views.ajaxSegments, name='ajaxSegments'),
)