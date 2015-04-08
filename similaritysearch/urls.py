from django.conf.urls import patterns, url

from similaritysearch import views


urlpatterns = patterns('',
    url(r'^referenceselection', views.ReferenceSelection.as_view(), name='referenceselection'),
    url(r'^targetselection', views.TargetSelection.as_view(), name='targetselection'),
    url(r'^segmentselection', views.SegmentSelection.as_view(), name='segmentselection'),
    url(r'^render', views.render_alignment, name='render'),
)