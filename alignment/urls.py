from django.conf.urls import patterns, url

from alignment import views


urlpatterns = patterns('',
    url(r'^targetselection', views.TargetSelection.as_view(), name='targetselection'),
    url(r'^segmentselection', views.SegmentSelection.as_view(), name='segmentselection'),
    url(r'^render', views.render_alignment, name='render'),
    url(r'^fasta', views.render_fasta_alignment, name='fasta'),
)