from django.conf.urls import patterns, url
from django.views.generic import TemplateView

from alignment import views


urlpatterns = patterns('',
    url(r'^targetselection', views.TargetSelection.as_view(), name='targetselection'),
    url(r'^segmentselection', views.SegmentSelection.as_view(), name='segmentselection'),
    url(r'^render/(?P<slug>[^/]+)/$', views.render_family_alignment, name='render-family'),
    url(r'^render', views.render_alignment, name='render'),
    url(r'^fasta/(?P<slug>[^/]+)/$', views.render_fasta_family_alignment, name='fasta-family'),
    url(r'^fasta', views.render_fasta_alignment, name='fasta'),
    url(r'^csv', views.render_csv_alignment, name='csv'),
    url(r'blastsearch$', TemplateView.as_view(template_name='sequence/blast_search.html'), name='blastsearch'),
    url(r'blastsearchresults', views.BlastSearchResults.as_view(), name='blastsearch'),
)