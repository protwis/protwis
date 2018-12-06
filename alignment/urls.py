from django.conf.urls import url
from django.views.generic import TemplateView
from django.views.decorators.cache import cache_page

from alignment import views



urlpatterns = [
    url(r'^targetselection', (views.TargetSelection.as_view()), name='targetselection'),
    url(r'^gproteinselection', (views.TargetSelectionGprotein.as_view()), name='targetselectiongprot'),
    url(r'^arrestinselection', (views.TargetSelectionArrestin.as_view()), name='arrestinselectionprot'),
    url(r'^segmentselectiongprot', views.SegmentSelectionGprotein.as_view(), name='segmentselectiongprot'),
    url(r'^segmentselectionarrestin', views.SegmentSelectionArrestin.as_view(), name='segmentselectionarrestin'),
    url(r'^segmentselection', views.SegmentSelection.as_view(), name='segmentselection'),
    url(r'^render/(?P<slug>[^/]+)/$', views.render_family_alignment, name='render-family'),
    url(r'^render', views.render_alignment, name='render'),
    url(r'^fasta/(?P<slug>[^/]+)/$', views.render_fasta_family_alignment, name='fasta-family'),
    url(r'^fasta', views.render_fasta_alignment, name='fasta'),
    url(r'^csv', views.render_csv_alignment, name='csv'),
    url(r'blastsearch$', views.BlastSearchInput.as_view(), name='blastsearch'),
    url(r'blastsearchresults', views.BlastSearchResults.as_view(), name='blastsearch_results'),
]
