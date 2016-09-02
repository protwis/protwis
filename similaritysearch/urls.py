from django.conf.urls import url

from similaritysearch import views


urlpatterns = [
    url(r'^referenceselection', views.ReferenceSelection.as_view(), name='referenceselection'),
    url(r'^targetselection', views.TargetSelection.as_view(), name='targetselection'),
    url(r'^segmentselection', views.SegmentSelection.as_view(), name='segmentselection'),
    url(r'^render', views.render_alignment, name='render'),
    url(r'^fasta', views.render_fasta_alignment, name='fasta'),
    url(r'^csv', views.render_csv_alignment, name='csv'),
]