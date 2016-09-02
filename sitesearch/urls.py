from django.conf.urls import url

from sitesearch import views


urlpatterns = [
    url(r'^targetselectionpdb', views.TargetSelectionPdb.as_view(), name='targetselectionpdb'),
    url(r'^targetselection', views.TargetSelection.as_view(), name='targetselection'),
    url(r'^segmentselectionpdb', views.SegmentSelectionPdb.as_view(), name='segmentselectionpdb'),
    url(r'^segmentselection', views.SegmentSelection.as_view(), name='segmentselection'),
    url(r'^structureupload', views.StructureUpload.as_view(), name='structureupload'),
    url(r'^sitesave', views.site_download, name='sitesave'),
    url(r'^siteupload', views.site_upload, name='siteupload'),
    url(r'^render', views.render_alignment, name='render'),
    url(r'^fasta', views.render_fasta_alignment, name='fasta'),
    url(r'^csv', views.render_csv_alignment, name='csv'),
]