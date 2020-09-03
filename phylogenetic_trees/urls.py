from django.conf.urls import url
from phylogenetic_trees import views


urlpatterns = [
   # url(r'^referenceselection', views.ReferenceSelection.as_view(), name='referenceselection'),
    url(r'^targetselection', views.TargetSelection.as_view(), name='targetselection'),
    url(r'^segmentselection', views.SegmentSelection.as_view(), name='segmentselection'),
    url(r'^treesettings', views.TreeSettings.as_view(), name='treesettings'),
    url(r'^signatureselection', views.signature_selection, name='sig_selection'),
    url(r'^render_v3', views.render_tree_v3, name='render_v3'),
    # url(r'^render_v2', views.render_tree_v2, name='render_v2'),
    # url(r'^render', views.render_tree, name='render'),
    url(r'^showrings', views.modify_tree, name='render'),
    url(r'^get_buttons', views.get_buttons, name='render'),

]
