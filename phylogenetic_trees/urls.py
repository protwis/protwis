from django.conf.urls import url
from phylogenetic_trees import views


urlpatterns = [
   # url(r'^referenceselection', views.ReferenceSelection.as_view(), name='referenceselection'),
    url(r'^targetselection', views.TargetSelection.as_view(), name='targetselection'),
    url(r'^segmentselection', views.SegmentSelection.as_view(), name='segmentselection'),
    url(r'^treesettings', views.TreeSettings.as_view(), name='treesettings'),
    url(r'^render_new', views.render_tree_new, name='render_new'),
    url(r'^render', views.render_tree, name='render'),
    url(r'^showrings', views.modify_tree, name='render'),
    url(r'^get_buttons', views.get_buttons, name='render'),

]
