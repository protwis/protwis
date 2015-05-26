from django.conf.urls import patterns, url, include

from common import views


urlpatterns = patterns('',
    url(r'^addtoselection', views.AddToSelection, name='addtoselection'),
    url(r'^removefromselection', views.RemoveFromSelection, name='removefromselection'),
    url(r'^clearselection', views.ClearSelection, name='clearselection'),
    url(r'^togglefamilytreenode', views.ToggleFamilyTreeNode, name='togglefamilytreenode'),
    url(r'^selectionannotation', views.SelectionAnnotation, name='selectionannotation'),
    url(r'^selectionspeciespredefined', views.SelectionSpeciesPredefined, name='selectionspeciespredefined'),
    url(r'^selectionspeciestoggle', views.SelectionSpeciesToggle, name='selectionspeciestoggle'),
    url(r'^expandsegment', views.ExpandSegment, name='expandsegment'),
)