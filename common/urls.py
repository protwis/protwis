from django.conf.urls import patterns, url, include

from common import views


urlpatterns = patterns('',
    url(r'^addtoselection', views.AddToSelection, name='addtoselection'),
    url(r'^removefromselection', views.RemoveFromSelection, name='removefromselection'),
    url(r'^clearselection', views.ClearSelection, name='clearselection'),
    url(r'^togglefamilytreenode', views.ToggleFamilyTreeNode, name='togglefamilytreenode'),
    url(r'^autocomplete/', include('autocomplete_light.urls')),
)