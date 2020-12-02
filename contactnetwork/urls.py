from django.conf.urls import url

from contactnetwork import views
# from django.views.generic import TemplateView

urlpatterns = [
    url(r'^clusteringdata$', views.ClusteringData, name='clusteringdata'),
    url(r'^clustering$', views.Clustering, name='clustering'),
    url(r'^structure_clustering$', views.Clustering, name='clustering'),
    url(r'^distances', views.ShowDistances, name='distances'),
    url(r'^distancedatagroups', views.DistanceDataGroups, name='distancedatagroups'),
    url(r'^distancedata', views.DistanceData, name='distancedata'),
    url(r'^interactions[/]?$', views.Interactions, name='interactions'),
    url(r'^comparative_analysis[/]?$', views.Interactions, name='interactions'),
    url(r'^interactiondata', views.InteractionData, name='interactiondata'),
    url(r'^browser[/]?$', views.InteractionBrowser, name='interactionsbrowser'),
    url(r'^browserdata', views.InteractionBrowserData, name='interactionsbrowserdata'),
    url(r'^state_contacts[/]?$', views.StateContacts, name='statecontacts'),
    url(r'^pdbtreedata', views.PdbTreeData, name='pdbtreedata'),
    url(r'^pdbtabledata', views.PdbTableData, name='pdbtabledata'),
    url(r'^pdb/(?P<pdbname>\w+)$', views.ServePDB, name='serve_pdb'),

]
