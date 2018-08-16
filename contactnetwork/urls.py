from django.conf.urls import url

from contactnetwork import views
# from django.views.generic import TemplateView

urlpatterns = [
    url(r'^interactions', views.Interactions, name='interactions'),
    url(r'^interactiondata', views.InteractionData, name='interactiondata'),
    url(r'^pdbtreedata', views.PdbTreeData, name='pdbtreedata'),
    url(r'^pdbtabledata', views.PdbTableData, name='pdbtreedata'),
    url(r'^pdb/(?P<pdbname>\w+)$', views.ServePDB, name='serve_pdb'),

]