from django.conf.urls import url

from interaction import views
# from django.views.generic import TemplateView

urlpatterns = [
    url(r'^$', views.InteractionSelection.as_view(), name='selection'),
    url(r'^list', views.list_structures, name='list'),
    # url(r'^calculate', views.showcalculation, name='showcalculation'),
    url(r'^sitesearch_calculate', views.calculate, {'redirect': '/sitesearch/segmentselectionpdb'},
        name='sitesearch_calculate'),
    url(r'^download', views.download, name='download'),
    url(r'^pdbfragment', views.pdbfragment, name='pdbfragment'),
    url(r'^pdb', views.pdb, name='pdb'),
    url(r'^complex/complexpdb', views.pdb, name='complexpdb'),
    url(r'^view', views.view, name='view'),
    url(r'^crystal', views.crystal, name='crystal'),
    url(r'^ligand', views.ligand, name='ligand'),
    url(r'^fragment', views.fragment, name='fragment'),
    url(r'^excel/(?P<session>session)/(?P<slug>[-\w]+)/$', views.excel, name='excel'),
    url(r'^excel/(?P<slug>[-\w]+)/$', views.excel, name='excel'),
    url(r'^ajax/(?P<slug>[-\w]+)/$', views.ajax, name='ajax'),
    url(r'^ajaxLigand/(?P<slug>[-\w]+)/(?P<ligand>.+)$', views.ajaxLigand, name='ajax'),
    url(r'^(?P<pdbname>\w+)$', views.StructureDetails, name='structure_details'),
    # url(r'^complex/(?P<pdbname>\w+)$', views.ComplexDetails, name='structure_details'),
]
