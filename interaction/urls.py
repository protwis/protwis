from django.conf.urls import patterns, url

from interaction import views


urlpatterns = patterns('',
    url(r'^$', views.index, name='index'),
    url(r'^list', views.list, name='list'),
    url(r'^calculate', views.calculate, name='calculate'),
    url(r'^download', views.download, name='download'),
    url(r'^pdb', views.pdb, name='pdb'),
    url(r'^view', views.view, name='view'),
    url(r'^ligand', views.ligand, name='ligand'),
)