from django.conf.urls import patterns, url

from documentation import views

urlpatterns = patterns('',
    url(r'^$', views.index, name='index'),
    url(r'^crystals', views.crystals, name='crystals'),
    url(r'^diagrams', views.diagrams, name='diagrams'),
    url(r'^numbering', views.numbering, name='numbering'),
    url(r'^pharmacophore', views.pharmacophore, name='update'),
	url(r'^sequences', views.sequences, name='sequences'),
    url(r'^similarities', views.similarities, name='similarities'),
)