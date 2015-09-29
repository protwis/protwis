from django.conf.urls import patterns, url

from pages import views

urlpatterns = patterns('',
    url(r'^$', views.index, name='index'),
	url(r'^releasenotes', views.releasenotes, name='releasenotes'),
    url(r'^contribute', views.contribute, name='contribute'),
    url(r'^contact', views.contact, name='contact'),
    url(r'^citing', views.citing, name='citing'),
    url(r'^meetings', views.meetings, name='meetings'),
    url(r'^servers', views.servers, name='servers'),
    url(r'^about', views.about, name='about'),
    url(r'^acknowledgements', views.acknowledgements, name='acknowledgements'),
    url(r'^legalnotice', views.legalnotice, name='legalnotice'),
    url(r'^linking', views.linking, name='linking'),
)