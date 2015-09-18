from django.conf.urls import patterns, url

from pages import views

urlpatterns = patterns('',
    url(r'^$', views.index, name='index'),
    url(r'^committees', views.committees, name='committees'),
    url(r'^contributors', views.contributors, name='contributors'),
    url(r'^contribute', views.contribute, name='contribute'),
    url(r'^contact', views.contact, name='contact'),
	url(r'^citing', views.citing, name='citing'),
	url(r'^poster', views.poster, name='poster'),
	url(r'^meetings', views.meetings, name='meetings'),
	url(r'^servers', views.servers, name='servers'),
	url(r'^update', views.update, name='update'),
)