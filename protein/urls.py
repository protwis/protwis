from django.conf.urls import patterns, url

from protein import views


urlpatterns = patterns('',
    url(r'^$', views.IndexView.as_view(), name='index'),
    url(r'^(?P<slug>[-\w]+)/$', views.DetailView.as_view(), name='detail'),
)