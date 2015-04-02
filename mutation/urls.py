from django.conf.urls import patterns, url

from mutation import views


urlpatterns = patterns('',
    url(r'^$', views.index, name='index'),
    url(r'^import', views.importmutation, name='import'),
)