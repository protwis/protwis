from django.conf.urls import patterns, url

from protein import views


urlpatterns = [
    url(r'^$', views.BrowseSelection.as_view(), name='index'),
    url(r'^autocomplete', views.SelectionAutocomplete, name='autocomplete'),
    url(r'^(?P<slug>[-\w]+)/$', views.detail, name='detail'),
]