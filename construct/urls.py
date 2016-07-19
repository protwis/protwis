from django.conf.urls import patterns, url

from construct import views


urlpatterns = [
    url(r'^$', views.ConstructBrowser.as_view(), name='browse'),
    url(r'^(?P<slug>[-\w]+)/$', views.detail, name='detail'),
]

