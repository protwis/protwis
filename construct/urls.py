from django.conf.urls import patterns, url

from construct import views


urlpatterns = [
    url(r'^$', views.ConstructBrowser.as_view(), name='browse'),
    url(r'^auto_webform/(?P<slug>[-\w]+)/$', views.fetch_pdb_for_webform, name='fetch'),
    url(r'^auto/(?P<slug>[-\w]+)/$', views.fetch_pdb, name='fetch'),
    url(r'^auto_all$', views.fetch_all_pdb, name='fetch'),
    url(r'^(?P<slug>[-\w]+)/$', views.detail, name='detail'),
]

