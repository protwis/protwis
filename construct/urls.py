from django.conf.urls import url
from django.views.decorators.cache import cache_page

from construct import views


urlpatterns = [
    url(r'^$', cache_page(60*60*24*7)(views.ConstructBrowser.as_view()), name='browse'),
    # url(r'^$', views.ConstructBrowser.as_view(), name='browse'), #no cache, for dev
    url(r'^auto_webform/(?P<slug>[-\w]+)/$', views.fetch_pdb_for_webform, name='fetch'),
    url(r'^auto/(?P<slug>[-\w]+)/$', views.fetch_pdb, name='fetch'),
    url(r'^auto_all$', views.fetch_all_pdb, name='fetch'),
    url(r'^(?P<slug>[-\w]+)/$', views.detail, name='detail'),
    url(r'^align$', views.align, name='align'),
]
