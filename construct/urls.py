
from django.conf.urls import patterns, url

from construct import views


urlpatterns = patterns('',
    url(r'^(?P<slug>[-\w]+)/$', views.detail, name='detail'),

)

