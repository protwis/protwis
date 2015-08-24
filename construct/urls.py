
from django.conf.urls import patterns, url

from construct import views


urlpatterns = patterns('',
   # url(r'^constructs', views.constructs, name='constructs'),
    url(r'^(?P<slug>[-\w]+)/$', views.detail, name='detail'),

#  url(r'^targetselection', views.TargetSelection.as_view(), name='targetselection'),
)

