
from django.conf.urls import patterns, url

from construct import views


urlpatterns = patterns('',
    url(r'^constructs', views.constructs, name='constructs'),

#  url(r'^targetselection', views.TargetSelection.as_view(), name='targetselection'),
)

