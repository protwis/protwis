from django.conf.urls import patterns, url

from residue import views


urlpatterns = patterns('',
    url(r'^targetselection', views.TargetSelection.as_view(), name='targetselection'),
)