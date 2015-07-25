from django.conf.urls import patterns, url

from residue import views


urlpatterns = patterns('',
    url(r'^targetselection', views.TargetSelection.as_view(), name='targetselection'),
    url(r'^residuetable', views.ResidueTablesSelection.as_view(), name='residuetable'),
)