from django.conf.urls import patterns, url

from common import views


urlpatterns = patterns('',
    url(r'^addtoselection', views.AddToSelection, name='addtoselection'),
    url(r'^clearselection', views.ClearSelection, name='clearselection'),
)