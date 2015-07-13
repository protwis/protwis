
from django.conf.urls import patterns, url

from structure import views


urlpatterns = patterns('',
    url(r'^constructs', views.constructs, name='constructs'),
)

