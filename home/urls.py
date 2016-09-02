from django.conf.urls import url
#from home import views
from django.views.generic import TemplateView
from django.conf import settings

from home import views

urlpatterns = [
    url(r'^$', views.index, name='index'), 
]