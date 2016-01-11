from django.conf.urls import patterns, url

from pages import views

urlpatterns = [
	url(r'^releasenotes', views.releasenotes, name='releasenotes'), 
]