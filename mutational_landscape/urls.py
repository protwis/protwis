from django.conf.urls import url

from mutational_landscape import views

urlpatterns = [
    url(r'^ajax/NaturalMutation/(?P<slug>[-\w]+)/$', views.ajaxNaturalMutation, name='ajaxNaturalMutation'),
]
