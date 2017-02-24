from django.conf.urls import url

from mutational_landscape import views

urlpatterns = [
    url(r'^ajax/NaturalMutation/(?P<slug>[-\w]+)/$', views.ajaxNaturalMutation, name='ajaxNaturalMutation'),
    url(r'^ajax/CancerMutation/(?P<slug>[-\w]+)/$', views.ajaxCancerMutation, name='ajaxCancerMutation'),
    url(r'^ajax/DiseaseMutation/(?P<slug>[-\w]+)/$', views.ajaxDiseaseMutation, name='ajaxDiseaseMutation'),
    url(r'^ajax/mutant_extract', views.mutant_extract, name='mutant_extract')
]
