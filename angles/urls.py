from django.conf.urls import url

from angles import views

urlpatterns = [
    url(r'^angleanalysis$', views.angleAnalysis, name='angleanalysis'),
    url(r'^angledat$', views.get_angles, name='anglejson')
]
