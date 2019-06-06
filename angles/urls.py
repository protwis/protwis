from django.conf.urls import url

from angles import views

urlpatterns = [
    url(r'^angleanalysis$', views.angleAnalyses, name='angleanalysis'),
    url(r'^angleanalyses$', views.angleAnalyses, name='angleanalyses'),
    url(r'^structurecheck$', views.structureCheck, name='structurecheck'),
    url(r'^angledata$', views.get_angles, name='anglejson')
]
