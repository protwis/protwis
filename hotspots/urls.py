from django.conf.urls import url

from hotspots import views

urlpatterns = [
    url(r'^hotspotsview$', views.hotspotsView, name='hotspotsView'),
    url(r'^hotspotsdata$', views.getHotspots, name='hotspotsJSON')
]
