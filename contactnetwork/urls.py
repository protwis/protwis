from django.conf.urls import url

from contactnetwork import views
# from django.views.generic import TemplateView

urlpatterns = [
    url(r'^heatmap$', views.HeatMap, name='heatmap'),
    url(r'^heatmapdata$', views.HeatMapDataJson, name='heatmapjsondata'),
]