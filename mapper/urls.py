from django.conf.urls import url
from django.views.decorators.cache import cache_page
from . import views

urlpatterns = [
    url(r'^MapperLandingPage/$', views.MapperLandingPageView.as_view(), name='MapperLandingPage'),
    url(r'^MapperGPCRomeWheel/$', views.MapperGPCRomeWheelView.as_view(), name='MapperGPCRomeWheel'),
    url(r'^MapperTree/$', views.MapperTreeView.as_view(), name='MapperTree'),
    url(r'^MapperHeatmap/$', views.MapperHeatmapView.as_view(), name='MapperHeatmap'),
    url(r'^MapperList/$', views.MapperListView.as_view(), name='MapperList'),
    url(r'^MapperCluster/$', views.MapperClusterView.as_view(), name='MapperCluster'),
    url(r'^Cluster', views.ClusterRender.as_view(), name='DataMapperCluster'),
]
