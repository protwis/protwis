from django.conf.urls import url
from django.views.decorators.cache import cache_page
from . import views

urlpatterns = [
    url(r'^DataMapperHome(?P<page>\w+)/$', views.DataMapperHome.as_view(), name='DataMapperHome'),
    url(r'^GPCRome', views.GPCRomeRender.as_view(), name='DataMapperGPCRome'),
    url(r'^Tree', views.TreeRender.as_view(), name='DataMapperTree'),
    url(r'^Cluster', views.ClusterRender.as_view(), name='DataMapperCluster'),
    url(r'^List', views.ListRender.as_view(), name='DataMapperList'),
    url(r'^Heatmap', views.HeatmapRender.as_view(), name='DataMapperHeatmap')
]
