from django.conf.urls import url
from angles import views
from django.conf import settings
from django.views.decorators.cache import cache_page

urlpatterns = [
    url(r'^angleanalysis$', views.angleanalysis, name='angleanalysis'),
    url(r'^angledat$', views.get_angles, name='anglejson'),
    url(r'^pdbtabledata$', views.PdbTableData, name='pdbtreedata'),
]
