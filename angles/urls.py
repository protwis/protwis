from django.conf.urls import url
from angles import views
from django.conf import settings
from django.views.decorators.cache import cache_page

urlpatterns = [
    url(r'^ttest$', views.testTemplate.as_view(), name='xxx'),
    url(r'^browsertest$', views.browsertest, name='browsertest'),
    url(r'^angledat$', views.get_angles, name='anglejson'),
    url(r'^pdb/(?P<pdbname>\w+)$', views.ServePDB, name='serve_pdb'),
    url(r'^pdbtabledata$', views.PdbTableData, name='pdbtreedata'),
    url(r'^pdbtabledata2$', views.PdbTableData2, name='pdbtreedata2'),
]
