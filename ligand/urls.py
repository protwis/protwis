from django.conf.urls import url
from django.views.decorators.cache import cache_page
from django.views.generic import TemplateView
from ligand.views import *

urlpatterns = [
    url(r'^$', LigandBrowser, name='ligand_browser'),
    url(r'^p/(?P<slug>[-\w]+)/$',p_detail, name='p_detail'),
    url(r'^l/(?P<ligand__id>[-\w]+)/$',l_detail, name='l_detail'),
]
