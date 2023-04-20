from django.conf.urls import url
from django.views.decorators.cache import cache_page

from drugs import views

urlpatterns = [
    url(r'^drugbrowser', cache_page(60 * 60 * 24 * 28)(views.DrugBrowser.as_view()), name='drugbrowser'),
    url(r'^drugstatistics',  views.drugstatistics, name='drugstatistics'),
    url(r'^drugmapping',  views.drugmapping, name='drugmapping'),
    url(r'^nhs/section/(?P<slug>[\w|\W]+)/$',  views.nhs_section, name='nhs_section'),
    url(r'^nhs/(?P<slug>[\w|\W]+)/$',  views.nhs_drug, name='nhs_drug'),
]
