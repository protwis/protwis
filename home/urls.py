from django.conf.urls import url
from django.views.generic import TemplateView
from django.conf import settings
from django.contrib.staticfiles.storage import staticfiles_storage
from django.views.generic.base import RedirectView
from home import views

urlpatterns = [
    url(r'^$', views.index, name='index'),
    # url(
    #     r'^favicon.ico$',
    #     RedirectView.as_view(
    #         url=staticfiles_storage.url('home/images/favicon_dev.ico'),
    #         permanent=False),
    #     name="favicon"
    # ),
    url(r'^citations', views.citations_json, name='citation'),
    url(r'^cite_gpcrdb', views.citeGPCRdb.as_view(), name='citation'),
    url(r'^cite_gproteindb', views.citeGproteinDb.as_view(), name='citation')
]
