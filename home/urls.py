from django.conf.urls import url
from django.views.generic import TemplateView
from django.conf import settings
from django.contrib.staticfiles.storage import staticfiles_storage
from django.views.generic.base import RedirectView
from home import views

from common.definitions import re_string_cite_us_page



urlpatterns = [
    url(r'^$', views.index, name='index'),
    # url(
    #     r'^favicon.ico$',
    #     RedirectView.as_view(
    #         url=staticfiles_storage.url('home/images/favicon_dev.ico'),
    #         permanent=False),
    #     name="favicon"
    # ),
    url(r'^citations/((?P<output_type>[-\w]+)/)?$', views.citations_json, name='citation'),
    url(r'^citation_by_url/((?P<output_type>[-\w]+)/)?$', views.citation_json_by_url, name='citation_by_url'),
    url(re_string_cite_us_page, views.cite_us, name='citation'),
]
