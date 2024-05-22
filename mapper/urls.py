from django.conf.urls import url
from django.urls import path
from django.views.decorators.cache import cache_page

from . import views


urlpatterns = [
    url(r'^$', views.LandingPage.as_view(), name='landing_page')
]