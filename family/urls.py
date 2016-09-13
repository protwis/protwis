from django.conf.urls import url

from family import views


urlpatterns = [
    url(r'^(?P<slug>[-\w]+)/$', views.detail, name='detail'),
]