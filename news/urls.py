from django.conf.urls import url

from news import views

urlpatterns = [
    url(r'^$', views.index, name='index'),
]