from django.conf.urls import url

from drugs import views

urlpatterns = [
    url(r'^drugbrowser',  views.drugbrowser, name='drugbrowser'),
    url(r'^drugstatistics',  views.drugstatistics, name='drugstatistics'),
    url(r'^drugmapping',  views.drugmapping, name='drugmapping'),
]