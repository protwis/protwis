from django.conf.urls import url

from signprot import views

urlpatterns = [
	url(r'^$', views.BrowseSelection.as_view(), name='index'),
    url(r'^(?P<slug>[-\w]+)/$', views.signprotdetail, name='signprotdetail'),
]
