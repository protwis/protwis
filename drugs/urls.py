from django.conf.urls import url

from drugs import views

urlpatterns = [
    url(r'^drugbrowser',  views.drugbrowser, name='drugbrowser'),
    url(r'^drugstatistics',  views.drugstatistics, name='drugstatistics'),
    url(r'^drugmapping',  views.drugmapping, name='drugmapping'),
    url(r'^nhs/section/(?P<slug>[\w|\W]+)/$',  views.nhs_section, name='nhs_section'),
    url(r'^nhs/(?P<slug>[\w|\W]+)/$',  views.nhs_drug, name='nhs_drug'),
    url(r'^NewDrugsBrowser', views.NewDrugsBrowser.as_view(), name='NewDrugsBrowser'),
]
