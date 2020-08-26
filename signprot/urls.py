from contactnetwork.views import PdbTableData
from django.conf.urls import url
from django.views.decorators.cache import cache_page
from signprot import views
from signprot.views import *

urlpatterns = [
    url(r'^$', views.BrowseSelection.as_view(), name='index'),
    url(r'^arrestin$', views.ArrestinSelection.as_view(), name='arrestin'),
    url(r'^statistics/(?P<dataset>[^/]*?)/$',  views.GProtein, name='gprotein'),
    url(r'^statistics',  views.GProtein, name='gprotein'),
#    url(r'^couplings',  views.couplings, name='couplings'),
    url(r'^couplings$', (CouplingBrowser.as_view()), name='coupling_browser'),
    url(r'^couplingsbrowser$', (CouplingBrowser.as_view()), name='coupling_browser'),
    url(r'^ginterface/(?P<protein>[^/]*?)/$', views.Ginterface, name='render'),
    url(r'^ginterface[/]?$', views.TargetSelection.as_view(), name='targetselection'),
    url(r'^ajax/barcode/(?P<slug>[-\w]+)/(?P<cutoff>\d+\.\d{0,2})/$', views.ajaxBarcode, name='ajaxBarcode'),
    url(r'^ajax/interface/(?P<slug>[-\w]+)/$', views.ajaxInterface, name='ajaxInterface'),
    url(r'^structure/(?P<pdbname>[-\w]+)/$', views.StructureInfo, name='StructureInfo'),
    url(r'^family/(?P<slug>[-\w]+)/$', views.familyDetail, name='familyDetail'),
    url(r'^matrix[/]?$', views.InteractionMatrix, name='InteractionMatrix'),
    url(r'^matrix/seqsig/', views.IMSequenceSignature, name='SequenceSignature'),
    url(r'^matrix/sigmat/', views.IMSignatureMatch, name='SignatureMatch'),
    url(r'^matrix/render_sigmat/$', views.render_IMSigMat, name='renderSignatureMatch'),
    url(r'^pdbtabledata', PdbTableData, name='pdbtreedata'),
    url(r'^(?P<slug>[-\w]+)/$', views.signprotdetail, name='signprotdetail'),
]
