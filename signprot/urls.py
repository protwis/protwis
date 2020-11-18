from . import views
from django.urls import path, re_path
from django.views.decorators.cache import cache_page
from signprot.views import *
from contactnetwork.views import PdbTableData

urlpatterns = [
    path('', views.BrowseSelection.as_view(), name='index'),
    re_path(r'arrestin$', views.ArrestinSelection.as_view(), name='arrestin'),
    re_path(r'statistics/(?P<dataset>[^/]*?)/$',  views.GProtein, name='gprotein'),
    path('statistics_venn',  views.GProteinVenn, name='gprotein'),
    path('statistics_tree',  views.GProteinTree, name='gprotein'),
    path('statistics',  views.GProtein, name='gprotein'),
#    path('couplings',  views.couplings, name='couplings'),
    re_path(r'couplings$', (CouplingBrowser.as_view()), name='coupling_browser'),
    re_path(r'ginterface/(?P<protein>[^/]*?)/$', views.Ginterface, name='render'),
    re_path(r'ginterface[/]?$', views.TargetSelection.as_view(), name='targetselection'),
    re_path(r'ajax/barcode/(?P<slug>[-\w]+)/(?P<cutoff>\d+\.\d{0,2})/$', views.ajaxBarcode, name='ajaxBarcode'),
    re_path(r'ajax/interface/(?P<slug>[-\w]+)/$', views.ajaxInterface, name='ajaxInterface'),
    re_path(r'structure/(?P<pdbname>[-\w]+)/$', views.StructureInfo, name='StructureInfo'),
    re_path(r'family/(?P<slug>[-\w]+)/$', views.familyDetail, name='familyDetail'),
    re_path(r'matrix[/]?$', views.InteractionMatrix, name='InteractionMatrix'),
    path('matrix/seqsig/', views.IMSequenceSignature, name='SequenceSignature'),
    path('matrix/sigmat/', views.IMSignatureMatch, name='SignatureMatch'),
    re_path(r'matrix/render_sigmat/$', views.render_IMSigMat, name='renderSignatureMatch'),
    path('pdbtabledata', PdbTableData, name='pdbtreedata'),
    re_path(r'(?P<slug>[-\w]+)/$', views.signprotdetail, name='signprotdetail'),
]
