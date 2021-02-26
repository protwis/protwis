from . import views
from django.urls import path
from django.views.decorators.cache import cache_page
from signprot.views import CouplingBrowser, CouplingBrowser2
from contactnetwork.views import PdbTableData

urlpatterns = [
    path('', views.BrowseSelection.as_view(), name='index'),
    path('arrestin', views.ArrestinSelection.as_view(), name='arrestin'),
    path('statistics/<dataset>/',  views.GProtein, name='gprotein'),
    path('statistics_venn',  views.GProteinVenn, name='gprotein'),
    path('statistics_tree',  views.GProteinTree, name='gprotein'),
    #path('statistics',  views.GProtein, name='gprotein'),
    path('statistics',  views.CouplingProfiles, name='coupling_profiles'),
    path('couplings', views.couplings, name='coupling_browser'),
    path('couplings1', cache_page(60*60*24*7)(CouplingBrowser.as_view()), name='coupling_browser1'),
#    path('couplings2', (CouplingBrowser2.as_view()), name='coupling_browser2'),
    path('couplings2', cache_page(60*60*24*7)(CouplingBrowser2.as_view()), name='coupling_browser2'),
    path('ginterface/<protein>/', views.Ginterface, name='render'),
    path('ginterface/', views.TargetSelection.as_view(), name='targetselection'),
    path('ajax/barcode/<slug>/<cutoff>/', views.ajaxBarcode, name='ajaxBarcode'),
    path('ajax/interface/<slug>/', views.ajaxInterface, name='ajaxInterface'),
    path('structure/<pdbname>/', views.StructureInfo, name='StructureInfo'),
    path('family/<slug>/', views.familyDetail, name='familyDetail'),
    path('matrix/', views.InteractionMatrix, name='InteractionMatrix'),
    path('matrix/seqsig/', views.IMSequenceSignature, name='SequenceSignature'),
    path('matrix/sigmat/', views.IMSignatureMatch, name='SignatureMatch'),
    path('matrix/render_sigmat/', views.render_IMSigMat, name='renderSignatureMatch'),
    path('pdbtabledata', PdbTableData, name='pdbtreedata'),
    path('<slug>/', views.signprotdetail, name='signprotdetail'),
]
