from . import views
from django.urls import path
from django.views.decorators.cache import cache_page
from django.views.generic.base import RedirectView
from signprot.views import CouplingBrowser, CouplingBrowser_deprecated
from contactnetwork.views import PdbTableData

urlpatterns = [
    path('', RedirectView.as_view(url='gprotein', permanent=False), name='index'),
    path('gprotein', views.BrowseSelection.as_view(), name='index'),
    path('arrestin', views.ArrestinSelection.as_view(), name='arrestin'),
    # path('arrestincouplings', cache_page(60*60*24*7)(ArrestinCoupling.as_view()), name='arrestin_coupling'),
    path('arrestincouplings', cache_page(60*60*24*7)(CouplingBrowser_deprecated.as_view(subunit_filter = "200_000_001", families = ["Beta"], page='arrestin')), name='arrestin_coupling'),
    # path('statistics/<dataset>/',  views.GProtein, name='gprotein'),
    path('statistics_venn',  views.GProteinVenn, name='gprotein'),
    path('statistics_tree',  views.GProteinTree, name='gprotein'),
    path('arrestin_venn',  views.ArrestinVenn, name='gprotein'),
    path('arrestin_tree',  views.ArrestinTree, name='gprotein'),
    #path('statistics',  views.GProtein, name='gprotein'),
    path('statistics',  views.CouplingProfiles, name='coupling_profiles'),
    path('coupling_datasets', views.CouplingDatasets, name='coupling_datasets'),
    path('coupling_biosensors', views.CouplingBiosensors, name='coupling_biosensors'),
    path('couplings', cache_page(60*60*24*7)(CouplingBrowser.as_view()), name='coupling_browser'),
    path('ginterface/<protein>/', views.Ginterface, name='render'),
    path('ginterface/', views.TargetSelection.as_view(), name='targetselection'),
    path('ajax/barcode/<slug>/<cutoff>/', views.ajaxBarcode, name='ajaxBarcode'),
    path('ajax/interface/<slug>/', views.ajaxInterface, name='ajaxInterface'),
    path('structure/<pdbname>/', views.StructureInfo, name='StructureInfo'),
    path('family/<slug>/', views.familyDetail, name='familyDetail'),
    #Have two set of URLs with two prefixes: GP, Arr
    path('matrix/', views.GProteinInteractionMatrix, name='InteractionMatrix'),
    path('arr_matrix/', views.ArrestinInteractionMatrix, name='InteractionMatrix'),
    path('matrix/seqsig/', views.IMSequenceSignature, name='SequenceSignature'),
    path('matrix/sigmat/', views.IMSignatureMatch, name='SignatureMatch'),
    path('matrix/render_sigmat/', views.render_IMSigMat, name='renderSignatureMatch'),
    path('pdbtabledata', PdbTableData, name='pdbtreedata'),
    path('matrix/AJAX_Interactions/', views.AJAX_Interactions, name='ajaxInteractions'),
    path('<slug>/', views.signprotdetail, name='signprotdetail'),
]
