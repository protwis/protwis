from django.conf.urls import url
from django.urls import path
from django.views.decorators.cache import cache_page

from drugs import views

urlpatterns = [
    url(r'^indication/(?P<code>[\w|\W]+)/$',  views.indication_detail, name='indication'),
    ############################################################################################
    url(r'^drug_selection', views.DrugSectionSelection.as_view(title='Agent / drug search', page='Drugs'), name='DrugSearch'),
    url(r'^target_selection', views.DrugSectionSelection.as_view(title='Target search', page='Targets'), name='TargetSearch'),
    url(r'^indication_selection', views.DrugSectionSelection.as_view(title='Disease search', page='Indications'), name='IndicationSearch'),
    ############################################################################################
    url(r'^overview', cache_page(60*60*24*7)(views.DruggedGPCRome.as_view()), name='DruggedGPCRome'),
    url(r'^disease_overview', views.DiseaseOverview.as_view(), name='DiseaseOverview'),
    ############################################################################################
    url(r'^drugs_venn$', views.Venn, {'origin': 'drugs'}, name='DrugsVenn'),
    url(r'^targets_venn$', views.Venn, {'origin': 'targets'}, name='TargetVenn'),
    ############################################################################################
    url(r'^TargetSelectionTool', cache_page(60*60*24*7)(views.TargetSelectionTool.as_view()), name='TargetSelectionTool'),
    path('get_sankey_data/', views.fetch_sankey_data_view, name='get_sankey_data'),
    path('get_sankey_indi_data/', views.fetch_sankey_indi_data_view, name='get_sankey_indi_data'),
]