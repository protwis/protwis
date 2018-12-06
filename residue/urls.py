from django.conf.urls import url

from residue import views


urlpatterns = [
    url(r'^targetselection', views.TargetSelection.as_view(), name='targetselection'),
    url(r'^residuetable$', views.ResidueTablesSelection.as_view(), name='residuetable'),
    url(r'^residuetabledisplay', views.ResidueTablesDisplay.as_view(), name='residuetable'),
    url(r'^residuefunctionbrowser$', views.ResidueFunctionBrowser.as_view(), name='residue_function_browser'),
]
