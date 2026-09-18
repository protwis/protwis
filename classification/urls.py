from django.conf.urls import url
from django.views.generic import RedirectView

from classification import views


urlpatterns = [
    # Redirect /classification/ to the main Classification page
    url(
        r'^$',
        RedirectView.as_view(
            url='/classification/Classification',
            permanent=False,
        ),
        name='classification-index',
    ),

    # Canonical URLs
    # GPCRBrowser was merged into Classification as its first tab ("GPCR list");
    # keep the old address working as a redirect for existing bookmarks/links.
    url(
        r'^GPCRBrowser[/]?$',
        RedirectView.as_view(
            url='/classification/Classification',
            permanent=False,
        ),
        name='classification-gpcrbrowser',
    ),
    url(r'^overview[/]?$', views.Classification.as_view(), name='classification-classification'),
    url(r'^visualizations[/]?$', views.ClassificationVisualizationsLanding.as_view(), name='classification-visualizations'),
    url(
        r'^visualizations/class/(?P<class_key>[A-Za-z0-9]+)[/]?$',
        views.ClassificationVisualizationDetail.as_view(),
        name='classification-visualizations-class',
    ),
    url(
        r'^visualizations/tree[/]?$',
        views.ClassificationTreeVisualizationDetail.as_view(),
        name='classification-visualizations-tree',
    ),
    url(
        r'^visualizations/superfamily[/]?$',
        views.GPCRSuperfamilyVisualizationDetail.as_view(),
        name='classification-visualizations-superfamily',
    ),
    url(
        r'^visualizations/family/(?P<family_key>[-a-z0-9]+)[/]?$',
        views.ReceptorFamilyVisualizationDetail.as_view(),
        name='classification-visualizations-family',
    ),
    url(r'^StructureSim[/]?$', views.StructureSim.as_view(), name='classification-structuresim'),
    # JSON-only API — consumed by the Cluster tab inlined into GPCRSuperfamilyVisualizationDetail.
    url(r'^NewClassClusterTree[/]?$', views.NewClassClusterTree.as_view(), name='classification-newclassclustertree'),
]

