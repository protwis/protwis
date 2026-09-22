from django.conf.urls import url
from django.views.generic import RedirectView

from classification import views


urlpatterns = [
    # Redirect /classification/ to the main Classification page
    url(
        r'^$',
        RedirectView.as_view(
            url='/classification/levels_and_terms',
            permanent=False,
        ),
        name='classification-index',
    ),

    # Canonical URLs
    # GPCRBrowser used to be its own page, then became the "GPCR list" tab on this
    # page, and now lives at /allgpcrs; keep the old address working as a redirect
    # for existing bookmarks/links.
    url(
        r'^GPCRBrowser[/]?$',
        RedirectView.as_view(
            url='/allgpcrs',
            permanent=False,
        ),
        name='classification-gpcrbrowser',
    ),
    url(r'^levels_and_terms[/]?$', views.Classification.as_view(), name='classification-classification'),
    url(r'^figures[/]?$', views.ClassificationVisualizationsLanding.as_view(), name='classification-visualizations'),
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

