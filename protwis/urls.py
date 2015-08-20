from django.conf.urls import patterns, include, url
from django.contrib import admin
from django.conf import settings


urlpatterns = patterns('',
    url(r'^$', include('home.urls')),
    url(r'^api/', include('api_' + settings.SITE_NAME + '.urls')),
    url(r'^admin/', include(admin.site.urls)),
    url(r'^common/', include('common.urls')),
    url(r'^protein/', include('protein.urls')),
    url(r'^family/', include('family.urls')),
    url(r'^mutations/', include('mutation.urls')),
    url(r'^interaction/', include('interaction.urls')),
    url(r'^residue/', include('residue.urls')),
    url(r'^alignment/', include('alignment.urls')),
    url(r'^similaritysearch/', include('similaritysearch.urls')),
    url(r'^phylogenetic_trees/', include('phylogenetic_trees.urls')),
    url(r'^similaritymatrix/', include('similaritymatrix.urls')),
    url(r'^structure/',include('structure.urls')),
)

if settings.DEBUG:
    import debug_toolbar
    urlpatterns += patterns('',
        url(r'^__debug__/', include(debug_toolbar.urls)),
    )
