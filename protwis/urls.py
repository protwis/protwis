from django.conf.urls import include, url,handler404, handler500
from django.contrib import admin
from django.conf import settings
from protwis import views


urlpatterns = [
    url(r'^', include('home.urls')),
    url(r'^services/', include('api_' + settings.SITE_NAME + '.urls')),
    url(r'^admin/', admin.site.urls),
    url(r'^common/', include('common.urls')),
    url(r'^protein/', include('protein.urls')),
    url(r'^family/', include('family.urls')),
    url(r'^mutations/', include('mutation.urls')),
    url(r'^news/', include('news.urls')),
    url(r'^interaction/', include('interaction.urls')),
    url(r'^residue/', include('residue.urls')),
    url(r'^alignment/', include('alignment.urls')),
    url(r'^similaritysearch/', include('similaritysearch.urls')),
    url(r'^pages/', include('pages.urls')),
    url(r'^phylogenetic_trees/', include('phylogenetic_trees.urls')),
    url(r'^similaritymatrix/', include('similaritymatrix.urls')),
    url(r'^structure/',include('structure.urls')),
    url(r'^construct/',include('construct.urls')),
    url(r'^sitesearch/',include('sitesearch.urls')),
    url(r'^drugs/',include('drugs.urls')),
    url(r'^signprot/',include('signprot.urls')),
    url(r'^ligand/',include('ligand.urls')),
    url(r'^mutational_landscape/',include('mutational_landscape.urls')),
    url(r'^contactnetwork/', include('contactnetwork.urls')),
    url(r'^structure_comparison/', include('contactnetwork.urls')), # Remapping this url to the contact network
    url(r'^seqsign/', include('seqsign.urls')),
    url(r'^angles/', include('angles.urls')),
    url(r'^hotspots/', include('hotspots.urls')),
]

handler404 = views.error404
handler500 = views.error500

if settings.DEBUG:
    import debug_toolbar
    urlpatterns.append( url(r'^__debug__/', include(debug_toolbar.urls)) )
