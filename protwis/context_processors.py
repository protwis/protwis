from django.conf import settings
from django.core.cache import cache

from common.models import Citation

domain_switches = {"gpcrdb.org" : "gpcr", "gproteindb.org" : "gprotein", "arrestindb.org": "arrestin", "biasedsignalingatlas.org": "biasedsignalingatlas"}
inverse_domain_switches = {v : k for k,v in domain_switches.items()}

def get_current_site(domain, return_domain=False):
    if not domain in domain_switches:
        if return_domain:
            return inverse_domain_switches[settings.DEFAULT_SITE]
        else:
            return settings.DEFAULT_SITE
    else:
        if return_domain:
            return domain
        else:
            return domain_switches[domain]

def current_site(request):
    domain = request.get_host().lower()
    return {
    'current_site': get_current_site(domain)
    }

def site_title(request):
    domain = current_site(request)["current_site"]
    domain_titles = {"gpcr": "GPCRdb", "gprotein": "GproteinDb", "arrestin": "ArrestinDb", "biasedsignalingatlas": "Biased Signaling Atlas"}

    if not domain in domain_titles:
        domain = settings.DEFAULT_SITE

    return {
       'site_title': domain_titles[domain]
     }

def canonical_tag(request):
    citation_dict = cache.get("citation_urls")
    if citation_dict == None:
        citation_dict = {}
        citation_urls = Citation.objects.all().values_list("url", flat = True)
        for url in citation_urls:
            path = url.split(".org")[1]
            citation_dict[path] = url
        cache.set("citation_urls", citation_dict, 60*60*24*7)

    if request.path in citation_dict:
        return {
           'canonical_tag': citation_dict[request.path]
         }
    elif request.path == "" or request.path == "/":
        return {
           'canonical_tag': "https://" + request.get_host()
         }
    else:
        return {
           'canonical_tag': "https://gpcrdb.org" + request.path
         }

def documentation_url(request):
    return {
        'documentation_url': settings.DOCUMENTATION_URL
    }

def google_analytics(request):
    """
    Use the variables returned in this function to
    render your Google Analytics tracking code template.
    """
    if settings.GOOGLE_ANALYTICS_KEY:
        return {
            'google_analytics': settings.GOOGLE_ANALYTICS_KEY
        }
    return {}
