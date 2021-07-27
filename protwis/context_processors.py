from django.conf import settings

def current_site(request):
    domain_switches = {"gpcrdb.org" : "gpcr", "gproteindb.org" : "gprotein", "arrestindb.org": "arrestin"}
    domain = request.get_host().lower()

    if not domain in domain_switches:
        domain = "gpcrdb.org"

    return {
       'current_site': domain_switches[domain]
     }

def site_title(request):
    domain = current_site(request)["current_site"]
    domain_titles = {"gpcr": "GPCRdb", "gprotein": "GproteinDb", "arrestin": "ArrestinDb"}

    if not domain in domain_titles:
        domain = "gpcr"

    return {
       'site_title': domain_titles[domain]
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
