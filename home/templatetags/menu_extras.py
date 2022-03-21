from django import template
from django.conf import settings

register = template.Library()

@register.simple_tag
def get_hostnames():
    if settings.DEBUG or not settings.HUB_ENABLED:
        return {
            'gpcr': "",
            'gprotein': "",
            'arrestin': "",
        }
    else:
        return {
            'gpcr': "https://gpcrdb.org",
            'gprotein': "https://gproteindb.org",
            'arrestin': "https://arrestindb.org",
        }
