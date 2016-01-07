from django import template
from django.conf import settings


register = template.Library()

def site_title():
    return {
        'site_title': settings.SITE_TITLE,
    }

register.inclusion_tag('home/title.html')(site_title)