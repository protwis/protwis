from django import template
from django.conf import settings

register = template.Library()

def mainmenu(domain):
    return {
        'site_title': settings.SITE_TITLE,
        'menu_template': 'home/mainmenu_' + settings.SITE_NAME + '.html',
        'logo_path': 'home/logo/' + settings.SITE_NAME + '/main.png',
        'documentation_url': settings.DOCUMENTATION_URL,
    }

register.inclusion_tag('home/mainmenu.html')(mainmenu)
