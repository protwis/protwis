from django.shortcuts import render
from django.conf import settings


def index(request):
    context = {
        'site_title': settings.SITE_TITLE,
        'menu_template': 'home/mainmenu_' + settings.SITE_NAME + '.html',
    }
    return render(request, 'home/index_' + settings.SITE_NAME + '.html', context)