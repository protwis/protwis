from django.shortcuts import render
from django.conf import settings

from news.models import News
from common.models import ReleaseNotes, ReleaseStatistics


def index(request):
    context = {}
    
    # title of the page
    context['site_title'] = settings.SITE_TITLE

    # get news 
    context['news'] = News.objects.order_by('-date').all()[:3]

    # get release notes
    context['release_notes'] = ReleaseNotes.objects.all()[0]

    # stats
    context['release_statistics'] = ReleaseStatistics.objects.filter(release=context['release_notes'])

    return render(request, 'home/index_{}.html'.format(settings.SITE_NAME), context)