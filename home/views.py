from django.shortcuts import render
from django.conf import settings

from news.models import News
from common.models import ReleaseNotes, ReleaseStatistics


def index(request):
    request.session.flush()

    context = {}
    
    # title of the page
    context['site_title'] = settings.SITE_TITLE

    # analytics
    context['google_analytics_key'] = settings.GOOGLE_ANALYTICS_KEY

    # get news
    context['news'] = News.objects.order_by('-date').all()[:3]

    # get release notes
    try:
        context['release_notes'] = ReleaseNotes.objects.all()[0]
        context['release_statistics'] = ReleaseStatistics.objects.filter(release=context['release_notes'])
    except IndexError:
        context['release_notes'] = ''
        context['release_statistics'] = []

    return render(request, 'home/index_{}.html'.format(settings.SITE_NAME), context)