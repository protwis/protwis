from django.shortcuts import render
from django.conf import settings
from django.views.decorators.cache import cache_page

from news.models import News
from common.models import ReleaseNotes, ReleaseStatistics

@cache_page(60 * 60 * 24)
def index(request):
    request.session.flush()

    context = {}

    # title of the page
    context['site_title'] = settings.SITE_TITLE
    context['documentation_url'] = settings.DOCUMENTATION_URL

    # analytics
    context['google_analytics_key'] = settings.GOOGLE_ANALYTICS_KEY

    if settings.GOOGLE_ANALYTICS_API:
         # Based on https://developers.google.com/analytics/devguides/reporting/core/v3/quickstart/service-py
        from apiclient.discovery import build
        from oauth2client.service_account import ServiceAccountCredentials
        # Define the auth scopes to request.
        scope = 'https://www.googleapis.com/auth/analytics.readonly'
        key_file_location = settings.GOOGLE_ANALYTICS_API

        # Fetched from API -- look at original code to re-fetch if changes.
        profile_id = '77082434' 

        # Authenticate and construct service.
        credentials = ServiceAccountCredentials.from_json_keyfile_name(
                key_file_location, scopes=[scope])
        # Build the service object.
        service = build('analytics', 'v3', credentials=credentials)

        users_year = service.data().ga().get(ids='ga:' + profile_id,start_date='365daysAgo',end_date='today',metrics='ga:users').execute().get('rows')[0][0]
        users_month = service.data().ga().get(ids='ga:' + profile_id,start_date='30daysAgo',end_date='today',metrics='ga:users').execute().get('rows')[0][0]

        context['users'] = "GPCRdb had {:,} users since this date last year and {:,} users in the last 30 days (<a href='https://analytics.google.com'>Google Analytics</a>).".format(int(users_year),int(users_month))

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
