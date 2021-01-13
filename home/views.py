from django.shortcuts import render
from django.conf import settings
from django.views.decorators.cache import cache_page
from django.http import JsonResponse


from news.models import News
from common.models import ReleaseNotes, ReleaseStatistics, Citation
from googleapiclient.discovery import build
from oauth2client.service_account import ServiceAccountCredentials



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
        # from googleapiclient.discovery import build
        # from oauth2client.service_account import ServiceAccountCredentials
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

        users_year = service.data().ga().get(ids='ga:' + profile_id, start_date='365daysAgo', end_date='today', metrics='ga:users').execute().get('rows')[0][0]
        users_month = service.data().ga().get(ids='ga:' + profile_id, start_date='30daysAgo', end_date='today', metrics='ga:users').execute().get('rows')[0][0]

        context['users'] = "GPCRdb had {:,} different users since this date last year and {:,} users in the last 30 days (<a href='https://analytics.google.com'>Google Analytics</a>).".format(int(users_year), int(users_month))

    # get news
    context['news'] = News.objects.order_by('-date').all()[:3]

    # get release notes
    try:
        context['release_notes'] = ReleaseNotes.objects.all()[0]
        rel_stats = list(ReleaseStatistics.objects.filter(release=context['release_notes'])\
                    .values_list("statistics_type__name", "value"))

        # Create dictionary and process part of the results
        stats = {}
        context['release_statistics'] = []
        for entry in rel_stats:
            stats[entry[0]] = entry[1]
            if "Exp." not in entry[0] and "models" not in entry[0]:
                context['release_statistics'].append({"statistics_type": entry[0], "value": entry[1]})


        # Adjusted formatting for release notes

        # To be extended wtih G proteins and Arrestins
        context['release_statistics'].insert(2, {"statistics_type": "Structures", "value": "GPCRs: {}, G proteins: {}".format(stats["Exp. GPCR structures"], stats["Exp. Gprotein structures"]) })
        #context['release_statistics'].insert(2, {"statistics_type": "Structures", "value": "GPCRs: {}, G proteins: {}, Arrestins: {}".format(stats["Exp. GPCR structures"], stats["Exp. Gprotein structures"], stats["Exp. Arrestin structures"]) })
        context['release_statistics'].insert(3, {"statistics_type": "Structure models", "value": "GPCRs: {}, GPCR-G protein complexes: {}".format(stats["GPCR structure models"], stats["GPCR-G protein structure models"]) })


    except IndexError:
        context['release_notes'] = ''
        context['release_statistics'] = []

    return render(request, 'home/index_{}.html'.format(settings.SITE_NAME), context)

# @cache_page(60 * 60 * 24)
def citations_json(request):
    context = {}
    citations_q = Citation.objects.all().values_list("url", "video", "docs", "main", "page_name", "publication__title", "publication__authors", "publication__year", "publication__reference",
                                                     "publication__journal__name", "publication__web_link__index").order_by("-publication__year", "page_name")
    response = JsonResponse(list(citations_q), safe=False)
    return response
