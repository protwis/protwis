from django.shortcuts import render
from django.conf import settings
from django.views.decorators.cache import cache_page
from django.http import JsonResponse
from django.db.models import F, Q
from django.views.generic import TemplateView

from protwis.context_processors import site_title
from news.models import News
from common.models import ReleaseNotes, ReleaseStatistics, Citation
from protein.models import Protein, ProteinCouplings
from structure.models import StructureComplexModel
from ligand.models import BiasedData, BiasedPathwaysAssay, Endogenous_GTP, BalancedLigands
from contactnetwork.models import InteractingResiduePair
from signprot.models import SignprotComplex, SignprotStructure
from googleapiclient.discovery import build
from oauth2client.service_account import ServiceAccountCredentials


# @cache_page(60 * 60 * 24)
def index(request):
    request.session.flush()

    context = {}

    # title of the page
    context["site_title"] = site_title(request)["site_title"]  # settings.SITE_TITLE
    context["documentation_url"] = settings.DOCUMENTATION_URL

    # development/production
    context["debug"] = settings.DEBUG

    # analytics
    context["google_analytics_key"] = settings.GOOGLE_ANALYTICS_KEY

    if settings.GOOGLE_ANALYTICS_API:
        # Based on https://developers.google.com/analytics/devguides/reporting/core/v3/quickstart/service-py
        # from googleapiclient.discovery import build
        # from oauth2client.service_account import ServiceAccountCredentials
        # Define the auth scopes to request.
        scope = "https://www.googleapis.com/auth/analytics.readonly"
        key_file_location = settings.GOOGLE_ANALYTICS_API

        # Fetched from API -- look at original code to re-fetch if changes.
        profile_id = "77082434"

        # Authenticate and construct service.
        credentials = ServiceAccountCredentials.from_json_keyfile_name(key_file_location, scopes=[scope])
        # Build the service object.
        service = build("analytics", "v3", credentials=credentials)

        users_year = service.data().ga().get(ids="ga:" + profile_id, start_date="365daysAgo", end_date="today", metrics="ga:users").execute().get("rows")[0][0]
        users_month = service.data().ga().get(ids="ga:" + profile_id, start_date="30daysAgo", end_date="today", metrics="ga:users").execute().get("rows")[0][0]

        context["users"] = "Together, they have served {:,}/".format(int(users_month)) +\
                           "{:,} users in the last month/year (<a href='https://analytics.google.com'>Google Analytics</a>)".format(int(users_year))

    # get news
    context["news"] = News.objects.order_by("-date").all()[:3]
    # Setting the headers for data boxes
    headers={'GproteinDb' : {0: {"statistics_type": '<span class="stats_title"><b>Sequences</b></span>', "value": ''},
                             2: {"statistics_type": '<span class="stats_title"><b>Couplings</b></span>', "value": ''},
                             3: {"statistics_type": '<span class="stats_title"><b>Structures</b></span>', "value": ''},
                             6: {"statistics_type": '<span class="stats_title"><b>Structure models</b></span>', "value": ''},
                             8: {"statistics_type": '<span class="stats_title"><b>Structure interactions</b></span>', "value": '',},
                             9: {"statistics_type": '<span class="stats_title"><b>Mutations</b></span>', "value": ''}},
             'ArrestinDb': {0: {"statistics_type": '<span class="stats_title"><b>Sequences</b></span>', "value": ''},
                            2: {"statistics_type": '<span class="stats_title"><b>Couplings</b></span>', "value": ''},
                            3: {"statistics_type": '<span class="stats_title"><b>Structures</b></span>', "value": ''},
                            6: {"statistics_type": '<span class="stats_title"><b>Structure interactions</b></span>', "value": '',},
                            7: {"statistics_type": '<span class="stats_title"><b>Mutations</b></span>', "value": ''}},
             'Biased Signaling Atlas': {0: {"statistics_type": '<span class="stats_title"><b>Biased ligands</b></span>', "value": ''},
                                      3: {"statistics_type": '<span class="stats_title"><b>Pathways</b></span>', "value": ''},
                                      4: {"statistics_type": '<span class="stats_title"><b>Pathway-preferring ligands</b></span>', "value": '',},
                                      5: {"statistics_type": '<span class="stats_title"><b>Reference ligands</b></span>', "value": ''}},
             'GPCRdb': {}}
    # get release notes
    try:
        context["release_notes"] = ReleaseNotes.objects.all()[0]
        ### DB specific release notes
        # context["release_notes"] = ReleaseNotes.objects.filter(database=context["site_title"])[0]
        rel_stats = list(ReleaseStatistics.objects.filter(release=context["release_notes"], database=context["site_title"]).values_list("statistics_type__name", "value"))
        # Create dictionary and process part of the results
        context["release_statistics"] = []
        count = 0
        for entry in rel_stats:
            if count in headers[context["site_title"]].keys():
                context["release_statistics"].append(headers[context["site_title"]][count])
            context["release_statistics"].append(
                {
                    "statistics_type": '<span class="stats_entry">' + entry[0].split(' '+context["site_title"])[0] + "</span>",
                    "value": '<span  class="stats_value">' + "{:,}".format(entry[1]) + "</span>",
                }
            )
            count+=1
    except IndexError:
        context["release_notes"] = ""
        context["release_statistics"] = []

    return render(request, "home/index.html", context)


@cache_page(60 * 60 * 24 * 7)
def citations_json(request):
    citations_q = (
        Citation.objects.all()
        .values_list(
            "url",
            "video",
            "docs",
            "main",
            "page_name",
            "publication__title",
            "publication__authors",
            "publication__year",
            "publication__reference",
            "publication__journal__name",
            "publication__web_link__index",
        )
        .order_by("-publication__year", "page_name")
    )
    response = JsonResponse(list(citations_q), safe=False)
    return response


class citeGPCRdb(TemplateView):
    template_name = 'home/cite_gpcrdb.html'

    def get_context_data(self, **kwargs):

        context = super().get_context_data(**kwargs)
        return context

class citeGproteinDb(TemplateView):
    template_name = 'home/cite_gproteindb.html'

    def get_context_data(self, **kwargs):

        context = super().get_context_data(**kwargs)
        return context
