from django.shortcuts import render
from django.conf import settings
from django.views.decorators.cache import cache_page
from django.http import JsonResponse
from django.db.models import F

from protwis.context_processors import site_title
from news.models import News
from common.models import ReleaseNotes, ReleaseStatistics, Citation
from protein.models import Protein, ProteinGProteinPair
from structure.models import StructureComplexModel
from contactnetwork.models import InteractingResiduePair
from signprot.models import SignprotComplex, SignprotStructure
from googleapiclient.discovery import build
from oauth2client.service_account import ServiceAccountCredentials



@cache_page(60 * 60 * 24)
def index(request):
    request.session.flush()

    context = {}

    # title of the page
    context['site_title'] = site_title(request)['site_title']#settings.SITE_TITLE
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
        context['release_statistics'] = []
        if context['site_title']=='GproteinDb':
            context['release_statistics'].append({"statistics_type": "<span class=\"stats_entry\">" + "Human G proteins" + "</span>", "value": "<span  class=\"stats_value\">" + "{:,}".format(Protein.objects.filter(family__parent__parent__name='Alpha', species__common_name='Human', accession__isnull=False).count()) + "</span>"})
            context['release_statistics'].append({"statistics_type": "<span class=\"stats_entry\">" + "Species orthologs" + "</span>", "value": "<span  class=\"stats_value\">" + "{:,}".format(Protein.objects.filter(family__parent__parent__name='Alpha', accession__isnull=False).count()) + "</span>"})

            context['release_statistics'].append({"statistics_type": "<span class=\"stats_entry stats_title\"><i>Experimental structures</i></span>", "value" : "<span  class=\"stats_value\"></span>"})
            signcomp = SignprotComplex.objects.all()
            context['release_statistics'].append({"statistics_type": "<span class=\"stats_entry\">" + "G proteins" + "</span>", "value": "<span  class=\"stats_value\">" + "{:,}".format(signcomp.count()+SignprotStructure.objects.all().count()) + "</span>"})
            context['release_statistics'].append({"statistics_type": "<span class=\"stats_entry\">" + "GPCR-G protein complexes" + "</span>", "value": "<span  class=\"stats_value\">" + "{:,}".format(SignprotComplex.objects.all().count()) + "</span>"})

            context['release_statistics'].append({"statistics_type": "<span class=\"stats_entry stats_title\"><i>Structure models</i></span>", "value" : "<span  class=\"stats_value\"></span>"})
            context['release_statistics'].append({"statistics_type": "<span class=\"stats_entry\">" + "GPCR-G protein complexes" + "</span>", "value": "<span  class=\"stats_value\">" + "{:,}".format(StructureComplexModel.objects.all().count()-SignprotComplex.objects.filter(structure__refined=True).count()) + "</span>"})
            context['release_statistics'].append({"statistics_type": "<span class=\"stats_entry\">" + "Refined complex structures" + "</span>", "value": "<span  class=\"stats_value\">" + "{:,}".format(signcomp.filter(structure__refined=True).count()) + "</span>"})

            context['release_statistics'].append({"statistics_type": "<span class=\"stats_entry stats_title\"><i>Structure interactions</i></span>", "value" : "<span  class=\"stats_value\"></span>"})
            interface_interactions_count = InteractingResiduePair.objects.filter(referenced_structure__in=signcomp.values_list('structure', flat=True)).exclude(res1__protein_conformation_id=F('res2__protein_conformation_id')).count()
            context['release_statistics'].append({"statistics_type": "<span class=\"stats_entry\">" + "GPCR-G protein interface" + "</span>", "value": "<span  class=\"stats_value\">" + "{:,}".format(interface_interactions_count) + "</span>"})

            context['release_statistics'].append({"statistics_type": "<span class=\"stats_entry stats_title\"><i>Couplings</i></span>", "value" : "<span  class=\"stats_value\"></span>"})
            context['release_statistics'].append({"statistics_type": "<span class=\"stats_entry\">" + "GPCR-G protein coupling" + "</span>", "value": "<span  class=\"stats_value\">" + "{:,}".format(ProteinGProteinPair.objects.all().count()) + "</span>"})

            context['release_statistics'].append({"statistics_type": "<span class=\"stats_entry stats_title\"><i>Mutations</i></span>", "value" : "<span  class=\"stats_value\"></span>"})
            context['release_statistics'].append({"statistics_type": "<span class=\"stats_entry\">" + "Interface mutations" + "</span>", "value": "<span  class=\"stats_value\">" + "{:,}".format(54) + "</span>"})
        else:
            rename_dictionary = {"Exp. GPCR structures" : "GPCRs", "Exp. Gprotein structures" : "G proteins", "GPCR structure models": "GPCRs", "GPCR-G protein structure models": "GPCR-G protein complexes", "Refined GPCR structures": "Refined GPCR structures"}
            first_struct = -1
            first_model = -1
            count = 0
            for entry in rel_stats:
                if first_struct < 0 and "Exp." in entry[0]:
                    first_struct = count
                elif first_model < 0 and "model" in entry[0]:
                    first_model = count

                if entry[0] in rename_dictionary:
                    context['release_statistics'].append({"statistics_type": "<span class=\"stats_entry\">" + rename_dictionary[entry[0]] + "</span>", "value": "<span class=\"stats_value\">" + "{:,}".format(entry[1]) + "</span>"})
                else:
                    context['release_statistics'].append({"statistics_type": "<span class=\"stats_entry\">" + entry[0] + "</span>", "value": "<span  class=\"stats_value\">" + "{:,}".format(entry[1]) + "</span>"})
                count += 1

            # Adjusted formatting for release notes
            context['release_statistics'].insert(first_model, {"statistics_type": "<span class=\"stats_entry stats_title\"><i>Structure models</i></span>", "value" : "<span  class=\"stats_value\"></span>"})
            context['release_statistics'].insert(first_struct, {"statistics_type": "<span class=\"stats_entry stats_title\"><i>Experimental structures</i></span>", "value" : "<span  class=\"stats_value\"></span>"})


    except IndexError:
        context['release_notes'] = ''
        context['release_statistics'] = []

    return render(request, 'home/index.html', context)

# @cache_page(60 * 60 * 24)
def citations_json(request):
    context = {}
    citations_q = Citation.objects.all().values_list("url", "video", "docs", "main", "page_name", "publication__title", "publication__authors", "publication__year", "publication__reference",
                                                     "publication__journal__name", "publication__web_link__index").order_by("-publication__year", "page_name")
    response = JsonResponse(list(citations_q), safe=False)
    return response
