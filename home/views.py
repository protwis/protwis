from django.shortcuts import render
from django.conf import settings
from django.views.decorators.cache import cache_page
from django.http import JsonResponse
from django.db.models import F, Q
from django.views.generic import TemplateView
from django.core.cache import cache
from django.utils.html import escape

from django.contrib.postgres.aggregates import ArrayAgg

from protwis.context_processors import site_title
from news.models import News
from common.models import ReleaseNotes, ReleaseStatistics, Citation, DefaultCitation
from protein.models import Protein, ProteinCouplings
from structure.models import StructureComplexModel
from ligand.models import BiasedData, BiasedPathwaysAssay, Endogenous_GTP, BalancedLigands
from contactnetwork.models import InteractingResiduePair
from signprot.models import SignprotComplex, SignprotStructure
from googleapiclient.discovery import build
from oauth2client.service_account import ServiceAccountCredentials

from common.definitions import re_string_cite_us_page
from protwis.context_processors import get_current_site


from collections import OrderedDict
from urllib.parse import urlparse, parse_qs
import os
import re

re_cite_us_page = re.compile(re_string_cite_us_page)

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
def citations_json(request, output_type='list'):
    PUBLICATION_KEY = 'publication'

    citation_fields = [
        "url",
        "video",
        "docs",
        "main",
        "page_name",
    ]
    publication_fields = [
        "title",
        "authors",
        "year",
        "reference", 
        "journal__name",
        "web_link__index",
    ]

    publication_fields_aliases = {
        "title": "title",
        "authors": "authors",
        "year": "year",
        "reference": "reference",
        "journal__name": "journal_name",
        "web_link__index": "doi",
    }

    aliased_publication_fields = [publication_fields_aliases[f] for f in publication_fields]

    all_fields = citation_fields + ['publication__' + f for f in publication_fields]

    citations_q = Citation.objects.all().prefetch_related('publication')
    citations_q = citations_q.values("id",'publication__id',*all_fields)
    citations_q = citations_q.order_by("id")

    # Agregate publications
    citations_dict = OrderedDict()
    for citation in citations_q:
        citation_id = citation['id']
        publication_id = citation['publication__id']
        if citation_id not in citations_dict:
            new_cit = {}
            citations_dict[citation_id] = new_cit
            for f in citation_fields:
                new_cit[f] = citation[f]
            pubs = {}
            new_cit[PUBLICATION_KEY] = pubs
        else:
            new_cit = citations_dict[citation_id]
            pubs = new_cit[PUBLICATION_KEY]
        if publication_id is not None:
            pubs[publication_id] = {publication_fields_aliases[f]: citation['publication__' + f] for f in publication_fields}

    # get order of publications from citation_publication_through
    citation_publication_through = Citation.publication.through
    qcitpub = citation_publication_through.objects.all()
    qcitpub = qcitpub.values_list('id', 'citation_id', 'publication_id')
    qcitpub = qcitpub.order_by('id')

    citation_pub_order_dict = {}
    for id, citation_id, publication_id in qcitpub:
        if citation_id not in citation_pub_order_dict:
            citation_pub_order_dict[citation_id] = {'count':0,'index':{}}
            c_citation_pub_order_dict = citation_pub_order_dict[citation_id]
        c_citation_pub_order_dict['index'][publication_id] = c_citation_pub_order_dict['count']
        c_citation_pub_order_dict['count'] += 1

    for citation_id, citation in citations_dict.items():
        if citation_id not in citation_pub_order_dict:
            pubs = None
        else:
            c_citation_pub_order_dict = citation_pub_order_dict[citation_id]
            pubs = [v 
                    for k,v in sorted(citation[PUBLICATION_KEY].items(),
                                    key=lambda x: c_citation_pub_order_dict['index'][x[0]])
                ]
        citation[PUBLICATION_KEY] = pubs

    if output_type == 'dict' or output_type == 'object':
        response = JsonResponse((list(citations_dict.values())), safe=False)
    else:
        citations_list = []
        empty_pubs = [{f: None for f in aliased_publication_fields}]
        for citation in citations_dict.values():
            citation_l = []
            for f in citation_fields:
                v = citation[f]
                if f == PUBLICATION_KEY:
                    continue
                citation_l.append(v)

            pubs = citation[PUBLICATION_KEY]
            if pubs is None:

                pubs = empty_pubs
            elif len(pubs) == 0:

                pubs = empty_pubs
                
            for pub in pubs:
                citation_pub_l = list(citation_l)
                for f in aliased_publication_fields:
                    citation_pub_l.append(pub[f])
                citations_list.append(citation_pub_l)
        response = JsonResponse((citations_list), safe=False)
   
    return response

def escape_but_none(value):
    if value is not None:
        return escape(value)
    return value

def get_default_citation(main, output_type='list'):
    PUBLICATION_KEY = 'publication'

    citation_fields = [
        "main",

    ]
    publication_fields = [
        "title",
        "authors",
        "year",
        "reference", 
        "journal__name",
        "web_link__index",
    ]

    publication_fields_aliases = {
        "title": "title",
        "authors": "authors",
        "year": "year",
        "reference": "reference",
        "journal__name": "journal_name",
        "web_link__index": "doi",
    }

    aliased_publication_fields = [publication_fields_aliases[f] for f in publication_fields]

    all_fields = citation_fields + ['publication__' + f for f in publication_fields]

    citations_q = DefaultCitation.objects.filter(main=main).prefetch_related('publication')
    citations_q = citations_q.values("id",'publication__id',*all_fields)
    citations_q = citations_q.order_by("id")
    citations_q = list(citations_q)

    # Agregate publications
    citations_dict = OrderedDict()
    for citation in citations_q:
        citation_id = citation['id']
        publication_id = citation['publication__id']
        if citation_id not in citations_dict:
            new_cit = {}
            citations_dict[citation_id] = new_cit
            for f in citation_fields:
                new_cit[f] = citation[f]
            pubs = {}
            new_cit[PUBLICATION_KEY] = pubs
        else:
            new_cit = citations_dict[citation_id]
            pubs = new_cit[PUBLICATION_KEY]
        if publication_id is not None:
            pubs[publication_id] = {publication_fields_aliases[f]: escape_but_none(citation['publication__' + f]) for f in publication_fields}

    # get order of publications from *_citation_publication_through
    citation_publication_through = DefaultCitation.publication.through
    qcitpub = citation_publication_through.objects.all()
    qcitpub = qcitpub.values_list('id', 'defaultcitation_id', 'publication_id')
    qcitpub = qcitpub.order_by('id')

    citation_pub_order_dict = {}
    for id, citation_id, publication_id in qcitpub:
        if citation_id not in citation_pub_order_dict:
            citation_pub_order_dict[citation_id] = {'count':0,'index':{}}
            c_citation_pub_order_dict = citation_pub_order_dict[citation_id]
        c_citation_pub_order_dict['index'][publication_id] = c_citation_pub_order_dict['count']
        c_citation_pub_order_dict['count'] += 1

    for citation_id, citation in citations_dict.items():
        if citation_id not in citation_pub_order_dict:
            pubs = None
        else:
            c_citation_pub_order_dict = citation_pub_order_dict[citation_id]
            pubs = [v 
                    for k,v in sorted(citation[PUBLICATION_KEY].items(),
                                    key=lambda x: c_citation_pub_order_dict['index'][x[0]])
                ]
        citation[PUBLICATION_KEY] = pubs

    if output_type == 'dict' or output_type == 'object':
        return citations_dict
    else:
        citations_list = []
        empty_pubs = [{f: None for f in aliased_publication_fields}]
        for citation in citations_dict.values():
            citation_l = []
            for f in citation_fields:
                v = citation[f]
                if f == PUBLICATION_KEY:
                    continue
                citation_l.append(v)

            pubs = citation[PUBLICATION_KEY]
            if pubs is None:

                pubs = empty_pubs
            elif len(pubs) == 0:

                pubs = empty_pubs
                
            for pub in pubs:
                citation_pub_l = list(citation_l)
                for f in aliased_publication_fields:
                    citation_pub_l.append(pub[f])
                citations_list.append(citation_pub_l)
        return citations_list

def citation_json_by_url(request, output_type='list'):
    cache_flag = True
    replace_hostname = False #True for debuging ONLY_DEFAULT_CITATION_DOMAINS

    PUBLICATION_KEY = 'publication'
    ONLY_DEFAULT_CITATION_DBS = ['arrestin','bias','gprotein']
    ONLY_DEFAULT_CITATION_DOMAINS_TAGS = {'gproteindb.org':'gprotein',
                                          'arrestindb.org':'arrestin',
                                          'biasedsignalingatlas.org':'bias'}
    

    citation_fields = [
        "url",
        "video",
        "docs",
        "main",
        "page_name",
    ]
    publication_fields = [
        "title",
        "authors",
        "year",
        "reference", 
        "journal__name",
        "web_link__index",
    ]

    publication_fields_aliases = {
        "title": "title",
        "authors": "authors",
        "year": "year",
        "reference": "reference",
        "journal__name": "journal_name",
        "web_link__index": "doi",
    }

    aliased_publication_fields = [publication_fields_aliases[f] for f in publication_fields]

    all_fields = citation_fields + ['publication__' + f for f in publication_fields]

    input_url = request.GET.get('url')

    # Break URL into components
    parsed_input_url  = urlparse(input_url)

    found_citation = False
    only_default_citation = False
    no_citation = False
    main = None

    parsed_input_url_path = parsed_input_url.path
    parsed_input_url_hostname = parsed_input_url.hostname

    if replace_hostname:
        parsed_input_url_hostname = get_current_site(parsed_input_url_hostname,return_domain=True)


    parsed_input_url_path_no_leading_slash = parsed_input_url_path.lstrip('/')
    parsed_input_url_path = '/'+parsed_input_url_path_no_leading_slash

    # Pages that never show citations
    if parsed_input_url_path in {'/','',
                                 '/biased_signalling/bias_guidelines/',
                                 '/biased_signalling/bias_guidelines',
                                 '/biased_signalling/reference_selection',
                                 '/biased_signalling/reference_selection/'}:
        citations_q = []
        no_citation = True
    elif re_cite_us_page.match(parsed_input_url_path_no_leading_slash):
        no_citation = True
        citations_q = []
    else:
        hostname = None                       
        if replace_hostname and settings.DEBUG:
            hostname = parsed_input_url_hostname
        elif not settings.DEBUG:
            hostname = parsed_input_url.hostname
        if hostname is not None and hostname in ONLY_DEFAULT_CITATION_DOMAINS_TAGS:
            
            citations_q = []
            only_default_citation = True
            db_citation_dict_data = {}
            if cache_flag:
                db_citation_dict_data = cache.get("db_main_citation_dict", {})
            
            if hostname not in db_citation_dict_data:
                dcit_q = DefaultCitation.objects.get(main__icontains=ONLY_DEFAULT_CITATION_DOMAINS_TAGS[hostname])
                main = dcit_q.main
                db_citation_dict_data[hostname] = main
                cache.set("db_main_citation_dict", db_citation_dict_data, timeout=60 * 60 * 24 * 7)
            else:
                main = db_citation_dict_data[hostname]
        else:

            if parsed_input_url_path in {'/construct/analysis','/construct/analysis/'}:
                parsed_input_url_path = '/construct/analysis'+'#'+parsed_input_url.fragment


            if parsed_input_url.path.startswith('/biased_signalling/') or parsed_input_url.path == '/biased_signalling':
                parsed_input_url_path = '/biased_signalling/'
            elif parsed_input_url.path == '/drugs/targets_venn':
                parsed_input_url_path = '/drugs/target_venn/'
            cannon_url = 'https://'+parsed_input_url_hostname + parsed_input_url_path
            cannon_url_star = 'https://'+parsed_input_url_hostname + os.path.join(parsed_input_url_path,'*')
            if parsed_input_url_path in {'/', '','/*','*'}:
                cannon_url_alt_slash = None
            elif parsed_input_url_path.endswith('/'):
                cannon_url_alt_slash = 'https://'+parsed_input_url_hostname + parsed_input_url_path[:-1]
            else:
                cannon_url_alt_slash = 'https://'+parsed_input_url_hostname + parsed_input_url_path + '/'
            for u in [cannon_url, cannon_url_star,cannon_url_alt_slash]:
                if u is None:
                    continue
                citations_q = Citation.objects.filter(url=u).prefetch_related('publication')
                citations_q = citations_q.values("id",'publication__id',*all_fields)
                citations_q = citations_q.order_by("id")
                citations_q = list(citations_q)
                if len(citations_q) > 0:
                    found_citation = True
                    break
            if not found_citation:
                db_citation_dict_data = cache.get("db_citation_dict", None)
                if db_citation_dict_data is None or not cache_flag:
                    db_citation_dict_data = {}
                    citations_arrestindb_q = Citation.objects.filter(main__icontains='arrestin')
                    citations_bsa_q = Citation.objects.filter(main__icontains='bias')
                    citations_gproteindb_q = Citation.objects.filter(main__icontains='gprotein')
                    citations_grpcrdb_q = Citation.objects.filter(main__icontains='gpcr')

                    db_citation_q_list = [citations_arrestindb_q, citations_bsa_q, citations_gproteindb_q, citations_grpcrdb_q]
                    db_citation_q_list = [db_citation_q.values('id','url','main').order_by('id') for db_citation_q in db_citation_q_list]

                    db_citation_q_list = [list(db_citation_q) for db_citation_q in db_citation_q_list]
                    main_names = [db_citation_q[0]['main'] for db_citation_q in db_citation_q_list]
                    db_citation_dict_data['main_names'] = main_names
                    url_path_dict = {}
                    for db_citation_q in db_citation_q_list:
                        for item in db_citation_q:
                            parsed_item_url = urlparse(item['url'])
                            if parsed_item_url.path in {'/', '','/*','*'}:
                                continue
                            v = (item['id'],item['main'])
                            parsed_item_url_path = parsed_item_url.path
                            if parsed_item_url_path in {'/construct/analysis','/construct/analysis/'}:
                                parsed_item_url_path = '/construct/analysis'+'#'+parsed_item_url.fragment
                            url_path_dict[parsed_item_url_path] = v
    
                            if parsed_item_url_path.endswith('*'):
                                url_path_dict[parsed_item_url_path[:-1]] = v
                            elif parsed_item_url_path.endswith('/'):
                                url_path_dict[parsed_item_url_path[:-1]] = v
                            else:
                                url_path_dict[parsed_item_url_path + '/'] = v
                    db_citation_dict_data['url_path_dict'] = url_path_dict
                    cache.set("db_citation_dict", db_citation_dict_data, timeout=60 * 60 * 24 * 7)
                url_path_dict = db_citation_dict_data['url_path_dict']
                citation_id,main = url_path_dict.get(parsed_input_url_path, (None,None))
                if citation_id is not None:
                    citations_q = Citation.objects.filter(id=citation_id).prefetch_related('publication')
                    citations_q = citations_q.values("id",'publication__id',*all_fields)
                    citations_q = citations_q.order_by("id")
                    citations_q = list(citations_q)
                else:
                    if db_citation_dict_data is None:
                        db_citation_dict_data = {}

                    if 'glob_path_dict' not in db_citation_dict_data or not cache_flag:
                        glob_path_dict = {}
                        url_path_dict = db_citation_dict_data['url_path_dict']
                        for k,v in url_path_dict.items():
                            if '*' in k:
                                k_re = '^'
                                for k_p in k.split('*'):
                                    k_re += re.escape(k_p)
                                    k_re += '.*'
                                k_re += '$'
                                glob_path_dict[re.compile(k_re)] = v
                        db_citation_dict_data['glob_path_dict'] = glob_path_dict
                        cache.set("db_citation_dict", db_citation_dict_data, timeout=60 * 60 * 24 * 7)
                    for k_re, v in db_citation_dict_data['glob_path_dict'].items():
                        if re.match(k_re, parsed_input_url_path):
                            citation_id,main = v
                            if citation_id is not None:
                                citations_q = Citation.objects.filter(id=citation_id).prefetch_related('publication')
                                citations_q = citations_q.values("id",'publication__id',*all_fields)
                                citations_q = citations_q.order_by("id")
                                citations_q = list(citations_q)
                            else:
                                citations_q = []
                            break

    # Agregate publications
    citations_dict = OrderedDict()
    for citation in citations_q:
        citation_id = citation['id']
        publication_id = citation['publication__id']
        if citation_id not in citations_dict:
            new_cit = {}
            citations_dict[citation_id] = new_cit
            for f in citation_fields:
                new_cit[f] = citation[f]
            pubs = {}
            new_cit[PUBLICATION_KEY] = pubs
        else:
            new_cit = citations_dict[citation_id]
            pubs = new_cit[PUBLICATION_KEY]
        if publication_id is not None:
            pubs[publication_id] = {publication_fields_aliases[f]: escape_but_none(citation['publication__' + f]) for f in publication_fields}

    if main is None:
        for cit in citations_dict.values():
            main = cit['main']
            break
    
    if not only_default_citation:
        for db in ONLY_DEFAULT_CITATION_DBS:
            if main is not None:
                if db.lower() in main.lower():
                    only_default_citation = True
                    break

    # get order of publications from citation_publication_through
    citation_publication_through = Citation.publication.through
    qcitpub = citation_publication_through.objects.all()
    qcitpub = qcitpub.values_list('id', 'citation_id', 'publication_id')
    qcitpub = qcitpub.order_by('id')

    citation_pub_order_dict = {}
    for id, citation_id, publication_id in qcitpub:
        if citation_id not in citation_pub_order_dict:
            citation_pub_order_dict[citation_id] = {'count':0,'index':{}}
            c_citation_pub_order_dict = citation_pub_order_dict[citation_id]
        c_citation_pub_order_dict['index'][publication_id] = c_citation_pub_order_dict['count']
        c_citation_pub_order_dict['count'] += 1

    for citation_id, citation in citations_dict.items():
        if citation_id not in citation_pub_order_dict:
            pubs = None
        else:
            c_citation_pub_order_dict = citation_pub_order_dict[citation_id]
            pubs = [v 
                    for k,v in sorted(citation[PUBLICATION_KEY].items(),
                                    key=lambda x: c_citation_pub_order_dict['index'][x[0]])
                ]
        citation[PUBLICATION_KEY] = pubs
    
    if not no_citation:
        default_citation = get_default_citation(main,output_type=output_type)

    if output_type == 'dict' or output_type == 'object':
        if no_citation:
            default_citation = {}
        response = JsonResponse((list(citations_dict.values()),main,only_default_citation,no_citation,list(default_citation.values())), safe=False)
    else:
        if no_citation:
            default_citation = []
        citations_list = []
        empty_pubs = [{f: None for f in aliased_publication_fields}]
        for citation in citations_dict.values():
            citation_l = []
            for f in citation_fields:
                v = citation[f]
                if f == PUBLICATION_KEY:
                    continue
                citation_l.append(v)

            pubs = citation[PUBLICATION_KEY]
            if pubs is None:

                pubs = empty_pubs
            elif len(pubs) == 0:

                pubs = empty_pubs
                
            for pub in pubs:
                citation_pub_l = list(citation_l)
                for f in aliased_publication_fields:
                    citation_pub_l.append(pub[f])
                citations_list.append(citation_pub_l)
        response = JsonResponse((citations_list,main,only_default_citation,no_citation,default_citation), safe=False)
   
    return response

def cite_us(request, site):
    context = {'site': site}
    return render(request, 'home/cite_us.html', context)
 