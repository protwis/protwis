from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.conf import settings
from django.db.models import Count
from django.core.cache import cache
from django.views.decorators.cache import cache_page

from drugs.models import Drugs
from protein.models import Protein

import re
import json
import numpy as np

def get_spaced_colors(n):
    max_value = 16581375 #255**3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]

    return ['#'+i for i in colors] # HEX colors
    # return [(int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors] # RGB colors

def striphtml(data):
    p = re.compile(r'<.*?>')
    return p.sub('', data)

# @cache_page(60*5) #  5 min
def drugstatistics(request):

    # ===== drugtargets =====
    drugtargets_raw = Protein.objects.filter(drugs__status='approved').values('entry_name').annotate(value=Count('drugs__name', distinct = True)).order_by('-value').order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugtargets_raw))
    drugtargets = []
    for i, drugtarget in enumerate(drugtargets_raw):
        drugtarget['label'] = drugtarget['entry_name'].replace("_human","").upper()
        # drugtarget['color'] = str(list_of_hec_colors[i])
        del drugtarget['entry_name']
        drugtargets.append(drugtarget)

    # ===== drugfamilies =====
    drugfamilies_raw = Protein.objects.filter(drugs__status='approved').values('family_id__parent__name').annotate(value=Count('drugs__name', distinct = True)).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugfamilies_raw))
    drugfamilies = []
    for i, drugfamily in enumerate(drugfamilies_raw):
        drugfamily['label'] = striphtml(drugfamily['family_id__parent__name']).replace(" receptors","")
        drugfamily['color'] = str(list_of_hec_colors[i])
        del drugfamily['family_id__parent__name']
        drugfamilies.append(drugfamily)

    # ===== drugclas =====
    drugclasses_raw = Protein.objects.filter(drugs__status='approved').values('family_id__parent__parent__parent__name').annotate(value=Count('drugs__name', distinct = True)).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugclasses_raw)+1)
    drugclasses = []
    for i, drugclas in enumerate(drugclasses_raw):
        drugclas['label'] = drugclas['family_id__parent__parent__parent__name']
        drugclas['color'] = str(list_of_hec_colors[i+1])
        del drugclas['family_id__parent__parent__parent__name']
        drugclasses.append(drugclas)

    # ===== drugtypes =====
    drugtypes_raw = Drugs.objects.values('drugtype').filter(status='approved').annotate(value=Count('name', distinct = True)).order_by('value')

    list_of_hec_colors = get_spaced_colors(len(drugtypes_raw)+5)
    drugtypes = []
    for i, drugtype in enumerate(drugtypes_raw):
        drugtype['label'] = drugtype['drugtype']
        drugtype['color'] = str(list_of_hec_colors[i])
        del drugtype['drugtype']
        drugtypes.append(drugtype)

    # ===== drugindications =====
    drugindications_raw = Drugs.objects.values('indication').filter(status='approved').annotate(value=Count('name', distinct = True)).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugindications_raw))
    drugindications = []
    for i, drugindication in enumerate(drugindications_raw):
        drugindication['label'] = drugindication['indication']
        drugindication['color'] = str(list_of_hec_colors[i])
        del drugindication['indication']
        drugindications.append(drugindication)

    # ===== drugtimes =====
    drugtime_raw = Drugs.objects.values('approval').filter(status='approved').annotate(y=Count('name', distinct = True)).order_by('approval')

    drugtimes = []
    running_total = 0

    for i, time in enumerate(range(1942,2017,1)):
        if str(time) in [i['approval'] for i in drugtime_raw]:
            y = [i['y'] for i in drugtime_raw if i['approval']==str(time)][0] + running_total
            x = time
            running_total = y
        else:
            x = time
            y = running_total
        if time % 2 == 0:
            drugtimes.append({'x':x,'y':y})

    drugs_over_time = [{"values": drugtimes, "yAxis": "1", "key": "GPCRs"}, {'values': [{'y': 2, 'x': '1942'}, {'x': '1944', 'y': 2}, {'y': 6, 'x': '1946'}, {'y': 9, 'x': '1948'}, {'y': 18, 'x': '1950'}, {'y': 30, 'x': '1952'}, {'y': 55, 'x': '1954'}, {'y': 72, 'x': '1956'}, {'y': 98, 'x': '1958'}, {'y': 131, 'x': '1960'}, {'y': 153, 'x': '1962'}, {'y': 171, 'x': '1964'}, {'y': 188, 'x': '1966'}, {'y': 205, 'x': '1968'}, {'y': 224, 'x': '1970'}, {'y': 242, 'x': '1972'}, {'y': 265, 'x': '1974'}, {'y': 300, 'x': '1976'}, {'y': 340, 'x': '1978'}, {'y': 361, 'x': '1980'}, {'y': 410, 'x': '1982'}, {'y': 442, 'x': '1984'}, {'y': 499, 'x': '1986'}, {'y': 542, 'x': '1988'}, {'y': 583, 'x': '1990'}, {'y': 639, 'x': '1992'}, {'y': 686, 'x': '1994'}, {'y': 779, 'x': '1996'}, {'y': 847, 'x': '1998'}, {'y': 909, 'x': '2000'}, {'y': 948, 'x': '2002'}, {'y': 1003, 'x': '2004'}, {'y': 1041, 'x': '2006'}, {'y': 1078, 'x': '2008'}, {'y': 1115, 'x': '2010'}, {'y': 1177, 'x': '2012'}, {'y': 1239, 'x': '2014'}, {'y': 1286, 'x': '2016'}], 'key': 'All FDA drugs', 'yAxis': '1'}]


    return render(request, 'drugstatistics.html', {'drugtypes':drugtypes, 'drugindications':drugindications, 'drugtargets':drugtargets, 'drugfamilies':drugfamilies, 'drugclasses':drugclasses, 'drugs_over_time':drugs_over_time})

@cache_page(60*5) #  5 min
def drugbrowser(request):
    # Get drugdata from here somehow

    name_of_cache = 'drug_browser2'

    context = cache.get(name_of_cache)

    if context==None:
        context = list()

        drugs = Drugs.objects.all().prefetch_related('target__family__parent__parent__parent')

        for drug in drugs:
            drugname = drug.name
            drugtype = drug.drugtype
            status = drug.status
            approval = drug.approval
            if approval==0:
                approval = '-'
            indication = drug.indication
            novelty = drug.novelty


            target_list = drug.target.all()
            targets = []
            for protein in target_list:
                # targets.append(str(protein))
                # jsondata = {'name':drugname, 'target': str(protein), 'approval': approval, 'indication': indication, 'status':status, 'drugtype':drugtype, 'novelty': novelty}
                
                clas = str(protein.family.parent.parent.parent.name)
                family = str(protein.family.parent.name)

                jsondata = {'name':drugname, 'target': str(protein), 'approval': approval, 'class':clas, 'family':family, 'indication': indication, 'status':status, 'drugtype':drugtype, 'novelty': novelty}
                context.append(jsondata)

            # jsondata = {'name':drugname, 'target': ', '.join(set(targets)), 'approval': approval, 'indication': indication, 'status':status, 'drugtype':drugtype, 'novelty': novelty}
            # context.append(jsondata)
        cache.set(name_of_cache, context, 60*60*24*1) # two days timeout on cache

    return render(request, 'drugbrowser.html', {'drugdata':context})

@cache_page(60*5) #  5 min
def drugmapping(request):
    context = dict()

    with open('/protwis/sites/protwis/drugs/flare.json') as data_file:    
        drugdata = json.load(data_file)
    
    context["drugdata"] = drugdata

    return render(request, 'drugmapping.html', {'drugdata':context})