from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.conf import settings
from django.db.models import Count
from django.core.cache import cache

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


def drugstatistics(request):

    # ===== drugtargets =====
    drugtargets_raw = Protein.objects.filter(drugs__status='approved').values('entry_name').annotate(value=Count('entry_name')).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugtargets_raw))
    drugtargets = []
    for i, drugtarget in enumerate(drugtargets_raw):
        drugtarget['label'] = drugtarget['entry_name'].replace("_human","").upper()
        # drugtarget['color'] = str(list_of_hec_colors[i])
        del drugtarget['entry_name']
        drugtargets.append(drugtarget)

    # ===== drugfamilies =====
    drugfamilies_raw = Protein.objects.filter().filter(drugs__status='approved').values('family_id__parent__name').annotate(value=Count('family_id__parent__name')).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugfamilies_raw))
    drugfamilies = []
    for i, drugfamily in enumerate(drugfamilies_raw):
        drugfamily['label'] = striphtml(drugfamily['family_id__parent__name']).replace(" receptors","")
        drugfamily['color'] = str(list_of_hec_colors[i])
        del drugfamily['family_id__parent__name']
        drugfamilies.append(drugfamily)

    # ===== drugclas =====
    drugclasses_raw = Protein.objects.filter().filter(drugs__status='approved').values('family_id__parent__parent__parent__name').annotate(value=Count('family_id__parent__parent__parent__name')).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugclasses_raw)+1)
    drugclasses = []
    for i, drugclas in enumerate(drugclasses_raw):
        drugclas['label'] = drugclas['family_id__parent__parent__parent__name']
        drugclas['color'] = str(list_of_hec_colors[i+1])
        del drugclas['family_id__parent__parent__parent__name']
        drugclasses.append(drugclas)

    # ===== drugtypes =====
    drugtypes_raw = Drugs.objects.values('drugtype').filter(status='approved').annotate(value=Count('drugtype')).order_by('value')

    list_of_hec_colors = get_spaced_colors(len(drugtypes_raw))
    drugtypes = []
    for i, drugtype in enumerate(drugtypes_raw):
        drugtype['label'] = drugtype['drugtype']
        drugtype['color'] = str(list_of_hec_colors[i])
        del drugtype['drugtype']
        drugtypes.append(drugtype)

    # ===== drugindications =====
    drugindications_raw = Drugs.objects.values('indication').filter(status='approved').annotate(value=Count('indication')).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugindications_raw))
    drugindications = []
    for i, drugindication in enumerate(drugindications_raw):
        drugindication['label'] = drugindication['indication']
        drugindication['color'] = str(list_of_hec_colors[i])
        del drugindication['indication']
        drugindications.append(drugindication)

    # ===== drugtimes =====
    drugtime_raw = Drugs.objects.values('approval').filter(status='approved').annotate(y=Count('approval')).order_by('approval')
    drugtimes = []
    running_total = 0
    for i, time in enumerate(drugtime_raw):
        if time['approval']!='-':
            time['x'] = int(time['approval'])
            time['y'] = int(time['y']) + running_total
            del time['approval']
            running_total = int(time['y'])

            if time['x'] % 2 == 0:
	            drugtimes.append(time)

    print(drugtimes)
    drugs_over_time = [{"values": drugtimes, "yAxis": "1", "key": "GPCRs"}]

    return render(request, 'drugstatistics.html', {'drugtypes':drugtypes, 'drugindications':drugindications, 'drugtargets':drugtargets, 'drugfamilies':drugfamilies, 'drugclasses':drugclasses, 'drugs_over_time':drugs_over_time})

def drugbrowser(request):
    # Get drugdata from here somehow

    name_of_cache = 'drug_browser'

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
        cache.set(name_of_cache, context, 60*60*24*2) # two days timeout on cache

    return render(request, 'drugbrowser.html', {'drugdata':context})

def drugmapping(request):
    context = dict()

    with open('/protwis/sites/protwis/drugs/flare.json') as data_file:    
        drugdata = json.load(data_file)
    
    context["drugdata"] = drugdata

    return render(request, 'drugmapping.html', {'drugdata':context})