from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.conf import settings
from django.db.models import Count

from drugs.models import Drugs
from protein.models import Protein

import re
import json
import numpy as np
import colorsys
import codecs

def get_spaced_colors(n):
    max_value = 16581375 #255**3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]

    return ['#'+i for i in colors] # HEX colors
    # return [(int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors] # RGB colors

def drugstatistics(request):
    # Get drugdata from here somehow
    drugtargets = Drugs.objects.all().filter(status='approved')

    # Query all proteins 

    drugtypes_raw = Drugs.objects.values('drugtype').filter(status='approved').annotate(value=Count('drugtype')).order_by('value')

    list_of_hec_colors = get_spaced_colors(len(drugtypes_raw))
    drugtypes = []
    for i, drugtype in enumerate(drugtypes_raw):
        drugtype['label'] = drugtype['drugtype']
        drugtype['color'] = str(list_of_hec_colors[i])
        del drugtype['drugtype']
        drugtypes.append(drugtype)

    drugindications_raw = Drugs.objects.values('indication').filter(status='approved').annotate(value=Count('indication')).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugindications_raw))
    drugindications = []
    for i, drugindication in enumerate(drugindications_raw):
        drugindication['label'] = drugindication['indication']
        drugindication['color'] = str(list_of_hec_colors[i])
        del drugindication['indication']
        drugindications.append(drugindication)
        if i >= 20:
        	break

    return render(request, 'drugstatistics.html', {'drugtypes':drugtypes, 'drugindications':drugindications})

def drugbrowser(request):
    # Get drugdata from here somehow
    context = list()

    drugs = Drugs.objects.all().prefetch_related()

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

    return render(request, 'drugbrowser.html', {'drugdata':context})

def drugmapping(request):
    context = dict()

    with open('/protwis/sites/protwis/drugs/flare.json') as data_file:    
        drugdata = json.load(data_file)
    
    context["drugdata"] = drugdata

    return render(request, 'drugmapping.html', {'drugdata':context})