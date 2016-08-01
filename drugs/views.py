from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.conf import settings

from drugs.models import Drugs
from protein.models import Protein

import re
import json

# Create your views here.

def drugstatistics(request):
    # Get drugdata from here somehow
    
    context = [{'name': 'Blah', 'target': 'adrb1_human'},{'name': 'Artemisin', 'target': 'adrb2_human'},{'name': '3', 'target': 'glp1_human'},{'name': '4', 'target': 'mtlnr_human'}]
    return render(request, 'drugstatistics.html', {'drugdata':context})

def drugbrowser(request):
    # Get drugdata from here somehow
    context = list()

    drugs = Drugs.objects.all()

    for drug in drugs:
    	target_list = drug.target.all()
    	for protein in target_list:
	    	drugname =  str(protein.drugs_set.all()[0])
	    	drugtype = str(protein.drugs_set.all()[0].drugtype)
	    	status = str(protein.drugs_set.all()[0].status)
	    	approval = protein.drugs_set.all()[0].approval
	    	if approval==0:
	    		approval = '-'
	    	indication = str(protein.drugs_set.all()[0].indication)
	    	novelty = str(protein.drugs_set.all()[0].novelty)

    		jsondata = {'name':drugname, 'target': str(protein), 'approval': approval, 'target': str(protein), 'indication': indication, 'status':status, 'drugtype':drugtype, 'novelty': novelty}
    		context.append(jsondata)

    return render(request, 'drugbrowser.html', {'drugdata':context})

def drugmapping(request):
    context = dict()

    with open('/protwis/sites/protwis/drugs/flare.json') as data_file:    
        drugdata = json.load(data_file)
    
    context["drugdata"] = drugdata

    return render(request, 'drugmapping.html', {'drugdata':context})