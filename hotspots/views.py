from django.conf import settings
from django.shortcuts import render
from django.db.models import Count, Avg, Min, Max, Q
from django.http import HttpResponse, JsonResponse, HttpResponseRedirect
from django.views.decorators.cache import cache_page
from django.views.generic import TemplateView, View

import contactnetwork.pdb as pdb
from structure.models import Structure
from residue.models import Residue

import copy
import io
from collections import OrderedDict
import numpy as np

def hotspotsView(request):
    """
    Show hotspots viewer page
    """
    return render(request, 'hotspots/hotspotsView.html')


def getHotspots(request):
    data = {'error': 0}

    return JsonResponse(data)
