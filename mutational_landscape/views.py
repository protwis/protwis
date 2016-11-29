from django.shortcuts import get_object_or_404, render
from django.http import HttpResponse
from django.core.cache import cache
from django.views.decorators.cache import cache_page

from protein.models import Protein, ProteinConformation, ProteinAlias, ProteinFamily, Gene, ProteinGProtein, ProteinGProteinPair
from residue.models import Residue, ResiduePositionSet
from mutational_landscape.models import NaturalMutations, CancerMutations, DiseaseMutations

from common import definitions
from collections import OrderedDict
from common.views import AbsTargetSelection

import json
# Create your views here.

def ajaxNaturalMutation(request, slug, **response_kwargs):

    name_of_cache = 'ajaxNaturalMutation_'+slug

    jsondata = cache.get(name_of_cache)

    if jsondata == None:
        jsondata = {}

        NMs = NaturalMutations.objects.filter(protein__entry_name=slug).prefetch_related('residue')

        for NM in NMs:

            SN = NM.residue.sequence_number
            jsondata[SN] = [NM.amino_acid, NM.allele_frequency, NM.allele_count, NM.allele_number, NM.number_homozygotes]



        jsondata = json.dumps(jsondata)
        response_kwargs['content_type'] = 'application/json'

        cache.set(name_of_cache, jsondata, 60*60*24*2) #two days timeout on cache

    return HttpResponse(jsondata, **response_kwargs)

def ajaxCancerMutation(request, slug, **response_kwargs):

    name_of_cache = 'ajaxCancerMutation_'+slug

    jsondata = cache.get(name_of_cache)

    if jsondata == None:
        jsondata = {}

        CMs = CancerMutations.objects.filter(protein__entry_name=slug).prefetch_related('residue')

        for CM in CMs:
            SN = CM.residue.sequence_number
            jsondata[SN] = [CM.amino_acid]

        jsondata = json.dumps(jsondata)
        response_kwargs['content_type'] = 'application/json'

        cache.set(name_of_cache, jsondata, 60*60*24*2) #two days timeout on cache

    return HttpResponse(jsondata, **response_kwargs)


def ajaxDiseaseMutation(request, slug, **response_kwargs):

    name_of_cache = 'ajaxDiseaseMutation_'+slug

    jsondata = cache.get(name_of_cache)

    if jsondata == None:
        jsondata = {}

        DMs = DiseaseMutations.objects.filter(protein__entry_name=slug).prefetch_related('residue')

        for DM in DMs:
            SN = NM.residue.sequence_number
            jsondata[SN] = [DM.amino_acid]

        jsondata = json.dumps(jsondata)
        response_kwargs['content_type'] = 'application/json'

        cache.set(name_of_cache, jsondata, 60*60*24*2) #two days timeout on cache

    return HttpResponse(jsondata, **response_kwargs)