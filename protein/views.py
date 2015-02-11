from django.shortcuts import get_object_or_404, render
from django.views import generic
from django.http import HttpResponse
from django.db.models import Q

from protein.models import Protein
from protein.models import ProteinAlias
from protein.models import Gene
from protein.models import ProteinFamily
from common.views import AbsTargetSelection
from common.views import AbsSegmentSelection

import json


class IndexView(generic.ListView):
    model = Protein

    def get_queryset(self):
        return Protein.objects.all()


class DetailView(generic.DetailView):
    model = Protein
    slug_field = 'entry_name'


def SelectionAutocomplete(request):
    if request.is_ajax():
        q = request.GET.get('term')
        type_of_selection = request.GET.get('type_of_selection')
        results = []

        if type_of_selection == 'targets':
            # find protein families
            pfs = ProteinFamily.objects.filter(name__icontains=q).exclude(slug='000')[:10]
            for pf in pfs:
                pf_json = {}
                pf_json['id'] = pf.id
                pf_json['label'] = pf.name
                pf_json['type'] = 'family'
                pf_json['category'] = 'Target families'
                results.append(pf_json)
        
        # find proteins
        ps = Protein.objects.filter(Q(name__icontains=q) | Q(entry_name__icontains=q) | Q(family__name__icontains=q))[:10]
        for p in ps:
            p_json = {}
            p_json['id'] = p.id
            p_json['label'] = p.name + " [" + p.species.common_name + "]"
            p_json['type'] = 'protein'
            p_json['category'] = 'Targets'
            results.append(p_json)

        # find protein aliases
        pas = ProteinAlias.objects.select_related('protein').filter(name__icontains=q)[:10]
        for pa in pas:
            pa_json = {}
            pa_json['id'] = pa.protein.id
            pa_json['label'] = pa.protein.name  + " [" + pa.protein.species.common_name + "]"
            pa_json['type'] = 'protein'
            pa_json['category'] = 'Targets'
            if pa_json not in results:
                results.append(pa_json)
        
        data = json.dumps(results)
    else:
        data = 'fail'
    mimetype = 'application/json'
    return HttpResponse(data, mimetype)