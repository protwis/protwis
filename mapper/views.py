from django.shortcuts import get_object_or_404, render, redirect
from django.views import generic
from django.http import JsonResponse, HttpResponse
from django.db.models import Q, F, Func, Value, Prefetch
from django.core.cache import cache
from django.views.decorators.cache import cache_page
from django.urls import reverse
from django.shortcuts import render
from django.conf import settings
from django.views.decorators.cache import cache_page
from django.views.generic import TemplateView


from protein.models import Protein, ProteinConformation, ProteinAlias, ProteinFamily, Gene, ProteinSegment
from residue.models import Residue
from structure.models import Structure, StructureModel, StructureExtraProteins
from interaction.models import ResidueFragmentInteraction,StructureLigandInteraction
from common.selection import Selection
from common.views import AbsBrowseSelection
from ligand.models import Ligand, LigandID

import json
from copy import deepcopy
from collections import OrderedDict


class LandingPage(TemplateView):
    template_name = 'mapper/data_mapper_landing.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        return context
        
#
# def LandingPage(request):
#     return render(request, 'mapper/data_mapper_landing.html')

def UploadFile(request):
    context = {
        'url' : 'mapper/data_mapper_landing.html',
        # 'view_name': request.resolver_match.view_name,
        'message_success' : 'File uploaded successfully'
    }
    return render(request, 'mapper/data_mapper_landing.html', context)