from django.shortcuts import render
import json

from contactnetwork.models import InteractingResiduePair
from structure.models import Structure

from django.http import JsonResponse

def HeatMap(request):
    """
    Show interaction heatmap
    """
    return render(request, 'contactnetwork/heatmap.html')

def HeatMapDataJson(request):
    pdbs = request.GET.getlist('pdbs')

    pairs = InteractingResiduePair.objects.filter(referenced_structure__protein_conformation__protein__entry_name__in=pdbs)

    return JsonResponse(pairs, safe=False)

