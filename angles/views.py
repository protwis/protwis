from django.shortcuts import render
from django.conf import settings
from django.views.generic import TemplateView, View
from django.http import HttpResponse, JsonResponse, HttpResponseRedirect
from django.views.decorators.cache import cache_page

import contactnetwork.pdb as pdb
from structure.models import Structure
from residue.models import Residue
from angles.models import ResidueAngle as Angle

import Bio.PDB
import copy
import io
from collections import OrderedDict
import numpy as np
from sklearn.decomposition import PCA
from numpy.core.umath_tests import inner1d
import freesasa
import scipy.stats as stats

def angleAnalysis(request):
    """
    Show angle analysis site
    """
    return render(request, 'angles/angleanalysis.html')

def get_angles(request):
    data = {'error': 0}

    # Request selection
    try:
        pdbs = request.GET.getlist('pdbs[]')
        pdbs = [pdb.upper() for pdb in pdbs]

        # Grab PDB data for (first?) PDB
        query = Angle.objects.filter(structure__pdb_code__index=pdbs[0]).prefetch_related("residue__generic_number")

        # Return prep data
        data['data'] = [[q.residue.generic_number.label,q.residue.sequence_number, q.a_angle, q.b_angle, q.hse, q.sasa, q.phi, q.psi, q.theta, q.tau ] for q in query]
    except IndexError:
        data['error'] = 1
        data['errorMessage'] = "No PDB(s) selection provided"

    return JsonResponse(data)

def ServePDB(request, pdbname):
    structure=Structure.objects.filter(pdb_code__index=pdbname.upper())
    if structure.exists():
        structure=structure.get()
    else:
        quit()

    if structure.pdb_data is None:
        quit()

    only_gns = list(structure.protein_conformation.residue_set.exclude(generic_number=None).values_list('protein_segment__slug','sequence_number','generic_number__label').all())
    only_gn = []
    gn_map = []
    segments = {}
    for gn in only_gns:
        only_gn.append(gn[1])
        gn_map.append(gn[2])
        if gn[0] not in segments:
            segments[gn[0]] = []
        segments[gn[0]].append(gn[1])
    data = {}
    data['pdb'] = structure.pdb_data.pdb
    data['only_gn'] = only_gn
    data['gn_map'] = gn_map
    data['segments'] = segments
    data['chain'] = structure.preferred_chain

    return JsonResponse(data)
