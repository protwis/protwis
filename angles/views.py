from django.conf import settings
from django.shortcuts import render
from django.db.models import Count, Avg, Min, Max
from django.http import HttpResponse, JsonResponse, HttpResponseRedirect
from django.views.decorators.cache import cache_page
from django.views.generic import TemplateView, View

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

        # Grab PDB data
        if len(pdbs)==1:
            query = Angle.objects.filter(structure__pdb_code__index=pdbs[0]).prefetch_related("residue__generic_number")
            # Prep data
            data['data'] = [[q.residue.generic_number.label,q.residue.sequence_number, q.a_angle, q.b_angle, q.hse, q.sasa, q.phi, q.psi, q.theta, q.tau ] for q in query]
        else:
            print(pdbs)
            query = Angle.objects.filter(structure__pdb_code__index__in=pdbs).prefetch_related("residue__generic_number") \
                    .values("residue__generic_number__label") \
                    .annotate(min_aangle = Min('a_angle'), avg_aangle=Avg('a_angle'), max_aangle = Max('a_angle'), \
                        min_bangle = Min('b_angle'), avg_bangle=Avg('b_angle'), max_bangle = Max('b_angle'), \
                        min_hse = Min('hse'), avg_hse=Avg('hse'), max_hse = Max('hse'), \
                        min_sasa = Min('sasa'), avg_sasa=Avg('sasa'), max_sasa = Max('sasa'), \
                        min_phi = Min('phi'), avg_phi=Avg('phi'), max_phi = Max('phi'), \
                        min_psi = Min('psi'), avg_psi=Avg('psi'), max_psi = Max('psi'), \
                        min_theta = Min('theta'), avg_theta=Avg('theta'), max_theta = Max('theta'), \
                        min_tau = Min('tau'), avg_tau=Avg('tau'), max_tau = Max('tau'))
            # Prep data
            data['data'] = [[q["residue__generic_number__label"], " ", \
                            [q["min_aangle"], q["avg_aangle"], q["max_aangle"]], \
                            [q["min_bangle"], q["avg_bangle"], q["max_bangle"]], \
                            [q["min_hse"], q["avg_hse"], q["max_hse"]], \
                            [q["min_sasa"], q["avg_sasa"], q["max_sasa"]], \
                            [q["min_phi"], q["avg_phi"], q["max_phi"]], \
                            [q["min_psi"], q["avg_psi"], q["max_psi"]], \
                            [q["min_theta"], q["avg_theta"], q["max_theta"]], \
                            [q["min_tau"], q["avg_tau"], q["max_tau"]], \
                            ] for q in query]

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
