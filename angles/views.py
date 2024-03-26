from django.conf import settings
from django.contrib.postgres.aggregates import ArrayAgg
from django.shortcuts import render
from django.db.models import Count, Avg, Min, Max, Q
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
import math
import cmath
from collections import OrderedDict
import numpy as np
# from sklearn.decomposition import PCA
from numpy.core.umath_tests import inner1d
import freesasa
import scipy.stats as stats

def angleAnalysis(request):
    """
    Show angle analysis page
    """
    return render(request, 'angles/angleanalysis.html')


def angleAnalyses(request):
    """
    Show angle analyses page
    """
    return render(request, 'angles/angleanalyses.html')

def structureCheck(request):
    """
    Show structure annotation check page
    """
    return render(request, 'angles/structurecheck.html')

def get_angles(request):
    data = {'error': 0}

    # angle names for custom averaging
    angles = ['avg_aangle', 'avg_bangle', 'avg_outer', 'avg_phi', 'avg_psi', 'avg_theta', 'avg_tau']

    # Request selection
    try:
    #if True:
        pdbs = request.GET.getlist('pdbs[]')
        pdbs = set([pdb.upper() for pdb in pdbs])
        print(pdbs)

        pdbs2 = request.GET.getlist('pdbs2[]')
        pdbs2 = set([pdb.upper() for pdb in pdbs2])
        print(pdbs2)

        # Grab PDB data
        if len(pdbs)==1 and len(pdbs2)==0:
            pdbs = list(pdbs)
            query = Angle.objects.filter(structure__pdb_code__index=pdbs[0]).prefetch_related("residue__generic_number").order_by('residue__display_generic_number__label')

            # Prep data
            #data['data'] = [[q.residue.generic_number.label,q.residue.sequence_number, q.a_angle, q.b_angle, q.outer_angle, q.hse, q.sasa, q.rsa, q.phi, q.psi, q.theta, q.tau, q.core_distance, q.ss_dssp, q.ss_stride ] for q in query ]
            data['data'] = []
            for q in query:
                if q.residue.display_generic_number != None:
                    data['data'].append([q.residue.short_display_generic_number(),q.residue.sequence_number, q.a_angle, q.b_angle, q.outer_angle, q.hse, q.sasa, q.rsa, q.phi, q.psi, q.theta, q.tau, q.core_distance, q.ss_dssp, q.ss_stride ])
                else:
                    data['data'].append(["-",q.residue.sequence_number, q.a_angle, q.b_angle, q.outer_angle, q.hse, q.sasa, q.rsa, q.phi, q.psi, q.theta, q.tau, q.core_distance, q.ss_dssp, q.ss_stride ])
            data['headers'] = [{"title" : "Value"}]
        else: # always a grouping or a comparison
            query = Angle.objects.filter(structure__pdb_code__index__in=pdbs).prefetch_related("residue__generic_number") \
                    .values("residue__generic_number__label") \
                    .order_by('residue__generic_number__label') \
                    .annotate(min_aangle = Min('a_angle'), avg_aangle=ArrayAgg('a_angle'), max_aangle = Max('a_angle'), \
                        min_bangle = Min('b_angle'), avg_bangle=ArrayAgg('b_angle'), max_bangle = Max('b_angle'), \
                        min_outer = Min('outer_angle'), avg_outer=ArrayAgg('outer_angle'), max_outer = Max('outer_angle'), \
                        min_hse = Min('hse'), avg_hse=Avg('hse'), max_hse = Max('hse'), \
                        min_sasa = Min('sasa'), avg_sasa=Avg('sasa'), max_sasa = Max('sasa'), \
                        min_rsa = Min('rsa'), avg_rsa=Avg('rsa'), max_rsa = Max('rsa'), \
                        min_phi = Min('phi'), avg_phi=ArrayAgg('phi'), max_phi = Max('phi'), \
                        min_psi = Min('psi'), avg_psi=ArrayAgg('psi'), max_psi = Max('psi'), \
                        min_theta = Min('theta'), avg_theta=ArrayAgg('theta'), max_theta = Max('theta'), \
                        min_tau = Min('tau'), avg_tau=ArrayAgg('tau'), max_tau = Max('tau'), \
                        min_distance = Min('core_distance'), avg_distance=Avg('core_distance'), max_distance = Max('core_distance'))

            # Process angle aggregates to angle averages
            for q in query:
                for angle in angles:
                    q[angle] = [ qa for qa in q[angle] if qa != None]
                    if angle in q and len(q[angle]) > 1:
                        # Sensible average for multiple angles (circular statistics: https://rosettacode.org/wiki/Averages/Mean_angle)
                        q[angle] = math.degrees(cmath.phase(sum(cmath.rect(1, math.radians(float(d))) for d in q[angle])/len(q[angle])))
                    elif len(q[angle]) == 1:
                        q[angle] = q[angle][0]

            # Prep data
            data['data'] = [ [q["residue__generic_number__label"], " ", \
                            [q["min_aangle"], q["avg_aangle"], q["max_aangle"]], \
                            [q["min_bangle"], q["avg_bangle"], q["max_bangle"]], \
                            [q["min_outer"], q["avg_outer"], q["max_outer"]], \
                            [q["min_hse"], q["avg_hse"], q["max_hse"]], \
                            [q["min_sasa"], q["avg_sasa"], q["max_sasa"]], \
                            [q["min_rsa"], q["avg_rsa"], q["max_rsa"]], \
                            [q["min_phi"], q["avg_phi"], q["max_phi"]], \
                            [q["min_psi"], q["avg_psi"], q["max_psi"]], \
                            [q["min_theta"], q["avg_theta"], q["max_theta"]], \
                            [q["min_tau"], q["avg_tau"], q["max_tau"]], \
                            [q["min_distance"], q["avg_distance"], q["max_distance"]], \
                            ] for q in query]

            print(data)
            if len(pdbs2)==0:
                data['headers'] = [{"title" : "Group<br/>Min"},{"title" : "Group<br/>Avg"},{"title" : "Group<br/>Max"}]
            else:
                data['headers'] = [{"title" : "Group 1<br/>Min"},{"title" : "Group 1<br/>Avg"},{"title" : "Group 1<br/>Max"}]

        # Select PDBs from same Class + same state
        data['headers2'] = [{"title" : "Group 2<br/>Min"},{"title" : "Group 2<br/>Avg"},{"title" : "Group 2<br/>Max"}]
        if len(pdbs2)==0:
            # select structure(s)
            structures = Structure.objects.filter(pdb_code__index__in=pdbs) \
                        .select_related('protein_conformation__protein__family','protein_conformation__state')

            # select PDBs
            states = set( structure.protein_conformation.state.slug for structure in structures )
            classes = set( structure.protein_conformation.protein.family.slug[:3] for structure in structures )

            query = Q()
            for classStart in classes:
                    query = query | Q(protein_conformation__protein__family__slug__startswith=classStart)
            set2 = Structure.objects.filter(protein_conformation__state__slug__in=states).exclude(structure_type__slug__startswith='af-').filter(query).values_list('pdb_code__index')

            pdbs2 = [ x[0] for x in set2 ]

            data['headers2'] = [{"title" : "Class<br/>Min"},{"title" : "Class<br/>Avg"},{"title" : "Class<br/>Max"}]

        query = Angle.objects.filter(structure__pdb_code__index__in=pdbs2).prefetch_related("residue__generic_number") \
                .values("residue__generic_number__label") \
                .annotate(min_aangle = Min('a_angle'), avg_aangle=ArrayAgg('a_angle'), max_aangle = Max('a_angle'), \
                    min_bangle = Min('b_angle'), avg_bangle=ArrayAgg('b_angle'), max_bangle = Max('b_angle'), \
                    min_outer = Min('outer_angle'), avg_outer=ArrayAgg('outer_angle'), max_outer = Max('outer_angle'), \
                    min_hse = Min('hse'), avg_hse=Avg('hse'), max_hse = Max('hse'), \
                    min_sasa = Min('sasa'), avg_sasa=Avg('sasa'), max_sasa = Max('sasa'), \
                    min_rsa = Min('rsa'), avg_rsa=Avg('rsa'), max_rsa = Max('rsa'), \
                    min_phi = Min('phi'), avg_phi=ArrayAgg('phi'), max_phi = Max('phi'), \
                    min_psi = Min('psi'), avg_psi=ArrayAgg('psi'), max_psi = Max('psi'), \
                    min_theta = Min('theta'), avg_theta=ArrayAgg('theta'), max_theta = Max('theta'), \
                    min_tau = Min('tau'), avg_tau=ArrayAgg('tau'), max_tau = Max('tau'), \
                    min_distance = Min('core_distance'), avg_distance=Avg('core_distance'), max_distance = Max('core_distance'))

        # Process angle aggregates to angle averages
        for q in query:
            for angle in angles:
                q[angle] = [ q for q in q[angle] if q != None]
                if angle in q and len(q[angle]) > 1:
                    # Sensible average for multiple angles (circular statistics: https://rosettacode.org/wiki/Averages/Mean_angle)
                    q[angle] = math.degrees(cmath.phase(sum(cmath.rect(1, math.radians(float(d))) for d in q[angle])/len(q[angle])))
                elif len(q[angle]) == 1:
                    q[angle] = q[angle][0]

        # Prep data
        data['data2'] = { q["residue__generic_number__label"]: [q["residue__generic_number__label"], " ", \
                        [q["min_aangle"], q["avg_aangle"], q["max_aangle"]], \
                        [q["min_bangle"], q["avg_bangle"], q["max_bangle"]], \
                        [q["min_outer"], q["avg_outer"], q["max_outer"]], \
                        [q["min_hse"], q["avg_hse"], q["max_hse"]], \
                        [q["min_sasa"], q["avg_sasa"], q["max_sasa"]], \
                        [q["min_rsa"], q["avg_rsa"], q["max_rsa"]], \
                        [q["min_phi"], q["avg_phi"], q["max_phi"]], \
                        [q["min_psi"], q["avg_psi"], q["max_psi"]], \
                        [q["min_theta"], q["avg_theta"], q["max_theta"]], \
                        [q["min_tau"], q["avg_tau"], q["max_tau"]], \
                        [q["min_distance"], q["avg_distance"], q["max_distance"]], \
                        ] for q in query}



    except IndexError:
    #else:
        data['error'] = 1
        data['errorMessage'] = "No PDB(s) selection provided"

    return JsonResponse(data)

def ServePDB(request, pdbname):
    # query = Angle.objects.filter(residue__protein_segment__slug__in=['TM1','TM2','TM3','TM4','TM5','TM6','TM7','H8']).prefetch_related("residue__generic_number") \
    #         .aggregate(total=Count('ss_stride'), \
    #         total2=Count('ss_dssp'))
    # print(query)
    #
    # query = Angle.objects.filter(residue__protein_segment__slug__in=['TM1','TM2','TM3','TM4','TM5','TM6','TM7','H8']).prefetch_related("residue__generic_number") \
    #         .values("ss_stride") \
    #         .annotate(total=Count('ss_stride')) \
    #         .order_by('ss_stride')
    # print(query)
    #
    # query = Angle.objects.filter(residue__protein_segment__slug__in=['TM1','TM2','TM3','TM4','TM5','TM6','TM7','H8']).prefetch_related("residue__generic_number") \
    #         .values("ss_dssp") \
    #         .annotate(total=Count('ss_dssp')) \
    #         .order_by('ss_dssp')
    # print(query)

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
