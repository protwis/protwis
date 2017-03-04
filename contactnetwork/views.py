from django.shortcuts import render
from django.db.models import Q

import json
import functools

from contactnetwork.models import *
from structure.models import Structure

from django.http import JsonResponse


def Interactions(request):
    """
    Show interaction heatmap
    """
    return render(request, 'contactnetwork/interactions.html')


def InteractionData(request):

    def gpcrdb_number_comparator(e1, e2):
            t1 = e1.split('x')
            t2 = e2.split('x')

            if e1 == e2:
                return 0

            if t1[0] == t2[0]:
                if t1[1] < t2[1]:
                    return -1
                else:
                    return 1

            if t1[0] < t2[0]:
                return -1
            else:
                return 1

    # PDB files
    try:
        pdbs = request.GET.getlist('pdbs[]')
    except IndexError:
        pdbs = []

    pdbs = [pdb.lower() for pdb in pdbs]

    # Segment filters
    try:
        segments = request.GET.getlist('segments[]')
    except IndexError:
        segments = []

    # Interaction types
    try:
        i_types = request.GET.getlist('interaction_types[]')
    except IndexError:
        i_types = []

    # Use generic numbers? Defaults to True.
    generic = True
    try:
        generic_string = request.GET.get('generic')
        if generic_string in ['false', 'False', 'FALSE', '0']:
            generic = False
    except IndexError:
        pass

    segment_filter_res1 = Q()
    segment_filter_res2 = Q()

    if segments:
        segment_filter_res1 |= Q(interacting_pair__res1__protein_segment__slug__in=segments)
        segment_filter_res2 |= Q(interacting_pair__res2__protein_segment__slug__in=segments)

    i_types_filter = Q()

    if i_types:
        i_types_filter |= Q(polymorphic_ctype__model__in=i_types)

    # Get the relevant interactions
    interactions = Interaction.objects.filter(
        interacting_pair__referenced_structure__protein_conformation__protein__entry_name__in=pdbs
    ).values(
        'interacting_pair__referenced_structure__protein_conformation__protein__entry_name',
        'interacting_pair__res1__sequence_number',
        'interacting_pair__res1__generic_number__label',
        'interacting_pair__res1__protein_segment__slug',
        'interacting_pair__res2__sequence_number',
        'interacting_pair__res2__generic_number__label',
        'interacting_pair__res2__protein_segment__slug',
        'polymorphic_ctype__model',
    ).filter(
        segment_filter_res1 & segment_filter_res2 & i_types_filter
    )

    # Initialize response dictionary
    data = {}
    data['interactions'] = {}
    data['pdbs'] = set()
    data['generic'] = generic
    data['segments'] = set()
    data['segment_map'] = {}

    # Dict to keep track of which residue numbers are in use
    number_dict = set()

    itypes = set()

    for i in interactions:
        pdb_name = i['interacting_pair__referenced_structure__protein_conformation__protein__entry_name']
        res1_seq = i['interacting_pair__res1__sequence_number']
        res2_seq = i['interacting_pair__res2__sequence_number']
        res1_gen = i['interacting_pair__res1__generic_number__label']
        res2_gen = i['interacting_pair__res2__generic_number__label']
        res1_seg = i['interacting_pair__res1__protein_segment__slug']
        res2_seg = i['interacting_pair__res2__protein_segment__slug']
        model = i['polymorphic_ctype__model']

        itypes |= {model}

        if generic and (not res1_gen or not res2_gen):
            continue

        # List PDB files that were found in dataset.
        data['pdbs'] |= {pdb_name}

        # Numbering convention
        res1 = res1_seq
        res2 = res2_seq

        if generic:
            res1 = res1_gen
            res2 = res2_gen

        # List which segments are available.
        data['segment_map'][res1] = res1_seg
        data['segment_map'][res2] = res2_seg
        data['segments'] |= {res1_seg} | {res2_seg}

        number_dict |= {res1, res2}

        if res1 < res2:
            coord = str(res1) + ',' + str(res2)
        else:
            coord = str(res2) + ',' + str(res1)

        if coord not in data['interactions']:
            data['interactions'][coord] = {}

        if pdb_name not in data['interactions'][coord]:
            data['interactions'][coord][pdb_name] = []

        data['interactions'][coord][pdb_name].append(model)

    if (generic):
        data['sequence_numbers'] = sorted(number_dict, key=functools.cmp_to_key(gpcrdb_number_comparator))
    else:
        data['sequence_numbers'] = sorted(number_dict)

    data['segments'] = list(data['segments'])
    data['pdbs'] = list(data['pdbs'])

    print(itypes)

    return JsonResponse(data)

