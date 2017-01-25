from django.shortcuts import render
from django.db.models import Q

import json

from contactnetwork.models import *
from structure.models import Structure

from django.http import JsonResponse


def Interactions(request):
    """
    Show interaction heatmap
    """
    return render(request, 'contactnetwork/interactions.html')


def InteractionData(request):
    # PDB files
    try:
        pdbs = request.GET.getlist('pdbs')
    except IndexError:
        pdbs = []

    # TM filters
    try:
        tms = request.GET.getlist('tms')
    except IndexError:
        tms = []


    # Use generic numbers? Defaults to True.
    generic = True
    try:
        generic_string = request.GET.get('generic')
        if generic_string in ['false', 'False', 'FALSE', '0']:
            generic = False
    except IndexError:
        pass

    tm_filter_res1 = Q()
    tm_filter_res2 = Q()

    if generic:
        for tm in tms:
            tm_filter_res1 = tm_filter_res1 | Q(interacting_pair__res1__generic_number__label__startswith=str(tm) + 'x')

        for tm in tms:
            tm_filter_res2 = tm_filter_res2 | Q(interacting_pair__res1__generic_number__label__startswith=str(tm) + 'x')

    tm_filter = tm_filter_res1 & tm_filter_res2

    # Get the relevant interactions
    interactions = Interaction.objects.filter(
        interacting_pair__referenced_structure__protein_conformation__protein__entry_name__in=pdbs
    ).values(
        'interacting_pair__referenced_structure__protein_conformation__protein__entry_name',
        'interacting_pair__res1__sequence_number',
        'interacting_pair__res1__generic_number__label',
        'interacting_pair__res2__sequence_number',
        'interacting_pair__res2__generic_number__label',
        'polymorphic_ctype__model',
    ).filter(
        tm_filter
    )

    # Initialize response dictionary
    data = {}
    data['interactions'] = {}
    data['pdbs'] = []
    data['generic'] = generic

    if generic:
        data['tms'] = []

    for i in interactions:
        pdb_name = i['interacting_pair__referenced_structure__protein_conformation__protein__entry_name']
        res1_seq = i['interacting_pair__res1__sequence_number']
        res2_seq = i['interacting_pair__res2__sequence_number']
        res1_gen = i['interacting_pair__res1__generic_number__label']
        res2_gen = i['interacting_pair__res2__generic_number__label']
        model = i['polymorphic_ctype__model']

        if generic and (not res1_gen or not res2_gen):
            continue

        # List PDB files that were found in dataset.
        if pdb_name not in data['pdbs']:
            data['pdbs'].append(pdb_name)

        # List which TMs are available.
        if generic:
            tm_res1 = res1_gen.split('x')[0]
            tm_res2 = res2_gen.split('x')[0]

            if tm_res1 not in data['tms']:
                data['tms'].append(tm_res1)

            if tm_res2 not in data['tms']:
                data['tms'].append(tm_res2)

        # Numbering convention
        res1 = res1_seq
        res2 = res2_seq

        if generic:
            res1 = res1_gen
            res2 = res2_gen

        if res1 < res2:
            coord = str(res1) + ',' + str(res2)
        else:
            coord = str(res2) + ',' + str(res1)

        if coord not in data['interactions']:
            data['interactions'][coord] = {}

        if pdb_name not in data['interactions'][coord]:
            data['interactions'][coord][pdb_name] = []

        data['interactions'][coord][pdb_name].append(model)

    return JsonResponse(data)

