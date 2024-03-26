from django.shortcuts import render
from django.db.models import Q, F, Prefetch, Avg, StdDev, IntegerField, Sum, Case, When, Min, Max
from django.db.models.functions import Concat
from django.views.decorators.cache import cache_page
from django.views.decorators.csrf import csrf_exempt
from django.core.cache import cache
from django.db import connection

from collections import defaultdict
from django.conf import settings

import json
import functools
import hashlib
import copy

from contactnetwork.models import *
from contactnetwork.distances import *
from contactnetwork.functions import *
from structure.models import Structure, StructureVectors, StructureExtraProteins
from structure.templatetags.structure_extras import *
from construct.models import Construct
from protein.models import Protein, ProteinSegment, ProteinCouplings, ProteinConformation
from residue.models import Residue, ResidueGenericNumber
from signprot.models import SignprotComplex
from interaction.models import StructureLigandInteraction, ResidueFragmentInteraction
from angles.models import ResidueAngle, get_angle_averages, get_all_angles
from mutation.models import MutationExperiment

Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')

from django.http import JsonResponse, HttpResponse
from collections import OrderedDict

import math, statistics
import cmath
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd
import time
import hashlib
import operator


def get_hash(data):
    # create unique hash key for alignment combo
    hash_key = ""
    for d in data:
        if isinstance(d, list):
            d_sorted = sorted(set(d))
            hash_key += "|" + "-".join(d_sorted)
        else:
            hash_key += "|" + str(d)

    return hashlib.md5(hash_key.encode('utf-8')).hexdigest()

def Clustering(request):
    """
    Show clustering page
    """
    return render(request, 'contactnetwork/clustering.html')

def Interactions(request):
    """
    Show interaction heatmap
    """

    template_data = {}
    # check for preselections in POST data
    if request.POST and (request.POST.get("pdbs1") != None or request.POST.get("pdbs2") != None):
        # handle post
        pdbs1 = request.POST.get("pdbs1")
        pdbs2 = request.POST.get("pdbs2")
        if pdbs1 == "":
            pdbs1 = None
        if pdbs2 == "":
            pdbs2 = None

        # create switch
        if pdbs1 != None and pdbs2 != None:
            template_data["pdbs1"] = '["' + '", "'.join(pdbs1.split("\r\n")) + '"]'
            template_data["pdbs2"] = '["' + '", "'.join(pdbs2.split("\r\n")) + '"]'
        else:
            if pdbs1 == None:
                template_data["pdbs"] = '["' + '", "'.join(pdbs2.split("\r\n")) + '"]'
            else:
                template_data["pdbs"] = '["' + '", "'.join(pdbs1.split("\r\n")) + '"]'

    return render(request, 'contactnetwork/interactions.html', template_data)

def InteractionBrowser(request):
    """
    Show interaction heatmap
    """

    template_data = {}
    # check for preselections in POST data
    if request.POST and (request.POST.get("pdbs1") != None or request.POST.get("pdbs2") != None):
        # handle post
        pdbs1 = request.POST.get("pdbs1")
        pdbs2 = request.POST.get("pdbs2")
        if pdbs1 == "":
            pdbs1 = None
        if pdbs2 == "":
            pdbs2 = None

        # create switch
        if pdbs1 != None and pdbs2 != None:
            template_data["pdbs1"] = '["' + '", "'.join(pdbs1.split("\r\n")) + '"]'
            template_data["pdbs2"] = '["' + '", "'.join(pdbs2.split("\r\n")) + '"]'
        else:
            if pdbs1 == None:
                template_data["pdbs"] = '["' + '", "'.join(pdbs2.split("\r\n")) + '"]'
            else:
                template_data["pdbs"] = '["' + '", "'.join(pdbs1.split("\r\n")) + '"]'

    return render(request, 'contactnetwork/browser.html', template_data)

def ShowDistances(request):
    """
    Show distances heatmap
    """

    template_data = {}

    if request.POST and (request.POST.get("pdbs1") != None or request.POST.get("pdbs2") != None):
        # check for preselections in POST data
        pdbs1 = request.POST.get("pdbs1")
        pdbs2 = request.POST.get("pdbs2")
        if pdbs1 == "":
            pdbs1 = None
        if pdbs2 == "":
            pdbs2 = None

        # create switch
        if pdbs1 != None and pdbs2 != None:
            template_data["pdbs1"] = '["' + '", "'.join(pdbs1.split("\r\n")) + '"]'
            template_data["pdbs2"] = '["' + '", "'.join(pdbs2.split("\r\n")) + '"]'
        else:
            if pdbs1 == None:
                template_data["pdbs"] = '["' + '", "'.join(pdbs2.split("\r\n")) + '"]'
            else:
                template_data["pdbs"] = '["' + '", "'.join(pdbs1.split("\r\n")) + '"]'

    return render(request, 'contactnetwork/distances.html', template_data)

def PdbTreeData(request):
    data = Structure.objects.exclude(structure_type__slug__startswith='af-').values(
        'representative',
        'pdb_code__index',
        'protein_conformation__protein__parent__family__parent__parent__parent__name',
        'protein_conformation__protein__parent__family__parent__parent__name',
        'protein_conformation__protein__parent__family__parent__name',
        'protein_conformation__protein__parent__family__name',
        'protein_conformation__protein__parent__family__parent__parent__parent__slug',
        'protein_conformation__protein__parent__family__parent__parent__slug',
        'protein_conformation__protein__parent__family__parent__slug',
        'protein_conformation__protein__parent__family__slug',
        'protein_conformation__state__slug'
        )

    # TODO: Use ordereddict
    l = lambda:defaultdict(l)
    data_dict = l()

    for d in data:
        pdb = d['pdb_code__index']
        rep = d['representative']
        l3 = d['protein_conformation__protein__parent__family__name']
        l2 = d['protein_conformation__protein__parent__family__parent__name']
        l1 = d['protein_conformation__protein__parent__family__parent__parent__name']
        l0 = d['protein_conformation__protein__parent__family__parent__parent__parent__name']
        s3 = d['protein_conformation__protein__parent__family__slug']
        s2 = d['protein_conformation__protein__parent__family__parent__slug']
        s1 = d['protein_conformation__protein__parent__family__parent__parent__slug']
        s0 = d['protein_conformation__protein__parent__family__parent__parent__parent__slug']
        state = d['protein_conformation__state__slug']

        if rep:
            rep = 'R'
        else:
            rep = 'N'

        if state == 'active':
            state = 'A'
        else:
            state = 'I'

        if not data_dict[s0 + ',' + l0][s1 + ',' + l1][s2 + ',' + l2][s3 + ',' + l3]:
            data_dict[s0 + ',' + l0][s1 + ',' + l1][s2 + ',' + l2][s3 + ',' + l3] = []

        data_dict[s0 + ',' + l0][s1 + ',' + l1][s2 + ',' + l2][s3 + ',' + l3].append(pdb + ' (' + state + ')'  + '(' + rep + ')')

    return JsonResponse(data_dict)

@cache_page(60*60*24*7)
def PdbTableData(request):
    exclude_non_interacting = True if request.GET.get('exclude_non_interacting') == 'true' else False
    effector = request.GET.get('effector') if request.GET.get('effector') != 'false' else False
    # interaction_protein_class = request.GET.get('interaction_protein_class')
    #constructs = Construct.objects.defer('schematics','snakecache').all().prefetch_related('crystallization__crystal_method')
    #methods = {}
    #for c in constructs:
        # print(c.name)
    #    if c.crystallization and c.crystallization.crystal_method:
    #        method = c.crystallization.crystal_method.name
    #    else:
    #        method = "N/A"
    #    methods[c.name] = method
    
    # data = Structure.objects.all().exclude(structure_type__slug__startswith='af-').prefetch_related(
    #             "pdb_code",
    #             "state",
    #             "stabilizing_agents",
    #             "structureligandinteraction_set__ligand__ligand_type",
    #             "structureligandinteraction_set__ligand_role",
    #             "structure_type",
    #             "protein_conformation__protein__parent__parent__parent",
    #             "protein_conformation__protein__parent__family__parent",
    #             "protein_conformation__protein__parent__family__parent__parent__parent",
    #             "protein_conformation__protein__species",Prefetch("ligands", queryset=StructureLigandInteraction.objects.filter(
    #             annotated=True).exclude(structure__structure_type__slug__startswith='af-').prefetch_related('ligand__ligand_type', 'ligand_role')))
    # get best signalprotein/species/receptor
    # if effector is defined (as one letter), filter by that
    # 'G alpha' = G proteins (all G protein classes starts with G)
    # 'A' = Arrestin
    if effector:
        data = Structure.objects.all().prefetch_related(
                "pdb_code",
                "state",
                "stabilizing_agents",
                "structureligandinteraction_set__ligand__ligand_type",
                "structureligandinteraction_set__ligand_role",
                "structure_type",
                "protein_conformation__protein__parent__parent__parent",
                "protein_conformation__protein__parent__family__parent",
                "protein_conformation__protein__parent__family__parent__parent__parent",
                "protein_conformation__protein__parent",
                "protein_conformation__protein__parent__parent",
                "protein_conformation__protein__family__parent",
                "protein_conformation__protein__family__parent__parent__parent",
                "protein_conformation__protein__species",Prefetch("ligands", queryset=StructureLigandInteraction.objects.filter(
                annotated=True).prefetch_related('ligand', 'ligand__ligand_type', 'ligand_role')))
        data = data.filter(extra_proteins__category__startswith=effector).prefetch_related(
        'extra_proteins__protein_conformation','extra_proteins__wt_protein').order_by(
        'extra_proteins__protein_conformation__protein__parent','state').annotate(
        res_count = Sum(Case(When(extra_proteins__structure__protein_conformation__residue__generic_number=None, then=0), default=1, output_field=IntegerField())))
        signal_ps = StructureExtraProteins.objects.filter(category__startswith=effector).values(
                            'structure__pdb_code__index','structure__protein_conformation__protein','structure__protein_conformation__protein__parent','display_name',
                            'wt_coverage','wt_protein__family__parent__parent__name','wt_protein__family__parent__name','category','note'
                            ).order_by().annotate(coverage = Max('wt_coverage'))
        ep = {s['structure__pdb_code__index']:s for s in signal_ps}
    else:
        data = Structure.objects.all().exclude(structure_type__slug__startswith='af-').prefetch_related(
                "pdb_code",
                "state",
                "stabilizing_agents",
                "structureligandinteraction_set__ligand__ligand_type",
                "structureligandinteraction_set__ligand_role",
                "structure_type",
                "protein_conformation__protein__parent__parent__parent",
                "protein_conformation__protein__parent__family__parent",
                "protein_conformation__protein__parent__family__parent__parent__parent",
                "protein_conformation__protein__parent",
                "protein_conformation__protein__parent__parent",
                "protein_conformation__protein__family__parent",
                "protein_conformation__protein__family__parent__parent__parent",
                "protein_conformation__protein__species",Prefetch("ligands", queryset=StructureLigandInteraction.objects.filter(
                annotated=True).exclude(structure__structure_type__slug__startswith='af-').prefetch_related('ligand', 'ligand__ligand_type', 'ligand_role')))
        data = data.prefetch_related('extra_proteins__protein_conformation','extra_proteins__wt_protein').order_by(
        'extra_proteins__protein_conformation__protein__parent','state').annotate(
        res_count = Sum(Case(When(extra_proteins__protein_conformation__residue__generic_number=None, then=0), default=1, output_field=IntegerField())))
        signal_ps = StructureExtraProteins.objects.all().exclude(category__in=['G beta','G gamma']).values(
                            'structure__pdb_code__index','structure__protein_conformation__protein__parent','display_name',
                            'wt_coverage','wt_protein__family__parent__parent__name','wt_protein__family__parent__name','category','note'
                            ).order_by().annotate(coverage = Max('wt_coverage'))
        ep = {s['structure__pdb_code__index']:s for s in signal_ps}

    if exclude_non_interacting and effector == 'G alpha':
        complex_structure_ids = SignprotComplex.objects.values_list('structure', flat=True)
        data = data.filter(id__in=complex_structure_ids)

    # get a gn residue count for all WT proteins
    proteins_pks = Structure.objects.all().exclude(structure_type__slug__startswith='af-').values_list("protein_conformation__protein__parent__pk", flat=True).distinct()
    proteins_af_pks = Structure.objects.all().filter(structure_type__slug__startswith='af-').values_list("protein_conformation__protein__pk", flat=True).distinct()
    if effector:
        proteins_pks = list(proteins_pks) + list(proteins_af_pks)
    residue_counts = ProteinConformation.objects.filter(protein__pk__in=proteins_pks).values('protein__pk').annotate(res_count = Sum(Case(When(residue__generic_number=None, then=0), default=1, output_field=IntegerField())))
    rcs = {}
    for rc in residue_counts:
        rcs[rc['protein__pk']] = rc['res_count']

    # get minimum resolution for every receptor/state pair
    # resolutions = Structure.objects.all().exclude(structure_type__slug__startswith='af-').values('protein_conformation__protein__parent','state__name').order_by().annotate(res = Min('resolution'))
    if effector:
        resolutions = Structure.objects.all().values('protein_conformation__protein__parent','state__name').order_by().annotate(res = Min('resolution'))
    else:
        resolutions = Structure.objects.all().exclude(structure_type__slug__startswith='af-').values('protein_conformation__protein__parent','state__name').order_by().annotate(res = Min('resolution'))
    best_resolutions = {}
    for r in resolutions:
        key = '{}_{}'.format(r['protein_conformation__protein__parent'], r['state__name'])
        best_resolutions[key] = r['res']

    best_signal_p = {}
    for ps in signal_ps:
        if effector:
            key = '{}_{}_{}'.format(ps['structure__protein_conformation__protein'], ps['structure__protein_conformation__protein__parent'], ps['display_name'])
        else:
            key = '{}_{}'.format(ps['structure__protein_conformation__protein__parent'], ps['display_name'])
        best_signal_p[key] = ps['coverage']

    data_dict = OrderedDict()

    if effector=='G alpha':
        signalling_header = 'G protein'
    elif effector=='A':
        signalling_header = 'Arrestin'
    else:
        signalling_header = 'Signalling protein'
    data_table = "<table id2='structure_selection' border=0 class='structure_selection row-border text-center compact text-nowrap' width='100%'> \
        <thead><tr> \
            <th rowspan=2> <input class ='form-check-input check_all' type='checkbox' value='' onclick='check_all(this);'> </th> \
            <th colspan=5>Receptor</th> \
            <th colspan=3>Species</th> \
            <th colspan=4>Structure</th> \
            <th colspan=3>Receptor state <a href=\"https://docs.gpcrdb.org/structures.html#structure-descriptors\" target=\"_blank\"><span class=\"glyphicon glyphicon-info-sign\"></span></a></th> \
            <th colspan=4>{}</th> \
            <th colspan=2>Auxiliary protein</th> \
            <th colspan=2>Ligand</th> \
        </tr> \
        <tr><th></th> \
            <th></th> \
            <th></th> \
            <th></th> \
            <th>% of Seq</th> \
            <th id=species></th> \
            <th></th> \
            <th>Identity %<br>to Human</th> \
            <th></th> \
            <th></th> \
            <th></th> \
            <th></th> \
            <th></th> \
            <th>Degree active (%)</th> \
            <th>TM6 tilt</th>".format(signalling_header)
#            <th><a href=\"http://docs.gpcrdb.org/structures.html\" target=\"_blank\">Cytosolic</br> opening</a></th>"
#            <th><a href=\"http://docs.gpcrdb.org/structures.html\" target=\"_blank\">7TM Open IC (Å)</a></th> \
#            <th>TM6 tilt (%, inactive: 0-X, intermed: X-Y, active Y-Z)</th> \
    data_table += "<th></th> \
            <th></th> \
            <th>Note</th> \
            <th>% of Seq</th> \
            <th></th> \
            <th></th> \
            <th></th> \
            <th></th> \
        </tr> \
        <tr> \
            <th colspan=6></th> \
            <th colspan=1 id=best_species class='text-center'></th> \
            <th colspan=4></th> \
            <th colspan=1 id=best_res class='text-center'></th> \
            <th colspan=12></th> \
        </tr></thead><tbody>\n"

    identity_lookup = {}
    for s in data:
        pdb_id = s.pdb_code.index
        if pdb_id in data_dict:
            continue

        r = {}
        #Setting a short path that addresses af models and regular structures
        #so we don't have to add infinite try/excepts
        if not s.protein_conformation.protein.parent:
            shorted = s.protein_conformation.protein
        else:
            shorted = s.protein_conformation.protein.parent

        r['protein'] = shorted.entry_short()
        r['protein_long'] = shorted.short()
        r['protein_family'] = shorted.family.parent.short()
        r['class'] = shorted.family.parent.parent.parent.shorter()
        r['species'] = s.protein_conformation.protein.species.common_name
        # # r['date'] = s.publication_date
        r['state'] = s.state.name
        r['distance_representative'] = 'Yes' if s.distance_representative else 'No'
        r['contact_representative'] = 'Yes' if s.contact_representative else 'No'
        r['class_consensus_based_representative'] = 'Yes' if s.class_contact_representative else 'No'

        #r['contact_representative_score'] = "{:.0%}".format(s.contact_representative_score)

        #r['active_class_contacts_fraction'] = "{:.0%}".format(s.active_class_contacts_fraction)
        #r['inactive_class_contacts_fraction'] = "{:.0%}".format(s.inactive_class_contacts_fraction)
        #r['diff_class_contacts_fraction'] = "{:.0%}".format(s.inactive_class_contacts_fraction - s.active_class_contacts_fraction)

        r['mammal'] = 'Only show mammalian receptor structures (even if the non-mammalian is the only)' if s.mammal else ''
        r['closest_to_human'] = 'Only show structures from human or the closest species (for each receptor and state)' if s.closest_to_human else ''
        r['closest_to_human_raw'] = s.closest_to_human

        r['extra_filter'] = []
        if s.mammal and s.closest_to_human:
            r['extra_filter'] = ['*Only show mammalian structures and those from human or closest species',r['mammal'],r['closest_to_human']]
        elif s.mammal:
            r['extra_filter'] = [r['mammal']]
        elif s.closest_to_human:
            r['extra_filter'] = [r['closest_to_human']]

        r['identity_to_human'] = "100"
        if r['species'] != 'Human':
            key = 'identity_to_human_{}_{}'.format(shorted.family.slug,s.protein_conformation.protein.species.pk)
            if key in identity_lookup:
                  r['identity_to_human'] = identity_lookup[key]
            else:
                r['identity_to_human'] = cache.get(key)
                if r['identity_to_human'] == None:
                    try:
                        a = Alignment()
                        ref_p = Protein.objects.get(family = shorted.family, species__common_name = 'Human', sequence_type__slug = 'wt')
                        a.load_reference_protein(ref_p)
                        a.load_proteins([shorted])
                        a.load_segments(ProteinSegment.objects.filter(slug__in=['TM1', 'TM2', 'TM3', 'TM4','TM5','TM6', 'TM7']))
                        a.build_alignment()
                        a.calculate_similarity()
                        a.calculate_statistics()
                        p = a.proteins[1]
                        r['identity_to_human'] = int(p.identity)
                    except:
                        r['identity_to_human'] = 0
                    cache.set(key,r['identity_to_human'], 24*7*3600)
                identity_lookup[key] = r['identity_to_human']

        residues_wt = rcs[shorted.pk]
        residues_s = s.res_count
        # residues_s = residues_wt
        #print(pdb,"residues",protein,residues_wt,residues_s,residues_s/residues_wt)
        r['fraction_of_wt_seq'] = int(100*residues_s/residues_wt)

        a_list = []
        for a in s.stabilizing_agents.all():
            a_list.append(a)
        g_protein = only_gproteins(a_list)
        arrestin = only_arrestins(a_list)
        fusion = only_fusions(a_list)
        antibody = only_antibodies(a_list)

        r['signal_protein'] = ''
        r['signal_protein_subtype'] = ''
        r['signal_protein_note'] = ''
        r['signal_protein_seq_cons'] = ''
        r['signal_protein_seq_cons_color'] = ''

        ### StructureExtraProtein data parsing
        if pdb_id in ep:
            if effector:
                try:
                    key = '{}_{}_{}'.format(s.protein_conformation.protein.pk, s.protein_conformation.protein.parent.pk, ep[pdb_id]['display_name'])
                except:
                    key = '{}_{}_{}'.format(s.protein_conformation.protein.pk, s.protein_conformation.protein.parent, ep[pdb_id]['display_name'])
            else:
                key = '{}_{}'.format(s.protein_conformation.protein.parent.pk,ep[pdb_id]['display_name'])

            if best_signal_p[key] == ep[pdb_id]['wt_coverage']:
                # this is the best coverage
                r['signal_protein_seq_cons_color'] = 'green'
            else:
                r['signal_protein_seq_cons_color'] = 'red'
            if ep[pdb_id]['category'] == "Arrestin":
                 r['signal_protein'] = ep[pdb_id]['wt_protein__family__parent__parent__name']
            else:
                 r['signal_protein'] = ep[pdb_id]['wt_protein__family__parent__name']

            # Slight reformatting along the lines of the structure browser
            r['signal_protein_subtype'] = ep[pdb_id]['display_name']
            if ep[pdb_id]['category'] == "G alpha" and r['signal_protein_subtype'][0] == 'G':
                r['signal_protein_subtype'] = '&alpha;' + r['signal_protein_subtype'][1:]

            note = ep[pdb_id]['note']
            if note:
                if len(note) > 20:
                    r['signal_protein_note'] = "<span title='{}'>{}...</span>".format(note, note[:20])
                else:
                    r['signal_protein_note'] = note

            r['signal_protein_seq_cons'] = ep[pdb_id]['wt_coverage']

        #if pdb_id in methods:
        #    r['method'] = methods[pdb_id]
        #else:
        #    r['method'] = "N/A"
        r['method'] = s.structure_type.type_short()
        try:
            r['resolution'] = "{0:.2g}".format(s.resolution)
            r['resolution_best'] = s.resolution==best_resolutions['{}_{}'.format(shorted.pk, s.state.name)]
        except:
            r['resolution'] = None
            r['resolution_best'] = None


        r['7tm_distance'] = s.distance
        r['tm6_angle'] = str(round(s.tm6_angle)) if s.tm6_angle != None else ''
        r['gprot_bound_likeness'] = str(round(s.gprot_bound_likeness)) if s.gprot_bound_likeness != None else ''

        # DEBUGGING - overwrite with distance to 6x38
#        tm6_distance = ResidueAngle.objects.filter(structure__pdb_code__index=pdb_id.upper(), residue__generic_number__label="6x38")
#        if len(tm6_distance)>0:
#            tm6_distance = tm6_distance[0].core_distance
#        else:
#            tm6_distance = -1
#        r['7tm_distance'] = tm6_distance

#        r['tm6_angle'] = s.tm6_angle if s.tm6_angle != None else 0

        # DEBUGGING - overwrite with distance to 6x38-2x41
        #tm2_tm6_distance = Distance.objects.filter(structure__pdb_code__index=pdb_id.upper(), res2__generic_number__label="6x38", res1__generic_number__label="2x41")
        # DEBUGGING - overwrite with distance to 6x37-2x46
        # tm2_tm6_distance = Distance.objects.filter(structure__pdb_code__index=pdb_id.upper(), gns_pair="2x46_6x37")
        # DEBUGGING - overwrite with distance to 5x59-6x37
        # tm5_tm6_distance = Distance.objects.filter(structure__pdb_code__index=pdb_id.upper(), gns_pair="5x59_6x37")
        # if len(tm2_tm6_distance)>0 and len(tm5_tm6_distance)>0:
        #     r['7tm_distance'] = "{} {} {}".format(tm2_tm6_distance[0].distance/100, tm5_tm6_distance[0].distance/100, (tm2_tm6_distance[0].distance-tm5_tm6_distance[0].distance)/100)
        # else:
        #     r['7tm_distance'] = -1

        # DEBUGGING - overwrite with distance to 3x39-6x41
        #tm3_tm6_distance = Distance.objects.filter(structure__pdb_code__index=pdb_id.upper(), gns_pair="3x39_6x41")
        #if len(tm3_tm6_distance)>0:
        #    r['7tm_distance'] = tm3_tm6_distance[0].distance
        #else:
        #    r['7tm_distance'] = -1

        # DEBUGGING - overwrite with tm6 tilt angle
#       r['7tm_distance'] = s.tm6_angle if s.tm6_angle != None else 0

        r['g_protein'] = g_protein
        r['arrestin']  = arrestin
        r['fusion'] = fusion
        if len(antibody) > 20:
            antibody = "<span title='{}'>{}</span>".format(antibody, antibody[:20] + "..")
        r['antibody'] = antibody

        r['ligand'] = "-"
        r['ligand_function'] = "-"
        r['ligand_type'] = "-"

        for l in s.ligands.all():
            r['ligand'] = l.ligand.name
            if len(r['ligand'])>20:
                r['ligand'] = r['ligand'][:20] + ".."
            r['ligand_function'] = l.ligand_role.name
            if l.ligand.ligand_type != None:
                r['ligand_type'] = l.ligand.ligand_type.name

        # if pdb_id.startswith('AFM_'):
        #     if len(pdb_id)==8:
        #         pdb_id = pdb_id.split('_')[1]+'_refined'
        #     else:
        #         pdb_id = pdb_id.replace('_HUMAN', '')

        data_dict[pdb_id] = r
        data_table += "<tr> \
                        <td data-sort='0'><input class='form-check-input pdb_selected' type='checkbox' value='' onclick='thisPDB(this);' representative='{}' distance_representative='{}' class_consensus_based_representative='{}' long='{}'  id='{}'></td> \
                        <td>{}</td> \
                        <td><span>{}</span></td> \
                        <td>{}</td> \
                        <td>{}</td> \
                        <td>{}</td> \
                        <td><p class='no_margins' style='color:{}'>{}</td> \
                        <td>{}</td> \
                        <td>{}</td> \
                        <td>{}</td> \
                        <td class='shorten'>{}</td> \
                        <td><p class='no_margins' style='color:{}'>{}</p></td> \
                        <td>{}</td> \
                        <td>{}</td> \
                        <td>{}</td> \
                        <td>{}</td> \
                        <td>{}</td> \
                        <td>{}</td> \
                        <td>{}</td> \
                        <td><p class='no_margins' style='color:{}'>{}</p></td> \
                        <td>{}</td> \
                        <td>{}</td> \
                        <td>{}</td> \
                        <td>{}</td> \
                        </tr> \n".format(
                                        r['contact_representative'],
                                        r['distance_representative'],
                                        r['class_consensus_based_representative'],
                                        r['protein_long'],
                                        pdb_id,
                                        r['protein'],
                                        r['protein_long'],
                                        r['protein_family'],
                                        r['class'],
                                        r['fraction_of_wt_seq'],
                                        'green' if r['closest_to_human_raw'] else 'red',
                                        r['species'],
                                        'Best' if r['closest_to_human_raw'] else '',
                                        r['identity_to_human'],
                                        r['method'],
                                        pdb_id,
                                        'green' if r['resolution_best'] else 'red',
                                        r['resolution'],
                                        'Best' if r['resolution_best'] else '',
                                        r['state'],
                                        r['gprot_bound_likeness'],
                                        r['tm6_angle'],
                                        r['signal_protein'],
                                        r['signal_protein_subtype'],
                                        r['signal_protein_note'],
                                        r['signal_protein_seq_cons_color'],
                                        r['signal_protein_seq_cons'],
                                        r['fusion'],
                                        r['antibody'],
                                        r['ligand'],
                                        r['ligand_function'],
                                        )
    data_table += "</tbody></table>"
    return HttpResponse(data_table)

    # return render(request, 'contactnetwork/test.html', {'data_table':data_table})

@csrf_exempt
def InteractionBrowserData(request):
    start_time = time.time()
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

    mode = 'single'
    request_method = request.GET
    if request.POST and (request.POST.get("pdbs[]") or request.POST.get("pdbs1[]")):
        request_method = request.POST

    if request_method.get("pdbs[]"):
        pdbs = request_method.getlist('pdbs[]')
    elif request_method.get("pdbs1[]") and request_method.get("pdbs2[]"):
        pdbs1 = request_method.getlist('pdbs1[]')
        pdbs2 = request_method.getlist('pdbs2[]')
        mode = 'double'
    else:
        return "No selection"

    if mode == 'double':
        pdbs1 = [pdb.lower() for pdb in pdbs1]
        pdbs2 = [pdb.lower() for pdb in pdbs2]
        pdbs = pdbs1 + pdbs2
    else:
        pdbs = [pdb.lower() for pdb in pdbs]
        pdbs1 = pdbs
        pdbs2 = []

    pdbs_upper = [pdb.upper() for pdb in pdbs]

    # Deduce class
    gpcr_class = Structure.objects.filter(pdb_code__index__in=pdbs_upper
                ).values_list('protein_conformation__protein__parent__family__parent__parent__parent__slug', flat=True).distinct()
    if len(gpcr_class)>1:
        print('ERROR mix of classes!', gpcr_class)
        return JsonResponse({'error':list(gpcr_class)})
    else:
        gpcr_class = gpcr_class[0]

    # Segment filters
    try:
        segments = request_method.getlist('segments[]')
    except IndexError:
        segments = []

    # Interaction types
    try:
        i_types = [x.lower() for x in request_method.getlist('interaction_types[]')]
        # Add unknown type, so no interactions are returned, otherwise all are returned
        if len(i_types) == 0:
            i_types = ["DOES_NOT_EXIST"]
    except IndexError:
        i_types = []

    # Strict interaction settings
    try:
        strict_interactions = [x.lower() for x in request_method.getlist('strict_interactions[]')]
    except IndexError:
        strict_interactions = []

    i_types_filter = Q()
    if i_types:
        if len(strict_interactions) == 0:
            i_types_filter |= Q(interaction_type__in=i_types)
        else:
            # Merging the interaction filter with filters for the strict settings
            for int_type in i_types:
                if int_type in strict_interactions:
                    if int_type == 'polar' or int_type == 'aromatic':
                        i_types_filter = i_types_filter | (Q(interaction_type=int_type) & Q(interaction_level=0))
                    elif int_type == 'hydrophobic' or int_type == 'van-der-waals':
                        i_types_filter = i_types_filter | (Q(interaction_type=int_type) & Q(atompaircount__gte=4))
                else:
                    i_types_filter = i_types_filter | Q(interaction_type=int_type)

    # Options settings
    try:
        contact_options = [x.lower() for x in request_method.getlist('options[]')]
    except IndexError:
        contact_options = []

    i_options_filter = Q()
    # Filter out contact within the same helix
    #if contact_options and len(contact_options) > 0:
    # if not contact_options or "intrahelical" not in contact_options:
    #    i_options_filter = ~Q(interacting_pair__res1__protein_segment=F('interacting_pair__res2__protein_segment'))

    # Filter out contact with backbone atoms
    # backbone_atoms = ["C", "O", "N", "CA"]
    # if not contact_options or "backbone" not in contact_options:
    #    i_options_filter = i_options_filter & ~Q(atomname_residue1__in=backbone_atoms) & ~Q(atomname_residue2__in=backbone_atoms)

    # Filters for inter- and intrasegment contacts + BB/SC filters
    backbone_atoms = ["C", "O", "N", "CA"]
    pure_backbone_atoms = ["C", "O", "N"]

    # INTERsegment interactions
    inter_segments = ~Q(interacting_pair__res1__protein_segment=F('interacting_pair__res2__protein_segment'))

    # Blocking filter for empty settings
    inter_options_filter = (inter_segments & Q(atomname_residue1="FOO"))

    # BB-BB
    if contact_options and "inter_bbbb" in contact_options:
        bbbb = Q(atomname_residue1__in=backbone_atoms) & Q(atomname_residue2__in=backbone_atoms) & inter_segments
        inter_options_filter = inter_options_filter | bbbb

    # SC-BB
    if contact_options and "inter_scbb" in contact_options:
        scbb = ((Q(atomname_residue1__in=backbone_atoms) & ~Q(atomname_residue2__in=pure_backbone_atoms)) | (~Q(atomname_residue1__in=pure_backbone_atoms) & Q(atomname_residue2__in=backbone_atoms))) & inter_segments
        inter_options_filter = inter_options_filter | scbb

    # SC-SC
    if contact_options and "inter_scsc" in contact_options:
        scsc = ~Q(atomname_residue1__in=pure_backbone_atoms) & ~Q(atomname_residue2__in=pure_backbone_atoms) & inter_segments
        inter_options_filter = inter_options_filter | scsc

    # INTRAsegment interactions
    intra_segments = Q(interacting_pair__res1__protein_segment=F('interacting_pair__res2__protein_segment'))

    # Blocking filter for empty settings
    intra_options_filter = (intra_segments & Q(atomname_residue1="FOO"))

    # BB-BB
    if contact_options and "intra_bbbb" in contact_options:
        bbbb = Q(atomname_residue1__in=backbone_atoms) & Q(atomname_residue2__in=backbone_atoms) & intra_segments
        intra_options_filter = intra_options_filter | bbbb

    # SC-BB
    if contact_options and "intra_scbb" in contact_options:
        scbb = ((Q(atomname_residue1__in=backbone_atoms) & ~Q(atomname_residue2__in=pure_backbone_atoms)) | (~Q(atomname_residue1__in=pure_backbone_atoms) & Q(atomname_residue2__in=backbone_atoms))) & intra_segments
        intra_options_filter = intra_options_filter | scbb

    # SC-SC
    if contact_options and "intra_scsc" in contact_options:
        scsc = ~Q(atomname_residue1__in=pure_backbone_atoms) & ~Q(atomname_residue2__in=pure_backbone_atoms) & intra_segments
        intra_options_filter = intra_options_filter | scsc

    i_options_filter = inter_options_filter | intra_options_filter

    # DISCUSS: cache hash now takes the normalize along, is this necessary
    normalized = "normalize" in contact_options

    forced_class_a = "classa" in contact_options

    # Segment filters are now disabled
    segment_filter_res1 = Q()
    segment_filter_res2 = Q()

    # Cache
    hash_list = [pdbs1,pdbs2,i_types, strict_interactions, contact_options]
    hash_cache_key = 'interactionbrowserdata_{}'.format(get_hash(hash_list))
    data = cache.get(hash_cache_key)

    # data = None
    if data==None:
        cache_key = 'amino_acid_pair_conservation_{}_{}'.format(gpcr_class,forced_class_a)
        print('Before getting class cache',time.time()-start_time)
        class_pair_lookup = cache.get(cache_key)
        print('After getting class cache',time.time()-start_time)
        # class_pair_lookup=None
        if class_pair_lookup==None or len(class_pair_lookup)==0:
            # Class pair conservation
            sum_proteins = Protein.objects.filter(family__slug__startswith=gpcr_class,sequence_type__slug='wt',species__common_name='Human').count()
            residues = Residue.objects.filter(protein_conformation__protein__family__slug__startswith=gpcr_class,
                                              protein_conformation__protein__sequence_type__slug='wt',
                                              protein_conformation__protein__species__common_name='Human',

                        ).exclude(display_generic_number=None).values('pk','sequence_number','generic_number__label','amino_acid','protein_conformation__protein__entry_name','display_generic_number__label').all()
            r_pair_lookup = defaultdict(lambda: defaultdict(lambda: set()))
            for r in residues:
                # use the class specific generic number
                r['display_generic_number__label'] = re.sub(r'\.[\d]+', '', r['display_generic_number__label'])
                if forced_class_a:
                    r_pair_lookup[r['generic_number__label']][r['amino_acid']].add(r['protein_conformation__protein__entry_name'])
                else:
                    r_pair_lookup[r['display_generic_number__label']][r['amino_acid']].add(r['protein_conformation__protein__entry_name'])
            class_pair_lookup = {}

            gen_keys = sorted(r_pair_lookup.keys(), key=functools.cmp_to_key(gpcrdb_number_comparator))
            for i,gen1 in enumerate(gen_keys):
                v1 = r_pair_lookup[gen1]
                temp_score_dict = []
                for aa, protein in v1.items():
                    temp_score_dict.append([aa,len(protein)/sum_proteins])

                most_freq_aa = sorted(temp_score_dict.copy(), key = lambda x: -x[1])[0]
                class_pair_lookup[gen1] = most_freq_aa
                for gen2 in gen_keys[i:]:
                    if gen1 == gen2:
                        continue
                    pairs = {}
                    v2 = r_pair_lookup[gen2]
                    coord = '{},{}'.format(gen1,gen2)
                    for aa1 in v1.keys():
                        p1 = v1[aa1]
                        class_pair_lookup[gen1+aa1] = round(100*len(p1)/sum_proteins)
                        for aa2 in v2.keys():
                            pair = '{}{}'.format(aa1,aa2)
                            p2 = v2[aa2]
                            p = p1.intersection(p2)
                            if p:
                                class_pair_lookup[coord+pair] = round(100*len(p)/sum_proteins)
            cache.set(cache_key, class_pair_lookup, 3600 * 24 * 7)

        ### Fetch class ligand / G-protein interactions for snakeplot colouring
        cache_key = 'class_ligand_interactions_{}_{}'.format(gpcr_class,forced_class_a)
        class_ligand_interactions = cache.get(cache_key)
        # class_ligand_interactions=None
        if class_ligand_interactions==None or len(class_ligand_interactions)==0:
            class_interactions = ResidueFragmentInteraction.objects.filter(
                structure_ligand_pair__structure__protein_conformation__protein__family__slug__startswith=gpcr_class, structure_ligand_pair__annotated=True).prefetch_related(
                'rotamer__residue__generic_number',
                'rotamer__residue__display_generic_number',
                'structure_ligand_pair__structure__protein_conformation__protein__family')

            class_ligand_interactions = {}
            for i in class_interactions:
                p = i.structure_ligand_pair.structure.protein_conformation.protein.family.slug
                if i.rotamer.residue.generic_number:

                    display_gn = re.sub(r'\.[\d]+', '', i.rotamer.residue.generic_number.label)
                    if forced_class_a:
                        gn = i.rotamer.residue.generic_number.label
                    else:
                        gn = display_gn

                else:
                    continue
                if gn not in class_ligand_interactions.keys():
                    class_ligand_interactions[gn] = set()

                class_ligand_interactions[gn].add(p)
            class_ligand_interactions = {key: len(value) for key, value in class_ligand_interactions.items()}
            cache.set(cache_key, class_ligand_interactions, 3600 * 24 * 7)

        cache_key = 'class_complex_interactions_{}_{}'.format(gpcr_class,forced_class_a)
        class_complex_interactions = cache.get(cache_key)
        # class_ligand_interactions=None
        if class_complex_interactions == None or len(class_complex_interactions) == 0:

            interactions = Interaction.objects.filter(
                interacting_pair__referenced_structure__protein_conformation__protein__family__slug__startswith=gpcr_class
            ).exclude(
                interacting_pair__res1__protein_conformation_id=F('interacting_pair__res2__protein_conformation_id') # Filter interactions with other proteins
            ).exclude(interacting_pair__res1__generic_number__isnull=True
            ).distinct(
            ).exclude(
                specific_type='water-mediated'
            ).values(
                'interaction_type',
                'interacting_pair__referenced_structure__protein_conformation__protein__family__slug',
                'interacting_pair__res1__generic_number__label',
                'interacting_pair__res1__display_generic_number__label',
                'interacting_pair__res2__pk',
            )

            class_complex_interactions = {}
            for i in interactions:
                p = i['interacting_pair__referenced_structure__protein_conformation__protein__family__slug']

                display_gn = re.sub(r'\.[\d]+', '', i['interacting_pair__res1__display_generic_number__label'])
                if forced_class_a:
                    gn = i['interacting_pair__res1__generic_number__label']
                else:
                    gn = display_gn

                if gn not in class_complex_interactions.keys():
                    class_complex_interactions[gn] = set()

                class_complex_interactions[gn].add(p)


            class_complex_interactions = {key: len(value) for key, value in class_complex_interactions.items()}
            cache.set(cache_key, class_complex_interactions, 3600 * 24 * 7)

        cache_key = 'class_mutation_positions_{}_{}'.format(gpcr_class,forced_class_a)
        class_mutations = cache.get(cache_key)
        # class_ligand_interactions=None
        if class_mutations == None or len(class_mutations) == 0:

            class_mutations_q = MutationExperiment.objects.filter(protein__family__slug__startswith=gpcr_class).prefetch_related('protein__family','protein__parent__family','residue__generic_number','residue__display_generic_number').order_by('foldchange','exp_qual')

            class_mutations = {}
            for m in class_mutations_q:

                if m.residue.generic_number and abs(m.foldchange)>5:
                    p = m.protein.family.slug
                    foldchange = m.foldchange
                    d_gn = m.residue.display_generic_number.label
                    gn = m.residue.generic_number.label
                    display_gn = re.sub(r'\.[\d]+', '', d_gn)
                    if not forced_class_a:
                        gn = display_gn
                    if gn not in class_mutations.keys():
                        class_mutations[gn] = set()

                    class_mutations[gn].add(p)

            class_mutations = {key: len(value) for key, value in class_mutations.items()}
            cache.set(cache_key, class_mutations, 3600 * 24 * 7)

        # Get the relevant interactions
        # TODO MAKE SURE ITs only gpcr residues..
        interactions = Interaction.objects.filter(
            interacting_pair__referenced_structure__pdb_code__index__in=pdbs_upper
        ).filter(
            interacting_pair__res1__protein_conformation_id=F('interacting_pair__res2__protein_conformation_id') # Filter interactions with other proteins
        ).filter(
            interacting_pair__res1__pk__lt=F('interacting_pair__res2__pk')
        ).exclude(interacting_pair__res1__generic_number__isnull=True
        ).filter(
            segment_filter_res1 & segment_filter_res2
        ).values(
            'interaction_type',
            'interacting_pair__referenced_structure__pk',
            'interacting_pair__res1__pk',
            'interacting_pair__res2__pk',
        ).distinct(
        ).annotate(
             atompaircount=Count('interaction_type'),
             arr=ArrayAgg('pk')
        ).exclude(
            specific_type='water-mediated'
        ).filter(
            i_types_filter
        ).filter(
            i_options_filter
        ).order_by(
            'interaction_type',
            'interacting_pair__referenced_structure__pk',
            'interacting_pair__res1__pk',
            'interacting_pair__res2__pk'
        )

        # FOR DEBUGGING interaction + strict filters
        # print(interactions.query)
        interactions = list(interactions)

        # Grab unique interaction_IDs
        interaction_ids = []
        for entry in interactions:
            interaction_ids.extend(entry['arr'])

        # Interaction type sort - optimize by statically defining interaction type order
        order = ['ionic', 'polar', 'aromatic', 'hydrophobic', 'van-der-waals','None']
        interactions = sorted(interactions, key=lambda x: order.index(x['interaction_type']))

        data = {}
        data['class_ligand_interactions'] = class_ligand_interactions
        data['class_complex_interactions'] = class_complex_interactions
        data['class_mutations'] = class_mutations

        data['gpcr_class'] = gpcr_class

        # with open('{}_{}.txt'.format(gpcr_class,"gprotein"), 'w') as f:
        #     for key,d in class_complex_interactions.items():
        #         print(key,d)
        #         f.write("%s,%s\n"%(key,d))
        # with open('{}_{}.txt'.format(gpcr_class,"ligand"), 'w') as f:
        #     for key,d in class_ligand_interactions.items():
        #         f.write("%s,%s\n"%(key,d))
        # with open('{}_{}.txt'.format(gpcr_class,"mutation"), 'w') as f:
        #     for key,d in class_mutations.items():
        #         f.write("%s,%s\n"%(key,d))

        data['segments'] = set()
        data['segment_map'] = {}
        data['interactions'] = {}
        data['pdbs'] = set()
        data['proteins'] = set()
        data['pfs'] = set()
        data['pfs_lookup'] = defaultdict(lambda: [])
        data['tab3'] = {}
        data['tab4'] = {}
        data['aa_map'] = {}
        data['gn_map'] = OrderedDict()
        data['pos_map'] = OrderedDict()

        if mode == 'double':
            data['pdbs1'] = set()
            data['pdbs2'] = set()
            data['proteins1'] = set()
            data['proteins2'] = set()
            data['pfs1'] = set()
            data['pfs2'] = set()
            data['pfs1_lookup'] = defaultdict(lambda: [])
            data['pfs2_lookup'] = defaultdict(lambda: [])

        structures = Structure.objects.filter(pdb_code__index__in=pdbs_upper
                     ).select_related('protein_conformation__protein'
                     ).values('pk','pdb_code__index',
                            'protein_conformation__protein__parent__entry_name',
                            'protein_conformation__protein__parent__family__slug',
                            'protein_conformation__protein__entry_name')
        s_lookup = {}
        pdb_lookup = {}
        for s in structures:
            protein, pdb_name,pf  = [s['protein_conformation__protein__parent__entry_name'],s['protein_conformation__protein__entry_name'],s['protein_conformation__protein__parent__family__slug']]
            s_lookup[s['pk']] = [protein, pdb_name,pf]
            pdb_lookup[pdb_name] = [protein, s['pk'],pf]
            data['pfs_lookup'][pf].append(pdb_name)
            # List PDB files that were found in dataset.
            data['pdbs'] |= {pdb_name}
            data['proteins'] |= {protein}
            data['pfs'] |= {pf}

            # Populate the two groups lists
            if mode == 'double':
                if pdb_name in pdbs1:
                    data['pdbs1'] |= {pdb_name}
                    data['proteins1'] |= {protein}
                    data['pfs1'] |= {pf}
                    data['pfs1_lookup'][pf].append(pdb_name)
                if pdb_name in pdbs2:
                    data['pdbs2'] |= {pdb_name}
                    data['proteins2'] |= {protein}
                    data['pfs2'] |= {pf}
                    data['pfs2_lookup'][pf].append(pdb_name)

        if mode == 'double':
            if normalized:
                data['set1_size'] = len(data['pfs1'])
                data['set2_size'] = len(data['pfs2'])
            else:
                data['set1_size'] = len(data['pdbs1'])
                data['set2_size'] = len(data['pdbs2'])
        else:
            if normalized:
                data['set_size'] = len(data['pfs'])
            else:
                data['set_size'] = len(data['pdbs'])



        # Get all unique GNS to populate all residue tables (tab4)
        # TODO, check if can be deleted... it is regenerated later with class_specific numbers
        # distinct_gns = list(Residue.objects.filter(protein_conformation__protein__entry_name__in=pdbs).exclude(generic_number=None).values_list('generic_number__label','protein_segment__slug').distinct().order_by())

        all_pdbs_pairs = cache.get("all_pdbs_aa_pairs")
        # all_pdbs_pairs = None
        if not all_pdbs_pairs:
            # To save less, first figure out all possible interaction pairs
            pos_interactions = list(Interaction.objects.all(
            ).exclude(interacting_pair__res1__generic_number__isnull=True
            ).values_list(
                'interacting_pair__res1__generic_number__label',
                'interacting_pair__res2__generic_number__label',
            ).filter(interacting_pair__res1__pk__lt=F('interacting_pair__res2__pk')).distinct())

            all_interaction_pairs = []
            all_interaction_residues = set()
            for i in pos_interactions:
                all_interaction_pairs.append('{},{}'.format(i[0],i[1]))
                all_interaction_residues.add(i[0])
                all_interaction_residues.add(i[1])
            all_interaction_residues = sorted(list(all_interaction_residues), key=functools.cmp_to_key(gpcrdb_number_comparator))

            all_pdbs = list(Structure.objects.all().exclude(structure_type__slug__startswith='af-').values_list('pdb_code__index', flat=True))
            all_pdbs = [x.lower() for x in all_pdbs]
            #generic_number__label__in=all_interaction_residues)
            residues = Residue.objects.filter(protein_conformation__protein__entry_name__in=all_pdbs).exclude(generic_number=None).values(
                        'pk','sequence_number','display_generic_number__label','generic_number__label','amino_acid','protein_conformation__protein__entry_name','protein_segment__slug').all()

            r_lookup = {}
            r_pair_lookup = defaultdict(lambda: defaultdict(lambda: []))
            segm_lookup = {}
            r_presence_lookup = defaultdict(lambda: [])
            r_class_translate = {}

            for r in residues:
                if r['generic_number__label'] not in all_interaction_residues:
                    continue
                r_lookup[r['pk']] = r
                r_pair_lookup[r['generic_number__label']][r['amino_acid']].append(r['protein_conformation__protein__entry_name'])
                r_presence_lookup[r['generic_number__label']].append(r['protein_conformation__protein__entry_name'])
                segm_lookup[r['generic_number__label']] = r['protein_segment__slug']
                r['display_generic_number__label'] = re.sub(r'\.[\d]+', '', r['display_generic_number__label'])
                r_class_translate[r['generic_number__label']] = r['display_generic_number__label']

            gen_keys = sorted(r_pair_lookup.keys(), key=functools.cmp_to_key(gpcrdb_number_comparator))
            all_pdbs_pairs = {}
            for i,gen1 in enumerate(all_interaction_residues):
                for gen2 in all_interaction_residues[i:]:
                    if gen1 == gen2:
                        continue
                    pairs = {}
                    v1 = r_pair_lookup[gen1]
                    v2 = r_pair_lookup[gen2]
                    coord = '{},{}'.format(gen1,gen2)
                    if coord not in all_interaction_pairs:
                        continue
                    for aa1 in v1.keys():
                        for aa2 in v2.keys():
                            pair = '{}{}'.format(aa1,aa2)
                            p1 = set(v1[aa1])
                            p2 = set(v2[aa2])
                            p = list(p1.intersection(p2))
                            if p:
                                if coord not in all_pdbs_pairs:
                                    all_pdbs_pairs[coord] = {}
                                all_pdbs_pairs[coord][pair] = p
            cache.set("all_pdbs_aa_pairs",all_pdbs_pairs,60*60*24*7) #Cache results
        residues = Residue.objects.filter(protein_conformation__protein__entry_name__in=pdbs
                ).exclude(generic_number=None).values('pk','sequence_number','generic_number__label','amino_acid','protein_conformation__protein__entry_name','protein_segment__slug','display_generic_number__label').all()
        r_lookup = {}
        r_pair_lookup = defaultdict(lambda: defaultdict(lambda: []))
        segm_lookup = {}
        r_presence_lookup = defaultdict(lambda: [])
        r_class_translate = {}
        r_class_translate_from_classA = {}

        distinct_gns = []

        for r in residues:

            # remove .50 number from the display number format (1.50x50), so only the GPCRdb number is left
            r['display_generic_number__label'] = re.sub(r'\.[\d]+', '', r['display_generic_number__label'])
            if forced_class_a:
                r_class_translate[r['generic_number__label']] = r['generic_number__label']
                r_class_translate_from_classA[r['generic_number__label']] = r['generic_number__label']
            else:
                # If not, then use the class relevant numbers
                r_class_translate[r['display_generic_number__label']] = r['generic_number__label']
                r_class_translate_from_classA[r['generic_number__label']] = r['display_generic_number__label']
                # change the generic number (class a) to the class specific one.
                r['generic_number__label'] = r['display_generic_number__label']

            r_lookup[r['pk']] = r
            r_pair_lookup[r['generic_number__label']][r['amino_acid']].append(r['protein_conformation__protein__entry_name'])
            segm_lookup[r['generic_number__label']] = r['protein_segment__slug']
            r_presence_lookup[r['generic_number__label']].append(r['protein_conformation__protein__entry_name'])
            data['segments'].add(r['protein_segment__slug'])

            # Generate all distinct (class-specific) GNs for tab4
            if [r['generic_number__label'],r['protein_segment__slug']] not in distinct_gns:
                distinct_gns.append([r['generic_number__label'],r['protein_segment__slug']])


        data['segment_map'] = segm_lookup

        updated_all_pdbs_pairs = {}
        for coord, d in all_pdbs_pairs.items():
            gen1 = coord.split(",")[0]
            gen2 = coord.split(",")[1]
            if gen1 in r_class_translate_from_classA and gen2 in r_class_translate_from_classA:
                # if gen1 and gen2 aren't in translate dictionary, then they're not going to be relevant later.
                coord_new = '{},{}'.format(r_class_translate_from_classA[gen1],r_class_translate_from_classA[gen2])
                updated_all_pdbs_pairs[coord_new] = d

        all_pdbs_pairs = updated_all_pdbs_pairs

        # Dict to keep track of which residue numbers are in use
        number_dict = set()

        print('Start going through interactions',time.time()-start_time)
        for i in interactions:
            s = i['interacting_pair__referenced_structure__pk']
            pdb_name = s_lookup[s][1]
            protein = s_lookup[s][0]
            pf = s_lookup[s][2]
            res1 = r_lookup[i['interacting_pair__res1__pk']]
            res2 = r_lookup[i['interacting_pair__res2__pk']]
            res1_seq = res1['sequence_number']
            res2_seq = res2['sequence_number']
            res1_aa = res1['amino_acid']
            res2_aa = res2['amino_acid']
            res1_gen = res1['generic_number__label']
            res2_gen = res2['generic_number__label']
            model = i['interaction_type']

            res1 = res1_gen
            res2 = res2_gen

            if res1 < res2 or res1_seq < res2_seq:
                coord = str(res1) + ',' + str(res2)
                classa_coord = str(r_class_translate[res1]) + ',' + str(r_class_translate[res2])
            else:
                coord = str(res2) + ',' + str(res1)
                classa_coord = str(r_class_translate[res2]) + ',' + str(r_class_translate[res1])
                res1_aa, res2_aa = res2_aa, res1_aa

            # Populate the AA map
            if pdb_name not in data['aa_map']:
                data['aa_map'][pdb_name] = {}

            data['aa_map'][res1] = res1_aa
            data['aa_map'][res2] = res2_aa

            number_dict |= {res1, res2}

            normalized_value = pdb_name
            if normalized: normalized_value = pf

            if mode == 'double':
                if res1 not in data['tab3']:
                    data['tab3'][res1] = {'set1':set(),'set2':set(), 'set1_count':set(), 'set2_count':set(), 'set1_aa': set(), 'set2_aa': set()}
                if res2 not in data['tab3']:
                    data['tab3'][res2] = {'set1':set(),'set2':set(), 'set1_count':set(), 'set2_count':set(), 'set1_aa': set(), 'set2_aa': set()}
                if coord not in data['interactions']:
                    data['interactions'][coord] = {'pdbs1':[], 'proteins1': [], 'pfs1':[], 'pdbs2':[], 'proteins2': [], 'pfs2':[], 'secondary1' : [], 'secondary2' : [], 'class_seq_cons' : [0,0], 'types' : [], 'types_count' : {}}
                    for i_type in ['ionic', 'polar', 'aromatic', 'hydrophobic', 'van-der-waals']:
                        data['interactions'][coord]['types_count'][i_type] = [{'pdbs':[],'pdb_freq':0,'pf_freq':0},{'pdbs':[],'pdb_freq':0,'pf_freq':0}] #set1, #set2

                if model in i_types:
                    if model not in data['interactions'][coord]['types']:
                        data['interactions'][coord]['types'].append(model)
                if pdb_name in pdbs1:
                    if model in i_types:
                        if pdb_name not in data['interactions'][coord]['pdbs1']:
                            data['interactions'][coord]['pdbs1'].append(pdb_name)
                        if protein not in data['interactions'][coord]['proteins1']:
                            data['interactions'][coord]['proteins1'].append(protein)
                        if pf not in data['interactions'][coord]['pfs1']:
                            data['interactions'][coord]['pfs1'].append(pf)

                    data['interactions'][coord]['secondary1'].append([model,res1_aa,res2_aa,pdb_name])
                    data['tab3'][res1]['set1'].add(res2)
                    data['tab3'][res2]['set1'].add(res1)
                    data['tab3'][res1]['set1_aa'].add(res1_aa)
                    data['tab3'][res2]['set1_aa'].add(res2_aa)
                    data['tab3'][res1]['set1_count'].add('{},{}'.format(pdb_name,res2))
                    data['tab3'][res2]['set1_count'].add('{},{}'.format(pdb_name,res1))

                    if normalized_value not in data['interactions'][coord]['types_count'][model][0]['pdbs']:
                            data['interactions'][coord]['types_count'][model][0]['pdbs'].append(pdb_name)

                if pdb_name in pdbs2:
                    if model in i_types:
                        if pdb_name not in data['interactions'][coord]['pdbs2']:
                            data['interactions'][coord]['pdbs2'].append(pdb_name)
                        if protein not in data['interactions'][coord]['proteins2']:
                            data['interactions'][coord]['proteins2'].append(protein)
                        if pf not in data['interactions'][coord]['pfs2']:
                            data['interactions'][coord]['pfs2'].append(pf)
                    data['interactions'][coord]['secondary2'].append([model,res1_aa,res2_aa,pdb_name])
                    data['tab3'][res1]['set2'].add(res2)
                    data['tab3'][res2]['set2'].add(res1)
                    data['tab3'][res1]['set2_aa'].add(res1_aa)
                    data['tab3'][res2]['set2_aa'].add(res2_aa)
                    data['tab3'][res1]['set2_count'].add('{},{}'.format(pdb_name,res2))
                    data['tab3'][res2]['set2_count'].add('{},{}'.format(pdb_name,res1))

                    if normalized_value not in data['interactions'][coord]['types_count'][model][1]['pdbs']:
                            data['interactions'][coord]['types_count'][model][1]['pdbs'].append(pdb_name)

                ## Presence lookup
                pdbs_with_res1 = r_presence_lookup[res1]
                pdbs_with_res2 = r_presence_lookup[res2]

                pdbs1_with_res1 = list(set(pdbs_with_res1).intersection(pdbs1))
                pdbs2_with_res1 = list(set(pdbs_with_res1).intersection(pdbs2))

                pdbs1_with_res2 = list(set(pdbs_with_res2).intersection(pdbs1))
                pdbs2_with_res2 = list(set(pdbs_with_res2).intersection(pdbs2))

                data['interactions'][coord]['pos1_presence'] = round(100*(len(pdbs1_with_res1) / len(pdbs1))-(100*len(pdbs2_with_res1) / len(pdbs2)))
                data['interactions'][coord]['pos2_presence'] = round(100*(len(pdbs1_with_res2) / len(pdbs1))-(100*len(pdbs2_with_res2) / len(pdbs2)))

            else:
                if res1 not in data['tab3']:
                    data['tab3'][res1] = {'unique':set(),'count':set(),'set_aa':set()}
                if res2 not in data['tab3']:
                    data['tab3'][res2] = {'unique':set(),'count':set(),'set_aa':set()}
                data['tab3'][res1]['unique'].add(res2)
                data['tab3'][res2]['unique'].add(res1)
                data['tab3'][res1]['set_aa'].add(res1_aa)
                data['tab3'][res2]['set_aa'].add(res2_aa)
                data['tab3'][res1]['count'].add('{},{}'.format(pdb_name,res2))
                data['tab3'][res2]['count'].add('{},{}'.format(pdb_name,res1))
                if coord not in data['interactions']:
                    data['interactions'][coord] = {'pdbs':[], 'proteins': [], 'pfs': [], 'secondary': [], 'class_seq_cons' : 0, 'types' : [], 'types_count' : {}, 'seq_pos':[res1_seq,res2_seq]}
                    for i_type in ['ionic', 'polar', 'aromatic', 'hydrophobic', 'van-der-waals']:
                        data['interactions'][coord]['types_count'][i_type] = {'pdbs':[],'pdb_freq':0,'pf_freq':0} #set

                if model in i_types or not i_types:
                    if model not in data['interactions'][coord]['types']:
                        data['interactions'][coord]['types'].append(model)
                if pdb_name not in data['interactions'][coord]['pdbs']:
                    data['interactions'][coord]['pdbs'].append(pdb_name)
                if protein not in data['interactions'][coord]['proteins']:
                    data['interactions'][coord]['proteins'].append(protein)
                if pf not in data['interactions'][coord]['pfs']:
                    data['interactions'][coord]['pfs'].append(pf)
                data['interactions'][coord]['secondary'].append([model,res1_aa,res2_aa,pdb_name])

                if normalized_value not in data['interactions'][coord]['types_count'][model]['pdbs']:
                    data['interactions'][coord]['types_count'][model]['pdbs'].append(pdb_name)

                ## Presence lookup
                pdbs_with_res1 = r_presence_lookup[res1]
                pdbs_with_res2 = r_presence_lookup[res2]

                pdbs1_with_res1 = list(set(pdbs_with_res1).intersection(pdbs1))
                pdbs1_with_res2 = list(set(pdbs_with_res2).intersection(pdbs1))

                data['interactions'][coord]['pos1_presence'] = round(100*len(pdbs1_with_res1) / len(pdbs1))
                data['interactions'][coord]['pos2_presence'] = round(100*len(pdbs1_with_res2) / len(pdbs1))
            data['interactions'][coord]['class_a_gns'] = classa_coord
        data['sequence_numbers'] = sorted(number_dict, key=functools.cmp_to_key(gpcrdb_number_comparator))

        ## MAKE TAB 4
        data['missing'] = {}
        for res in distinct_gns:
            if res[0] not in data['tab4']:
                data['tab4'][res[0]] = {'ps':res[1], 'angles':[], 'angles_set':[]}
            if res[0] not in data['missing']:
                data['missing'][res[0]] = {'present':set()}

            for pdb in pdbs:
                if pdb in r_presence_lookup[res[0]]:
                    if normalized:
                        data['missing'][res[0]]['present'].add(pdb_lookup[pdb][2])
                    else:
                        data['missing'][res[0]]['present'].add(pdb)
            data['missing'][res[0]]['present'] = list(data['missing'][res[0]]['present'])

        print('Do Secondary data',time.time()-start_time)
        data['secondary'] = {}
        secondary_dict = {'set1':0 , 'set2':0, 'aa_pairs':OrderedDict()}
        secondary_dict_single = {'set':0 , 'aa_pairs':OrderedDict()}
        aa_pairs_dict = {'set1':0 , 'set2':0, 'class':{}}
        aa_pairs_dict_single = {'set':0, 'class':{}}
        delete_coords = []
        # print(len(data['interactions']),'interactions2')
        for c,v in data['interactions'].items():
            if mode == 'double':
                if len(v["pdbs1"])+len(v["pdbs2"])==0:
                    #empty
                    delete_coords.append(c)
                    continue
                data['secondary'][c] = OrderedDict()
                current = {}
                current["set1"] = pdbs1.copy()
                current["set2"] = pdbs2.copy()
                v["pdbs_freq_1"] = len(v["pdbs1"]) / len(pdbs1)
                v["pdbs_freq_2"] = len(v["pdbs2"]) / len(pdbs2)

                #pf freq
                v["pf_freq_1"] = 0
                for pf, pf_pdbs in data['pfs1_lookup'].items():
                    pf_contribution = 0
                    for pdb in pf_pdbs:
                        # if pdb from pf has an interaction, add the fraction of the pf set
                        if pdb in v["pdbs1"]:
                            pf_contribution += 1 / len(pf_pdbs)
                    v["pf_freq_1"] += pf_contribution
                v["pf_freq_1"] /= len(data['pfs1'])

                v["pf_freq_2"] = 0
                for pf, pf_pdbs in data['pfs2_lookup'].items():
                    pf_contribution = 0
                    for pdb in pf_pdbs:
                        # if pdb from pf has an interaction, add the fraction of the pf set
                        if pdb in v["pdbs2"]:
                            pf_contribution += 1 / len(pf_pdbs)
                    v["pf_freq_2"] += pf_contribution
                v["pf_freq_2"] /= len(data['pfs2'])

                for i_t, vals in v['types_count'].items():
                    vals[0]['pdb_freq'] = len(vals[0]["pdbs"]) / len(pdbs1)
                    vals[1]['pdb_freq'] = len(vals[1]["pdbs"]) / len(pdbs2)

                    for pf, pf_pdbs in data['pfs1_lookup'].items():
                        pf_contribution = 0
                        for pdb in pf_pdbs:
                            # if pdb from pf has an interaction, add the fraction of the pf set
                            if pdb in vals[0]["pdbs"]:
                                pf_contribution += 1 / len(pf_pdbs)
                        vals[0]["pf_freq"] += pf_contribution
                    vals[0]["pf_freq"] /= len(data['pfs1'])

                    for pf, pf_pdbs in data['pfs2_lookup'].items():
                        pf_contribution = 0
                        for pdb in pf_pdbs:
                            # if pdb from pf has an interaction, add the fraction of the pf set
                            if pdb in vals[1]["pdbs"]:
                                pf_contribution += 1 / len(pf_pdbs)
                        vals[1]["pf_freq"] += pf_contribution
                    vals[1]["pf_freq"] /= len(data['pfs2'])

                for setname,iset in [['set1','secondary1'],['set2','secondary2']]:
                    distinct_aa_pairs = set()
                    for s in v[iset]:
                        i = s[0]
                        aa_pair = ''.join(s[1:3])
                        distinct_aa_pairs.add(aa_pair)
                        if s[3] in current[setname]:
                            #remove PDB from current set, to deduce those without an interaction
                            current[setname].remove(s[3])
                        if i not in data['secondary'][c]:
                            data['secondary'][c][i] = copy.deepcopy(secondary_dict)
                        data['secondary'][c][i][setname] += 1
                        if aa_pair not in data['secondary'][c][i]['aa_pairs']:
                            data['secondary'][c][i]['aa_pairs'][aa_pair] = copy.deepcopy(aa_pairs_dict)
                            # Count overall occurances in sets
                            aa1 = s[1]
                            aa2 = s[2]
                            gen1 = c.split(",")[0]
                            gen2 = c.split(",")[1]
                            pdbs_with_aa1 = r_pair_lookup[gen1][aa1]
                            pdbs_with_aa2 = r_pair_lookup[gen2][aa2]
                            pdbs_intersection = list(set(pdbs_with_aa1).intersection(pdbs_with_aa2))
                            pdbs1_with_pair = list(set(pdbs_intersection).intersection(pdbs1))
                            pdbs2_with_pair = list(set(pdbs_intersection).intersection(pdbs2))
                            data['secondary'][c][i]['aa_pairs'][aa_pair]['pair_set1'] = pdbs1_with_pair
                            data['secondary'][c][i]['aa_pairs'][aa_pair]['pair_set2'] = pdbs2_with_pair

                            if c+aa_pair in class_pair_lookup:
                                data['secondary'][c][i]['aa_pairs'][aa_pair]['class'] = class_pair_lookup[c+aa_pair]
                            else:
                                data['secondary'][c][i]['aa_pairs'][aa_pair]['class'] = "-"

                        data['secondary'][c][i]['aa_pairs'][aa_pair][setname] += 1

                    sum_cons = 0
                    for aa_pair in distinct_aa_pairs:
                        if c+aa_pair in class_pair_lookup:
                            sum_cons += class_pair_lookup[c+aa_pair]

                    i = 0 if setname=='set1' else 1
                    v['class_seq_cons'][i] = sum_cons

                i = 'None' ## Remember to also have this name in the "order" dict.
                data['secondary'][c][i] = copy.deepcopy(secondary_dict)
                for setname in ['set1','set2']:
                    data['secondary'][c][i][setname] += len(current[setname])
                    for aa_pair, pdbs in all_pdbs_pairs[c].items():
                        if aa_pair not in data['secondary'][c][i]['aa_pairs']:
                            data['secondary'][c][i]['aa_pairs'][aa_pair] = copy.deepcopy(aa_pairs_dict)
                            if c+aa_pair in class_pair_lookup:
                                data['secondary'][c][i]['aa_pairs'][aa_pair]['class'] = class_pair_lookup[c+aa_pair]
                            else:
                                data['secondary'][c][i]['aa_pairs'][aa_pair]['class'] = "-"

                            aa1 = aa_pair[0]
                            aa2 = aa_pair[1]
                            gen1 = c.split(",")[0]
                            gen2 = c.split(",")[1]
                            pdbs_with_aa1 = r_pair_lookup[gen1][aa1]
                            pdbs_with_aa2 = r_pair_lookup[gen2][aa2]
                            pdbs_intersection = list(set(pdbs_with_aa1).intersection(pdbs_with_aa2))
                            pdbs1_with_pair = list(set(pdbs_intersection).intersection(pdbs1))
                            pdbs2_with_pair = list(set(pdbs_intersection).intersection(pdbs2))
                            data['secondary'][c][i]['aa_pairs'][aa_pair]['pair_set1'] = pdbs1_with_pair
                            data['secondary'][c][i]['aa_pairs'][aa_pair]['pair_set2'] = pdbs2_with_pair


                        for pdb in current[setname]:
                            if pdb in pdbs:
                                # if pdb without interaction is in pdbs of aa_pair, add one.
                                data['secondary'][c][i]['aa_pairs'][aa_pair][setname] += 1
                        if setname == 'set2':
                            # if 2nd run
                            if data['secondary'][c][i]['aa_pairs'][aa_pair]["set1"] == 0 and data['secondary'][c][i]['aa_pairs'][aa_pair]["set2"] == 0:
                                del data['secondary'][c][i]['aa_pairs'][aa_pair]

                # Order based on AA counts
                for i in data['secondary'][c].keys():
                    data['secondary'][c][i]['aa_pairs'] = OrderedDict(sorted(data['secondary'][c][i]['aa_pairs'].items(), key=lambda x: x[1]["set1"]+x[1]["set2"], reverse = True))

                data['secondary'][c] = OrderedDict(sorted(data['secondary'][c].items(), key=lambda x: order.index(x[0])))
            elif mode =='single':
                # continue
                if len(v["pdbs"])==0:
                    #empty
                    delete_coords.append(c)
                    continue
                data['secondary'][c] = OrderedDict()
                current = {}
                current["set"] = pdbs1.copy()
                setname = "set"
                distinct_aa_pairs = set()

                #Calculate pdb_freqs
                v["pdbs_freq"] = len(v["pdbs"]) / len(pdbs1)

                #pf freq
                v["pf_freq"] = 0
                for pf, pf_pdbs in data['pfs_lookup'].items():
                    pf_contribution = 0
                    for pdb in pf_pdbs:
                        # if pdb from pf has an interaction, add the fraction of the pf set
                        if pdb in v["pdbs"]:
                            pf_contribution += 1 / len(pf_pdbs)
                    v["pf_freq"] += pf_contribution
                v["pf_freq"] /= len(data['pfs'])

                for i_t, vals in v['types_count'].items():
                    vals['pdb_freq'] = len(vals["pdbs"]) / len(pdbs1)
                    for pf, pf_pdbs in data['pfs_lookup'].items():
                        pf_contribution = 0
                        for pdb in pf_pdbs:
                            # if pdb from pf has an interaction, add the fraction of the pf set
                            if pdb in vals["pdbs"]:
                                pf_contribution += 1 / len(pf_pdbs)
                        vals["pf_freq"] += pf_contribution
                    vals["pf_freq"] /= len(data['pfs'])


                for s in v['secondary']:
                    i = s[0]
                    aa_pair = ''.join(s[1:3])
                    distinct_aa_pairs.add(aa_pair)
                    if s[3] in current[setname]:
                        #remove PDB from current set, to deduce those without an interaction
                        current[setname].remove(s[3])
                    if i not in data['secondary'][c]:
                        data['secondary'][c][i] = copy.deepcopy(secondary_dict_single)
                    data['secondary'][c][i][setname] += 1
                    if aa_pair not in data['secondary'][c][i]['aa_pairs']:
                        data['secondary'][c][i]['aa_pairs'][aa_pair] = copy.deepcopy(aa_pairs_dict_single)
                        # Count overall occurances in sets
                        aa1 = s[1]
                        aa2 = s[2]
                        gen1 = c.split(",")[0]
                        gen2 = c.split(",")[1]
                        pdbs_with_aa1 = r_pair_lookup[gen1][aa1]
                        pdbs_with_aa2 = r_pair_lookup[gen2][aa2]
                        pdbs_intersection = list(set(pdbs_with_aa1).intersection(pdbs_with_aa2))
                        pdbs_with_pair = list(set(pdbs_intersection).intersection(pdbs1))
                        data['secondary'][c][i]['aa_pairs'][aa_pair]['pair_set'] = pdbs_with_pair

                        if c+aa_pair in class_pair_lookup:
                            data['secondary'][c][i]['aa_pairs'][aa_pair]['class'] = class_pair_lookup[c+aa_pair]
                        else:
                            data['secondary'][c][i]['aa_pairs'][aa_pair]['class'] = "-"

                    data['secondary'][c][i]['aa_pairs'][aa_pair][setname] += 1

                sum_cons = 0
                for aa_pair in distinct_aa_pairs:
                    if c+aa_pair in class_pair_lookup:
                        sum_cons += class_pair_lookup[c+aa_pair]

                v['class_seq_cons'] = sum_cons
                i = 'None' ## Remember to also have this name in the "order" dict.
                data['secondary'][c][i] = copy.deepcopy(secondary_dict_single)
                data['secondary'][c][i][setname] += len(current[setname])
                for aa_pair, pdbs in all_pdbs_pairs[c].items():
                    if aa_pair not in data['secondary'][c][i]['aa_pairs']:
                        data['secondary'][c][i]['aa_pairs'][aa_pair] = copy.deepcopy(aa_pairs_dict_single)
                        if c+aa_pair in class_pair_lookup:
                            data['secondary'][c][i]['aa_pairs'][aa_pair]['class'] = class_pair_lookup[c+aa_pair]
                        else:
                            data['secondary'][c][i]['aa_pairs'][aa_pair]['class'] = "-"

                        aa1 = aa_pair[0]
                        aa2 = aa_pair[1]
                        gen1 = c.split(",")[0]
                        gen2 = c.split(",")[1]
                        pdbs_with_aa1 = r_pair_lookup[gen1][aa1]
                        pdbs_with_aa2 = r_pair_lookup[gen2][aa2]
                        pdbs_intersection = list(set(pdbs_with_aa1).intersection(pdbs_with_aa2))
                        pdbs_with_pair = list(set(pdbs_intersection).intersection(pdbs1))
                        data['secondary'][c][i]['aa_pairs'][aa_pair]['pair_set'] = pdbs_with_pair


                    for pdb in current[setname]:
                        if pdb in pdbs:
                            # if pdb without interaction is in pdbs of aa_pair, add one.
                            data['secondary'][c][i]['aa_pairs'][aa_pair][setname] += 1
                    if data['secondary'][c][i]['aa_pairs'][aa_pair]["set"] == 0:
                        del data['secondary'][c][i]['aa_pairs'][aa_pair]
        # print(delete_coords)
        for d in delete_coords:
            del data['interactions'][d]

        # del class_pair_lookup
        del r_lookup
        # del r_pair_lookup


        ## PREPARE ADDITIONAL DATA (INTERACTIONS AND ANGLES)
        print('Prepare distance values for',mode,'mode',time.time()-start_time)
        interaction_keys = [k.replace(",","_") for k in data['interactions'].keys()]
        interaction_keys = [v['class_a_gns'].replace(",","_") for k,v in data['interactions'].items()]
        if mode == "double":

            group_1_distances = get_distance_averages(data['pdbs1'],s_lookup, interaction_keys,normalized, standard_deviation = False)
            group_2_distances = get_distance_averages(data['pdbs2'],s_lookup, interaction_keys,normalized, standard_deviation = False)

            print('got distance values for',mode,'mode',time.time()-start_time)
            for coord in data['interactions']:
                distance_coord = coord.replace(",", "_")
                # Replace coord to ensure using classA as distances are indexed with those
                distance_coord = data['interactions'][coord]['class_a_gns'].replace(",", "_")
                if distance_coord in group_1_distances and distance_coord in group_2_distances:
                    distance_diff = round(group_1_distances[distance_coord]-group_2_distances[distance_coord],0)
                else:
                    distance_diff = ""
                data['interactions'][coord]['distance'] = distance_diff
            print('Done merging distance values for',mode,'mode',time.time()-start_time)
        else:
            group_distances = get_distance_averages(data['pdbs'],s_lookup, interaction_keys,normalized, standard_deviation = True)
            for coord in data['interactions']:
                distance_coord = coord.replace(",","_")
                # Replace coord to ensure using classA as distances are indexed with those
                distance_coord = data['interactions'][coord]['class_a_gns'].replace(",", "_")
                if distance_coord in group_distances:
                    distance = round(group_distances[distance_coord],0)
                else:
                    distance = ""
                data['interactions'][coord]['distance'] = distance

        # del class_pair_lookup
        # del r_pair_lookup
        print('Prepare all angles values for',mode,'mode',time.time()-start_time)
        data['all_angles'] = get_all_angles(pdbs_upper,data['pfs'],normalized, forced_class_a = forced_class_a)
        print('Prepare angles values for',mode,'mode',time.time()-start_time)

        if mode == "double":

            group_1_angles = get_angle_averages(data['pdbs1'],s_lookup, normalized, forced_class_a = forced_class_a)
            group_2_angles = get_angle_averages(data['pdbs2'],s_lookup, normalized, forced_class_a = forced_class_a)
            data['all_angles_set1'] = get_all_angles(data['pdbs1'],data['pfs1'],normalized, forced_class_a = forced_class_a)
            data['all_angles_set2'] = get_all_angles(data['pdbs2'],data['pfs2'],normalized, forced_class_a = forced_class_a)

            print('got angles values for',mode,'mode',time.time()-start_time)
            custom_angles = ['a_angle', 'outer_angle', 'phi', 'psi', 'theta', 'tau']
            index_names = {0:'core_distance',1:'a_angle',2:'outer_angle',3:'tau',4:'phi',5:'psi',6: 'sasa',7: 'rsa',8:'theta',9:'hse',10:'tau_angle', 11:'tau', 12:'rotation_angle'}

            for coord in data['interactions']:
                gn1 = coord.split(",")[0]
                gn2 = coord.split(",")[1]

                gn1_values = [['','','']] * 11
                if gn1 in group_1_angles and gn1 in group_2_angles:
                    gn1_values = []
                    for i,v in enumerate(group_1_angles[gn1]):
                        try:
                            if index_names[i] in custom_angles:
                                diff = abs(v-group_2_angles[gn1][i])
                                diff = round(min(360-diff,diff))
                            else:
                                diff = round(v-group_2_angles[gn1][i],0)
                            gn1_values.append([diff,v,group_2_angles[gn1][i]])
                        except:
                            # Fails if there is a None (like gly doesnt have outer angle?)
                            gn1_values.append(['','',''])

                gn2_values = [['','','']] * 11
                if gn2 in group_1_angles and gn2 in group_2_angles:
                    gn2_values = []
                    for i,v in enumerate(group_1_angles[gn2]):
                        try:
                            if index_names[i] in custom_angles:
                                diff = abs(v-group_2_angles[gn2][i])
                                diff = round(min(360-diff,diff))
                            else:
                                diff = round(v-group_2_angles[gn2][i],0)
                            gn2_values.append([diff,v,group_2_angles[gn2][i]])
                        except:
                            # Fails if there is a None (like gly doesnt have outer angle?)
                            gn2_values.append(['','',''])
                data['interactions'][coord]['angles'] = [gn1_values,gn2_values]

            for gn in data['tab4'].keys():
                gn_values = [['','','']] * 11
                if gn in group_1_angles and gn in group_2_angles:
                    gn_values = []
                    for i,v in enumerate(group_1_angles[gn]):
                        try:
                            if index_names[i] in custom_angles:
                                diff = abs(v-group_2_angles[gn][i])
                                diff = round(min(360-diff,diff))
                            else:
                                diff = round(v-group_2_angles[gn][i],0)
                            gn_values.append([diff,v,group_2_angles[gn][i]])
                        except:
                            # Fails if there is a None (like gly doesnt have outer angle?)
                            gn_values.append(['','',''])
                data['tab4'][gn]['angles_set'] = gn_values
                data['tab4'][gn]['angles'] = gn_values

                if gn in group_1_angles:
                    data['tab4'][gn]['angles_set1'] = [ round(elem) if isinstance(elem, float) else '' for elem in group_1_angles[gn] ]
                else:
                    data['tab4'][gn]['angles_set1'] = [''] * 11
                if gn in group_2_angles:
                    data['tab4'][gn]['angles_set2'] = [ round(elem) if isinstance(elem, float) else '' for elem in group_2_angles[gn] ]
                else:
                    data['tab4'][gn]['angles_set2'] = [''] * 11


            print('Done combining data',mode,'mode',time.time()-start_time)
        else:

            # get_angle_averages gets "mean" in case of single pdb
            group_angles = get_angle_averages(data['pdbs'],s_lookup, normalized, standard_deviation=True, forced_class_a = forced_class_a)
            for coord in data['interactions']:
                gn1 = coord.split(",")[0]
                gn2 = coord.split(",")[1]

                gn1_values = [''] * 11
                if gn1 in group_angles:
                    gn1_values = []
                    for i,v in enumerate(group_angles[gn1]):
                        try:
                            gn1_values.append(round(v))
                        except:
                            # Fails if there is a None (like gly doesnt have outer angle?)
                            gn1_values.append("")

                gn2_values = [''] * 11
                if gn2 in group_angles:
                    gn2_values = []
                    for i,v in enumerate(group_angles[gn2]):
                        try:
                            gn2_values.append(round(v))
                        except:
                            # Fails if there is a None (like gly doesnt have outer angle?)
                            gn2_values.append("")
                data['interactions'][coord]['angles'] = [gn1_values,gn2_values]
                data['tab3'][gn1]['angles_set'] = gn1_values
                data['tab3'][gn2]['angles_set'] = gn2_values
                data['tab3'][gn1]['angles'] = gn1_values
                data['tab3'][gn2]['angles'] = gn2_values
                for gn in data['tab4'].keys():
                    gn_values = [''] * 11
                    if gn in group_angles:
                        gn_values = []
                        for i,v in enumerate(group_angles[gn]):
                            try:
                                gn_values.append(round(v))
                            except:
                                # Fails if there is a None (like gly doesnt have outer angle?)
                                gn_values.append("")
                    data['tab4'][gn]['angles_set'] = gn_values
                    data['tab4'][gn]['angles'] = gn_values

        # Tab 2 data generation
        # Get the relevant interactions
        data['tab2'] = {}
        # del class_pair_lookup
        # del r_pair_lookup
        if mode == "double":

            set_id = 'set1'
            aa_pair_data = data['tab2']
            interactions = list(Interaction.objects.filter(
                    interacting_pair__referenced_structure__pdb_code__index__in=[ pdb.upper() for pdb in data['pdbs1']]
                ).filter(
                    id__in=interaction_ids
                ).exclude(
                    interacting_pair__res1__generic_number=None,
                    interacting_pair__res2__generic_number=None
                ).annotate(
                    gn1=F('interacting_pair__res1__generic_number__label'),
                    gn2=F('interacting_pair__res2__generic_number__label'),
                    aa1=F('interacting_pair__res1__amino_acid'),
                    aa2=F('interacting_pair__res2__amino_acid'),
                ).values(
                    'gn1',
                    'gn2',
                    'aa1',
                    'aa2',
                ).distinct().annotate(
                    i_types=ArrayAgg('interaction_type'),
                    structures=ArrayAgg('interacting_pair__referenced_structure__pdb_code__index'),
                    pfs=ArrayAgg('interacting_pair__referenced_structure__protein_conformation__protein__parent__family__slug'),
                    structuresC=Count('interacting_pair__referenced_structure',distinct=True),
                    pfsC=Count('interacting_pair__referenced_structure__protein_conformation__protein__parent__family__name',distinct=True)
                ))
            for i in interactions:
                key = '{},{}{}{}'.format(r_class_translate_from_classA[i['gn1']],r_class_translate_from_classA[i['gn2']],i['aa1'],i['aa2'])
                if key not in aa_pair_data:
                    aa_pair_data[key] = {'classA':'{},{}'.format(i['gn1'],i['gn2']),'set1':{'interaction_freq':0,'interaction_freq_pf':0, 'types_count':defaultdict(set)}, 'set2':{'interaction_freq':0,'interaction_freq_pf':0, 'types_count':defaultdict(set)}, 'types':[]}
                aa_pair_data[key]['types'] += i['i_types']
                d = aa_pair_data[key][set_id]
                d['interaction_freq'] = round(100*i['structuresC'] / len(data['pdbs1']),0)
                d['interaction_freq_pf'] = round(100*i['pfsC'] / len(data['pfs1']),0)
                if normalized:
                    merged_types_structures = list(zip(i['i_types'],i['pfs']))
                else:
                    merged_types_structures = list(zip(i['i_types'],i['structures']))
                d['types_count'] = defaultdict(set)
                for key, val in merged_types_structures:
                    d['types_count'][key].add(val)
            print('Gotten first set occurance calcs',time.time()-start_time)

            set_id = 'set2'
            interactions = list(Interaction.objects.filter(
                    interacting_pair__referenced_structure__pdb_code__index__in=[ pdb.upper() for pdb in data['pdbs2']]
                ).filter(
                    id__in=interaction_ids
                ).exclude(
                    interacting_pair__res1__generic_number=None,
                    interacting_pair__res2__generic_number=None
                ).annotate(
                    gn1=F('interacting_pair__res1__generic_number__label'),
                    gn2=F('interacting_pair__res2__generic_number__label'),
                    aa1=F('interacting_pair__res1__amino_acid'),
                    aa2=F('interacting_pair__res2__amino_acid'),
                ).values(
                    'gn1',
                    'gn2',
                    'aa1',
                    'aa2',
                ).distinct().annotate(
                    i_types=ArrayAgg('interaction_type'),
                    structures=ArrayAgg('interacting_pair__referenced_structure__pdb_code__index'),
                    pfs=ArrayAgg('interacting_pair__referenced_structure__protein_conformation__protein__parent__family__slug'),
                    structuresC=Count('interacting_pair__referenced_structure',distinct=True),
                    pfsC=Count('interacting_pair__referenced_structure__protein_conformation__protein__parent__family__name',distinct=True)
                ))

            for i in interactions:
                key = '{},{}{}{}'.format(r_class_translate_from_classA[i['gn1']],r_class_translate_from_classA[i['gn2']],i['aa1'],i['aa2'])
                if key not in aa_pair_data:
                    aa_pair_data[key] = {'classA':'{},{}'.format(i['gn1'],i['gn2']),'set1':{'interaction_freq':0,'interaction_freq_pf':0, 'types_count':defaultdict(set)}, 'set2':{'interaction_freq':0,'interaction_freq_pf':0, 'types_count':defaultdict(set)}, 'types':[]}
                aa_pair_data[key]['types'] += i['i_types']
                d = aa_pair_data[key][set_id]
                d['interaction_freq'] = round(100*i['structuresC'] / len(data['pdbs2']),0)
                d['interaction_freq_pf'] = round(100*i['pfsC'] / len(data['pfs2']),0)
                if normalized:
                    merged_types_structures = list(zip(i['i_types'],i['pfs']))
                else:
                    merged_types_structures = list(zip(i['i_types'],i['structures']))

                for key, val in merged_types_structures:
                    d['types_count'][key].add(val)
            print('Gotten second set occurance calcs',time.time()-start_time)

            ## Fill in remaining data
            pdbs1 = data['pdbs1']
            pdbs2 = data['pdbs2']
            for key,d in aa_pair_data.items():

                gen1 = key.split(',')[0]
                gen2 = key.split(',')[1][:-2]
                d['pos_key'] = '{},{}'.format(gen1,gen2)
                aa1 = key[-2]
                aa2 = key[-1]
                d['aa1'] = aa1
                d['aa2'] = aa2

                d['types'] = list(set(d['types']))

                order = ['ionic', 'polar', 'aromatic', 'hydrophobic', 'van-der-waals','None']
                d['types'] = sorted(d['types'], key=lambda x: order.index(x))

                if key in class_pair_lookup:
                    d['class'] = class_pair_lookup[key]
                else:
                    d['class'] = ""

                #Find the individual keys for AA format 1x50,2x50A
                key_aa1 = gen1+aa1
                if key_aa1 in class_pair_lookup:
                    d['class_aa1'] = class_pair_lookup[key_aa1]
                else:
                    d['class_aa1'] = ""

                key_aa2 = gen2+aa2
                if key_aa2 in class_pair_lookup:
                    d['class_aa2'] = class_pair_lookup[key_aa2]
                else:
                    d['class_aa2'] = ""

                # Work out the occurance of interaction types in each set..
                d['types_freq'] = {}
                for i_type in ['ionic', 'polar', 'aromatic', 'hydrophobic', 'van-der-waals']:
                    set1_type_freq = round(100*len(d['set1']['types_count'][i_type])/data['set1_size'])
                    set2_type_freq = round(100*len(d['set2']['types_count'][i_type])/data['set2_size'])
                    d['types_freq'][i_type] = [set1_type_freq,
                                               set2_type_freq,
                                               set1_type_freq-set2_type_freq] #set1, #set2
                del d['set1']['types_count']
                del d['set2']['types_count']
                pdbs_with_aa1 = r_pair_lookup[gen1][aa1]
                pdbs_with_aa2 = r_pair_lookup[gen2][aa2]

                # SET 1
                pdbs1_with_aa1 = list(set(pdbs_with_aa1).intersection(pdbs1))
                pdbs1_with_aa2 = list(set(pdbs_with_aa2).intersection(pdbs1))
                # SET 2
                pdbs2_with_aa1 = list(set(pdbs_with_aa1).intersection(pdbs2))
                pdbs2_with_aa2 = list(set(pdbs_with_aa2).intersection(pdbs2))


                pdbs_intersection = list(set(pdbs_with_aa1).intersection(pdbs_with_aa2))
                pdbs1_with_pair = list(set(pdbs_intersection).intersection(pdbs1))
                pdbs2_with_pair = list(set(pdbs_intersection).intersection(pdbs2))


                if normalized:
                    pfs1_with_aa1 = list(set([ pdb_lookup[p][2] for p in pdbs1_with_aa1]))
                    pfs1_with_aa2 = list(set([ pdb_lookup[p][2] for p in pdbs1_with_aa2]))
                    pfs1_with_pair = list(set([ pdb_lookup[p][2] for p in pdbs1_with_pair]))


                    pfs2_with_aa1 = list(set([ pdb_lookup[p][2] for p in pdbs2_with_aa1]))
                    pfs2_with_aa2 = list(set([ pdb_lookup[p][2] for p in pdbs2_with_aa2]))
                    pfs2_with_pair = list(set([ pdb_lookup[p][2] for p in pdbs2_with_pair]))

                    d['set1']['occurance'] = {'aa1':pfs1_with_aa1,'aa2':pfs1_with_aa2,'pair':pfs1_with_pair}
                    d['set2']['occurance'] = {'aa1':pfs2_with_aa1,'aa2':pfs2_with_aa2,'pair':pfs2_with_pair}
                else:
                    d['set1']['occurance'] = {'aa1':pdbs1_with_aa1,'aa2':pdbs1_with_aa2,'pair':pdbs1_with_pair}
                    d['set2']['occurance'] = {'aa1':pdbs2_with_aa1,'aa2':pdbs2_with_aa2,'pair':pdbs2_with_pair}
        else:
            # Single set!
            # TODO: fix the interaction filter subselection
            aa_pair_data = data['tab2']
            interactions = list(Interaction.objects.filter(
                    id__in=interaction_ids
                ).exclude(
                    interacting_pair__res1__generic_number=None,
                    interacting_pair__res2__generic_number=None
                ).annotate(
                    gn1=F('interacting_pair__res1__generic_number__label'),
                    gn2=F('interacting_pair__res2__generic_number__label'),
                    aa1=F('interacting_pair__res1__amino_acid'),
                    aa2=F('interacting_pair__res2__amino_acid'),
                ).values(
                    'gn1',
                    'gn2',
                    'aa1',
                    'aa2',
                ).distinct().annotate(
                    i_types=ArrayAgg('interaction_type'),
                    structures=ArrayAgg('interacting_pair__referenced_structure__pdb_code__index'),
                    pfs=ArrayAgg('interacting_pair__referenced_structure__protein_conformation__protein__parent__family__slug'),
                    structuresC=Count('interacting_pair__referenced_structure',distinct=True),
                    pfsC=Count('interacting_pair__referenced_structure__protein_conformation__protein__parent__family__name',distinct=True)
                ))

            for i in interactions:
                key = '{},{}{}{}'.format(r_class_translate_from_classA[i['gn1']],r_class_translate_from_classA[i['gn2']],i['aa1'],i['aa2'])
                if key not in aa_pair_data:
                    aa_pair_data[key] = {'classA':'{},{}'.format(i['gn1'],i['gn2']),'set':{'interaction_freq':0,'interaction_freq_pf':0, 'types_count':defaultdict(set)}, 'types':[]}
                aa_pair_data[key]['types'] += i['i_types']
                d = aa_pair_data[key]['set']
                d['interaction_freq'] = round(100*i['structuresC'] / len(data['pdbs']),0)
                d['interaction_freq_pf'] = round(100*i['pfsC'] / len(data['pfs']),0)
                d['structures'] = i['structures']
                if normalized:
                    merged_types_structures = list(zip(i['i_types'],i['pfs']))
                else:
                    merged_types_structures = list(zip(i['i_types'],i['structures']))

                for key, val in merged_types_structures:
                    d['types_count'][key].add(val)

            ## Fill in remaining data
            pdbs1 = data['pdbs']

            for key,d in aa_pair_data.items():
                gen1 = key.split(',')[0]
                gen2 = key.split(',')[1][:-2]
                d['pos_key'] = '{},{}'.format(gen1,gen2)
                aa1 = key[-2]
                aa2 = key[-1]
                d['aa1'] = aa1
                d['aa2'] = aa2

                d['types'] = list(set(d['types']))

                order = ['ionic', 'polar', 'aromatic', 'hydrophobic', 'van-der-waals','None']
                d['types'] = sorted(d['types'], key=lambda x: order.index(x))

                if key in class_pair_lookup:
                    d['class'] = class_pair_lookup[key]
                else:
                    d['class'] = ""

                #Find the individual keys for AA format 1x50,2x50A
                key_aa1 = gen1+aa1
                if key_aa1 in class_pair_lookup:
                    d['class_aa1'] = class_pair_lookup[key_aa1]
                else:
                    d['class_aa1'] = ""

                key_aa2 = gen2+aa2
                if key_aa2 in class_pair_lookup:
                    d['class_aa2'] = class_pair_lookup[key_aa2]
                else:
                    d['class_aa2'] = ""

                # Work out the occurance of interaction types in each set..
                d['types_freq'] = {}
                for i_type in ['ionic', 'polar', 'aromatic', 'hydrophobic', 'van-der-waals']:
                    set_type_freq = round(100*len(d['set']['types_count'][i_type])/data['set_size'])
                    d['types_freq'][i_type] = set_type_freq #set
                del d['set']['types_count']

                pdbs_with_aa1 = r_pair_lookup[gen1][aa1]
                pdbs_with_aa2 = r_pair_lookup[gen2][aa2]

                # SET 1
                pdbs1_with_aa1 = list(set(pdbs_with_aa1).intersection(pdbs1))
                pdbs1_with_aa2 = list(set(pdbs_with_aa2).intersection(pdbs1))
                pdbs_intersection = list(set(pdbs_with_aa1).intersection(pdbs_with_aa2))
                pdbs1_with_pair = list(set(pdbs_intersection).intersection(pdbs1))

                if normalized:
                    pfs1_with_aa1 = list(set([ pdb_lookup[p][2] for p in pdbs1_with_aa1]))
                    pfs1_with_aa2 = list(set([ pdb_lookup[p][2] for p in pdbs1_with_aa2]))
                    pfs1_with_pair = list(set([ pdb_lookup[p][2] for p in pdbs1_with_pair]))

                    d['set']['occurance'] = {'aa1':pfs1_with_aa1,'aa2':pfs1_with_aa2,'pair':pfs1_with_pair}
                else:
                    d['set']['occurance'] = {'aa1':pdbs1_with_aa1,'aa2':pdbs1_with_aa2,'pair':pdbs1_with_pair}

        print('Prepare distance values for aa/gen for',mode,'mode',time.time()-start_time)
        interaction_keys = [k.replace(",","_") for k in data['interactions'].keys()]
        interaction_keys = [v['class_a_gns'].replace(",","_") for k,v in data['interactions'].items()]
        if mode == "double":


            group_1_distances = get_distance_averages(data['pdbs1'],s_lookup, interaction_keys,normalized, standard_deviation = False, split_by_amino_acid = True)

            group_2_distances = get_distance_averages(data['pdbs2'],s_lookup, interaction_keys,normalized, standard_deviation = False, split_by_amino_acid = True)

            print('got distance values for',mode,'mode',time.time()-start_time)
            for key, d in data['tab2'].items():
                class_a_key = '{}{}'.format(d['classA'], key[-2:])
                if class_a_key in group_1_distances and class_a_key in group_2_distances:
                    distance_diff = round(group_1_distances[class_a_key]-group_2_distances[class_a_key],2)
                else:
                    distance_diff = ""
                d['distance'] = distance_diff
            print('Done merging distance values for',mode,'mode',time.time()-start_time)
        else:
            group_distances = get_distance_averages(data['pdbs'],s_lookup, interaction_keys,normalized, standard_deviation = False, split_by_amino_acid = True)
            for key, d in data['tab2'].items():
                class_a_key = '{}{}'.format(d['classA'], key[-2:])
                if class_a_key in group_distances:
                    distance_diff = round(group_distances[class_a_key],2)
                else:
                    distance_diff = ""
                d['distance'] = distance_diff
        # del class_pair_lookup
        # del r_pair_lookup
        print('calculate angles per gen/aa',time.time()-start_time)
        if mode == "double":

            group_1_angles_aa = get_angle_averages(data['pdbs1'],s_lookup, normalized, standard_deviation = False, split_by_amino_acid = True, forced_class_a = forced_class_a)
            group_2_angles_aa = get_angle_averages(data['pdbs2'],s_lookup, normalized, standard_deviation = False, split_by_amino_acid = True, forced_class_a = forced_class_a)

            for key, d in data['tab2'].items():
                gen1 = key.split(',')[0]
                gen2 = key.split(',')[1][:-2]
                aa1 = key[-2]
                aa2 = key[-1]

                gn1 = '{},{}'.format(gen1,aa1)
                gn2 = '{},{}'.format(gen2,aa2)

                gn1_values = [['','','']] * 11
                if gn1 in group_1_angles_aa and gn1 in group_2_angles_aa:
                    gn1_values = []
                    for i,v in enumerate(group_1_angles_aa[gn1]):
                        try:
                            if index_names[i] in custom_angles:
                                diff = abs(v-group_2_angles_aa[gn1][i])
                                diff = round(min(360-diff,diff))
                            else:
                                diff = round(v-group_2_angles_aa[gn1][i],0)
                            gn1_values.append([diff,v,group_2_angles_aa[gn1][i]])
                        except:
                            # Fails if there is a None (like gly doesnt have outer angle?)
                            gn1_values.append(['','',''])

                gn2_values = [['','','']] * 11
                if gn2 in group_1_angles_aa and gn2 in group_2_angles_aa:
                    gn2_values = []
                    for i,v in enumerate(group_1_angles_aa[gn2]):
                        try:
                            if index_names[i] in custom_angles:
                                diff = abs(v-group_2_angles_aa[gn2][i])
                                diff = round(min(360-diff,diff))
                            else:
                                diff = round(v-group_2_angles_aa[gn2][i],0)
                            gn2_values.append([diff,v,group_2_angles_aa[gn2][i]])
                        except:
                            # Fails if there is a None (like gly doesnt have outer angle?)
                            gn2_values.append(['','',''])

                # print(key,gn1_values,gn2_values)
                d['angles'] = [gn1_values,gn2_values]

            del group_1_angles_aa
            del group_2_angles_aa
        else:

            group_angles_aa = get_angle_averages(data['pdbs'],s_lookup, normalized, standard_deviation = True, split_by_amino_acid = True, forced_class_a = forced_class_a)

            for key, d in data['tab2'].items():
                gen1 = key.split(',')[0]
                gen2 = key.split(',')[1][:-2]
                aa1 = key[-2]
                aa2 = key[-1]

                gn1 = '{},{}'.format(gen1,aa1)
                gn2 = '{},{}'.format(gen2,aa2)

                gn1_values = [''] * 11
                if gn1 in group_angles_aa:
                    gn1_values = []
                    for i,v in enumerate(group_angles_aa[gn1]):
                        try:
                            gn1_values.append(round(v,1))
                        except:
                            # Fails if there is a None (like gly doesnt have outer angle?)
                            gn1_values.append("")

                gn2_values = [''] * 11
                if gn2 in group_angles_aa:
                    gn2_values = []
                    for i,v in enumerate(group_angles_aa[gn2]):
                        try:
                            gn2_values.append(round(v,1))
                        except:
                            # Fails if there is a None (like gly doesnt have outer angle?)
                            gn2_values.append("")
                d['angles'] = [gn1_values,gn2_values]



        #print(data['tab3'])
        print('calculate tab4',time.time()-start_time)
        del aa_pair_data
        del all_pdbs_pairs
        # del ds
        for res1, d in data['tab4'].items():
            for key, values in d.items():
                if type(values) is set:
                    data['tab4'][res1][key] = list(values)

            if mode == "double":
                # set_1_avg_freq = 0
                # set_2_avg_freq = 0
                # for s in ['1','2']:
                #     # Need to figure all the individual frequencies by looking in main dictionary
                #     count = len(d['set{}'.format(s)])
                #     running_sum = 0
                #     for res2 in d['set{}'.format(s)]:
                #         pair = '{},{}'.format(res1,res2)
                #         pair_reverse = '{},{}'.format(res2,res1)
                #         if pair in data['interactions']:
                #             len_pdbs = len(data['interactions'][pair]['pdbs{}'.format(s)])
                #         if pair_reverse in data['interactions']:
                #             len_pdbs = len(data['interactions'][pair_reverse]['pdbs{}'.format(s)])
                #         running_sum += len_pdbs/len(data['pdbs{}'.format(s)])
                #     if count:
                #         avg_freq = running_sum/count
                #         if s == '1':
                #             set_1_avg_freq = avg_freq
                #         else:
                #             set_2_avg_freq = avg_freq

                # absolute_diff = abs(set_1_avg_freq-set_2_avg_freq)
                # # print(res1,set_1_avg_freq,set_2_avg_freq,absolute_diff)
                # data['tab4'][res1]['avg_freq_diff_sets'] = absolute_diff


                # Figure out most frequent AA and % in set
                cons_aa_set1 = ['','']
                cons_aa_set2 = ['','']
                aa_at_pos = r_pair_lookup[res1]
                temp_score_dict = []
                for aa, pdbs in aa_at_pos.items():

                    pdbs1_with_aa = list(set(pdbs).intersection(pdbs1))
                    pdbs2_with_aa = list(set(pdbs).intersection(pdbs2))
                    temp_score_dict.append([aa,len(pdbs1_with_aa),len(pdbs2_with_aa), len(pdbs)])

                most_freq_set1 = sorted(temp_score_dict.copy(), key = lambda x: -x[1])
                most_freq_set2 = sorted(temp_score_dict.copy(), key = lambda x: -x[2])
                most_freq_combi = sorted(temp_score_dict.copy(), key = lambda x: -x[3])
                data['tab4'][res1]['set1_seq_cons'] = most_freq_set1[0]
                data['tab4'][res1]['set2_seq_cons'] = most_freq_set2[0]
                data['tab4'][res1]['all_seq_cons'] = most_freq_combi[0]

                if res1 in group_1_angles:
                    data['tab4'][res1]['angles_set1'] = [ round(elem) if isinstance(elem, float) else '' for elem in group_1_angles[res1] ]
                else:
                    data['tab4'][res1]['angles_set1'] = [''] * 11
                if res1 in group_2_angles:
                    data['tab4'][res1]['angles_set2'] = [ round(elem) if isinstance(elem, float) else '' for elem in group_2_angles[res1] ]
                else:
                    data['tab4'][res1]['angles_set2'] = [''] * 11
                # Get angle data for res1
                res1_values = [['','','']] * 11
                if res1 in group_1_angles and res1 in group_2_angles:
                    res1_values = []
                    for i,v in enumerate(group_1_angles[res1]):
                        try:
                            if index_names[i] in custom_angles:
                                diff = abs(v-group_2_angles[res1][i])
                                diff = round(min(360-diff,diff))
                            else:
                                diff = round(v-group_2_angles[res1][i],0)
                            res1_values.append([diff,v,group_2_angles[res1][i]])
                        except:
                            # Fails if there is a None (like gly doesnt have outer angle?)
                            res1_values.append(['','',''])
                data['tab4'][res1]['angles'] = res1_values
            else:
                #TODO SINGLE SET
                aa_at_pos = r_pair_lookup[res1]
                temp_score_dict = []
                for aa, pdbs in aa_at_pos.items():

                    pdbs1_with_aa = list(set(pdbs).intersection(pdbs1))
                    temp_score_dict.append([aa,len(pdbs1_with_aa)])

                most_freq_set = sorted(temp_score_dict.copy(), key = lambda x: -x[1])
                data['tab4'][res1]['set_seq_cons'] = most_freq_set[0]
            # Common for all modes
            if res1 in class_pair_lookup:
                # print(res1,class_pair_lookup[res1])
                data['tab4'][res1]['class_cons'] = class_pair_lookup[res1]
            else:
                data['tab4'][res1]['class_cons'] = ['','']
                print('no res1',res1,'in class lookup')

        # calculate information for 2D helical displacement plot
        if mode == "double":
            pdbs1_upper = [pdb.upper() for pdb in pdbs1]
            pdbs2_upper = [pdb.upper() for pdb in pdbs2]
            helical_time = time.time()
            print("Start helical movements")
            data['tm_movement_2D'] = {}
            # data['tm_movement_2D']["classA_ligands"] = tm_movement_2D(pdbs1_upper, pdbs2_upper, 2, data, r_class_translate_from_classA)
            data['tm_movement_2D']["membrane_mid"] = tm_movement_2D(pdbs1_upper, pdbs2_upper, 3, data, r_class_translate_from_classA)
            data['tm_movement_2D']["intracellular"] = tm_movement_2D(pdbs1_upper, pdbs2_upper, 1, data, r_class_translate_from_classA)
            data['tm_movement_2D']["extracellular"] = tm_movement_2D(pdbs1_upper, pdbs2_upper, 0, data, r_class_translate_from_classA)

            # viewbox
            diff_x = 0
            diff_y = 0
            for x in data['tm_movement_2D']: # 2D set
                setx = [z["x"] for y in ["coordinates_set1", "coordinates_set2"] for z in data['tm_movement_2D'][x][y]]
                sety = [z["y"] for y in ["coordinates_set1", "coordinates_set2"] for z in data['tm_movement_2D'][x][y]]

                if diff_x < (max(setx) - min(setx)):
                    diff_x = max(setx) - min(setx)

                if diff_y < (max(sety) - min(sety)):
                    diff_y = max(sety) - min(sety)

            data['tm_movement_2D']["viewbox_size"] = {"diff_x": diff_x, "diff_y" : diff_y}

            print("Helical movement calculations", time.time() - helical_time)

        # calculate distance movements
        if mode == "double":
            pdbs1_upper = [pdb.upper() for pdb in pdbs1]
            pdbs2_upper = [pdb.upper() for pdb in pdbs2]
            print('Start 1')
            start = time.time()
            dis1 = Distances()
            dis1.load_pdbs(pdbs1_upper)
            dis1.fetch_and_calculate(with_arr = True)
            print('done fetching set 1',time.time()-start)
            # dis1.calculate_window(list_of_gns)
        #    dis1.calculate()

            start = time.time()
            dis2 = Distances()
            dis2.load_pdbs(pdbs2_upper)
            dis2.fetch_and_calculate(with_arr = True)
            # dis2.calculate_window(list_of_gns)
            #dis2.calculate()
            print('done fetching set 2',time.time()-start)

            diff = OrderedDict()
            from math import sqrt
            from scipy import stats

            # for d1 in dis1.stats_window_reduced:
            # for d1 in dis1.stats_window:

            start = time.time()
            total = {}
            common_labels = list(set(dis1.stats_key.keys()).intersection(dis2.stats_key.keys()))
            for label in common_labels:

                # Get variables
                d1, d2 = dis1.stats_key[label],dis2.stats_key[label]
                # Correct decimal
                mean1, mean2 = d1[1],d2[1]
                std1,std2 = d1[2],d2[2]
                var1,var2 = std1**2,std2**2
                n1, n2 = d1[4],d2[4]

                mean_diff = mean1-mean2

                # Save values for NGL calcs
                gn1,gn2 = label.split("_")
                if gn1 not in total:
                    total[gn1] = {}
                if gn2 not in total:
                    total[gn2] = {}
                total[gn1][gn2] = total[gn2][gn1] = round(mean_diff,1)
                # # Make easier readable output
                # individual_pdbs_1 = dict(zip(d1[6], d1[5]))
                # individual_pdbs_2 = dict(zip(d2[6], d2[5]))

                # if n1>1 and n2>1 and var1>0 and var2>0:
                #     ## T test to assess seperation of data (only if N>1 and there is variance)
                #     t_stat_welch = abs(mean1-mean2)/(sqrt( (var1/n1) + (var2/n2)  ))
                #     df = n1+n2 - 2
                #     p = 1 - stats.t.cdf(t_stat_welch,df=df)
                # else:
                #     p = 0

                # diff[label] = [round(mean_diff,1),[std1,std2],[mean1,mean2],[n1,n2],p,[individual_pdbs_1,individual_pdbs_2]]

            # diff =  OrderedDict(sorted(diff.items(), key=lambda t: -abs(t[1][0])))
            # print(diff)
            print('done diff',time.time()-start)
            start = time.time()
            ngl_max_diff = 0
            for gn1 in total.keys():
                vals = []
                for gn,val in total[gn1].items():
                    if gn[0]!=gn1[0]:
                        #if not same segment
                        vals.append(val)
                total[gn1]['avg'] = round(float(sum(vals))/max(len(vals),1),1)
                if abs(total[gn1]['avg'])>ngl_max_diff:
                    ngl_max_diff = round(abs(total[gn1]['avg']),1)

            print('done ngl', time.time() - start)
            data['distances'] = total
            data['ngl_max_diff_distance'] = ngl_max_diff

        data['tab3'] = {}
        data['pdbs'] = list(data['pdbs'])
        data['proteins'] = list(data['proteins'])
        data['pfs'] = list(data['pfs'])
        data['pfs_lookup'] = dict(data['pfs_lookup'])

        data['segm_lookup'] = segm_lookup
        data['segm_lookup_ordered'] = sorted(segm_lookup, key=functools.cmp_to_key(gpcrdb_number_comparator))
        data['segments'] = list(data['segments'])
        data['normalized'] = normalized
        data['forced_class_a'] = forced_class_a
        data['residue_table'] = r_class_translate


        excluded_segment = ['C-term','N-term'] #'ICL1','ECL1','ECL2','ICL2'
        segments = ProteinSegment.objects.all().exclude(slug__in = excluded_segment)
        proteins =  Protein.objects.filter(entry_name__in=list(data['proteins'])).distinct().all()
        a = Alignment()
        a.ignore_alternative_residue_numbering_schemes = True
        a.load_proteins(proteins)
        a.load_segments(segments) #get all segments to make correct diagrams
        # build the alignment data matrix
        a.build_alignment()
        # calculate consensus sequence + amino acid and feature frequency
        a.calculate_statistics()
        consensus = a.full_consensus
        data['snakeplot_lookup'] = {}
        data['snakeplot_lookup_aa'] = {}
        data['snakeplot_lookup_aa_cons'] = {}
        for a in consensus:
            if a.display_generic_number:
                # Be sure to use the correct GN
                gn = re.sub(r'\.[\d]+', '', a.display_generic_number.label)
                if forced_class_a:
                    gn = a.family_generic_number
                a.display_generic_number.label = gn
                data['snakeplot_lookup'][gn] = a.sequence_number
                data['snakeplot_lookup_aa'][gn] = a.amino_acid
                gen_aa = gn + a.amino_acid
                if gen_aa in class_pair_lookup:
                    data['snakeplot_lookup_aa_cons'][gn] = class_pair_lookup[gen_aa]
                else:
                    data['snakeplot_lookup_aa_cons'][gn] = 0
        from common.diagrams_gpcr import DrawSnakePlot
        snakeplot = DrawSnakePlot(consensus, 'Class A', 'family_diagram_preloaded_data',nobuttons = True)
        data['snakeplot'] = str(snakeplot)
        if mode == 'double':
            data['pdbs1'] = list(data['pdbs1'])
            data['pdbs2'] = list(data['pdbs2'])
            data['proteins1'] = list(data['proteins1'])
            data['proteins2'] = list(data['proteins2'])
            data['pfs1'] = list(data['pfs1'])
            data['pfs2'] = list(data['pfs2'])
            data['pfs1_lookup'] = dict(data['pfs1_lookup'])
            data['pfs2_lookup'] = dict(data['pfs2_lookup'])
        else:
            data['pdbs'] = list(data['pdbs'])
        cache.set(hash_cache_key,data,3600*24)
    print('Done',time.time()-start_time)

    return JsonResponse(data)

def DistanceDataGroups(request):
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
        pdbs1 = request.GET.getlist('pdbs1[]')
        pdbs2 = request.GET.getlist('pdbs2[]')
    except IndexError:
        pdbs1 = []



    pdbs1 = [pdb.upper() for pdb in pdbs1]
    pdbs2 = [pdb.upper() for pdb in pdbs2]
    pdbs1_lower = [pdb.lower() for pdb in pdbs1]
    pdbs2_lower = [pdb.lower() for pdb in pdbs2]


    cache_key = ",".join(sorted(pdbs1)) + "_" + ",".join(sorted(pdbs2))
    cache_key = hashlib.md5(cache_key.encode('utf-8')).hexdigest()

    data = cache.get(cache_key)
    # data = None
    if data!=None:
        print('Result cached')
        return JsonResponse(data)

    # Segment filters
    try:
        segments = request.GET.getlist('segments[]')
    except IndexError:
        segments = []

    # Use generic numbers? Defaults to True.
    generic = True

    # Initialize response dictionary
    data = {}
    data['interactions'] = OrderedDict()
    data['pdbs'] = set()
    data['generic'] = generic
    data['segments'] = set()
    data['segment_map'] = {}
    # For Max schematics TODO -- make smarter.
    data['segment_map'] = {}
    data['aa_map'] = {}


    data['gn_map'] = OrderedDict()
    data['pos_map'] = OrderedDict()
    data['segment_map_full'] = OrderedDict()
    data['segment_map_full_gn'] = OrderedDict()
    data['generic_map_full'] = OrderedDict()


    list_of_gns = []

    excluded_segment = ['C-term','N-term','ICL1','ECL1','ECL2','ICL2','H8']
    segments = ProteinSegment.objects.all().exclude(slug__in = excluded_segment)
    proteins =  Protein.objects.filter(protein__entry_name__in=pdbs1_lower+pdbs2_lower).distinct().all()
    if len(proteins)>1:
        a = Alignment()
        a.ignore_alternative_residue_numbering_schemes = True;
        a.load_proteins(proteins)
        a.load_segments(segments) #get all segments to make correct diagrams
        # build the alignment data matrix
        a.build_alignment()
        # calculate consensus sequence + amino acid and feature frequency
        a.calculate_statistics()
        consensus = a.full_consensus

        for aa in consensus:
            if 'x' in aa.family_generic_number:
                list_of_gns.append(aa.family_generic_number)
                data['gn_map'][aa.family_generic_number] = aa.amino_acid
                data['pos_map'][aa.sequence_number] = aa.amino_acid
                data['segment_map_full_gn'][aa.family_generic_number] = aa.segment_slug
    else:
        rs = Residue.objects.filter(protein_conformation__protein=proteins[0]).prefetch_related('protein_segment','display_generic_number','generic_number')
        for r in rs:
            if (not generic):
                data['pos_map'][r.sequence_number] = r.amino_acid
                data['segment_map_full'][r.sequence_number] = r.protein_segment.slug
                if r.display_generic_number:
                    data['generic_map_full'][r.sequence_number] = r.short_display_generic_number()
            else:
                if r.generic_number:
                    list_of_gns.append(r.generic_number.label)
                    data['gn_map'][r.generic_number.label] = r.amino_acid
                    data['pos_map'][r.sequence_number] = r.amino_acid
                    data['segment_map_full_gn'][r.generic_number.label] = r.protein_segment.slug

    print('Start 1')
    start = time.time()
    dis1 = Distances()
    dis1.load_pdbs(pdbs1)
    dis1.fetch_and_calculate(with_arr = True)
    print('done fetching set 1',time.time()-start)
    # dis1.calculate_window(list_of_gns)
#    dis1.calculate()

    start = time.time()
    dis2 = Distances()
    dis2.load_pdbs(pdbs2)
    dis2.fetch_and_calculate(with_arr = True)
    # dis2.calculate_window(list_of_gns)
    #dis2.calculate()
    print('done fetching set 2',time.time()-start)

    diff = OrderedDict()
    from math import sqrt
    from scipy import stats

    # for d1 in dis1.stats_window_reduced:
    # for d1 in dis1.stats_window:

    start = time.time()
    total = {}
    common_labels = list(set(dis1.stats_key.keys()).intersection(dis2.stats_key.keys()))
    for label in common_labels:

        # Get variables
        d1, d2 = dis1.stats_key[label],dis2.stats_key[label]
        # Correct decimal
        mean1, mean2 = d1[1]/100,d2[1]/100
        std1,std2 = d1[2]/100,d2[2]/100
        var1,var2 = std1**2,std2**2
        n1, n2 = d1[4],d2[4]

        mean_diff = mean2-mean1

        # Save values for NGL calcs
        gn1,gn2 = label.split("_")
        if gn1 not in total:
            total[gn1] = {}
        if gn2 not in total:
            total[gn2] = {}
        total[gn1][gn2] = total[gn2][gn1] = round(mean_diff,1)
        # Make easier readable output
        d1[5] = [x / 100 for x in d1[5]]
        d2[5] = [x / 100 for x in d2[5]]
        individual_pdbs_1 = dict(zip(d1[6], d1[5]))
        individual_pdbs_2 = dict(zip(d2[6], d2[5]))

        if n1>1 and n2>1 and var1>0 and var2>0:
            ## T test to assess seperation of data (only if N>1 and there is variance)
            t_stat_welch = abs(mean1-mean2)/(sqrt( (var1/n1) + (var2/n2)  ))
            df = n1+n2 - 2
            p = 1 - stats.t.cdf(t_stat_welch,df=df)
        else:
            p = 0

        diff[label] = [round(mean_diff,1),[std1,std2],[mean1,mean2],[n1,n2],p,[individual_pdbs_1,individual_pdbs_2]]

    diff =  OrderedDict(sorted(diff.items(), key=lambda t: -abs(t[1][0])))

    print('done diff',time.time()-start)
    compared_stats = {}

    # Remove differences that seem statistically irrelevant
    for label,d in diff.items():
        if d[4]>0.05:
            #print(label,d)
            d[0] = 0
    print('done')

    # Dict to keep track of which residue numbers are in use
    number_dict = set()
    max_diff = 0
    #for d in dis.stats:
    for key,d in diff.items():

        # Ignore low means diff.
        if d[0]<0.5 and d[0]>-0.5:
            continue
        res1 = key.split("_")[0]
        res2 = key.split("_")[1]
        res1_seg = res1.split("x")[0]
        res2_seg = res2.split("x")[0]
        data['segment_map'][res1] = res1_seg
        data['segment_map'][res2] = res2_seg
        data['segments'] |= {res1_seg} | {res2_seg}

        # Populate the AA map
        if 1==2:
            #When showing pdbs
            if pdb_name not in data['aa_map']:
                data['aa_map'][pdb_name] = {}

        number_dict |= {res1, res2}

        if res1 < res2:
            coord = str(res1) + ',' + str(res2)
        else:
            coord = str(res2) + ',' + str(res1)

        # if pdb_name not in data['interactions'][coord]:
        #     data['interactions'][coord][pdb_name] = []

        if d:
            if len(data['interactions'])<10000:

                data['interactions'][coord] = d
            else:
                break

        if abs(d[0])>max_diff:
            max_diff = round(abs(d[0]))
        # data['sequence_numbers'] = sorted(number_dict)
    if (generic):
        data['sequence_numbers'] = sorted(number_dict, key=functools.cmp_to_key(gpcrdb_number_comparator))
    # else:
    #     data['sequence_numbers'] = sorted(number_dict)
        # break

    data['segments'] = list(data['segments'])
    data['pdbs'] = list(pdbs1+pdbs2)
    data['pdbs1'] = list(pdbs1)
    data['pdbs2'] = list(pdbs2)
    data['max_diff'] = max_diff
    print(len(data['interactions']),'len data')
    #total = {}
    start = time.time()
    ngl_max_diff = 0
    for gn1 in total.keys():
        vals = []
        for gn,val in total[gn1].items():
            if gn[0]!=gn1[0]:
                #if not same segment
                vals.append(val)
        total[gn1]['avg'] = round(float(sum(vals))/max(len(vals),1),1)
        if abs(total[gn1]['avg'])>ngl_max_diff:
            ngl_max_diff = round(abs(total[gn1]['avg']),1)

    print('done ngl',time.time()-start)
    data['ngl_data'] = total
    data['ngl_max_diff'] = ngl_max_diff
    print('send json')
    cache.set(cache_key,data,3600*24*7)
    return JsonResponse(data)


def originMatrix(pdbs):
    # select all TM7 distances to core
    ds = list(ResidueAngle.objects.filter(structure__pdb_code__index__in=pdbs)\
                        .exclude(residue__generic_number=None) \
                        .exclude(core_distance=None) \
                        .values('structure__pdb_code__index', 'residue__generic_number__label', 'mid_distance'))

    # create dictionary of all structures and all distances
    core_distances = {}
    for i,d in enumerate(ds):
        if not d['structure__pdb_code__index'] in core_distances:
            core_distances[d['structure__pdb_code__index']] = {}
        core_distances[d['structure__pdb_code__index']][d['residue__generic_number__label']] = d['mid_distance']

    pdbs = list(core_distances.keys())

    distance_matrix = np.full((len(pdbs), len(pdbs)), 0.0)
    for i, pdb1 in enumerate(pdbs):
        for j in range(i+1, len(pdbs)):
            pdb2 = pdbs[j]

            # Get common GNs between two PDBs
            common_between_pdbs = sorted(list(set(dict.keys(core_distances[pdb1])).intersection(core_distances[pdb2])))
            # Get distance between cells that have both GNs.
            distance = np.sum([ np.abs(core_distances[pdb1][key] - core_distances[pdb2][key]) for key in common_between_pdbs])
            # normalize
            distance_matrix[i, j] = pow(distance,2)/(len(common_between_pdbs)*len(common_between_pdbs))
            distance_matrix[j, i] = distance_matrix[i, j]

    return [distance_matrix, pdbs]


def coreMatrix(pdbs, core = True, middle = False):
    # select all TM7 distances to core
    ds = list(ResidueAngle.objects.filter(structure__pdb_code__index__in=pdbs)\
                        .exclude(residue__generic_number=None) \
                        .exclude(core_distance=None) \
                        .values('structure__pdb_code__index', 'residue__generic_number__label', 'core_distance', 'midplane_distance'))

    # create dictionary of all structures and all distances
    core_distances = {}
    for i,d in enumerate(ds):
        if not d['structure__pdb_code__index'] in core_distances:
            core_distances[d['structure__pdb_code__index']] = {}
        if core:
            core_distances[d['structure__pdb_code__index']][d['residue__generic_number__label']+"_core"] = d['core_distance']
        if middle:
            core_distances[d['structure__pdb_code__index']][d['residue__generic_number__label']+"_mid"] = d['midplane_distance']

    # IN some cases the 7TM distances are missing e.g. due to missing (structure or annotation) of a TM bundle
    #pdbs_present = list(ResidueAngle.objects.filter(structure__pdb_code__index__in=pdbs)\
    #                    .exclude(residue__generic_number=None) \
    #                    .exclude(core_distance=None) \
    #                    .distinct('structure__pdb_code__index').values('structure__pdb_code__index'))

    #pdbs = [pdb['structure__pdb_code__index'] for pdb in pdbs_present]

    pdbs = list(core_distances.keys())

    distance_matrix = np.full((len(pdbs), len(pdbs)), 0.0)
    for i, pdb1 in enumerate(pdbs):
        for j in range(i+1, len(pdbs)):
            pdb2 = pdbs[j]

            # Get common GNs between two PDBs
            common_between_pdbs = sorted(list(set(dict.keys(core_distances[pdb1])).intersection(core_distances[pdb2])))
            # Get distance between cells that have both GNs.
            distance = np.sum([ np.abs(core_distances[pdb1][key] - core_distances[pdb2][key]) for key in common_between_pdbs])
            # normalize
            distance_matrix[i, j] = pow(distance,2)/(len(common_between_pdbs)*len(common_between_pdbs))
            distance_matrix[j, i] = distance_matrix[i, j]

    return [distance_matrix, pdbs]


def stableResMatrix(pdbs):
    # all classes
    # classes = ['001', '002', '003', '004', '006']
    # results = []
    # for selclass in classes:
    #     numStructs = Distance.objects.filter(structure__protein_conformation__protein__family__slug__startswith=selclass).values('structure_id').distinct().count()
    #     conformations = Distance.objects.filter(structure__protein_conformation__protein__family__slug__startswith=selclass).values('structure__protein_conformation')
    #
    #     # All GNs with at least 90% presence in all structures of this class
    #     structure_gn = Residue.objects.filter(protein_conformation__in=conformations) \
    #         .exclude(generic_number=None) \
    #         .exclude(generic_number__label__startswith='8x') \
    #         .exclude(generic_number__label__startswith='12x') \
    #         .exclude(generic_number__label__startswith='23x') \
    #         .exclude(generic_number__label__startswith='34x') \
    #         .exclude(generic_number__label__startswith='45x') \
    #         .values('generic_number__label') \
    #         .annotate(c = Count('protein_conformation', distinct=True)) \
    #         .order_by('generic_number__label') \
    #         .filter(c__gte=int(numStructs*0.9))
    #
    #     common_gn = [ entry["generic_number__label"] for entry in structure_gn ]
    #
    #     structure_gn = Residue.objects.filter(protein_conformation__in=conformations) \
    #         .exclude(generic_number=None) \
    #         .exclude(generic_number__label__startswith='8x') \
    #         .exclude(generic_number__label__startswith='12x') \
    #         .exclude(generic_number__label__startswith='23x') \
    #         .exclude(generic_number__label__startswith='34x') \
    #         .exclude(generic_number__label__startswith='45x') \
    #         .values('generic_number__label') \
    #         .annotate(c = Count('protein_conformation', distinct=True)) \
    #         .order_by('generic_number__label') \
    #         .filter(c__gte=int(numStructs))
    #
    #     always_gn = [ entry["generic_number__label"] for entry in structure_gn ]
    #
    #     print("Class {} with {} structures".format(selclass, numStructs))
    #     print(common_gn)
    #     print(always_gn)
    #
    #     for gn in always_gn:
    #         # collect all distances merging pairs => stdev + average
    ##         ds = list(Distance.objects.filter(structure__protein_conformation__protein__family__slug__startswith=selclass) \
    ##                         .filter(gns_pair__contains=gn)
    ##                         .filter(gn1__in=common_gn).filter(gn2__in=common_gn) \
    ##                         .values('gns_pair') \
    ##                         .annotate(mean = Avg('distance'), std = StdDev('distance'), c = Count('distance')))
    ##
    ##         # calculate average variation
    ##         totalnorm = sum([ entry['std']/entry['mean'] for entry in ds])/len(ds)
    ##         totalstd = sum([ entry['std'] for entry in ds])/len(ds)
    ##         print("{} has a variation of {} - {}".format(gn, totalnorm, totalstd))
    ##         results.append([selclass, gn, totalnorm, totalstd])
    #
    # print(results)


    # Most stable residues from each class
    # stable_residues = {'001':'3x53', '002':'1x44', '003':'', '004':'5x42', '006':'1x29'}
    stable_residues = {'001':'4x50', '002':'1x44', '003':'', '004':'5x42', '006':'1x29'}

    stable_distances = {}
    pdb_classes = {}
    for selclass in ['001', '002', '003', '004', '006']:
        # select all distances to selected residue
        reference = stable_residues[selclass]
        ds = list(Distance.objects.filter(structure__pdb_code__index__in=pdbs) \
                                .filter(structure__protein_conformation__protein__family__slug__startswith=selclass) \
                                .filter(Q(gn1=reference) | Q(gn2=reference)) \
                                .values('structure__pdb_code__index', 'gns_pair', 'distance'))

        # create dictionary of all structures and all distances
        for i,d in enumerate(ds):
            if not d['structure__pdb_code__index'] in stable_distances:
                stable_distances[d['structure__pdb_code__index']] = {}
                pdb_classes[d['structure__pdb_code__index']] = selclass

            gn_label = d['gns_pair'].replace(reference, "").replace("_", "")
            stable_distances[d['structure__pdb_code__index']][gn_label] = d['distance']


    pdbs = list(stable_distances.keys())

    distance_matrix = np.full((len(pdbs), len(pdbs)), 0.0)
    for i, pdb1 in enumerate(pdbs):
        for j in range(i+1, len(pdbs)):
            pdb2 = pdbs[j]

            if pdb_classes[pdb1] == pdb_classes[pdb2]:
                # Get common GNs between two PDBs
                common_between_pdbs = sorted(list(set(dict.keys(stable_distances[pdb1])).intersection(stable_distances[pdb2])))
                # Get distance between cells that have both GNs.
                distance = np.sum([ np.abs(stable_distances[pdb1][key] - stable_distances[pdb2][key]) for key in common_between_pdbs])
                # normalize
                distance_matrix[i, j] = pow(distance,2)/(len(common_between_pdbs)*len(common_between_pdbs))
                distance_matrix[j, i] = distance_matrix[i, j]
            else:
                # Comparison accross classes - set very high distance
                distance_matrix[i, j] = 100000
                distance_matrix[j, i] = distance_matrix[i, j]

    return [distance_matrix, pdbs]

def ClusteringData(request):
    # PDB files
    try:
        pdbs = request.GET.get('pdbs').split(',')
    except IndexError:
        pdbs = []

    if len(pdbs) == 0:
        quit()

    # Exclude structures with incomplete GN annotation
    excluded_pdbs = [ "5LWE", "5ZKP" ]
    pdbs = [pdb.upper() for pdb in pdbs if pdb.upper() not in excluded_pdbs ]

    cluster_method = 0
    if 'cluster-method' in request.GET:
        cluster_method = request.GET.get('cluster-method')

    cache_key = ",".join(sorted(pdbs))
    cache_key = "structure_clustering_" + cluster_method + "_" + hashlib.md5(cache_key.encode('utf-8')).hexdigest()

    data = cache.get(cache_key)
    if data == None:
        # output dictionary
        data = {}

        # load all
        # DEBUG set clustering method hardcoded:
        # cluster_method = '7'
        if cluster_method == '1':
            [distance_matrix, pdbs] = coreMatrix(pdbs)
        elif cluster_method == '2':
            [distance_matrix, pdbs] = stableResMatrix(pdbs) # replace with distance to most stable residue
        elif cluster_method == '3':
            dis = Distances()
            dis.load_pdbs(pdbs)
            distance_matrix = dis.get_distance_matrix(normalize = False)

            # pdbs have been reordered -> map back to be consistent with the distance matrix
            pdbs = dis.pdbs
        elif cluster_method == '4': # distance to membrane mid
            [distance_matrix, pdbs] = coreMatrix(pdbs, middle = True, core = False)
        elif cluster_method == '5': # distance to membrane mid and 7TM axis
            [distance_matrix, pdbs] = coreMatrix(pdbs, middle = True)
        elif cluster_method == '6': # distance to origin
            [distance_matrix, pdbs] = originMatrix(pdbs)
        elif cluster_method == '7': # distance to origin
            dis = Distances()
            dis.filtered_gns = True
            print("Print setting lower only")
            dis.load_pdbs(pdbs)
            distance_matrix = dis.get_distance_matrix(normalize = False)

            # pdbs have been reordered -> map back to be consistent with the distance matrix
            pdbs = dis.pdbs
        else:
            dis = Distances()
            dis.load_pdbs(pdbs)
            distance_matrix = dis.get_distance_matrix()

            # pdbs have been reordered -> map back to be consistent with the distance matrix
            pdbs = dis.pdbs

        # Collect structure annotations
        pdb_annotations = {}

        # Grab all annotations and all the ligand role when present in aggregates
        # NOTE: we can probably remove the parent step and go directly via family in the query
        annotations = Structure.objects.filter(pdb_code__index__in=pdbs) \
                        .values_list('pdb_code__index','state__slug','protein_conformation__protein__parent__entry_name','protein_conformation__protein__parent__name','protein_conformation__protein__parent__family__parent__name', \
                        'protein_conformation__protein__parent__family__parent__parent__name', 'protein_conformation__protein__parent__family__parent__parent__parent__name', 'structure_type__name', 'protein_conformation__protein__family__slug', 'tm6_angle', 'gprot_bound_likeness')\
                        .annotate(arr=ArrayAgg('structureligandinteraction__ligand_role__slug', filter=Q(structureligandinteraction__annotated=True))) \

        # Adding signaling protein data on top
        signaling_proteins = {}
        common_signaling = {'Gi1': "Gi/o", 'Gi2': "Gi/o", 'Go': "Gi/o", 'Gq': "Gq/11", 'G11': "Gq/11", 'Gs': "Gs", 'G12': "G12/13", 'G13': "G12/13", 'Gt1': "Gt", 'Gt3': "Gt",
                            'GPa1': "Gpa1", 'Beta-arrestin-1': "Arrestin", 'S-arrestin':  "Arrestin"}
        signal_ps = StructureExtraProteins.objects.filter(structure__pdb_code__index__in=pdbs, category__in=["G alpha", "Arrestin"]).values('structure__pdb_code__index','display_name', 'wt_protein__family__name').order_by().annotate(coverage = Max('wt_coverage'))
        for ps in signal_ps:
            if not ps["structure__pdb_code__index"] in signaling_proteins:
                if ps["display_name"] in common_signaling:
                    signaling_proteins[ps["structure__pdb_code__index"]] = common_signaling[ps['display_name']]
                else:
                    signaling_proteins[ps["structure__pdb_code__index"]] = ps['display_name']

        # Check for GRK complexes
        grk_complexes = list(Structure.objects.filter(stabilizing_agents__name__contains="GRK").exclude(structure_type__slug__startswith='af-').values_list("pdb_code__index", flat = True))
        for pdb in grk_complexes:
            if not pdb in signaling_proteins:
                signaling_proteins[pdb] = "GRK"

        protein_slugs = set()
        for an in annotations:
            pdb_annotations[an[0]] = list(an[1:])

            # add slug to lists
            slug = pdb_annotations[an[0]][7]
            protein_slugs.add(slug)

            # UGLY needs CLEANUP in data - replace agonist-partial with partial-agonist ()
            pdb_annotations[an[0]][10] = ["partial-agonist" if x=="agonist-partial" else x for x in pdb_annotations[an[0]][10]]

            # Cleanup the aggregates as None values are introduced
            pdb_annotations[an[0]][7] = list(filter(None.__ne__, pdb_annotations[an[0]][10]))

            # SUPERUGLY - replace
            holder = pdb_annotations[an[0]][8]
            holder2 = pdb_annotations[an[0]][9]
            pdb_annotations[an[0]][8] = slug
            pdb_annotations[an[0]][9] = holder
            pdb_annotations[an[0]][10] = holder2

            if an[0] in signaling_proteins:
                pdb_annotations[an[0]].append(signaling_proteins[an[0]])
            else:
                pdb_annotations[an[0]].append("")

        data['annotations'] = pdb_annotations

        # Grab G-protein coupling profile for all receptors covered by the selection
        # TODO: make general cache(s) for these kinds of data
        selectivitydata = {}
        coupling = ProteinCouplings.objects.filter(protein__family__slug__in=protein_slugs, source="GuideToPharma").values_list('protein__family__slug', 'transduction').annotate(arr=ArrayAgg('g_protein__name'))

        for pairing in coupling:
            if pairing[0] not in selectivitydata:
                selectivitydata[pairing[0]] = {}
            selectivitydata[pairing[0]][pairing[1]] = pairing[2]

        data['Gprot_coupling'] = selectivitydata

        # hierarchical clustering
        hclust = sch.linkage(ssd.squareform(distance_matrix), method='average')
        tree = sch.to_tree(hclust, False)

        #inconsistency = sch.inconsistent(hclust)
        #inconsistency = sch.maxinconsts(hclust, inconsistency)
        silhouette_coefficient = {}
        getSilhouetteIndex(tree, distance_matrix, silhouette_coefficient)
        data['tree'] = getNewick(tree, "", tree.dist, pdbs, silhouette_coefficient)

        # Order distance_matrix by hclust
        N = len(distance_matrix)
        res_order = seriation(hclust, N, N + N-2)
        seriated_dist = np.zeros((N,N))
        a,b = np.triu_indices(N,k=1)
        seriated_dist[a,b] = distance_matrix[ [res_order[i] for i in a], [res_order[j] for j in b]]
        seriated_dist[b,a] = seriated_dist[a,b]

        data['distance_matrix'] = seriated_dist.tolist()
        data['dm_labels'] = [pdbs[i] for i in res_order]

    cache.set(cache_key, data, 60*60*24*7)
    return JsonResponse(data)

# For reordering matrix based on h-tree
# Borrowed from https://gmarti.gitlab.io/ml/2017/09/07/how-to-sort-distance-matrix.html
def seriation(Z,N,cur_index):
    '''
        input:
            - Z is a hierarchical tree (dendrogram)
            - N is the number of points given to the clustering process
            - cur_index is the position in the tree for the recursive traversal
        output:
            - order implied by the hierarchical tree Z

        seriation computes the order implied by a hierarchical tree (dendrogram)
    '''
    if cur_index < N:
        return [cur_index]
    else:
        left = int(Z[cur_index-N,0])
        right = int(Z[cur_index-N,1])
        return (seriation(Z,N,left) + seriation(Z,N,right))

def getNewick(node, newick, parentdist, leaf_names, silhouette_coefficient):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        si_node = silhouette_coefficient[node.id]
        if len(newick) > 0:
            newick = ")%.2f:%.2f%s" % (si_node, parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names, silhouette_coefficient)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names, silhouette_coefficient)
        newick = "(%s" % (newick)
        return newick

def getSilhouetteIndex(node, distance_matrix, results):
    # set rootnode (DEBUG purposes)
    if node.id not in results:
        results[node.id] = 0

    if not node.is_leaf():
        # get list of indices cluster left (A)
        a = node.get_left().pre_order(lambda x: x.id)

        # get list of indices cluster right (B)
        b = node.get_right().pre_order(lambda x: x.id)

        if len(a) > 1:
            # calculate average Si - cluster A
            si_a = calculateSilhouetteIndex(distance_matrix, a, b)
            results[node.get_left().id] = si_a

            getSilhouetteIndex(node.get_left(), distance_matrix, results)

        if len(b) > 1:
            # calculate average Si - cluster B
            si_b = calculateSilhouetteIndex(distance_matrix, b, a)
            results[node.get_right().id] = si_b

            getSilhouetteIndex(node.get_right(), distance_matrix, results)

# Implementation based on Rousseeuw, P.J. J. Comput. Appl. Math. 20 (1987): 53-65
def calculateSilhouetteIndex(distance_matrix, a, b):
    si = 0
    for i in a:
        # calculate ai - avg distance within cluster
        ai = 0
        for j in a:
            if i != j:
                ai += distance_matrix[i,j]/(len(a)-1)

        # calculate bi - avg distance to closest cluster
        bi = 0
        for j in b:
            bi += distance_matrix[i,j]/len(b)

        # silhouette index (averaged)
        si += (bi-ai)/max(ai,bi)/len(a)

    return si


def DistanceData(request):
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

    pdbs = [pdb.upper() for pdb in pdbs]
    pdbs_lower = [pdb.lower() for pdb in pdbs]

    cache_key = ",".join(sorted(pdbs))
    cache_key = hashlib.md5(cache_key.encode('utf-8')).hexdigest()

    data = cache.get(cache_key)
    # data = None
    if data!=None:
        print('Result cached')
        return JsonResponse(data)

    # Segment filters
    try:
        segments = request.GET.getlist('segments[]')
    except IndexError:
        segments = []

    # Use generic numbers? Defaults to True.
    generic = True

    # Initialize response dictionary
    data = {}
    data['interactions'] = OrderedDict()
    data['pdbs'] = set()
    data['generic'] = generic
    data['segments'] = set()
    data['segment_map'] = {}
    # For Max schematics TODO -- make smarter.
    data['segment_map'] = {}
    data['aa_map'] = {}


    data['gn_map'] = OrderedDict()
    data['pos_map'] = OrderedDict()
    data['segment_map_full'] = OrderedDict()
    data['segment_map_full_gn'] = OrderedDict()
    data['generic_map_full'] = OrderedDict()

    dis = Distances()
    dis.load_pdbs(pdbs)
    start = time.time()
    dis.fetch_and_calculate(with_arr = True)
    print('done fetching',time.time()-start)
    # dis.calculate_window()

    excluded_segment = ['C-term','N-term']
    segments = ProteinSegment.objects.all().exclude(slug__in = excluded_segment)
    proteins =  Protein.objects.filter(protein__entry_name__in=pdbs_lower).distinct().all()

    start = time.time()
    list_of_gns = []
    if len(proteins)>1:
        a = Alignment()
        a.ignore_alternative_residue_numbering_schemes = True;
        a.load_proteins(proteins)
        a.load_segments(segments) #get all segments to make correct diagrams
        # build the alignment data matrix
        a.build_alignment()
        # calculate consensus sequence + amino acid and feature frequency
        a.calculate_statistics()
        consensus = a.full_consensus

        for aa in consensus:
            if 'x' in aa.family_generic_number:
                list_of_gns.append(aa.family_generic_number)
                data['gn_map'][aa.family_generic_number] = aa.amino_acid
                data['pos_map'][aa.sequence_number] = aa.amino_acid
                data['segment_map_full_gn'][aa.family_generic_number] = aa.segment_slug
    else:
        rs = Residue.objects.filter(protein_conformation__protein=proteins[0]).prefetch_related('protein_segment','display_generic_number','generic_number')
        for r in rs:
            if (not generic):
                data['pos_map'][r.sequence_number] = r.amino_acid
                data['segment_map_full'][r.sequence_number] = r.protein_segment.slug
                if r.display_generic_number:
                    data['generic_map_full'][r.sequence_number] = r.short_display_generic_number()
            else:
                if r.generic_number:
                    list_of_gns.append(r.generic_number.label)
                    data['gn_map'][r.generic_number.label] = r.amino_acid
                    data['pos_map'][r.sequence_number] = r.amino_acid
                    data['segment_map_full_gn'][r.generic_number.label] = r.protein_segment.slug

    print('done alignment',time.time()-start)
    # Dict to keep track of which residue numbers are in use
    number_dict = set()
    max_dispersion = 0
    start = time.time()
    for d in dis.stats:
        # print(d)
        res1 = d[0].split("_")[0]
        res2 = d[0].split("_")[1]
        res1_seg = res1.split("x")[0]
        res2_seg = res2.split("x")[0]
        data['segment_map'][res1] = res1_seg
        data['segment_map'][res2] = res2_seg
        data['segments'] |= {res1_seg} | {res2_seg}

        # Populate the AA map
        if 1==2:
            #When showing pdbs
            if pdb_name not in data['aa_map']:
                data['aa_map'][pdb_name] = {}

        number_dict |= {res1, res2}

        if res1 < res2:
            coord = str(res1) + ',' + str(res2)
        else:
            coord = str(res2) + ',' + str(res1)


        # if pdb_name not in data['interactions'][coord]:
        #     data['interactions'][coord][pdb_name] = []
        if len(pdbs) > 1:
            if d[3]:
                # correct data decimal
                d[5] = [x / 100 for x in d[5]]
                # Make easier readable output
                individual_pdbs = dict(zip(d[6], d[5]))

                if len(data['interactions'])<2000:
                    data['interactions'][coord] = [round(d[1])/100,round(d[3],3),individual_pdbs]
                else:
                    break

                if d[3]>max_dispersion:
                    max_dispersion = round(d[3],3)
        else:
            if d[1]:
                if len(data['interactions'])<2000:
                    data['interactions'][coord] = [round(d[1])/100,round(d[1],3)/100,d[-1]]
                else:
                    break

                if d[1]>max_dispersion*100:
                    max_dispersion = round(d[1],3)/100
        # data['sequence_numbers'] = sorted(number_dict)
    if (generic):
        data['sequence_numbers'] = sorted(number_dict, key=functools.cmp_to_key(gpcrdb_number_comparator))
    # else:
    #     data['sequence_numbers'] = sorted(number_dict)
        # break

    data['segments'] = list(data['segments'])
    data['pdbs'] = list(pdbs)
    data['max_dispersion'] = max_dispersion

    print('done data prep',time.time()-start)
    total = {}
    ngl_max_diff = 0
    for i,gn1 in enumerate(list_of_gns):
        if gn1 not in total:
            total[gn1] = {}
        for gn2 in list_of_gns[i:]:
            if gn2 not in total:
                total[gn2] = {}
            label = "{}_{}".format(gn1,gn2)
            if label in dis.stats_key:
                value = dis.stats_key[label][3]
                total[gn1][gn2] =  value
                total[gn2][gn1] =  value
        vals = []
        for gn,val in total[gn1].items():
            if gn[0]!=gn1[0]:
                #if not same segment
                vals.append(val)
        total[gn1]['avg'] = round(float(sum(vals))/max(len(vals),1),3)**2
        if abs(total[gn1]['avg'])>ngl_max_diff:
            ngl_max_diff = round(abs(total[gn1]['avg']),3)

    data['ngl_data'] = total
    data['ngl_max_diff'] = ngl_max_diff

    # Cache result 7 days
    cache.set(cache_key,data,3600*24*7)
    return JsonResponse(data)

# DEPRECATED FUNCTION?
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
        i_types_filter |= Q(interaction_type__in=i_types)


    hash_list = [pdbs,i_types,generic]
    hash_cache_key = 'interactiondata_{}'.format(get_hash(hash_list))
    data = cache.get(hash_cache_key)
    if data==None:

        # Get the relevant interactions
        interactions = Interaction.objects.filter(
            interacting_pair__referenced_structure__protein_conformation__protein__entry_name__in=pdbs
        ).values(
            'interacting_pair__referenced_structure__protein_conformation__protein__entry_name',
            'interacting_pair__res1__amino_acid',
            'interacting_pair__res2__amino_acid',
            'interacting_pair__res1__sequence_number',
            'interacting_pair__res1__generic_number__label',
            'interacting_pair__res1__protein_segment__slug',
            'interacting_pair__res2__sequence_number',
            'interacting_pair__res2__generic_number__label',
            'interacting_pair__res2__protein_segment__slug',
            'interaction_type',
        ).filter(
            segment_filter_res1 & segment_filter_res2 & i_types_filter
        )

        # Interaction type sort - optimize by statically defining interaction type order
        order = ['ionic', 'polar', 'aromatic', 'hydrophobic', 'van-der-waals']
        interactions = sorted(interactions, key=lambda x: order.index(x['interaction_type']))

        # Initialize response dictionary
        data = {}
        data['interactions'] = {}
        data['pdbs'] = set()
        data['generic'] = generic
        data['segments'] = set()
        data['segment_map'] = {}
        # For Max schematics TODO -- make smarter.
        data['segment_map'] = {}
        data['aa_map'] = {}

        # Create a consensus sequence.

        excluded_segment = ['C-term','N-term']
        segments = ProteinSegment.objects.all().filter(proteinfamily='GPCR').exclude(slug__in = excluded_segment)
        proteins =  Protein.objects.filter(protein__entry_name__in=pdbs).all()

        data['gn_map'] = OrderedDict()
        data['pos_map'] = OrderedDict()
        data['segment_map_full'] = OrderedDict()
        data['segment_map_full_gn'] = OrderedDict()
        data['generic_map_full'] = OrderedDict()

        if len(proteins)>1:
            a = Alignment()
            a.ignore_alternative_residue_numbering_schemes = True;
            a.load_proteins(proteins)
            a.load_segments(segments) #get all segments to make correct diagrams
            # build the alignment data matrix
            a.build_alignment()
            # calculate consensus sequence + amino acid and feature frequency
            a.calculate_statistics()
            consensus = a.full_consensus

            for aa in consensus:
                if 'x' in aa.family_generic_number:
                    data['gn_map'][aa.family_generic_number] = aa.amino_acid
                    data['pos_map'][aa.sequence_number] = aa.amino_acid
                    data['segment_map_full_gn'][aa.family_generic_number] = aa.segment_slug
        else:
            rs = Residue.objects.filter(protein_conformation__protein=proteins[0]).prefetch_related('protein_segment','display_generic_number','generic_number')
            for r in rs:
                if (not generic):
                    data['pos_map'][r.sequence_number] = r.amino_acid
                    data['segment_map_full'][r.sequence_number] = r.protein_segment.slug
                    if r.display_generic_number:
                        data['generic_map_full'][r.sequence_number] = r.short_display_generic_number()

        for i in interactions:
            pdb_name = i['interacting_pair__referenced_structure__protein_conformation__protein__entry_name']
            if not pdb_name in data['pdbs']:
                data['pdbs'].add(pdb_name)

        # Map from ordinary residue numbers to generic where available
        if (not generic):
            data['generic_map'] = {}

        # Dict to keep track of which residue numbers are in use
        number_dict = set()

        for i in interactions:
            pdb_name = i['interacting_pair__referenced_structure__protein_conformation__protein__entry_name']
            res1_seq = i['interacting_pair__res1__sequence_number']
            res2_seq = i['interacting_pair__res2__sequence_number']
            res1_gen = i['interacting_pair__res1__generic_number__label']
            res2_gen = i['interacting_pair__res2__generic_number__label']
            res1_seg = i['interacting_pair__res1__protein_segment__slug']
            res2_seg = i['interacting_pair__res2__protein_segment__slug']
            res1_aa = i['interacting_pair__res1__amino_acid']
            res2_aa = i['interacting_pair__res2__amino_acid']
            model = i['interaction_type']

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

            if not generic and res1_gen:
                data['generic_map'][res1] = res1_gen

            if not generic and res2_gen:
                data['generic_map'][res2] = res2_gen

            # List which segments are available.
            data['segment_map'][res1] = res1_seg
            data['segment_map'][res2] = res2_seg
            data['segments'] |= {res1_seg} | {res2_seg}

            # Populate the AA map
            if pdb_name not in data['aa_map']:
                data['aa_map'][pdb_name] = {}

            data['aa_map'][pdb_name][res1] = res1_aa
            data['aa_map'][pdb_name][res2] = res2_aa

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
        cache.set(hash_cache_key, data, 3600*24)
    return JsonResponse(data)

def ServePDB(request, pdbname):
    structure=Structure.objects.filter(pdb_code__index=pdbname.upper())
    if structure.exists():
        structure=structure.get()
    else:
        quit() #quit!

    if structure.pdb_data is None:
        quit()

    only_gns = list(structure.protein_conformation.residue_set.exclude(generic_number=None).values_list('protein_segment__slug','sequence_number','generic_number__label','display_generic_number__label').all())
    only_gn = []
    gn_map = []
    gn_map_classa = []
    segments = {}
    for gn in only_gns:
        only_gn.append(gn[1])
        # Use and format the display generic number to get the shorthand class specific number.
        gn_map.append(re.sub(r'\.[\d]+', '', gn[3]))
        gn_map_classa.append(gn[2])
        if gn[0] not in segments:
            segments[gn[0]] = []
        segments[gn[0]].append(gn[1])

    data = {}
    data['pdb'] = structure.pdb_data.pdb
    data['only_gn'] = only_gn
    data['gn_map'] = gn_map
    data['gn_map_classa'] = gn_map_classa
    data['segments'] = segments
    data['chain'] = structure.preferred_chain
    # positioning data
    sv = StructureVectors.objects.filter(structure=structure)
    if sv.exists():
        sv = sv.get()
        data['translation'] = sv.translation
        data['center_axis'] = sv.center_axis

    return JsonResponse(data)

def StateContacts(request):
    contacts = ConsensusInteraction.objects.filter(state_specific = True).prefetch_related('gn1','gn2','state','protein_class').all()
    context = {}
    context['contacts'] = contacts
    return render(request, 'contactnetwork/state_contacts.html', context)
