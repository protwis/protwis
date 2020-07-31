from django.shortcuts import get_object_or_404, render
from django.http import HttpResponse
from django.core.cache import cache
from django.core.cache import caches
try:
    cache_variation = caches['variation']
except:
    cache_variation = cache

from django.db.models import Count, Min, Sum, Avg, Q
from django.views.decorators.cache import cache_page

import hashlib

from protein.models import Protein, ProteinConformation, ProteinAlias, ProteinFamily, Gene, ProteinGProtein, ProteinGProteinPair
from residue.models import Residue, ResiduePositionSet, ResidueSet
from mutational_landscape.models import NaturalMutations, CancerMutations, DiseaseMutations, PTMs, NHSPrescribings

from common.diagrams_gpcr import DrawHelixBox, DrawSnakePlot

from drugs.models import Drugs

from mutation.functions import *
from mutation.models import *

from interaction.models import *
from interaction.views import ajax #import x-tal interactions

from common import definitions
from collections import OrderedDict
from common.views import AbsTargetSelection
from common.views import AbsSegmentSelection
from family.views import linear_gradient, color_dict, RGB_to_hex, hex_to_RGB

import re
import json
import numpy as np
from collections import OrderedDict
from copy import deepcopy

from io import BytesIO
import re
import math
import unicodedata
import urllib
import xlsxwriter #sudo pip3 install XlsxWriter
import operator
import string

class TargetSelection(AbsTargetSelection):
    step = 1
    number_of_steps = 1
    filters = False
    psets = False
    # docs = 'mutations.html#mutation-browser'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Show missense variants',
            'url': '/mutational_landscape/render',
            'color': 'success',
        },
    }
    default_species = False

#@cache_page(60*60*24*21)
def render_variants(request, protein=None, family=None, download=None, receptor_class=None, gn=None, aa=None, **response_kwargs):

    simple_selection = request.session.get('selection', False)
    proteins = []
    target_type = 'protein'
    if protein:  # if protein static page
        proteins.append(Protein.objects.get(entry_name=protein.lower()))

    # flatten the selection into individual proteins
    elif simple_selection:
        for target in simple_selection.targets:
            if target.type == 'protein':
                proteins.append(target.item)
            elif target.type == 'family':
                target_type = 'family'
                familyname = target.item
                # species filter
                species_list = []
                for species in simple_selection.species:
                    species_list.append(species.item)

                # annotation filter
                protein_source_list = []
                for protein_source in simple_selection.annotation:
                    protein_source_list.append(protein_source.item)

                if species_list:
                    family_proteins = Protein.objects.filter(family__slug__startswith=target.item.slug,
                        species__in=(species_list),
                        source__in=(protein_source_list)).select_related('residue_numbering_scheme', 'species')
                else:
                    family_proteins = Protein.objects.filter(family__slug__startswith=target.item.slug,
                        source__in=(protein_source_list)).select_related('residue_numbering_scheme', 'species')

                for fp in family_proteins:
                    proteins.append(fp)

    # Caching results for unique protein sets
    cache_key = "VARIATION_"+hashlib.md5(str(proteins).encode('utf-8')).hexdigest()
    if not cache_variation.has_key(cache_key):
        NMs = NaturalMutations.objects.filter(Q(protein__in=proteins)).prefetch_related('residue__generic_number','residue__display_generic_number','residue__protein_segment','protein')
        ptms = PTMs.objects.filter(Q(protein__in=proteins)).prefetch_related('residue')
        ptms_dict = {}

        ## MICROSWITCHES
        micro_switches_rset = ResiduePositionSet.objects.get(name="State (micro-)switches")
        ms_label = []
        for residue in micro_switches_rset.residue_position.all():
            ms_label.append(residue.label)

        ms_object = Residue.objects.filter(protein_conformation__protein=proteins[0], generic_number__label__in=ms_label)
        ms_sequence_numbers = []
        for ms in ms_object:
            ms_sequence_numbers.append(ms.sequence_number)

        ## SODIUM POCKET
        sodium_pocket_rset = ResiduePositionSet.objects.get(name="Sodium ion pocket")
        sp_label = []
        for residue in sodium_pocket_rset.residue_position.all():
            sp_label.append(residue.label)

        sp_object = Residue.objects.filter(protein_conformation__protein=proteins[0], generic_number__label__in=ms_label)
        sp_sequence_numbers = []
        for sp in sp_object:
            sp_sequence_numbers.append(sp.sequence_number)

        for ptm in ptms:
            ptms_dict[ptm.residue.sequence_number] = ptm.modification

        ## G PROTEIN INTERACTION POSITIONS
        # THIS SHOULD BE CLASS SPECIFIC (different set)
        rset = ResiduePositionSet.objects.get(name='G-protein interface')
        gprotein_generic_set = []
        for residue in rset.residue_position.all():
            gprotein_generic_set.append(residue.label)

        ### GET LB INTERACTION DATA
        # get also ortholog proteins, which might have been crystallised to extract
        # interaction data also from those
        if protein:
            orthologs = Protein.objects.filter(family__slug=proteins[0].family.slug, sequence_type__slug='wt')
        else:
            orthologs = Protein.objects.filter(family__slug__startswith=proteins[0].family.slug, sequence_type__slug='wt')

        interactions = ResidueFragmentInteraction.objects.filter(
            structure_ligand_pair__structure__protein_conformation__protein__parent__in=orthologs, structure_ligand_pair__annotated=True).exclude(interaction_type__type ='hidden').all()
        interaction_data = {}
        for interaction in interactions:
            if interaction.rotamer.residue.generic_number:
                sequence_number = interaction.rotamer.residue.sequence_number
                # sequence_number = lookup[interaction.rotamer.residue.generic_number.label]
                label = interaction.rotamer.residue.generic_number.label
                aa = interaction.rotamer.residue.amino_acid
                interactiontype = interaction.interaction_type.name
                if sequence_number not in interaction_data:
                    interaction_data[sequence_number] = []
                if interactiontype not in interaction_data[sequence_number]:
                    interaction_data[sequence_number].append(interactiontype)

        # Fixes fatal error - in case of receptor family selection (e.g. H1 receptors)
        if target_type == 'family' and len(proteins[0].family.slug) < 15:
            pc = ProteinConformation.objects.get(protein__family__name=familyname, protein__sequence_type__slug='consensus')
            residuelist = Residue.objects.filter(protein_conformation=pc).order_by('sequence_number').prefetch_related('protein_segment', 'generic_number', 'display_generic_number')
        else:
            residuelist = Residue.objects.filter(protein_conformation__protein=proteins[0]).prefetch_related('protein_segment', 'display_generic_number', 'generic_number')

        jsondata = {}
        for NM in NMs:
            functional_annotation = ''
            SN = NM.residue.sequence_number
            if NM.residue.generic_number:
                GN = NM.residue.generic_number.label
            else:
                GN = ''
            if SN in sp_sequence_numbers:
                functional_annotation += 'SodiumPocket '
            if SN in ms_sequence_numbers:
                functional_annotation += 'MicroSwitch '
            if SN in ptms_dict:
                functional_annotation += 'PTM (' + ptms_dict[SN] + ') '
            if SN in interaction_data:
                functional_annotation += 'LB (' + ', '.join(interaction_data[SN]) + ') '
            if GN in gprotein_generic_set:
                functional_annotation += 'GP (contact) '

            ms_type = NM.type
            if ms_type == 'missense':
                effect = 'deleterious' if NM.sift_score <= 0.05 or NM.polyphen_score >= 0.1 else 'tolerated'
                color = '#e30e0e' if NM.sift_score <= 0.05 or NM.polyphen_score >= 0.1 else '#70c070'
            else:
                effect = 'deleterious'
                color = '#575c9d'
            # account for multiple mutations at this position!
            NM.functional_annotation = functional_annotation
            # print(NM.functional_annotation)
            jsondata[SN] = [NM.amino_acid, NM.allele_frequency, NM.allele_count, NM.allele_number, NM.number_homozygotes, NM.type, effect, color, functional_annotation]


        natural_mutation_list = {}
        max_snp_pos = 1
        for NM in NMs:
            if NM.residue.generic_number:
                if NM.residue.generic_number.label in natural_mutation_list:
                    natural_mutation_list[NM.residue.generic_number.label]['val'] += 1
                    if not str(NM.amino_acid) in natural_mutation_list[NM.residue.generic_number.label]['AA']:
                        natural_mutation_list[NM.residue.generic_number.label]['AA'] = natural_mutation_list[NM.residue.generic_number.label]['AA'] + str(NM.amino_acid) + ' '

                    if natural_mutation_list[NM.residue.generic_number.label]['val'] > max_snp_pos:
                        max_snp_pos = natural_mutation_list[NM.residue.generic_number.label]['val']
                else:
                    natural_mutation_list[NM.residue.generic_number.label] = {'val':1, 'AA': NM.amino_acid + ' '}

        jsondata_natural_mutations = {}

        for r in residuelist:
            if r.generic_number:
                if r.generic_number.label in natural_mutation_list:
                    jsondata_natural_mutations[r.sequence_number] = natural_mutation_list[r.generic_number.label]

        jsondata_natural_mutations['color'] = linear_gradient(start_hex="#c79494", finish_hex="#c40100", n=max_snp_pos)
        # jsondata_cancer_mutations['color'] = linear_gradient(start_hex="#d8baff", finish_hex="#422d65", n=max_cancer_pos)
        # jsondata_disease_mutations['color'] = linear_gradient(start_hex="#ffa1b1", finish_hex="#6e000b", n=max_disease_pos)
        #
        SnakePlot = DrawSnakePlot(residuelist, "Class A", protein, nobuttons=1)
        HelixBox = DrawHelixBox(residuelist, 'Class A', protein, nobuttons=1)

        cache_data = {'mutations': NMs, 'type': target_type, 'HelixBox': HelixBox, 'SnakePlot': SnakePlot, 'receptor': str(proteins[0].entry_name), 'mutations_pos_list': json.dumps(jsondata), 'natural_mutations_pos_list': json.dumps(jsondata_natural_mutations)}
        cache_variation.set(cache_key, cache_data, 60*60*24*21)
    else:
        cache_data = cache_variation.get(cache_key)

    # EXCEL TABLE EXPORT
    if download:
        data = []
        for r in cache_data['mutations']:
            values = r.__dict__
            complete = {'protein': r.protein.name, 'sequence_number': r.residue.sequence_number, 'orig_amino_acid': r.residue.amino_acid}
            complete.update(values)
            data.append(complete)
        headers = ['protein', 'sequence_number', 'orig_amino_acid', 'type', 'amino_acid', 'allele_count', 'allele_number', 'allele_frequency', 'polyphen_score', 'sift_score', 'number_homozygotes', 'functional_annotation']

        # EXCEL SOLUTION
        output = BytesIO()
        workbook = xlsxwriter.Workbook(output)
        worksheet = workbook.add_worksheet()

        col = 0
        for h in headers:
            worksheet.write(0, col, h)
            col += 1
        row = 1
        for d in data:
            col = 0
            for h in headers:
                worksheet.write(row, col, str(d[h]))
                col += 1
            row += 1
        workbook.close()
        output.seek(0)
        xlsx_data = output.read()

        response = HttpResponse(xlsx_data, content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
        if target_type == 'family':
            response['Content-Disposition'] = 'attachment; filename=' + clean_filename('GPCRdb_' + str(familyname) + '_variant_data.xlsx')  # % 'mutations'
        else:
            response['Content-Disposition'] = 'attachment; filename=GPCRdb_' + proteins[0].entry_name + '_variant_data.xlsx'  # % 'mutations'
        return response

    return render(request, 'browser.html', cache_data)

def ajaxNaturalMutation(request, slug, **response_kwargs):

    name_of_cache = 'ajaxNaturalMutation_'+slug

    ptms = PTMs.objects.filter(protein__entry_name=slug).prefetch_related('residue')
    ptms_dict = {}

    for ptm in ptms:
        ptms_dict[ptm.residue.sequence_number] = ptm.modification

    ## MICROSWITCHES
    micro_switches_rset = ResiduePositionSet.objects.get(name="State (micro-)switches")
    ms_label = []
    for residue in micro_switches_rset.residue_position.all():
        ms_label.append(residue.label)

    ms_object = Residue.objects.filter(protein_conformation__protein__entry_name=slug, generic_number__label__in=ms_label)
    ms_sequence_numbers = []
    for ms in ms_object:
        ms_sequence_numbers.append(ms.sequence_number)

    ## SODIUM POCKET
    sodium_pocket_rset = ResiduePositionSet.objects.get(name="Sodium ion pocket")
    sp_label = []
    for residue in sodium_pocket_rset.residue_position.all():
        sp_label.append(residue.label)

    sp_object = Residue.objects.filter(protein_conformation__protein__entry_name=slug, generic_number__label__in=ms_label)
    sp_sequence_numbers = []
    for sp in sp_object:
        sp_sequence_numbers.append(sp.sequence_number)

    ## G PROTEIN INTERACTION POSITIONS
    # THIS SHOULD BE CLASS SPECIFIC (different set)
    rset = ResiduePositionSet.objects.get(name='G-protein interface')
    gprotein_generic_set = []
    for residue in rset.residue_position.all():
        gprotein_generic_set.append(residue.label)

    ### GET LB INTERACTION DATA
    # get also ortholog proteins, which might have been crystallised to extract
    # interaction data also from those
    p = Protein.objects.get(entry_name=slug)
    orthologs = Protein.objects.filter(family__slug__startswith=p.family.slug, sequence_type__slug='wt')

    interactions = ResidueFragmentInteraction.objects.filter(
        structure_ligand_pair__structure__protein_conformation__protein__parent__in=orthologs, structure_ligand_pair__annotated=True).exclude(interaction_type__type ='hidden').order_by('rotamer__residue__sequence_number')
    interaction_data = {}
    for interaction in interactions:
        if interaction.rotamer.residue.generic_number:
            sequence_number = interaction.rotamer.residue.sequence_number
            # sequence_number = lookup[interaction.rotamer.residue.generic_number.label]
            label = interaction.rotamer.residue.generic_number.label
            aa = interaction.rotamer.residue.amino_acid
            interactiontype = interaction.interaction_type.name
            if sequence_number not in interaction_data:
                interaction_data[sequence_number] = []
            if interactiontype not in interaction_data[sequence_number]:
                interaction_data[sequence_number].append(interactiontype)

    jsondata = cache.get(name_of_cache)

    if jsondata == None:
        jsondata = {}

        NMs = NaturalMutations.objects.filter(protein__entry_name=slug).prefetch_related('residue')

        for NM in NMs:

            SN = NM.residue.sequence_number
            type = NM.type

            if type == 'missense':
                if NM.sift_score != None and NM.polyphen_score != None:
                    effect = 'deleterious' if NM.sift_score <= 0.05 or NM.polyphen_score >= 0.1 else 'tolerated'
                    color = '#e30e0e' if NM.sift_score <= 0.05 or NM.polyphen_score >= 0.1 else '#70c070'
                else:
                    effect = 'unknown'
                    color = '#818181'
            else:
                effect = 'deleterious'
                color = '#65368e'

            functional_annotation = ''
            SN = NM.residue.sequence_number
            if NM.residue.generic_number:
                GN = NM.residue.generic_number.label
            else:
                GN = ''
            if SN in sp_sequence_numbers:
                functional_annotation += 'SodiumPocket '
            if SN in ms_sequence_numbers:
                functional_annotation += 'MicroSwitch '
            if SN in ptms_dict:
                functional_annotation += 'PTM (' + ptms_dict[SN] + ') '
            if SN in interaction_data:
                functional_annotation += 'LB (' + ', '.join(interaction_data[SN]) + ') '
            if GN in gprotein_generic_set:
                functional_annotation += 'GP (contact) '

            if functional_annotation == '':
                functional_annotation = '-'
            # account for multiple mutations at this position!
            jsondata[SN] = [NM.amino_acid, NM.allele_frequency, NM.allele_count, NM.allele_number, NM.number_homozygotes, NM.type, effect, color, functional_annotation]

        jsondata = json.dumps(jsondata)
        response_kwargs['content_type'] = 'application/json'

        cache.set(name_of_cache, jsondata, 20) # 60*60*24*2 two days timeout on cache

    return HttpResponse(jsondata, **response_kwargs)

def ajaxPTMs(request, slug, **response_kwargs):

    name_of_cache = 'ajaxPTMs_'+slug

    jsondata = cache.get(name_of_cache)

    if jsondata == None:
        jsondata = {}

        NMs = PTMs.objects.filter(protein__entry_name=slug).prefetch_related('residue')

        for NM in NMs:

            SN = NM.residue.sequence_number
            mod = NM.modification

            jsondata[SN] = [mod]

        jsondata = json.dumps(jsondata)
        response_kwargs['content_type'] = 'application/json'

        cache.set(name_of_cache, jsondata, 20) # 60*60*24*2 two days timeout on cache

    return HttpResponse(jsondata, **response_kwargs)

# def ajaxCancerMutation(request, slug, **response_kwargs):
#
#     name_of_cache = 'ajaxCancerMutation_'+slug
#
#     jsondata = cache.get(name_of_cache)
#
#     if jsondata == None:
#         jsondata = {}
#
#         CMs = CancerMutations.objects.filter(protein__entry_name=slug).prefetch_related('residue')
#
#         for CM in CMs:
#             SN = CM.residue.sequence_number
#             jsondata[SN] = [CM.amino_acid]
#
#         jsondata = json.dumps(jsondata)
#         response_kwargs['content_type'] = 'application/json'
#
#         cache.set(name_of_cache, jsondata, 20) #two days timeout on cache
#
#     return HttpResponse(jsondata, **response_kwargs)
#
# def ajaxDiseaseMutation(request, slug, **response_kwargs):
#
#     name_of_cache = 'ajaxDiseaseMutation_'+slug
#
#     jsondata = cache.get(name_of_cache)
#
#     if jsondata == None:
#         jsondata = {}
#
#         DMs = DiseaseMutations.objects.filter(protein__entry_name=slug).prefetch_related('residue')
#
#         for DM in DMs:
#             SN = DM.residue.sequence_number
#             jsondata[SN] = [DM.amino_acid]
#
#         jsondata = json.dumps(jsondata)
#         response_kwargs['content_type'] = 'application/json'
#
#         cache.set(name_of_cache, jsondata, 20) #two days timeout on cache
#
#     return HttpResponse(jsondata, **response_kwargs)

def mutant_extract(request):
    import pandas as pd
    mutations = MutationExperiment.objects.all().prefetch_related('residue__display_generic_number','protein__family','exp_func','exp_type','ligand','ligand_role','refs','mutation')
    # mutations = MutationExperiment.objects.filter(protein__entry_name__startswith=slug_without_species).order_by('residue__sequence_number').prefetch_related('residue')

    temp = pd.DataFrame(columns=['EntryName','Family','LigandType','Class','SequenceNumber','GPCRdb','Segment','WTaa','Mutantaa','foldchange','Ligand','LigandRole','ExpQual','ExpWTValue','ExpWTVUnit','ExpMutantValue','ExpMutantSign','ExpType','ExpFunction'])
    row = 0
    for mutation in mutations:
        if mutation.ligand:
            ligand = mutation.ligand.name
        else:
            ligand = 'NaN'

        if mutation.exp_qual:
            qual = mutation.exp_qual.qual
        else:
            qual = 'NaN'

        if mutation.exp_func_id:
            func = mutation.exp_func.func
        else:
            func = 'NaN'

        if mutation.ligand_role_id:
            lrole = mutation.ligand_role.name
        else:
            lrole = 'NaN'

        if mutation.exp_type_id:
            etype = mutation.exp_type.type
        else:
            etype = 'NaN'

        if mutation.residue.display_generic_number:
            gpcrdb = mutation.residue.display_generic_number.label
        else:
            gpcrdb = 'NaN'

        if mutation.foldchange != 0:
            # print(mutation.protein.entry_name, mutation.residue.sequence_number, mutation.residue.amino_acid, mutation.mutation.amino_acid, mutation.foldchange,ligand, lrole,qual,mutation.wt_value, mutation.wt_unit, mutation.mu_value, mutation.mu_sign, etype, func)
            temp.loc[row] = pd.Series({'EntryName': mutation.protein.entry_name, 'Family': mutation.protein.family.parent.name,'LigandType': mutation.protein.family.parent.parent.name,'Class': mutation.protein.family.parent.parent.parent.name, 'SequenceNumber': int(mutation.residue.sequence_number), 'GPCRdb': gpcrdb, 'Segment': mutation.residue.protein_segment.slug,'WTaa': mutation.residue.amino_acid, 'Mutantaa': mutation.mutation.amino_acid, 'foldchange': mutation.foldchange, 'Ligand': ligand, 'LigandRole': lrole, 'ExpQual':  qual, 'ExpWTValue': mutation.wt_value, 'ExpWTVUnit': mutation.wt_unit, 'ExpMutantValue': mutation.mu_value, 'ExpMutantSign': mutation.mu_sign, 'ExpType': etype, 'ExpFunction': func})
            row += 1
        if row % 200 == 0 and row != 0:
            print(row)

    temp.to_csv('170125_GPCRdb_mutation.csv')
        # jsondata[mutation.residue.sequence_number].append([mutation.foldchange,ligand,qual])
    # print(jsondata)

@cache_page(60*60*24*21)
def statistics(request):

    context = dict()

    families = ProteinFamily.objects.all()
    lookup = {}
    for f in families:
        lookup[f.slug] = f.name.replace("receptors","").replace(" receptor","").replace(" hormone","").replace("/neuropeptide","/").replace(" (G protein-coupled)","").replace(" factor","").replace(" (LPA)","").replace(" (S1P)","").replace("GPR18, GPR55 and GPR119","GPR18/55/119").replace("-releasing","").replace(" peptide","").replace(" and oxytocin","/Oxytocin").replace("Adhesion class orphans","Adhesion orphans").replace("muscarinic","musc.").replace("-concentrating","-conc.")

    class_proteins = Protein.objects.filter(family__slug__startswith="00",source__name='SWISSPROT', species_id=1).prefetch_related('family').order_by('family__slug')

    temp = OrderedDict([
                    ('name',''),
                    ('number_of_variants', 0),
                    ('number_of_children', 0),
                    ('receptor_t',0),
                    ('density_of_variants', 0),
                    ('children', OrderedDict())
                    ])

    coverage = OrderedDict()

    # Make the scaffold
    for p in class_proteins:
        #print(p,p.family.slug)
        fid = p.family.slug.split("_")
        if fid[0] not in coverage:
            coverage[fid[0]] = deepcopy(temp)
            coverage[fid[0]]['name'] = lookup[fid[0]]
        if fid[1] not in coverage[fid[0]]['children']:
            coverage[fid[0]]['children'][fid[1]] = deepcopy(temp)
            coverage[fid[0]]['children'][fid[1]]['name'] = lookup[fid[0]+"_"+fid[1]]
        if fid[2] not in coverage[fid[0]]['children'][fid[1]]['children']:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]] = deepcopy(temp)
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['name'] = lookup[fid[0]+"_"+fid[1]+"_"+fid[2]][:28]
        if fid[3] not in coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children']:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]] = deepcopy(temp)
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['name'] = p.entry_name.split("_")[0] #[:10]
            coverage[fid[0]]['receptor_t'] += 1
            coverage[fid[0]]['children'][fid[1]]['receptor_t'] += 1
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['receptor_t'] += 1
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['receptor_t'] = 1
    # # POULATE WITH DATA
    variants_target = Protein.objects.filter(family__slug__startswith="00", entry_name__icontains='_human').values('family_id__slug').annotate(value=Count('naturalmutations__residue_id', distinct = True))
    protein_lengths = Protein.objects.filter(family__slug__startswith="00", entry_name__icontains='_human').values('family_id__slug','sequence')
    protein_lengths_dict = {}
    for i in protein_lengths:
        protein_lengths_dict[i['family_id__slug']] =  i['sequence']
    for i in variants_target:
        # print(i)
        fid = i['family_id__slug'].split("_")
        coverage[fid[0]]['number_of_variants'] += i['value']
        coverage[fid[0]]['children'][fid[1]]['number_of_variants'] += i['value']
        coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['number_of_variants'] += i['value']
        coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['number_of_variants'] += i['value']
        density = float(i['value'])/len(protein_lengths_dict[i['family_id__slug']])
        coverage[fid[0]]['density_of_variants'] += round(density,2)
        coverage[fid[0]]['children'][fid[1]]['density_of_variants'] += round(density,2)
        coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['density_of_variants'] += round(density,2)
        coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['density_of_variants'] += round(density,2)
        coverage[fid[0]]['number_of_children'] += 1
        coverage[fid[0]]['children'][fid[1]]['number_of_children'] += 1
        coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['number_of_children'] += 1
        coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['number_of_children'] += 1

    # MAKE THE TREE
    tree = OrderedDict({'name':'GPCRs','children':[]})
    i = 0
    n = 0
    for c,c_v in coverage.items():
        c_v['name'] = c_v['name'].split("(")[0]
        if c_v['name'].strip() in ['Other GPCRs']:
            # i += 1
            continue
            # pass
        children = []
        for lt,lt_v in c_v['children'].items():
            if lt_v['name'].strip() == 'Orphan' and c_v['name'].strip()=="Class A":
                # $pass
                continue
            children_rf = []
            for rf,rf_v in lt_v['children'].items():
                rf_v['name'] = rf_v['name'].split("<")[0]
                children_r = []
                for r,r_v in rf_v['children'].items():
                    r_v['sort'] = n
                    children_r.append(r_v)
                    n += 1
                rf_v['children'] = children_r
                rf_v['sort'] = n
                children_rf.append(rf_v)
            lt_v['children'] = children_rf
            lt_v['sort'] = n
            children.append(lt_v)
        c_v['children'] = children
        c_v['sort'] = n
        tree['children'].append(c_v)
        #tree = c_v
        #break
        i += 1

    context['tree'] = json.dumps(tree)

    ## Overview statistics
    total_receptors = NaturalMutations.objects.filter(type='missense').values('protein_id').distinct().count()
    total_mv = len(NaturalMutations.objects.filter(type='missense'))
    total_lof = len(NaturalMutations.objects.exclude(type='missense'))
    if total_receptors > 0:
        total_av_rv = round(len(NaturalMutations.objects.filter(type='missense', allele_frequency__lt=0.001))/ total_receptors,1)
        total_av_cv = round(len(NaturalMutations.objects.filter(type='missense', allele_frequency__gte=0.001))/ total_receptors,1)
    else:
        total_av_rv = 0
        total_av_cv = 0

    context['stats'] = {'total_mv':total_mv,'total_lof':total_lof,'total_av_rv':total_av_rv, 'total_av_cv':total_av_cv}

    return render(request, 'variation_statistics.html', context)

def get_functional_sites(protein):

    ## PTMs
    ptms = list(PTMs.objects.filter(protein=protein).values_list('residue', flat=True).distinct())

    ## MICROSWITCHES
    micro_switches_rset = ResiduePositionSet.objects.get(name="State (micro-)switches")
    ms_label = []
    for residue in micro_switches_rset.residue_position.all():
        ms_label.append(residue.label)

    ms_object = list(Residue.objects.filter(protein_conformation__protein=protein, generic_number__label__in=ms_label).values_list('id', flat=True).distinct())

    ## SODIUM POCKET
    sodium_pocket_rset = ResiduePositionSet.objects.get(name="Sodium ion pocket")
    sp_label = []
    for residue in sodium_pocket_rset.residue_position.all():
        sp_label.append(residue.label)
    sp_object = list(Residue.objects.filter(protein_conformation__protein=protein, generic_number__label__in=ms_label).values_list('id', flat=True).distinct())


    ## G PROTEIN INTERACTION POSITIONS
    # THIS SHOULD BE CLASS SPECIFIC (different set)
    rset = ResiduePositionSet.objects.get(name='G-protein interface')
    gprotein_generic_set = []
    for residue in rset.residue_position.all():
        gprotein_generic_set.append(residue.label)
    GP_object = list(Residue.objects.filter(protein_conformation__protein=protein, generic_number__label__in=gprotein_generic_set).values_list('id', flat=True).distinct())

    ### GET LB INTERACTION DATA
    ## get also ortholog proteins, which might have been crystallised to extract
    ## interaction data also from those
    orthologs = Protein.objects.filter(family__slug__startswith=protein.family.slug, sequence_type__slug='wt').prefetch_related('protein__family')
    interaction_residues = ResidueFragmentInteraction.objects.filter(
        structure_ligand_pair__structure__protein_conformation__protein__parent__in=orthologs, structure_ligand_pair__annotated=True).exclude(interaction_type__type ='hidden').values_list('rotamer__residue_id', flat=True).distinct()

    ## Get variants of these known residues:
    known_function_sites = set(x for l in [GP_object,sp_object,ms_object,ptms,interaction_residues] for x in l)
    NMs = NaturalMutations.objects.filter(residue_id__in=known_function_sites)
    return len(NMs)

# Based on https://gist.github.com/wassname/1393c4a57cfcbf03641dbc31886123b8
def clean_filename(filename, replace=' '):
    char_limit = 255

    # replace spaces
    filename.replace(' ','_')

    # keep only valid ascii chars
    cleaned_filename = unicodedata.normalize('NFKD', filename).encode('ASCII', 'ignore').decode()

    # keep only whitelisted chars
    whitelist = "-_.() %s%s" % (string.ascii_letters, string.digits)
    cleaned_filename = ''.join(c for c in cleaned_filename if c in whitelist)
    if len(cleaned_filename)>char_limit:
        print("Warning, filename truncated because it was over {}. Filenames may no longer be unique".format(char_limit))
    return cleaned_filename[:char_limit]

@cache_page(60*60*24*21)
def economicburden(request):
    economic_data = [{'values': [{'y': 0, 'x': 'known-homozygous'}, {'y': 0, 'x': 'known-all variants'}, {'y': 29574708, 'x': 'putative-homozygous'}, {'y': 186577951, 'x': 'putative-all variants'}], 'key': 'Analgesics'}, {'values': [{'y': 0, 'x': 'known-homozygous'}, {'y': 0, 'x': 'known-all variants'}, {'y': 0, 'x': 'putative-homozygous'}, {'y': 14101883, 'x': 'putative-all variants'}], 'key': 'Antidepressant Drugs'}, {'values': [{'y': 0, 'x': 'known-homozygous'}, {'y': 0, 'x': 'known-all variants'}, {'y': 0, 'x': 'putative-homozygous'}, {'y': 10637449, 'x': 'putative-all variants'}], 'key': 'Antihist, Hyposensit & Allergic Emergen'}, {'values': [{'y': 0, 'x': 'known-homozygous'}, {'y': 0, 'x': 'known-all variants'}, {'y': 0, 'x': 'putative-homozygous'}, {'y': 6633692, 'x': 'putative-all variants'}], 'key': 'Antispasmod.&Other Drgs Alt.Gut Motility'}, {'values': [{'y': 0, 'x': 'known-homozygous'}, {'y': 0, 'x': 'known-all variants'}, {'y': 8575714, 'x': 'putative-homozygous'}, {'y': 27008513, 'x': 'putative-all variants'}], 'key': 'Beta-Adrenoceptor Blocking Drugs'}, {'values': [{'y': 0, 'x': 'known-homozygous'}, {'y': 10108322, 'x': 'known-all variants'}, {'y': 25187489, 'x': 'putative-homozygous'}, {'y': 89224667, 'x': 'putative-all variants'}], 'key': 'Bronchodilators'}, {'values': [{'y': 0, 'x': 'known-homozygous'}, {'y': 5466184, 'x': 'known-all variants'}, {'y': 0, 'x': 'putative-homozygous'}, {'y': 10313279, 'x': 'putative-all variants'}], 'key': 'Drugs For Genito-Urinary Disorders'}, {'values': [{'y': 13015487, 'x': 'known-homozygous'}, {'y': 44334808, 'x': 'known-all variants'}, {'y': 13015487, 'x': 'putative-homozygous'}, {'y': 45130626, 'x': 'putative-all variants'}], 'key': 'Drugs Used In Diabetes'}, {'values': [{'y': 0, 'x': 'known-homozygous'}, {'y': 0, 'x': 'known-all variants'}, {'y': 0, 'x': 'putative-homozygous'}, {'y': 12168533, 'x': 'putative-all variants'}], 'key': "Drugs Used In Park'ism/Related Disorders"}, {'values': [{'y': 0, 'x': 'known-homozygous'}, {'y': 0, 'x': 'known-all variants'}, {'y': 0, 'x': 'putative-homozygous'}, {'y': 28670250, 'x': 'putative-all variants'}], 'key': 'Drugs Used In Psychoses & Rel.Disorders'}, {'values': [{'y': 0, 'x': 'known-homozygous'}, {'y': 0, 'x': 'known-all variants'}, {'y': 0, 'x': 'putative-homozygous'}, {'y': 11069531, 'x': 'putative-all variants'}], 'key': 'Drugs Used In Substance Dependence'}, {'values': [{'y': 0, 'x': 'known-homozygous'}, {'y': 0, 'x': 'known-all variants'}, {'y': 0, 'x': 'putative-homozygous'}, {'y': 8694786, 'x': 'putative-all variants'}], 'key': 'Hypothalamic&Pituitary Hormones&Antioest'}, {'values': [{'y': 0, 'x': 'known-homozygous'}, {'y': 0, 'x': 'known-all variants'}, {'y': 0, 'x': 'putative-homozygous'}, {'y': 9855456, 'x': 'putative-all variants'}], 'key': 'Sex Hormones & Antag In Malig Disease'}, {'values': [{'y': 0, 'x': 'known-homozygous'}, {'y': 0, 'x': 'known-all variants'}, {'y': 7848808, 'x': 'putative-homozygous'}, {'y': 25446045, 'x': 'putative-all variants'}], 'key': 'Treatment Of Glaucoma'}, {'values': [{'y': 864112, 'x': 'known-homozygous'}, {'y': 6107013, 'x': 'known-all variants'}, {'y': 19047162, 'x': 'putative-homozygous'}, {'y': 15754588, 'x': 'putative-all variants'}], 'key': 'other'}]

    ### PER DRUG TABLE

    ## drug data
    nhs_sections = NHSPrescribings.objects.all().values("drugname__name", "bnf_section").distinct()
    section_dict = {}
    for drug in nhs_sections:
        if drug['drugname__name'] in section_dict:
            section_dict[drug['drugname__name']].append(drug['bnf_section'])
        else:
            section_dict[drug['drugname__name']] = [drug['bnf_section']]

    nhs_data = NHSPrescribings.objects.all().values('drugname__name').annotate(Avg('actual_cost'), Avg('items'), Avg('quantity'))

    drug_data = []
    temp = {}
    for i in nhs_data:
        ## druginformation
        drugname = i['drugname__name']
        average_cost = int(i['actual_cost__avg'])
        average_quantity = int(i['quantity__avg'])
        average_items = int(i['items__avg'])
        section = section_dict[drugname]

        if average_items > 0:
            item_cost= round(float(average_cost)/average_items,1)
        else:
            item_cost = 0

        ## get target information
        protein_targets = Protein.objects.filter(drugs__name=drugname).distinct()
        targets = [p.entry_name.split('_human')[0].upper() for p in list(protein_targets)]

        known_functional = 0
        for target in protein_targets:
            if target.entry_name in temp:
                known_functional += temp[target.entry_name]
            else:
                function_sites = get_functional_sites(target)
                known_functional += function_sites
                temp[target.entry_name] = function_sites

        putative_func = len(NaturalMutations.objects.filter(Q(protein__in=protein_targets), Q(sift_score__lte=0.05) | Q(polyphen_score__gte=0.1)).annotate(count_putative_func=Count('id')))

        jsondata = {'drugname':drugname, 'targets': targets, 'average_cost': average_cost, 'average_quantity': average_quantity, 'average_items':average_items, 'item_cost':item_cost, 'known_func': known_functional, 'putative_func':putative_func, 'section':section}
        drug_data.append(jsondata)


    return render(request, 'economicburden.html', {'data':economic_data, 'drug_data':drug_data})
