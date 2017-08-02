from django.shortcuts import get_object_or_404, render
from django.http import HttpResponse
from django.core.cache import cache
from django.db.models import Count, Min, Sum, Avg, Q
from django.core.cache import cache
from django.views.decorators.cache import cache_page

from protein.models import Protein, ProteinConformation, ProteinAlias, ProteinFamily, Gene, ProteinGProtein, ProteinGProteinPair
from residue.models import Residue, ResiduePositionSet, ResidueSet
from mutational_landscape.models import NaturalMutations, CancerMutations, DiseaseMutations, PTMs

from common.diagrams_gpcr import DrawHelixBox, DrawSnakePlot

from mutation.functions import *
from mutation.models import *

from interaction.models import *
from interaction.views import ajax #import x-tal interactions

from common import definitions
from collections import OrderedDict
from common.views import AbsTargetSelection
from common.views import AbsSegmentSelection

import re
import json
import numpy as np
from collections import OrderedDict
from copy import deepcopy

from io import BytesIO
import re
import math
import urllib
import xlsxwriter #sudo pip3 install XlsxWriter
import operator


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

def render_variants(request, protein = None, family = None, download = None, receptor_class = None, gn = None, aa = None, **response_kwargs):
    simple_selection = request.session.get('selection', False)
    proteins = []

    if protein: # if protein static page
        proteins.append(Protein.objects.get(entry_name = protein))

    # flatten the selection into individual proteins

    if simple_selection:
        for target in simple_selection.targets:
            if target.type == 'protein':
                proteins.append(target.item)
            elif target.type == 'family':
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

    NMs = NaturalMutations.objects.filter(Q(protein__in=proteins)).prefetch_related('residue')
    ptms = PTMs.objects.filter(Q(protein__in=proteins)).prefetch_related('residue')
    ptms_dict = {}

    for ptm in ptms:
        ptms_dict[ptm.residue.sequence_number] = ptm.modification

    ## G PROTEIN INTERACTION POSITIONS
    rset = ResiduePositionSet.objects.get(name='Signalling protein pocket')
    gprotein_generic_set = []
    for residue in rset.residue_position.all():
        gprotein_generic_set.append(residue.label)

    ### GET LB INTERACTION DATA
    # get also ortholog proteins, which might have been crystallised to extract
    # interaction data also from those
    orthologs = Protein.objects.filter(family__slug__startswith=proteins[0].family.slug, sequence_type__slug='wt')

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

    jsondata = {}
    for NM in NMs:
        functional_annotation = ''
        SN = NM.residue.sequence_number
        if NM.residue.generic_number:
            GN = NM.residue.generic_number.label
        else:
            GN = ''

        if SN in ptms_dict:
            functional_annotation +=  'PTM (' + ptms_dict[SN] + ')'
        if SN in interaction_data:
            functional_annotation +=  'LB (' + ', '.join(interaction_data[SN]) + ')'
        if GN in gprotein_generic_set:
            functional_annotation +=  'GP (contact)'

        type = NM.type
        if type == 'missense':
            effect = 'deleterious' if NM.sift_score <= 0.05 or NM.polyphen_score >= 0.1 else 'tolerated'
            color = '#e30e0e' if NM.sift_score <= 0.05 or NM.polyphen_score >= 0.1 else '#70c070'
        else:
            effect = 'deleterious'
            color = '#575c9d'
        # account for multiple mutations at this position!
        NM.functional_annotation = functional_annotation
        # print(NM.functional_annotation)
        jsondata[SN] = [NM.amino_acid, NM.allele_frequency, NM.allele_count, NM.allele_number, NM.number_homozygotes, NM.type, effect, color, functional_annotation]

    residuelist = Residue.objects.filter(protein_conformation__protein__entry_name=proteins[0].entry_name).prefetch_related('protein_segment','display_generic_number','generic_number')
    SnakePlot = DrawSnakePlot(
                residuelist, "Class A", protein, nobuttons=1)
    HelixBox = DrawHelixBox(residuelist,'Class A', protein, nobuttons = 1)

    if download:

        data = []
        for r in NMs:
            values = r.__dict__
            print(values)
            data.append(values)
        headers = ['amino_acid', 'allele_count','allele_number', 'allele_frequency', 'polyphen_score', 'sift_score', 'number_homozygotes']

        #EXCEL SOLUTION
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

        response = HttpResponse(xlsx_data,content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
        response['Content-Disposition'] = 'attachment; filename=GPCRdb_'+proteins[0].entry_name+'_variant_data.xlsx' #% 'mutations'
        return response
    return render(request, 'browser.html', {'mutations': NMs, 'HelixBox':HelixBox, 'SnakePlot':SnakePlot, 'receptor': str(proteins[0].entry_name),'mutations_pos_list' : json.dumps(jsondata)})

def ajaxNaturalMutation(request, slug, **response_kwargs):

    name_of_cache = 'ajaxNaturalMutation_'+slug

    jsondata = cache.get(name_of_cache)

    if jsondata == None:
        jsondata = {}

        NMs = NaturalMutations.objects.filter(protein__entry_name=slug).prefetch_related('residue')

        for NM in NMs:

            SN = NM.residue.sequence_number
            type = NM.type

            if type == 'missense':
                effect = 'deleterious' if NM.sift_score <= 0.05 or NM.polyphen_score >= 0.1 else 'tolerated'
                color = '#e30e0e' if NM.sift_score <= 0.05 or NM.polyphen_score >= 0.1 else '#70c070'
            else:
                effect = 'deleterious'
                color = '#575c9d'
            # account for multiple mutations at this position!
            jsondata[SN] = [NM.amino_acid, NM.allele_frequency, NM.allele_count, NM.allele_number, NM.number_homozygotes, NM.type, effect, color]

        jsondata = json.dumps(jsondata)
        response_kwargs['content_type'] = 'application/json'

        cache.set(name_of_cache, jsondata, 20) # 60*60*24*2 two days timeout on cache

    return HttpResponse(jsondata, **response_kwargs)

def ajaxCancerMutation(request, slug, **response_kwargs):

    name_of_cache = 'ajaxCancerMutation_'+slug

    jsondata = cache.get(name_of_cache)

    if jsondata == None:
        jsondata = {}

        CMs = CancerMutations.objects.filter(protein__entry_name=slug).prefetch_related('residue')

        for CM in CMs:
            SN = CM.residue.sequence_number
            jsondata[SN] = [CM.amino_acid]

        jsondata = json.dumps(jsondata)
        response_kwargs['content_type'] = 'application/json'

        cache.set(name_of_cache, jsondata, 20) #two days timeout on cache

    return HttpResponse(jsondata, **response_kwargs)

def ajaxDiseaseMutation(request, slug, **response_kwargs):

    name_of_cache = 'ajaxDiseaseMutation_'+slug

    jsondata = cache.get(name_of_cache)

    if jsondata == None:
        jsondata = {}

        DMs = DiseaseMutations.objects.filter(protein__entry_name=slug).prefetch_related('residue')

        for DM in DMs:
            SN = DM.residue.sequence_number
            jsondata[SN] = [DM.amino_acid]

        jsondata = json.dumps(jsondata)
        response_kwargs['content_type'] = 'application/json'

        cache.set(name_of_cache, jsondata, 20) #two days timeout on cache

    return HttpResponse(jsondata, **response_kwargs)

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

# @cache_page(60*60*24*2) #  2 days
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
    total_av_rv = round(len(NaturalMutations.objects.filter(type='missense', allele_frequency__lt=0.001))/ total_receptors,1)
    total_av_cv = round(len(NaturalMutations.objects.filter(type='missense', allele_frequency__gte=0.001))/ total_receptors,1)
    context['stats'] = {'total_mv':total_mv,'total_lof':total_lof,'total_av_rv':total_av_rv, 'total_av_cv':total_av_cv}

    return render(request, 'variation_statistics2.html', context)

def economicburden(request):
    data = [{'values': [{'y': 0.886, 'x': 1}], 'key': 'analgesics', 'yAxis': 'Scaling factor 1'},
    {'values': [{'y': 0.118, 'x': 1}], 'key': 'antidepressant drugs', 'yAxis': 'Scaling factor 1'},
    {'values': [{'y': 0.203, 'x': 1}], 'key': 'beta-adrenoceptor blocking drugs', 'yAxis': 'Scaling factor 1'},
    {'values': [{'y': 0.524, 'x': 1}], 'key': 'bronchodilators', 'yAxis': 'Scaling factor 1'},
    {'values': [{'y': 0.312, 'x': 1}], 'key': 'drugs used in diabetes', 'yAxis': 'Scaling factor 1'},
    {'values': [{'y': 0.080, 'x': 1}], 'key': 'drugs used in parkinson/related disorders', 'yAxis': 'Scaling factor 1'},
    {'values': [{'y': 0.205, 'x': 1}], 'key': 'drugs used in psychoses & relelated disorders', 'yAxis': 'Scaling factor 1'},
    {'values': [{'y': 0.102, 'x': 1}], 'key': 'drugs used in substance dependence', 'yAxis': 'Scaling factor 1'},
    {'values': [{'y': 0.080, 'x': 1}], 'key': 'hormones & antagonists in malignant disease', 'yAxis': 'Scaling factor 1'},
    {'values': [{'y': 0.072, 'x': 1}], 'key': 'hypothalamic & pituitary hormones', 'yAxis': 'Scaling factor 1'},
    {'values': [{'y': 0.257, 'x': 1}], 'key': 'other', 'yAxis': 'Scaling factor 1'}]

    return render(request, 'economicburden.html', {'data':data})
