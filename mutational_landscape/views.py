from django.shortcuts import get_object_or_404, render
from django.http import HttpResponse
from django.core.cache import cache
from django.views.decorators.cache import cache_page

from protein.models import Protein, ProteinConformation, ProteinAlias, ProteinFamily, Gene, ProteinGProtein, ProteinGProteinPair
from residue.models import Residue, ResiduePositionSet
from mutational_landscape.models import NaturalMutations, CancerMutations, DiseaseMutations

from mutation.functions import *
from mutation.models import *

from common import definitions
from collections import OrderedDict
from common.views import AbsTargetSelection

import json
# Create your views here.

def ajaxNaturalMutation(request, slug, **response_kwargs):

    name_of_cache = 'ajaxNaturalMutation_'+slug

    jsondata = cache.get(name_of_cache)

    if jsondata == None:
        jsondata = {}

        NMs = NaturalMutations.objects.filter(protein__entry_name=slug).prefetch_related('residue')

        for NM in NMs:

            SN = NM.residue.sequence_number
            # account for multiple mutations at this position!
            jsondata[SN] = [NM.amino_acid, NM.allele_frequency, NM.allele_count, NM.allele_number, NM.number_homozygotes]

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
