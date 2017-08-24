from django.shortcuts import render
from django.views.generic import TemplateView, View
from django.http import HttpResponse
from django.views.decorators.cache import cache_page
from django.views.decorators.csrf import csrf_exempt
from django.conf import settings


from common.diagrams_gpcr import DrawSnakePlot
from common.definitions import AA_PROPENSITY, HYDROPHOBICITY
from common.views import AbsTargetSelection
from common.definitions import AMINO_ACIDS, AMINO_ACID_GROUPS, STRUCTURAL_RULES
from construct.models import *
from construct.functions import *
from construct.tool import *
from protein.models import Protein, ProteinConformation, ProteinSegment
from structure.models import Structure
from mutation.models import Mutation
from residue.models import ResiduePositionSet



from datetime import datetime
import json
import copy
import re
from collections import OrderedDict

Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')


# Create your views here.
#@cache_page(60 * 60 * 24)
def detail(request, slug):

    # get constructs
    c = Construct.objects.filter(name=slug).prefetch_related(
                "crystal","mutations","purification","protein__family__parent__parent__parent", "insertions__insert_type","expression","solubilization", "modifications", "deletions",
                "crystallization__crystal_method", "crystallization__crystal_type",
                "crystallization__chemical_lists", "crystallization__chemical_lists__chemicals__chemical__chemical_type",
                "protein__species","structure__pdb_code","structure__publication__web_link", "contributor").all()[0]

    # get residues
    residues = Residue.objects.filter(protein_conformation__protein=c.protein).order_by('sequence_number').prefetch_related(
        'protein_segment', 'generic_number', 'display_generic_number')

    residues_lookup = {}
    for r in residues:
        residues_lookup[r.sequence_number] = r

    schematics = c.schematic()

    snake = c.snake()

    # snake = cache.get(c.name+'_snake')
    # snake= None
    # if snake==None:
    #     print(c.name+'_snake no cache')
    #     snake = cache.set(c.name+'_snake', DrawSnakePlot(residues,c.protein.get_protein_class(),str(c.protein),nobuttons = True), 60*60*24*2) #two days
    #     snake = cache.get(c.name+'_snake')
    # else:
    #     print(c.name+'_snake used cache')

    chunk_size = 10
    context = {'c':c, 'chunk_size': chunk_size, 'snake': snake, 'annotations': json.dumps(schematics['annotations']), 'schematics': schematics, 'residues_lookup': residues_lookup}
    return render(request,'construct/construct_detail.html',context)

class ConstructStatistics(TemplateView):
    """
    Fetching construct data for browser
    """

    template_name = "construct/statistics.html"

    def get_context_data (self, **kwargs):

        context = super(ConstructStatistics, self).get_context_data(**kwargs)
        cons = Construct.objects.all().order_by("protein__family__slug","protein__entry_name").prefetch_related(
            "crystal","mutations","purification","protein__family__parent__parent__parent", "insertions__insert_type", "modifications", "deletions", "crystallization__chemical_lists",
            "protein__species","structure__pdb_code","structure__publication__web_link", "contributor")

        #PREPARE DATA
        proteins = Construct.objects.all().values_list('protein', flat = True)
        pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).filter(residue__generic_number__label__in=['1x50','8x50','5x50','6x50','3x50','4x50']).values_list('protein__entry_name','residue__sequence_number','residue__generic_number__label')
        #print(pconfs)
        x50s = {}
        for pc in pconfs:
            if pc[0] not in x50s:
                x50s[pc[0]] = {}
            x50s[pc[0]][pc[2]] = pc[1]

        pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).filter(residue__protein_segment__slug__in=['TM3','TM4','TM5','TM6']).values('protein__entry_name','residue__protein_segment__slug').annotate(start=Min('residue__sequence_number'),GN=Max('residue__generic_number__label'),GN2=Min('residue__generic_number__label'),end=Max('residue__sequence_number'))
        # print(pconfs)
        # x50s = {}
        track_anamalities = {}
        for pc in pconfs:
            #print(pc)
            entry_name = pc['protein__entry_name']
            helix = pc['residue__protein_segment__slug'][-1]
            if entry_name not in track_anamalities:
                track_anamalities[entry_name] = {}
            if helix not in track_anamalities[entry_name]:
                track_anamalities[entry_name][helix] = [0,0]
            x50 = x50s[entry_name][helix+"x50"]
            gn_start = int(pc['GN2'][-2:])
            gn_end  = int(pc['GN'][-2:])
            seq_start = pc['start']
            seq_end = pc['end']
            seq_range_start = x50-seq_start
            seq_range_end = seq_end-x50
            gn_range_start = 50-gn_start
            gn_range_end = gn_end-50
            if seq_range_start!=gn_range_start:
                # print(entry_name,"Helix",helix, "has anamolity in start",gn_range_start-seq_range_start)
                track_anamalities[entry_name][helix][0] = gn_range_start-seq_range_start
            if seq_range_end!=gn_range_end:
                # print(entry_name,"Helix",helix, "has anamolity in end",gn_range_end-seq_range_end)
                track_anamalities[entry_name][helix][1] = gn_range_end-seq_range_end
            #print(pc,helix,x50,gn_start,gn_end,seq_start,seq_end,,x50-seq_start,50-gn_start,gn_end-50)
            #print(x50s[entry_name])
            # if pc[0] not in x50s:
            #     x50s[pc[0]] = {}
            # x50s[pc[0]][pc[2]] = pc[1]
        # print(track_anamalities)
        pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).prefetch_related('protein').filter(residue__protein_segment__slug='TM1').annotate(start=Min('residue__sequence_number'))
        #pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).filter(residue__generic_number__label__in=['1x50']).values_list('protein__entry_name','residue__sequence_number','residue__generic_number__label')
        tm1_start = {}
        for pc in pconfs:
            tm1_start[pc.protein.entry_name] = pc.start

        pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).prefetch_related('protein').filter(residue__protein_segment__slug='C-term').annotate(start=Min('residue__sequence_number'),end=Max('residue__sequence_number'))
        cterm_start = {}
        cterm_end = {}
        for pc in pconfs:
            cterm_start[pc.protein.entry_name] = pc.start
            cterm_end[pc.protein.entry_name] = pc.end

        #GRAB RESIDUES for mutations
        mutations = []
        positions = []
        proteins = []
        for c in cons:
            p = c.protein
            entry_name = p.entry_name
            p_class = p.family.slug.split('_')[0]
            pdb = c.crystal.pdb_code
            pdb = '' # do not count same mutation many times
            for mutation in c.mutations.filter(effects__slug='thermostabilising').all():
                if p.entry_name not in proteins:
                    proteins.append(entry_name)
                mutations.append((mutation,entry_name,pdb,p_class))
                if mutation.sequence_number not in positions:
                    positions.append(mutation.sequence_number)
        rs = Residue.objects.filter(protein_conformation__protein__entry_name__in=proteins, sequence_number__in=positions).prefetch_related('generic_number','protein_conformation__protein')

        rs_lookup = {}
        gns = []
        for r in rs:
            if not r.generic_number:
                continue #skip non GN
            entry_name = r.protein_conformation.protein.entry_name
            pos = r.sequence_number
            if entry_name not in rs_lookup:
                rs_lookup[entry_name] = {}
            if pos not in rs_lookup[entry_name]:
                rs_lookup[entry_name][pos] = r

        truncations = {}
        truncations_new = {}
        truncations['nterm'] = {}
        truncations['nterm_fusion'] = {}
        truncations_new['nterm'] = OrderedDict()
        truncations_new['cterm'] = OrderedDict()
        truncations_new['nterm_fusion'] = OrderedDict()
        truncations_new['icl3_fusion'] = OrderedDict()
        truncations_new['icl2_fusion'] = OrderedDict()
        track_fusions = OrderedDict()
        track_fusions2 = OrderedDict()
        track_without_fusions = OrderedDict()
        truncations_new_possibilties = {}
        truncations_new_sum = {}
        truncations['cterm'] = {}
        truncations['icl3'] = {}
        truncations['icl3_fusion'] = {}
        truncations['icl2'] = {}
        truncations['icl2_fusion'] = {}
        class_names = {}
        for c in cons:
            p = c.protein
            entry_name = p.entry_name
            state = c.structure.state.slug
            #print(c.structure.state.slug)
            p_class = p.family.slug.split('_')[0]
            if p_class not in class_names:
                class_names[p_class] =  re.sub(r'\([^)]*\)', '', p.family.parent.parent.parent.name)
            p_class_name = class_names[p_class]
            if state=='active':
                p_class_name += "_active"
            fusion_n = False
            fusion_icl3 = False

            fusion_position, fusions = c.fusion()
            found_nterm = False
            found_cterm = False

            if p_class_name not in track_fusions:
                track_fusions[p_class_name] = OrderedDict()

            if entry_name not in track_fusions[p_class_name]:
                track_fusions[p_class_name][entry_name] = {'found':[],'for_print':[]}

            if fusions:
                fusion_name = fusions[0][2]
                if fusion_name not in track_fusions2:
                    track_fusions2[fusion_name] = {'found':[],'for_print':[]}
            # if entry_name=='aa2ar_human':
            #     print(state,p_class_name)
            for deletion in c.deletions.all():
                # if entry_name=='aa2ar_human':
                #     print(entry_name,deletion.start,cterm_start[entry_name],c.name) # lpar1_human

                if deletion.end <= x50s[entry_name]['1x50']:
                    found_nterm = True
                    bw = "1."+str(50-x50s[entry_name]['1x50']+deletion.end)
                    #bw = bw + " " + str(x50s[entry_name]['1x50']-deletion.end)
                    from_tm1 = tm1_start[entry_name] - deletion.end
                    # if entry_name=='aa2ar_human':
                    #     print(from_tm1,entry_name)

                    position = 'nterm'
                    if fusion_position=='nterm':
                        position = 'nterm_fusion'
                        if from_tm1 not in track_fusions[p_class_name][entry_name]['found']:
                            track_fusions[p_class_name][entry_name]['found'].append(from_tm1)
                        if from_tm1 not in track_fusions2[fusion_name]['found']:
                            track_fusions2[fusion_name]['found'].append(from_tm1)

                    if p_class_name not in truncations[position]:
                        truncations[position][p_class_name] = {}
                    if bw not in truncations[position][p_class_name]:
                        truncations[position][p_class_name][bw] = []
                    if entry_name not in truncations[position][p_class_name][bw]:
                        truncations[position][p_class_name][bw].append(entry_name)

                    if position not in truncations_new_possibilties:
                        truncations_new_possibilties[position] = []
                    if from_tm1 not in truncations_new_possibilties[position]:
                        truncations_new_possibilties[position].append(from_tm1)
                        truncations_new_possibilties[position] = sorted(truncations_new_possibilties[position])

                    if position not in truncations_new_sum:
                        truncations_new_sum[position] = {}
                    if p_class_name not in truncations_new_sum[position]:
                        truncations_new_sum[position][p_class_name] = {}

                    if p_class_name not in truncations_new[position]:
                        truncations_new[position][p_class_name] = {'receptors':OrderedDict(),'no_cut':[], 'possiblities':[]}
                    if entry_name not in truncations_new[position][p_class_name]['receptors']:
                        truncations_new[position][p_class_name]['receptors'][entry_name] = [[],[],[tm1_start[entry_name]]]
                    if fusion_position!='nterm':
                        if from_tm1 not in truncations_new[position][p_class_name]['receptors'][entry_name][0]:
                            truncations_new[position][p_class_name]['receptors'][entry_name][0].append(from_tm1)
                            if from_tm1 not in truncations_new_sum[position][p_class_name]:
                                truncations_new_sum[position][p_class_name][from_tm1] = 0
                            truncations_new_sum[position][p_class_name][from_tm1] += 1
                    # if from_tm1 not in truncations_new[position][p_class_name]['possiblities']:
                    #     truncations_new[position][p_class_name]['possiblities'].append(from_tm1)
                    #     truncations_new[position][p_class_name]['possiblities'] = sorted(truncations_new[position][p_class_name]['possiblities'])
                    # if from_tm1==0:
                    #     print(state,entry_name,p_class_name,truncations_new[position][p_class_name]['receptors'][entry_name])

                if deletion.start >= x50s[entry_name]['8x50']:
                    found_cterm = True

                    bw = x50s[entry_name]['8x50']-deletion.start
                    bw = "8."+str(50-x50s[entry_name]['8x50']+deletion.start)

                    from_h8 = deletion.start - cterm_start[entry_name]

                    if p_class_name not in truncations['cterm']:
                        truncations['cterm'][p_class_name] = {}
                    if bw not in truncations['cterm'][p_class_name]:
                        truncations['cterm'][p_class_name][bw] = []
                    if entry_name not in truncations['cterm'][p_class_name][bw]:
                        truncations['cterm'][p_class_name][bw].append(entry_name)

                    position = 'cterm'
                    if deletion.start>1000:
                        #TODO there are some wrong ones, can be seen by having >1000 positions which are fusion
                        continue
                        print(deletion.start,from_h8,cterm_start[entry_name],c.crystal.pdb_code )

                    if position not in truncations_new_possibilties:
                        truncations_new_possibilties[position] = []
                    if from_h8 not in truncations_new_possibilties[position]:
                        truncations_new_possibilties[position].append(from_h8)
                        truncations_new_possibilties[position] = sorted(truncations_new_possibilties[position])

                    if position not in truncations_new_sum:
                        truncations_new_sum[position] = {}
                    if p_class_name not in truncations_new_sum[position]:
                        truncations_new_sum[position][p_class_name] = {}




                    if p_class_name not in truncations_new[position]:
                        truncations_new[position][p_class_name] = {'receptors':OrderedDict(),'no_cut':[], 'possiblities':[]}
                    if entry_name not in truncations_new[position][p_class_name]['receptors']:
                        truncations_new[position][p_class_name]['receptors'][entry_name] = [[],[],[cterm_end[entry_name]-cterm_start[entry_name]]]
                    if from_h8 not in truncations_new[position][p_class_name]['receptors'][entry_name][0]:
                        truncations_new[position][p_class_name]['receptors'][entry_name][0].append(from_h8)
                        if from_h8 not in truncations_new_sum[position][p_class_name]:
                            truncations_new_sum[position][p_class_name][from_h8] = 0
                        truncations_new_sum[position][p_class_name][from_h8] += 1

                if deletion.start > x50s[entry_name]['5x50'] and deletion.start < x50s[entry_name]['6x50']:
                    # if fusion_icl3:
                    #      print(entry_name,c.name,deletion.start,deletion.end,x50s[entry_name]['5x50'])
                    fusion_icl3 = True
                    bw = x50s[entry_name]['5x50']-deletion.start
                    bw = "5x"+str(50-x50s[entry_name]['5x50']+deletion.start+track_anamalities[entry_name]['5'][1])
                    bw2 = "6x"+str(50-x50s[entry_name]['6x50']+deletion.end+track_anamalities[entry_name]['6'][0])
                    bw_combine = bw+"-"+bw2
                    position = 'icl3'
                    if fusion_position=='icl3':
                        position = 'icl3_fusion'
                        if bw not in track_fusions2[fusion_name]['found']:
                            track_fusions2[fusion_name]['found'].append(bw)
                        if bw2 not in track_fusions2[fusion_name]['found']:
                            track_fusions2[fusion_name]['found'].append(bw2)

                    if fusion_position=='icl3':
                        #Track those with fusion
                        if bw not in track_fusions[p_class_name][entry_name]['found']:
                            track_fusions[p_class_name][entry_name]['found'].append(bw)
                        if bw2 not in track_fusions[p_class_name][entry_name]['found']:
                            track_fusions[p_class_name][entry_name]['found'].append(bw2)
                    else:
                        print('ICL3 CUT WITHOUT FUSION',bw_combine,entry_name,c.name)
                        if p_class_name not in track_without_fusions:
                            track_without_fusions[p_class_name] = OrderedDict()

                        if entry_name not in track_without_fusions[p_class_name]:
                            track_without_fusions[p_class_name][entry_name] = {'found':[],'for_print':[]}

                        #Track those without fusion
                        if bw not in track_without_fusions[p_class_name][entry_name]['found']:
                            track_without_fusions[p_class_name][entry_name]['found'].append(bw)
                        if bw2 not in track_without_fusions[p_class_name][entry_name]['found']:
                            track_without_fusions[p_class_name][entry_name]['found'].append(bw2)


                    if p_class_name not in truncations[position]:
                        truncations[position][p_class_name] = {}
                    if bw_combine not in truncations[position][p_class_name]:
                        truncations[position][p_class_name][bw_combine] = []
                    if entry_name not in truncations[position][p_class_name][bw_combine]:
                        truncations[position][p_class_name][bw_combine].append(entry_name)


                    if position+"_start" not in truncations_new_possibilties:
                        truncations_new_possibilties[position+"_start"] = []
                    if position+"_end" not in truncations_new_possibilties:
                        truncations_new_possibilties[position+"_end"] = []
                    if bw not in truncations_new_possibilties[position+"_start"]:
                        truncations_new_possibilties[position+"_start"].append(bw)
                        truncations_new_possibilties[position+"_start"] = sorted(truncations_new_possibilties[position+"_start"])
                    if bw2 not in truncations_new_possibilties[position+"_end"]:
                        truncations_new_possibilties[position+"_end"].append(bw2)
                        truncations_new_possibilties[position+"_end"] = sorted(truncations_new_possibilties[position+"_end"])
                    


                if deletion.start > x50s[entry_name]['3x50'] and deletion.start < x50s[entry_name]['4x50']:
                    # if fusion_icl3:
                    #      print(entry_name,c.name,deletion.start,deletion.end,x50s[entry_name]['5x50'])
                    fusion_icl3 = True
                    bw = x50s[entry_name]['5x50']-deletion.start
                    bw = "3x"+str(50-x50s[entry_name]['3x50']+deletion.start+track_anamalities[entry_name]['3'][1])
                    bw2 = "4x"+str(50-x50s[entry_name]['4x50']+deletion.end+track_anamalities[entry_name]['4'][0])
                    bw_combine = bw+"-"+bw2
                    position = 'icl2'
                    if fusion_position=='icl3':
                        position = 'icl2_fusion'
                        if bw not in track_fusions2[fusion_name]['found']:
                            track_fusions2[fusion_name]['found'].append(bw)
                        if bw2 not in track_fusions2[fusion_name]['found']:
                            track_fusions2[fusion_name]['found'].append(bw2)


                    if bw not in track_fusions[p_class_name][entry_name]['found']:
                            track_fusions[p_class_name][entry_name]['found'].append(bw)
                    if bw2 not in track_fusions[p_class_name][entry_name]['found']:
                            track_fusions[p_class_name][entry_name]['found'].append(bw2)

                    if p_class_name not in truncations[position]:
                        truncations[position][p_class_name] = {}
                    if bw_combine not in truncations[position][p_class_name]:
                        truncations[position][p_class_name][bw_combine] = []
                    if entry_name not in truncations[position][p_class_name][bw_combine]:
                        truncations[position][p_class_name][bw_combine].append(entry_name)

                    if position+"_start" not in truncations_new_possibilties:
                        truncations_new_possibilties[position+"_start"] = []
                    if position+"_end" not in truncations_new_possibilties:
                        truncations_new_possibilties[position+"_end"] = []
                    if bw not in truncations_new_possibilties[position+"_start"]:
                        truncations_new_possibilties[position+"_start"].append(bw)
                        truncations_new_possibilties[position+"_start"] = sorted(truncations_new_possibilties[position+"_start"])
                    if bw2 not in truncations_new_possibilties[position+"_end"]:
                        truncations_new_possibilties[position+"_end"].append(bw2)
                        truncations_new_possibilties[position+"_end"] = sorted(truncations_new_possibilties[position+"_end"])

            position = 'nterm'
            if fusion_position=='nterm':
                position = 'nterm_fusion'
            if p_class_name not in truncations_new[position]:
                truncations_new[position][p_class_name] = {'receptors':OrderedDict(),'no_cut':[], 'possiblities':[]}


            # if entry_name=='aa2ar_human':
            #     print(found_nterm,entry_name,position,p_class_name)

            from_tm1 = tm1_start[entry_name]
            if not found_nterm:
                if position not in truncations_new_sum:
                    truncations_new_sum[position] = {}
                if p_class_name not in truncations_new_sum[position]:
                    truncations_new_sum[position][p_class_name] = {}
                if from_tm1 not in truncations_new_sum[position][p_class_name]:
                        truncations_new_sum[position][p_class_name][from_tm1] = 0
                #if full receptor in xtal
                if entry_name not in truncations_new[position][p_class_name]['receptors']:
                    truncations_new[position][p_class_name]['receptors'][entry_name] = [[],[from_tm1],[from_tm1]]

                    # add one for this position if it is first time receptor is mentioned
                    truncations_new_sum[position][p_class_name][from_tm1] += 1
                else:
                    if from_tm1 not in truncations_new[position][p_class_name]['receptors'][entry_name][1]:
                        truncations_new[position][p_class_name]['receptors'][entry_name][1].append(from_tm1)
                        truncations_new_sum[position][p_class_name][from_tm1] += 1
            # else:
            #     #if full was found, fill in the max
            #     #print(entry_name,found_nterm)
            #     if from_tm1 not in truncations_new[position][p_class_name]['receptors'][entry_name][1]:
            #         truncations_new[position][p_class_name]['receptors'][entry_name][2].append(from_tm1)

            if position!='nterm_fusion' and from_tm1 not in truncations_new_possibilties[position]:
                truncations_new_possibilties[position].append(from_tm1)
                truncations_new_possibilties[position] = sorted(truncations_new_possibilties[position])

            position = 'cterm'
            if p_class_name not in truncations_new[position]:
                truncations_new[position][p_class_name] = {'receptors':OrderedDict(),'no_cut':[], 'possiblities':[]}
            if position not in truncations_new_possibilties:
                truncations_new_possibilties[position] = []


            from_h8 = cterm_end[entry_name] - cterm_start[entry_name]
            if not found_cterm:
                if position not in truncations_new_sum:
                    truncations_new_sum[position] = {}
                if p_class_name not in truncations_new_sum[position]:
                    truncations_new_sum[position][p_class_name] = {}
                if from_h8 not in truncations_new_sum[position][p_class_name]:
                        truncations_new_sum[position][p_class_name][from_h8] = 0
                if entry_name not in truncations_new[position][p_class_name]['receptors']:
                    truncations_new[position][p_class_name]['receptors'][entry_name] = [[],[from_h8],[from_h8]]
                    # add one for this position if it is first time receptor is mentioned
                    truncations_new_sum[position][p_class_name][from_h8] += 1
                else:
                    if from_h8 not in truncations_new[position][p_class_name]['receptors'][entry_name][1]:
                        truncations_new[position][p_class_name]['receptors'][entry_name][1].append(from_h8)
                        truncations_new_sum[position][p_class_name][from_h8] += 1
            # else:
            #     if from_h8 not in truncations_new[position][p_class_name]['receptors'][entry_name][1]:
            #         truncations_new[position][p_class_name]['receptors'][entry_name][1].append(from_h8)

            if from_h8 not in truncations_new_possibilties[position]:
                truncations_new_possibilties[position].append(from_h8)
            truncations_new_possibilties[position] = sorted(truncations_new_possibilties[position])


        for pos, p_vals in truncations_new_sum.items():
            for pclass, c_vals in p_vals.items():
                new_list = OrderedDict()
                for position in truncations_new_possibilties[pos]:
                    if position in c_vals:
                        new_list[position] = c_vals[position]
                    else:
                        new_list[position] = ''
                # print(pclass,c_vals,new_list)
                truncations_new_sum[pos][pclass] = new_list
        # print(truncations_new)
        #truncations = OrderedDict(truncations)
        ordered_truncations = OrderedDict()
        for segment, s_vals in sorted(truncations.items()):
            #print(segment)
            ordered_truncations[segment] = OrderedDict()
            for p_class, c_vals in sorted(s_vals.items()):
                #print(p_class)
                ordered_truncations[segment][p_class] = OrderedDict()
                for pos, p_vals in sorted(c_vals.items(),key=lambda x: (len(x[1]),x[0]), reverse=True):
                    #print(pos, len(p_vals))
                    ordered_truncations[segment][p_class][pos] = p_vals


        fusion_possibilities = truncations_new_possibilties['nterm_fusion'][::-1] + truncations_new_possibilties['icl2_fusion_start'] + truncations_new_possibilties['icl3_fusion_start'] + truncations_new_possibilties['icl2_fusion_end'] + truncations_new_possibilties['icl3_fusion_end']
        # fusion_possibilities = truncations_new_possibilties['nterm_fusion'][::-1]  + truncations_new_possibilties['icl3_start'] + truncations_new_possibilties['icl3_end']
        track_fusion_sums = OrderedDict()
        track_without_fusion_sums = OrderedDict()
        for pclass, receptors in track_fusions.items():
            track_fusion_sums[pclass] = OrderedDict()
            for p in fusion_possibilities:
                track_fusion_sums[pclass][p] = 0
            for receptor, vals in receptors.items():
                temp = []
                for p in fusion_possibilities:
                    if p in vals['found']:
                        temp.append(1)
                        track_fusion_sums[pclass][p] += 1
                    else:
                        temp.append(0)
                vals['for_print'] = temp
        for pclass, receptors in track_without_fusions.items():
            track_without_fusion_sums[pclass] = OrderedDict()
            for p in truncations_new_possibilties['icl3_start'] + truncations_new_possibilties['icl3_end']:
                track_without_fusion_sums[pclass][p] = 0
            for receptor, vals in receptors.items():
                temp = []
                for p in truncations_new_possibilties['icl3_start'] + truncations_new_possibilties['icl3_end']:
                    if p in vals['found']:
                        temp.append(1)
                        track_without_fusion_sums[pclass][p] += 1
                    else:
                        temp.append(0)
                vals['for_print'] = temp
        # print(track_fusion_sums)
        for fusion, vals in track_fusions2.items():
            temp = []
            for p in fusion_possibilities:
                if p in vals['found']:
                    temp.append(1)
                else:
                    temp.append(0)
            vals['for_print'] = temp
        # print(track_without_fusions)
        #truncations =  OrderedDict(sorted(truncations.items(), key=lambda x: x[1]['hits'],reverse=True))
        #print(ordered_truncations)
        # print(track_fusions2)
        context['truncations'] = ordered_truncations
        context['truncations_new'] = truncations_new
        context['truncations_new_possibilties'] = truncations_new_possibilties
        context['truncations_new_sum'] = truncations_new_sum
        context['fusion_possibilities'] = fusion_possibilities
        context['test'] = track_fusions
        context['test2'] = track_fusions2
        context['track_fusion_sums'] = track_fusion_sums
        context['track_without_fusions'] = track_without_fusions

        mutation_list = OrderedDict()
        mutation_type = OrderedDict()
        mutation_wt = OrderedDict()
        mutation_mut = OrderedDict()

        mutation_matrix = OrderedDict()
        mutation_track = []
        aa_list = list(AMINO_ACIDS.keys())[:20]
        mutation_matrix_sum_mut = OrderedDict()
        #print(aa_list)

        for i, mut in enumerate(AMINO_ACIDS):
            if i==20:
                break
            mutation_matrix[mut] = OrderedDict()
            for aa in aa_list:
                mutation_matrix[mut][aa] = [0,0]
            mutation_matrix[mut][mut] = [0,'-']
            mutation_matrix[mut]['sum'] = [0,0]
            mutation_matrix_sum_mut[mut] = [0,0]
        #print(mutation_matrix)
        for mutation in mutations:
            wt = mutation[0].wild_type_amino_acid
            mut = mutation[0].mutated_amino_acid
            entry_name = mutation[1]
            pos = mutation[0].sequence_number
            p_class = mutation[3]
            p_class = class_names[p_class]
            pdb = mutation[2]

            mut_uniq = entry_name+'_'+str(pos)+'_'+wt+'_'+mut

            if mut_uniq not in mutation_track:
                # print(mut_uniq)
                #do not count the same mutation (from different Xtals) multiple times
                mutation_track.append(mut_uniq)
                mutation_matrix[wt][mut][1] += 1
                mutation_matrix[wt][mut][0] = min(1,round(mutation_matrix[wt][mut][1]/30,2))
                mutation_matrix[wt]['sum'][1] += 1
                mutation_matrix[wt]['sum'][0] = min(1,round(mutation_matrix[wt]['sum'][1]/30,2))
                mutation_matrix_sum_mut[mut][1] += 1
                mutation_matrix_sum_mut[mut][0] = min(1,round(mutation_matrix_sum_mut[mut][1]/30,2))
            # print(entry_name,"\t", pdb,"\t", pos,"\t", wt,"\t", mut)


            if p_class not in mutation_type:
                mutation_type[p_class] = OrderedDict()
            if wt+"=>"+mut not in mutation_type[p_class]:
                mutation_type[p_class][wt+"=>"+mut] = {'hits':0, 'proteins':[]}
            if entry_name not in mutation_type[p_class][wt+"=>"+mut]['proteins']:
                mutation_type[p_class][wt+"=>"+mut]['proteins'].append(entry_name)
                mutation_type[p_class][wt+"=>"+mut]['hits'] += 1


            if p_class not in mutation_wt:
                mutation_wt[p_class] = OrderedDict()
            if wt not in mutation_wt[p_class]:
                mutation_wt[p_class][wt] = {'hits':0, 'proteins':[]}
            if entry_name not in mutation_wt[p_class][wt]['proteins']:
                mutation_wt[p_class][wt]['proteins'].append(entry_name)
                mutation_wt[p_class][wt]['hits'] += 1

            if p_class not in mutation_mut:
                mutation_mut[p_class] = OrderedDict()
            if mut not in mutation_mut[p_class]:
                mutation_mut[p_class][mut] = {'hits':0, 'proteins':[]}
            if entry_name not in mutation_mut[p_class][mut]['proteins']:
                mutation_mut[p_class][mut]['proteins'].append(entry_name)
                mutation_mut[p_class][mut]['hits'] += 1

            if entry_name not in rs_lookup:
                continue
            if pos not in rs_lookup[entry_name]:
                continue
            gn = rs_lookup[entry_name][pos].generic_number.label

            if p_class not in mutation_list:
                mutation_list[p_class] = OrderedDict()
            if gn not in mutation_list[p_class]:
                mutation_list[p_class][gn] = {'proteins':[], 'hits':0, 'mutation':[]}
            if entry_name not in mutation_list[p_class][gn]['proteins']:
                mutation_list[p_class][gn]['proteins'].append(entry_name)
                mutation_list[p_class][gn]['hits'] += 1
                mutation_list[p_class][gn]['mutation'].append((mutation[0].wild_type_amino_acid,mutation[0].mutated_amino_acid))


        #print(mutation_matrix)

        for p_class, values in mutation_list.items():
            for gn, vals in values.items():
                if vals['hits']<2:
                    pass
                    #values.pop(gn, None)
            mutation_list[p_class] = OrderedDict(sorted(values.items(), key=lambda x: x[1]['hits'],reverse=True))
            #mutation_list = OrderedDict(sorted(mutation_list.items(), key=lambda x: x[1]['hits'],reverse=True))

        for p_class, values in mutation_type.items():
            mutation_type[p_class] = OrderedDict(sorted(values.items(), key=lambda x: x[1]['hits'],reverse=True))

        for p_class, values in mutation_wt.items():
            mutation_wt[p_class] = OrderedDict(sorted(values.items(), key=lambda x: x[1]['hits'],reverse=True))

        for p_class, values in mutation_mut.items():
            mutation_mut[p_class] = OrderedDict(sorted(values.items(), key=lambda x: x[1]['hits'],reverse=True))

        context['mutation_list'] = mutation_list
        context['mutation_type'] = mutation_type
        context['mutation_wt'] = mutation_wt
        context['mutation_mut'] = mutation_mut
        context['mutation_matrix'] = mutation_matrix
        context['mutation_matrix_sum_mut'] = mutation_matrix_sum_mut

        for c in cons:
            pass

        return context


class ConstructTable(TemplateView):
    """
    Fetching construct data for browser
    """

    template_name = "construct/residuetable.html"

    def get_context_data (self, **kwargs):

        context = super(ConstructTable, self).get_context_data(**kwargs)
        cons = Construct.objects.all().prefetch_related(
            "crystal","mutations","purification","protein__family__parent__parent__parent", "insertions__insert_type", "modifications", "deletions", "crystallization__chemical_lists",
            "protein__species","structure__pdb_code","structure__publication__web_link", "contributor")

        #PREPARE DATA
        proteins = Construct.objects.all().values_list('protein', flat = True)

        #GRAB RESIDUES for mutations
        mutations = []
        positions = []
        proteins = []
        class_names = {}
        classes = []
        for c in cons:
            p = c.protein
            entry_name = p.entry_name
            p_class = p.family.slug.split('_')[0]
            if p_class not in classes:
                classes.append(p_class)
            pdb = c.crystal.pdb_code
            for mutation in c.mutations.all():
                if p.entry_name not in proteins:
                    proteins.append(entry_name)
                mutations.append((mutation,entry_name,pdb,p_class,c.name))
                if mutation.sequence_number not in positions:
                    positions.append(mutation.sequence_number)
        rs = Residue.objects.filter(protein_conformation__protein__entry_name__in=proteins, sequence_number__in=positions).prefetch_related('generic_number','protein_conformation__protein','protein_segment')

        #print(classes)

        excluded_segment = ['C-term','N-term']
        list(settings.REFERENCE_POSITIONS.keys())
        align_segments = ProteinSegment.objects.all().filter(slug__in = list(settings.REFERENCE_POSITIONS.keys())).prefetch_related()

        amino_acids_stats = {}
        amino_acids_groups_stats = {}
        accessible_in_class = {}
        for c in classes:

            #if c=='001':
            #    continue
            amino_acids_stats[c] = {}
            amino_acids_groups_stats[c] = {}
            accessible_in_class[c] = []

            if c =='001':
                residue_set_name = 'Class A binding pocket'
            elif c=='004':
                residue_set_name = 'Class C binding pocket'
            elif c=='005':
                residue_set_name = 'Class F binding pocket'
            else:
                residue_set_name = ''

            if residue_set_name:
                rset = ResiduePositionSet.objects.get(name=residue_set_name)
                for residue in rset.residue_position.all():
                    accessible_in_class[c].append(residue.label)
            #print(accessible_in_class)

            alignment_proteins = Protein.objects.filter(family__slug__startswith=c, species__common_name='Human', source__name='SWISSPROT')
            #print(c,len(alignment_proteins))

            a = Alignment()

            a.load_proteins(alignment_proteins)

            a.load_segments(align_segments) #get all segments to make correct diagrams

            # build the alignment data matrix
            a.build_alignment()

            # calculate consensus sequence + amino acid and feature frequency
            a.calculate_statistics()

            #print(a.amino_acid_stats)
                        # {% for ns, segments in a.generic_numbers.items %}
                        # <tr>
                        #     {% for s, num in segments.items %}
                        #         {% for n, dn in num.items %}
                        #             {% if 'Common G-alpha numbering scheme' in a.numbering_schemes.0 %}
                        #             <td class="ali-td-generic-num">{{ dn|make_list|slice:'2:'|join:''}}</td>
                        #             {% else %}
                        #             <td class="ali-td-generic-num">{{ dn|safe }}</td>
                        #             {% endif %}
                        #         {% endfor %}
                        #         <td class="ali-td">&nbsp;</td>
                        #     {% endfor %}
                        # </tr>
                        # {% endfor %}
            s_id = 0
            a_id = 0
            for ns, segments in a.generic_numbers.items():
                for s, num in segments.items():
                    for n, dn in num.items():
                        temp = []
                        temp2 = []
                        for i, aa in enumerate(AMINO_ACIDS):
                            temp.append(a.amino_acid_stats[i][s_id][a_id])

                        for i, aa in enumerate(AMINO_ACID_GROUPS):
                            temp2.append(a.feature_stats[i][s_id][a_id])
                        amino_acids_stats[c][n] = temp
                        amino_acids_groups_stats[c][n] = temp2
                    a_id += 1
                s_id += 1

        rs_lookup = {}
        gns = []
        for r in rs:
            entry_name = r.protein_conformation.protein.entry_name
            pos = r.sequence_number
            segment = r.protein_segment.slug
            if entry_name not in rs_lookup:
                rs_lookup[entry_name] = {}
            if pos not in rs_lookup[entry_name]:
                rs_lookup[entry_name][pos] = r

        count_per_gn = {}
        for mutation in mutations:
            entry_name = mutation[1]
            pos = mutation[0].sequence_number
            if rs_lookup[entry_name][pos].generic_number:
                gn = rs_lookup[entry_name][pos].generic_number.label
                if gn not in count_per_gn:
                    count_per_gn[gn] = {'hits': 0, 'proteins': []}
                if entry_name not in count_per_gn[gn]['proteins']:
                    count_per_gn[gn]['proteins'].append(entry_name)
                    count_per_gn[gn]['hits'] += 1

        #print(count_per_gn)

        mutation_list = []
        for mutation in mutations:
            wt = mutation[0].wild_type_amino_acid
            mut = mutation[0].mutated_amino_acid
            entry_name = mutation[1]
            pdb = mutation[2]
            cname = mutation[4]
            pos = mutation[0].sequence_number
            p_class = mutation[3]
            if p_class not in class_names:
                class_names[p_class] = p.family.parent.parent.parent.name
            p_class_name = class_names[p_class]
            p_class = class_names[p_class]



            if entry_name not in rs_lookup:
                continue
            if pos not in rs_lookup[entry_name]:
                continue
            segment = rs_lookup[entry_name][pos].protein_segment.slug
            if rs_lookup[entry_name][pos].generic_number:
                gn = rs_lookup[entry_name][pos].generic_number.label
                stats = amino_acids_stats[mutation[3]][gn]
                stats2 = amino_acids_groups_stats[mutation[3]][gn]
                if gn in accessible_in_class[mutation[3]]:
                    accessible = 'yes'
                else:
                    accessible = 'no'
                count = count_per_gn[gn]['hits']
            else:
                gn = ''
                stats = ''
                stats2 = ''
                accessible = 'N/A'



            mutation_list.append({'entry_name':entry_name,'pdb':pdb,'cname':cname, 'segment':segment,'pos': pos, 'gn': gn, 'wt': wt, 'mut': mut,'p_class': p_class, 'amino_acids':stats, 'amino_acids_groups':stats2, 'accessible': accessible, 'count': count})


        context['amino_acids'] = AMINO_ACIDS
        context['amino_groups'] = AMINO_ACID_GROUPS
        context['mutation_list'] = mutation_list
        #print(mutation_list)

        return context

class ConstructMutations(TemplateView):
    """
    Fetching construct data for browser
    """

    template_name = "construct/mutations.html"

    def get_context_data (self, **kwargs):

        context = super(ConstructMutations, self).get_context_data(**kwargs)
        cons = Construct.objects.all().prefetch_related(
            "crystal","mutations__effects","purification","protein__family__parent__parent__parent", "insertions__insert_type", "modifications", "deletions", "crystallization__chemical_lists",
            "protein__species","structure__pdb_code","structure__publication__web_link", "contributor")

        #PREPARE DATA
        proteins = Construct.objects.all().values_list('protein', flat = True)

        #GRAB RESIDUES for mutations
        mutations = []
        positions = []
        proteins = []
        class_names = {}
        for c in cons:
            # print(c)
            p = c.protein
            entry_name = p.entry_name
            p_class = p.family.slug.split('_')[0]
            pdb = c.crystal.pdb_code
            for mutation in c.mutations.all():
                if p.entry_name not in proteins:
                    proteins.append(entry_name)
                mutations.append((mutation,entry_name,pdb,p_class,c.name))
                if mutation.sequence_number not in positions:
                    positions.append(mutation.sequence_number)
        rs = Residue.objects.filter(protein_conformation__protein__entry_name__in=proteins, sequence_number__in=positions).prefetch_related('generic_number','protein_conformation__protein','protein_segment')

        rs_lookup = {}
        gns = []
        for r in rs:
            # print("r",r)
            entry_name = r.protein_conformation.protein.entry_name
            pos = r.sequence_number
            segment = r.protein_segment.slug
            if entry_name not in rs_lookup:
                rs_lookup[entry_name] = {}
            if pos not in rs_lookup[entry_name]:
                rs_lookup[entry_name][pos] = r

        mutation_list = []
        for mutation in mutations:
            # print("m",mutation)
            wt = mutation[0].wild_type_amino_acid
            mut = mutation[0].mutated_amino_acid

            mut_type = mutation[0].effects.all().values('slug')
            mut_types = []
            for eff in mut_type:
                mut_types.append(eff['slug'])
            mut_type = ",".join(mut_types)
            # # for
            # print(mut_type)
            # mut_type = ''

            entry_name = mutation[1]
            pdb = mutation[2]
            cname = mutation[4]
            pos = mutation[0].sequence_number
            p_class = mutation[3]
            if p_class not in class_names:
                class_names[p_class] = p.family.parent.parent.parent.name
            p_class_name = class_names[p_class]
            p_class = class_names[p_class]


            if entry_name not in rs_lookup:
                continue
            if pos not in rs_lookup[entry_name]:
                continue
            segment = rs_lookup[entry_name][pos].protein_segment.slug
            if rs_lookup[entry_name][pos].generic_number:
                gn = rs_lookup[entry_name][pos].generic_number.label
            else:
                gn = ''


            mutation_list.append({'entry_name':entry_name,'pdb':pdb,'cname':cname, 'segment':segment,'pos': pos, 'gn': gn, 'wt': wt, 'mut': mut,'p_class': p_class, 'type': mut_type})
            if (pdb == 'adrb1_melga'): print(gn)


        context['mutation_list'] = mutation_list
        print(mutation_list)

        return context



def stabilisation_browser(request):
    ''' View to display and summarise mutation data for thermostabilising mutational constructs. '''

    # Set up: Restructure the STRUCTURAL_RULES for the constructs into a crude-tree like structure to enable
    # quick and concise searching within the for loops below.
    structural_rule_tree = create_structural_rule_trees(STRUCTURAL_RULES)

    # Get a list of all constructs.
    constructs = Construct.objects.all()\
            .order_by().only(
                "protein__entry_name",
                # "mutations__sequence_number",
                # "mutations__residue__generic_number",
                # "mutations__residue__protein_segment__slug",
                # "mutations__mutated_amino_acid",
                # "mutations__wild_type_amino_acid",
                "protein__family__slug",
                "protein__family__parent__parent__parent__name",
                "structure__state__name",
                "crystal__pdb_code")\
            .prefetch_related(
                "structure__state",
                "mutations__residue__generic_number",
                "mutations__residue__protein_segment",
                "mutations__effects__bar",
                "protein__family__parent__parent__parent",
                "crystal")
                
    # Get a list of all relevant proteins and generic numbers
    conservation_proteins = constructs.values_list('protein__family__parent__parent__parent__name',
                                                   flat=True)\
                            .distinct()
    conservation_gen_nums = constructs.values_list('mutations__residue__generic_number__label', flat=True).distinct()

    # Calculate the conservation values for the mutations across their receptor families and protein classes.
    # Alignment performed using generic numbers.
    conservation = conservation_table(conservation_proteins, conservation_gen_nums)


    # For each analysis mode, define the information that is to be used as a unique identifier for grouping data.
    # I.e. for position_only, grouping is performed by class an position.  Hence each row will have a unique class &
    # position.  This is used as the unique identifier, or ID. -- recorded in 'include_in_id'
    # Each row has some calculated or 'unique' values, as well as the id. This is found below.  However, for example,
    # the wild type AA is not unique accross the pos_and_mut group, as so this must be removed from the row-info.
    # This is recorded in 'exclude_from_info'
    groupings = {
        "all":{"include_in_id":['class', 'gen_num', 'wild_type', 'mutant'], "exclude_from_info":['']},
        "pos_and_wt":{"include_in_id":['class', 'gen_num', 'wild_type'],
                      "exclude_from_info":['ala_leu_subset', 'ala_subset', 'mutant']},
        "pos_and_mut":{"include_in_id":['class', 'gen_num', 'mutant'],
                       "exclude_from_info":['ala_leu_subset', 'wild_type']},
        "position_only":{"include_in_id":['class', 'gen_num'],
                         "exclude_from_info":['ala_leu_subset', 'ala_subset', 'wild_type', 'mutant']}
        }

    # Set up dictionaries to record information.
    mutation_groups = {"position_only":{}, "all":{}, "pos_and_wt":{}, "pos_and_mut":{}}

    # Grab thermostabilising mutations
    mutations_thermo = ConstructMutation.objects.filter(effects__slug='thermostabilising').all()\
                .prefetch_related(
                    "construct__structure__state",
                    "residue__generic_number",
                    "residue__protein_segment",
                    "construct__protein__family__parent__parent__parent",
                    "construct__crystal")


    # For each construct, get the needed information, and add to the context dictionary called mutation_list.
    for mutant in mutations_thermo:

        # Get info for the construct
        struct_id = mutant.construct.structure_id
        state = mutant.construct.structure.state.name
        prot = mutant.construct.protein
        p_class = prot.family.parent.parent.parent.name
        p_ligand = prot.family.parent.parent.name
        p_receptor = prot.family.parent.name
        pdb = mutant.construct.crystal.pdb_code

        # Get the generic number and segment, if known.
        try:
            if mutant.residue.generic_number is None:
                generic_number = u'\u2014'
            else:
                generic_number = mutant.residue.generic_number.label
            segment = mutant.residue.protein_segment.slug
        except AttributeError:
            generic_number = u'\u2014'
            segment = u'\u2014'

        # Collect the mutation info needed to create a unique group id, and the info relevant to the full row.
        mutant_id = {'gen_num':generic_number, 'wild_type':mutant.wild_type_amino_acid,
                     'mutant':mutant.mutated_amino_acid, 'GPCR_count':0, 'segment':segment, 'class': p_class}
        mutant_info = {'pdb':pdb,
                       'ligand': p_ligand,
                       'receptor': p_receptor,
                       'wild_type':mutant_id["wild_type"],
                       'mutant':mutant_id['mutant'],
                       'state':state,
                       'struct_id':struct_id}

        # Check if the calculated columns have been calculated for the pos, wt & mut grouping.
        # If so, all groups already have the column calculations needed added.
        # If not, all other grouping info must be calculated anyway to retrieve the site info for the wt & mut
        # grouping.
        wt_mut_group_id = ",".join([str(val) for key, val in mutant_id.items()
                                    if key in groupings['all']['include_in_id']])

        if wt_mut_group_id not in mutation_groups['all']:
            # In here: insert the code to find the site info
            calced_cols = get_calculated_columns(structural_rule_tree,
                                                 mutant_id['mutant'],
                                                 mutant_id['wild_type'],
                                                 generic_number,
                                                 p_class,
                                                 p_receptor,
                                                 conservation)

        # For each group, add the required info.
        for group_name, attr in groupings.items():
            # Create a dictionary of information pertaining to the whole group to which the mutant belongs
            #
            group_info = {key:item for key, item in mutant_id.items() if key not in attr['exclude_from_info']}
            # Create a group ID (which will be unique for each grouping)
            group_id = ",".join([str(val) for key, val in mutant_id.items()
                                 if key in attr['include_in_id']])

            # Get the context dict entry for which the mutant should be added.
            # If none, create one with the group_info
            group = mutation_groups[group_name].setdefault(group_id,
                                                           [group_info, {}]
                                                          )

            # If the group is newly created, calculate the values for the Frequency and Conservation Cols
            if group[1] == {}:
                # Get propensity and hydrophobicity values.
                group[0]['propensity'],\
                group[0]['hydro'],\
                group[0]["class_cons"],\
                group[0]["receptor_fam_cons"],\
                group[0]["ionic_lock"],\
                group[0]["sodium_ion"],\
                group[0]["res_switch"]\
                     = calced_cols[group_name]

                # Add further information to group_info allow for fast mutation subset filtering.
                if group_name == "all":
                    if mutant_id['mutant'] == 'A':
                        in_ala_subset = 'ala_subset'
                    elif mutant_id['wild_type'] == 'A' and mutant_id['mutant'] == 'L':
                        in_ala_subset = 'ala_subset'
                    else:
                        in_ala_subset = 'no_subset'

                    group[0]['ala_subset'] = in_ala_subset


            # Count the number of construct mutations recorded in the row.
            group[0]['GPCR_count'] += 1

            # Remove unnecessary items from the mutant info
            info = {key:set((item,)) for key, item in mutant_info.items() if key not in attr['include_in_id']}

            if group[1] == {}:
                # Initialise the dict with the first mutant.
                group[1].update(info)
            else:
                 # Add the specific mutant info.
                for key, item in info.items():
                    group[1][key].update(item)
                # Remove receptor family conservation info if row refers to >1 receptor family
                if len(group[1]['receptor']) != 1:
                    group[0]["receptor_fam_cons"] = u'\u2014'

    # Send the context dictionary to the template to be rendered
    return render(request, "construct/stabilisation_browser.html",
                  {'pos_and_mut': mutation_groups['pos_and_mut'],
                   'pos_and_wt': mutation_groups['pos_and_wt'],
                   'all': mutation_groups['all'],
                   'position_only': mutation_groups["position_only"]})

def conservation_table(prot_classes, gen_nums):
    '''Calculate the conservation values needed for the thermostabilisation view'''
    table = {}

    # Collect residue counts for all residues in the protein classes and at the generic number positions within the
    # prot_classes and gen_nums set, grouped by amino acid, generic number, protein receptor family, and protein class.
    residues = Residue.objects.order_by()\
        .only(
            "amino_acid",
            "generic_number__label",
            "protein_conformation__protein__species_id",
            "protein_conformation__protein__source_id",
            "protein_conformation__protein__family__parent__parent__parent__name")\
        .prefetch_related(
            "protein_conformation__protein__family__parent__parent__parent",
            "protein_conformation__protein__species",
            "protein_conformation__protein__source",
            "generic_number")\
        .filter(
            protein_conformation__protein__family__parent__parent__parent__name__in=list(prot_classes),
            protein_conformation__protein__species_id="1", protein_conformation__protein__source_id="1",
            generic_number__label__in=list(gen_nums))\
        .values(
            'amino_acid',
            'protein_conformation__protein__family__parent__parent__parent__name',
            "protein_conformation__protein__family__parent__name",
            "generic_number__label")\
        .annotate(Count('amino_acid'))

    # Restructure the data into table format, where each row contains the count for an amino acid at generic number
    # position, for either a given protein class or receptor family.
    for dic in residues:
        prot_row = table.setdefault(
            (dic['protein_conformation__protein__family__parent__parent__parent__name'], dic['generic_number__label']),
            {'total':0})
        prot_row['total'] += dic['amino_acid__count']
        prot_row.setdefault(dic['amino_acid'], 0)
        prot_row[dic['amino_acid']] += dic['amino_acid__count']

        rec_row = table.setdefault(
            (dic['protein_conformation__protein__family__parent__name'], dic['generic_number__label']), {'total':0})
        rec_row['total'] += dic['amino_acid__count']
        rec_row.setdefault(dic['amino_acid'], 0)
        rec_row[dic['amino_acid']] += dic['amino_acid__count']

    # Divide each row by it's total to get the frequency of each amino acid across the row (rather than it's count).
    for _, row in table.items():
        for amino_acid, count in row.items():
            if amino_acid != 'total':
                row[amino_acid] = round(count/row['total'], 2)

    return table

def get_calculated_columns(rule_tree, mutant, wild_type, g_n, prot_class, rec_fam, conservation): # pylint: disable=too-many-arguments

    ''' Calculate the propensity, hydrophobicity and site info for the given mut & wt for each grouping'''

    # Get the conservation values for the protein class and receptor family
    class_cons = conservation.get((prot_class, g_n), {})
    fam_cons = conservation.get((rec_fam, g_n), {})

    # Get the part of the structural_rule_tree relevant to the position and generic number (& hence to all groupings).
    related_rules = {
        'ionic_lock_tree':rule_tree["ionic_lock_tree"].get(prot_class[6], {}).get(g_n, {}),
        'sodium_ion_tree':rule_tree["sodium_ion_tree"].get(prot_class[6], {}).get(g_n, {}),
        'residue_switch_tree':rule_tree["residue_switch_tree"].get(prot_class[6], {}).get(g_n, {}),
    }

    # Return a dictionary consisting of the data and site column entries for each grouping / data analysis mode.
    return {
        'position_only': get_data_pos_grouping(related_rules),
        'pos_and_mut':get_data_mut_grouping(related_rules, mutant, class_cons, fam_cons),
        'pos_and_wt':get_data_wt_grouping(related_rules, wild_type, class_cons, fam_cons),
        'all': get_data_all_grouping(related_rules, mutant, wild_type, class_cons, fam_cons)
        }



def get_data_pos_grouping(rules):
    '''
    Calculate the Data and Site columns in the browser view for the position only analysis mode
    '''
    # Note: an empty dictionary evaluates to False in an if statement,
    ionic_lock = 'Pos Match' if rules['ionic_lock_tree'] else u'\u2014'
    sodium_ion = 'Pos Match' if rules['sodium_ion_tree'] else u'\u2014'
    residue_switch = 'Pos Match' if rules['residue_switch_tree'] else u'\u2014'

    # There is no mutant or wild type info, so all data cols are returned as u'\u2014'
    return (u'\u2014', u'\u2014', u'\u2014', u'\u2014', ionic_lock, sodium_ion, residue_switch)

def get_data_mut_grouping(rules, mutant, class_cons, fam_cons):
    '''
    Calculate the Data and Site columns in the browser view for the pos & mut analysis mode
    '''

    # Note: an empty dictionary evaluates to False in an if statement,
    # Check that rules exist that apply to the class, position and gn.
    if rules['ionic_lock_tree']:
        #  If so, check if there is a rule relevant to the mutant
        if rules['ionic_lock_tree'].get(mutant, {}):
            ionic_lock = 'Pos & Mutant AA Match'
        else:
            ionic_lock = 'Pos Match (But Not Mutant AA)'
    else:
        ionic_lock = u'\u2014'


    if rules['sodium_ion_tree']:
        if rules['sodium_ion_tree'].get(mutant, {}):
            sodium_ion = 'Pos & AA Mutant Match'
        else:
            sodium_ion = 'Pos Match (But Not Mutant AA)'
    else:
        sodium_ion = u'\u2014'

    if rules['residue_switch_tree']:
        if rules['residue_switch_tree'].get(mutant, {}):
            residue_switch = 'Pos & Mutant AA Match'
        else:
            residue_switch = 'Pos Match (But Not Mutant AA)'
    else:
        residue_switch = u'\u2014'


    return (AA_PROPENSITY.get(mutant, u'\u2014'),
            HYDROPHOBICITY.get(mutant, u'\u2014'),
            class_cons.get(mutant, u'\u2014'),
            fam_cons.get(mutant, u'\u2014'),
            ionic_lock, sodium_ion, residue_switch)


def get_data_wt_grouping(rules, wild_type, class_cons, fam_cons):
    '''
    Calculate the Data and Site columns in the browser view for the pos & wt analysis mode
    '''
    # # Note: an empty dictionary evaluates to False in an if statement,
    if rules['ionic_lock_tree']:
        # Note:  This is the simpliest, but not the most concise code.
        # However, okay as code is VERY rarely used.
        ionic_lock_set = set()
        for _, wt_rule_dict  in rules['ionic_lock_tree'].items():
            for key in wt_rule_dict:
                ionic_lock_set.add(key)
        if wild_type in ionic_lock_set:
            ionic_lock = 'Pos & Wild Type AA Match'
        else:
            ionic_lock = 'Pos Match (But Not Wild Type AA)'
    else:
        ionic_lock = u'\u2014'

    # Check that rules exist that apply to the class, position and gn.
    if rules['sodium_ion_tree']:
        sodium_ion_set = set()
         # If so, check if there is a rule relevant to the wild type.  As the dictionary tree is constructed so that
         # the mutant is in the 3rd level, and the wold type in the 4th.  Hence each mutant branch must be checked
         # for the wild type.
        for _, wt_rule_dict  in rules['sodium_ion_tree'].items():
            for key in wt_rule_dict:
                sodium_ion_set.add(key)
        if wild_type in sodium_ion_set:
            sodium_ion = 'Pos & Wild Type AA Match'
        else:
            sodium_ion = 'Pos Match (But Not Wild Type AA)'
    else:
        sodium_ion = u'\u2014'

    if rules['residue_switch_tree']:
        residue_switch_set = set()
        for _, wt_rule_dict  in rules['residue_switch_tree'].items():
            for key in wt_rule_dict:
                residue_switch_set.add(key)
        if wild_type in residue_switch_set:
            residue_switch = 'Pos & Wild Type AA Match'
        else:
            residue_switch = 'Pos Match (But Not Wild Type AA)'
    else:
        residue_switch = u'\u2014'

    return (AA_PROPENSITY.get(wild_type, u'\u2014'),
            HYDROPHOBICITY.get(wild_type, u'\u2014'),
            class_cons.get(wild_type, u'\u2014'),
            fam_cons.get(wild_type, u'\u2014'),
            ionic_lock, sodium_ion, residue_switch)

def get_data_all_grouping(rules, mutant, wild_type, class_cons, fam_cons):
    '''
    Calculate the Data and Site columns in the browser view for the pos, mut & wt analysis mode
    '''
    # Get propensity fold change where possible
    mut = AA_PROPENSITY.get(mutant, u'\u2014')
    w_t = AA_PROPENSITY.get(wild_type, u'\u2014')
        #  Where possible, calculate the fold change
    prop = u'\u2014' if isinstance(mut, str) | isinstance(w_t, str) else str(round(mut-w_t, 2))
        # Append the mut and wt values to the end.
    prop = prop + ' (' + str(mut) + u'\u2212'+ str(w_t) +')'

    # Get hydrophobicity fold change where possible
    mut = HYDROPHOBICITY.get(mutant, u'\u2014')
    w_t = HYDROPHOBICITY.get(wild_type, u'\u2014')
    hydro = u'\u2014' if isinstance(mut, str) | isinstance(w_t, str) else str(round(mut-w_t, 2))
    hydro = hydro + ' (' + str(mut) + u'\u2212'+ str(w_t) +')'

    # Get the receptor family conservation fold change where possible
    mut = fam_cons.get(mutant, u'\u2014')
    w_t = fam_cons.get(wild_type, u'\u2014')
    rec_cons = u'\u2014' if isinstance(mut, str) | isinstance(w_t, str) else str(round(mut-w_t, 2))
    rec_cons += ' ('+str(mut)+u'\u2212'+str(w_t)+')'

    # Get the protein class conservation fold change where possible
    mut = class_cons.get(mutant, u'\u2014')
    w_t = class_cons.get(wild_type, u'\u2014')
    prot_cons = u'\u2014' if isinstance(mut, str) | isinstance(w_t, str) else str(round(mut-w_t, 2))
    prot_cons += ' ('+str(mut)+u'\u2212'+str(w_t)+')'

    # Get site info from the structural site rules
    ionic_lock = rules['ionic_lock_tree'].get(mutant, {}).get(wild_type, u'\u2014')
    sodium_ion = rules['sodium_ion_tree'].get(mutant, {}).get(wild_type, u'\u2014')
    residue_switch = rules['residue_switch_tree'].get(mutant, {}).get(wild_type, u'\u2014')

    return (prop,
            hydro,
            prot_cons,
            rec_cons,
            ionic_lock,
            sodium_ion,
            residue_switch)



def parse_rule_definition(rule_def):
    '''
        Take in a rule definition from the structural rules, and parse so that's it's suitable both for display and
        use in the rule dictionaries.

        Args:
        -  rule_def should be of the form:
                        Ionic / Sodium / Residue + ... ... + removal / contraction / addition

        Returns:
          site - meaning type of site the definiton refers to.  to be 'ionic_lock', 'sodium_ion', or 'residue_switch'
          definiton - the action at the site.  to be 'Removed', 'Contracted', or 'Added'
    '''
    # Get the type of action in the definition
    if rule_def[-7:] == 'removal':
        result = 'Removed'
    elif rule_def[-11:] == 'contraction':
        result = 'Contracted'
    else:
        result = 'Added'

    # Get action placement from the definition
    if rule_def[:5] == 'Ionic':
        site = 'ionic_lock'
    elif rule_def[:6] == 'Sodium':
        site = 'sodium_ion'
    elif rule_def[:7] == 'Residue':
        site = 'residue_switch'
    else: # Then there is no sensible way to understand this rule.
        site = 'other'
        result = rule_def # Override previous rule finding.

    return (site, result)


def create_structural_rule_trees(rule_dictionary):
    '''
     Restructure the structural rules from a list of dictionaries to a tree-like nested dictionary,
     so that they may be easily and quickly searched.

     I.e. each type of site gets its own tree/dictionary, as it has it's own column,
     This allows for simplier code when querying the rules.
    '''
    structural_rule_trees = {'ionic_lock_tree':{}, 'sodium_ion_tree':{}, 'residue_switch_tree':{}, 'other_tree':{}}

    # List of classes included by the 'All' class designation.
    classes = {'A', 'B', 'C', 'F'}
    # List of amino acids included by the 'X' amino acid designation.
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
                   'T', 'V', 'W', 'Y', 'B', 'Z', 'J']

    # For each tree, initiate the inner class dictionary, for each class.
    for _, tree in structural_rule_trees.items():
        for prot_class in classes:
            tree.setdefault(prot_class, {})

    # For each class type in the Structural rules list, iterate through the contained dictionaries.
    for item in {'A', 'B', 'C', 'All'}:
        for rule in rule_dictionary[item]:
            # Get the dictionary to which the rule pertains
            site, definition = parse_rule_definition(rule['Definition'])
            tree = structural_rule_trees[site+"_tree"]

            # Get a set of the classes and wild type aas that the rule affects
            rule_class = classes if rule['Class'] == 'All' else {rule['Class']}
            rule_wt = amino_acids if rule['Wt AA'] == 'X' else rule['Wt AA'].split('/')

            # Iterate through the keys in each rule, adding a 'branch' to the nested dictionary, as needed.
            for prot_class in rule_class:
                node = tree.setdefault(prot_class, {})\
                                           .setdefault(rule['Generic Position'], {})\
                                           .setdefault(rule['Mut AA'], {})
                for acid in rule_wt:
                    # If the rule definition is already stored, append the next definition to it.
                    # Otherwise, create a new entry, consisting of the rule definiton.
                    acid_node = node.get(acid, "")
                    if acid_node == "":
                        # Then no previous rules.
                        node[acid] = definition
                    else: # Add to the previous results
                        node[acid] = acid_node + ", " + definition

    return structural_rule_trees



def fetch_all_pdb(request):

    structures = Structure.objects.all()

    for s in structures:
        pdbname = str(s)
        print(pdbname)
        failed = []
        try:
            protein = Protein.objects.filter(entry_name=pdbname.lower()).get()
            d = fetch_pdb_info(pdbname,protein)

            #delete before adding new
            Construct.objects.filter(name=d['construct_crystal']['pdb_name']).delete()
            add_construct(d)
        except:
            print(pdbname,'failed')
            failed.append(pdbname)


    context = {'failed':failed}

    return render(request,'pdb_all.html',context)

def fetch_pdb(request, slug):

    protein = Protein.objects.filter(entry_name=slug.lower()).get()

    d = fetch_pdb_info(slug,protein)

    #delete before adding new
    print(d['construct_crystal']['pdb_name'])
    #Construct.objects.filter(name__iexact=d['construct_crystal']['pdb_name']).delete()
    #add_construct(d)

    #cache.delete(d['construct_crystal']['pdb_name']+'_schematics')
    #cache.delete(d['construct_crystal']['pdb_name']+'_snake')
    context = {'d':d}

    return render(request,'pdb_fetch.html',context)

def fetch_pdb_for_webform(request, slug, **response_kwargs):

    slug = slug.lower()
    protein = Protein.objects.filter(entry_name=slug).get()

    d = fetch_pdb_info(slug,protein)
    d = convert_ordered_to_disordered_annotation(d)

    jsondata = json.dumps(d)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)

class ConstructBrowser(TemplateView):
    """
    Fetching construct data for browser
    """

    template_name = "construct_browser.html"

    def get_context_data (self, **kwargs):

        context = super(ConstructBrowser, self).get_context_data(**kwargs)
        try:
            cons = Construct.objects.all().prefetch_related(
                "crystal","mutations","purification","protein__family__parent__parent__parent", "insertions__insert_type", "modifications", "deletions", "crystallization__chemical_lists",
                "protein__species","structure__pdb_code","structure__publication__web_link", "contributor")

            context['constructs'] = []
            for c in cons:
                #c.schematics = c.schematic()
                c.wt_schematic = c.wt_schematic()
                c.cons_schematic = c.cons_schematic()
                context['constructs'].append(c)

        except Construct.DoesNotExist as e:
            pass

        return context

class ExperimentBrowser(TemplateView):
    """
    Fetching construct data for browser
    """

    template_name = "experimental_browser.html"

    def get_context_data (self, **kwargs):

        context = super(ExperimentBrowser , self).get_context_data(**kwargs)
        try:
            cons = Construct.objects.all().prefetch_related(
                "crystal","mutations","purification","protein__family__parent__parent__parent", "insertions__insert_type","expression","solubilization", "modifications", "deletions",
                "crystallization__crystal_method", "crystallization__crystal_type",
                "crystallization__chemical_lists", "crystallization__chemical_lists__chemicals__chemical__chemical_type",
                "protein__species","structure__pdb_code","structure__publication__web_link", "contributor").annotate(pur_count = Count('purification__steps')).annotate(sub_count = Count('solubilization__chemical_list__chemicals'))

            #context['constructs'] = cache.get('construct_browser')
            #if context['constructs']==None:
            context['constructs'] = []
            context['schematics'] = []
            for c in cons:
                # c.schematic_cache = c.schematic()
                c.summary = c.chem_summary()
                context['constructs'].append(c)

            #cache.set('construct_browser', context['constructs'], 60*60*24*2) #two days
            # else:
            #     print('construct_browser used cache')

        except Construct.DoesNotExist as e:
            pass

        return context

class design(AbsTargetSelection):

    # Left panel
    step = 1
    number_of_steps = 1
    # docs = 'generic_numbering.html'  # FIXME

    # description = 'Select receptors to index by searching or browsing in the middle column. You can select entire' \
    #     + ' receptor families and/or individual receptors.\n\nSelected receptors will appear in the right column,' \
    #     + ' where you can edit the list.\n\nSelect which numbering schemes to use in the middle column.\n\nOnce you' \
    #     + ' have selected all your receptors, click the green button.'

    description = 'Get construct suggestions based on published constructs.'

    # Middle section
    numbering_schemes = False
    filters = False
    search = True
    title = "Select a receptor"

    template_name = 'designselection.html'
    type_of_selection = 'targets'
    selection_only_receptors = True

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])

    # Buttons
    buttons = {
        'continue': {
            'label': 'Show results',
            'onclick': 'submitupload()',
            'color': 'success',
            #'url': 'calculate/'
        }
    }

    redirect_on_select = False

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        return context

@csrf_exempt #jquery send post, so no csrf
def align(request):

    ids = json.loads(request.POST.get('ids'))

    c_ids = []
    s_ids = []
    for i in ids:
        if i.startswith('align'):
            s_ids.append(i.split('_')[1])
        else:
            c_ids.append(i)

    cons = Construct.objects.filter(pk__in=c_ids).prefetch_related(
                "crystal","mutations","purification","protein__family__parent__parent__parent", "insertions", "modifications", "deletions", "crystallization__chemical_lists",
                "protein__species","structure__pdb_code","structure__publication__web_link", "contributor")
    proteins = []
    constructs = OrderedDict()
    annotations = {}
    for c in cons:
        # print(c)
        proteins.append(c.protein)
        constructs[c.name] = c.protein.entry_name
        annotations[c.name] = c.schematic()['annotations']

    print(annotations)

    if len(s_ids):
        rs = Residue.objects.filter(protein_conformation__protein__in=proteins, protein_segment__slug__in=s_ids).prefetch_related(
        'protein_conformation__protein', 'protein_conformation__state', 'protein_segment',
        'generic_number__scheme', 'display_generic_number__scheme')
    else:
        s_ids = ['N-term','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','ICL4','H8','C-term']
        rs = Residue.objects.filter(protein_conformation__protein__in=proteins).prefetch_related(
        'protein_conformation__protein', 'protein_conformation__state', 'protein_segment',
        'generic_number__scheme', 'display_generic_number__scheme')

    print("residues",len(rs))

    distinct_gn = []
    ordered_gn = OrderedDict()
    distinct_segments = []
    overview = OrderedDict()
    segment_length = OrderedDict()
    for s in s_ids:
        overview[s] = OrderedDict()
        segment_length[s] = {'aligned':0, 'before':0,'after':0,'total':0}

    protein_lookup = {}
    print('build stuff')

    segment = ''
    protein = ''

    track_unaligned = {}

    #Find all possible generic numbers, to ensure gaps
    for r in rs.order_by('protein_conformation__id','sequence_number'):
        if segment!=r.protein_segment.slug or protein!=r.protein_conformation.protein.entry_name:
            no_encountered_gn = True
            length = 0
            length_before = 0
            length_after = 0

        segment = r.protein_segment.slug
        protein = r.protein_conformation.protein.entry_name

        if protein not in protein_lookup:
            protein_lookup[protein] = {}
            track_unaligned[protein] = {}

        if segment not in track_unaligned[protein]:
            track_unaligned[protein][segment] = {'before':[],'after':[]}

        if segment not in distinct_segments:
            distinct_segments.append(segment)
            overview[segment] = OrderedDict()

        if r.generic_number:
            no_encountered_gn = False
            gn = r.generic_number.label
            protein_lookup[protein][gn] = {'aa':r.amino_acid,'pos':r.sequence_number}
            gn_sort = gn.split('x')[1]
            gn_sort = float("0."+gn_sort)
            if gn not in distinct_gn:
                distinct_gn.append(gn)
                overview[segment][gn_sort] = [gn,{'aa':'-','pos':''}]
            length += 1
        else:
            if no_encountered_gn:
                track_unaligned[protein][segment]['before'].append({'aa':r.amino_acid,'pos':r.sequence_number})
                length_before += 1
            else:
                track_unaligned[protein][segment]['after'].append({'aa':r.amino_acid,'pos':r.sequence_number})
                length_after += 1

        if len(overview[segment])>segment_length[segment]['aligned']:
            segment_length[segment]['aligned'] = len(overview[segment])
        if length_before>segment_length[segment]['before']:
            segment_length[segment]['before'] = length_before
        if length_after>segment_length[segment]['after']:
            segment_length[segment]['after'] = length_after
        if segment_length[segment]['aligned']+segment_length[segment]['before']+segment_length[segment]['after']>segment_length[segment]['total']:
            segment_length[segment]['total'] = segment_length[segment]['aligned']+segment_length[segment]['before']+segment_length[segment]['after']

    # SORT generic residues to ensure correct order
    gn_list = ""
    ordered_summary = OrderedDict()
    for seg,gns in overview.items():
        ordered_summary[seg] = OrderedDict()
        #GN LIST
        gn_list += """<td class="ali-td ali-residue res-color-X">&nbsp;</td>"""
        if seg!='C-term':
            for _ in range(segment_length[seg]['before']):
                gn_list += """<td class="ali-td">&nbsp;</td>"""
        for gn in sorted(gns):
            ordered_summary[seg][gns[gn][0]] = {'aa':'-','pos':''}
            gn_list += """<td class="ali-td-generic-num">{}</td>""".format("x"+gns[gn][0].split("x")[1])

        if seg=='C-term':
            for _ in range(segment_length[seg]['before']):
                gn_list += """<td class="ali-td">&nbsp;</td>"""

        for _ in range(segment_length[seg]['after']):
            gn_list += """<td class="ali-td">&nbsp;</td>"""

    alignment = OrderedDict()
    alignment_print_sequence = ""
    for c,p in constructs.items():
        alignment[c] = copy.deepcopy(ordered_summary)
        alignment_print_sequence += '<tr>'
        for seg,gns in alignment[c].items():

            if p not in track_unaligned:
                track_unaligned[p] = {seg: {'before':[],'after':[]}}

            if p not in protein_lookup:
                protein_lookup[p] = {}

            if seg not in track_unaligned[p]:
                track_unaligned[p][seg] = {'before':[],'after':[]}

            alignment_print_sequence += """<td class="ali-td ali-residue res-color-_">&nbsp;</td>"""

            if seg!='C-term':
                for _ in range(segment_length[seg]['before']-len(track_unaligned[p][seg]['before'])):
                    alignment_print_sequence += """<td class="ali-td ali-residue res-color-">
                                                    -</td>"""

            for aa in track_unaligned[p][seg]['before']:
                if aa['pos'] in annotations[c]:
                    annotation = annotations[c][aa['pos']][0]
                    annotation_text = "<br>"+annotations[c][aa['pos']][1]
                else:
                    annotation = ''
                    annotation_text = ''

                alignment_print_sequence += """<td class="ali-td ali-residue res-color-{}">
                                                <div data-toggle="tooltip" data-placement="top" data-html="true"
                                                title="{}{}{}">{}</div></td>""".format(annotation,aa['aa'],aa['pos'],annotation_text,aa['aa'])
            for gn, aa in gns.items():
                if gn in protein_lookup[p]:
                    aa = protein_lookup[p][gn]
                    alignment[c][seg][gn] = aa

                if aa['pos'] in annotations[c]:
                    annotation = annotations[c][aa['pos']][0]
                    annotation_text = "<br>"+annotations[c][aa['pos']][1]
                else:
                    annotation = ''
                    annotation_text = ''
                alignment_print_sequence += """<td class="ali-td ali-residue res-color-{}">
                                                <div data-toggle="tooltip" data-placement="top" data-html="true"
                                                title="{}{}<br>SCHEME: {}{}">{}</div></td>""".format(annotation,aa['aa'],aa['pos'],gn,annotation_text,aa['aa'])

            for aa in track_unaligned[p][seg]['after']:
                if aa['pos'] in annotations[c]:
                    annotation = annotations[c][aa['pos']][0]
                    annotation_text = "<br>"+annotations[c][aa['pos']][1]
                else:
                    annotation = ''
                    annotation_text = ''
                alignment_print_sequence += """<td class="ali-td ali-residue res-color-{}">
                                                <div data-toggle="tooltip" data-placement="top" data-html="true"
                                                title="{}{}{}">{}</div></td>""".format(annotation,aa['aa'],aa['pos'],annotation_text,aa['aa'])

            for _ in range(segment_length[seg]['after']-len(track_unaligned[p][seg]['after'])):
                alignment_print_sequence += """<td class="ali-td ali-residue res-color-">
                                                -</td>"""
            if seg=='C-term':
                for _ in range(segment_length[seg]['before']-len(track_unaligned[p][seg]['before'])):
                    alignment_print_sequence += """<td class="ali-td ali-residue res-color-">
                                                    -</td>"""

        alignment_print_sequence += '</tr>'

    print('done',len(alignment_print_sequence))
    context = {'constructs': constructs,'alignment_print_sequence': alignment_print_sequence, 'segment_length' : segment_length, 'gn_list' : gn_list, 'segments': s_ids, 'c_ids': json.dumps(c_ids)} #, 'alignment_print_sequence': alignment_print_sequence

    return render(request,'align.html',context)
