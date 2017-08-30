from django.shortcuts import render
from django.http import HttpResponse
from django.db.models import Min, Count, Max
from django.conf import settings
from django.views.decorators.cache import cache_page
from django import forms

from construct.models import *
from structure.models import Structure
from protein.models import ProteinConformation, Protein, ProteinSegment, ProteinFamily
from alignment.models import AlignmentConsensus
from common.definitions import AMINO_ACIDS, AMINO_ACID_GROUPS, STRUCTURAL_RULES

import json
from collections import OrderedDict
import re
import xlrd
import yaml
import os
import time
import pickle

Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')

class FileUploadForm(forms.Form):
    file_source = forms.FileField()

def parse_excel(path):
    workbook = xlrd.open_workbook(path)
    worksheets = workbook.sheet_names()
    d = {}
    for worksheet_name in worksheets:
        if worksheet_name in d:
            print('Error, worksheet with this name already loaded')
            continue

        #d[worksheet_name] = OrderedDict()
        worksheet = workbook.sheet_by_name(worksheet_name)

        num_rows = worksheet.nrows - 1
        num_cells = worksheet.ncols
        curr_row = 0 #skip first, otherwise -1

        headers = []
        for i in range(num_cells):
            h = worksheet.cell_value(0, i)
            if h=="":
                h = "i_"+str(i)
            if h in headers:
                h += "_"+str(i)
            headers.append(worksheet.cell_value(0, i))
        for curr_row in range(1,num_rows+1):
            row = worksheet.row(curr_row)
            key = worksheet.cell_value(curr_row, 0)

            if key=='':
                continue
            if key not in d:
                d[key] = []
            temprow = OrderedDict()
            for curr_cell in range(num_cells):
                cell_value = worksheet.cell_value(curr_row, curr_cell)
                if headers[curr_cell] not in temprow:
                    temprow[headers[curr_cell]] = cell_value
            d[key].append(temprow)
    return d

def compare_family_slug(a,b):
    a = a.split("_")
    b = b.split("_")

    if a[0]!=b[0]:
        return 0,"Different Class"
    elif a[1]!=b[1]:
        return 1,"Class"
    elif a[2]!=b[2]:
        return 2,"Ligand Type"
    elif a[3]!=b[3]:
        return 3,"Receptor Family"
    else:
        return 4,"Receptor"

def new_tool(request):

    simple_selection = request.session.get('selection', False)
    proteins = []
    for target in simple_selection.targets:
        if target.type == 'protein':
            proteins.append(target.item)
    context = {}

    context['target'] = proteins[0]

    level = proteins[0].family.slug
    if level.split("_")[0]=='001':
        c_level = 'A'
    elif level.split("_")[0]=='002':
        c_level = 'B'
    elif level.split("_")[0]=='003':
        c_level = 'B'
    elif level.split("_")[0]=='004':
        c_level = 'C'
    elif level.split("_")[0]=='005':
        c_level = 'F'
    else:
        c_level = ''

    states = list(Structure.objects.filter(protein_conformation__protein__family__slug__startswith=level.split("_")[0]).all().values_list('state__slug', flat = True).distinct())
    if 'active' in states:
        active_xtals = True
    else:
        active_xtals = False


    rs = Residue.objects.filter(protein_conformation__protein=proteins[0]).prefetch_related('protein_segment','display_generic_number','generic_number')

    residues = {}
    residues_gn = {}
    residues_pos = {}
    for r in rs:
        segment = r.protein_segment.slug
        segment = segment.replace("-","")
        if segment not in residues:
            residues[segment] = []
        residues[segment].append(r)
        label = ''
        if r.generic_number:
            residues_gn[r.generic_number.label] = r
            label = r.display_generic_number.label

        residues_pos[r.sequence_number] = [r.amino_acid,r.protein_segment.slug,label]


    cons = Construct.objects.all().prefetch_related('crystal', 'protein__family','deletions','structure__state','insertions__insert_type')
    
    inserts = {}
    inserts['fusions'] = []
    inserts['other'] = {}
    for ins in ConstructInsertionType.objects.all().order_by('name','subtype'):
       # print(ins.name,ins.subtype,ins.sequence)
        if ins.name == 'fusion':
            inserts['fusions'].append(ins.subtype)
        else:
            if ins.name not in inserts['other']:
                inserts['other'][ins.name] = []
            if ins.subtype not in inserts['other'][ins.name]:
                inserts['other'][ins.name].append(ins.subtype)
        # fusion, f_results = c.fusion()
        # if fusion:
        #     f_protein = f_results[0][2]
        #     if f_protein not in inserts['fusions']:
        #         inserts['fusions'].append(f_protein)
        # else:
        #     for ins in c.insertions.all():
        #         print(ins)
    context['ICL_max'] = {'ICL2': residues['ICL2'][-1].sequence_number, 'ICL3': residues['ICL3'][-1].sequence_number}
    context['ICL_min'] = {'ICL2': residues['ICL2'][0].sequence_number,'ICL3': residues['ICL3'][0].sequence_number}
    context['residues'] = residues
    context['residues_gn'] = residues_gn
    context['residues_pos'] = residues_pos
    context['class'] = c_level
    context['active_xtals'] = active_xtals
    context['inserts'] = inserts
    context['form'] = FileUploadForm
    #print(residues)

    return render(request,'new_tool.html',context)

def tool(request):

    simple_selection = request.session.get('selection', False)
    proteins = []
    for target in simple_selection.targets:
        if target.type == 'protein':
            proteins.append(target.item)
    print(proteins)
    context = {}

    context['target'] = proteins[0]

    level = proteins[0].family.slug
    if level.split("_")[0]=='001':
        c_level = 'A'
    elif level.split("_")[0]=='002':
        c_level = 'B'
    elif level.split("_")[0]=='003':
        c_level = 'B'
    elif level.split("_")[0]=='004':
        c_level = 'C'
    elif level.split("_")[0]=='005':
        c_level = 'F'
    else:
        c_level = ''

    states = list(Structure.objects.filter(protein_conformation__protein__family__slug__startswith=level.split("_")[0]).all().values_list('state__slug', flat = True).distinct())
    if 'active' in states:
        active_xtals = True
    else:
        active_xtals = False


    rs = Residue.objects.filter(protein_conformation__protein=proteins[0]).prefetch_related('protein_segment','display_generic_number','generic_number')

    residues = {}
    residues_gn = {}
    residues_pos = {}
    for r in rs:
        segment = r.protein_segment.slug
        segment = segment.replace("-","")
        if segment not in residues:
            residues[segment] = []
        residues[segment].append(r)
        label = ''
        if r.generic_number:
            residues_gn[r.generic_number.label] = r
            label = r.display_generic_number.label

        residues_pos[r.sequence_number] = [r.amino_acid,r.protein_segment.slug,label]


    cons = Construct.objects.all().prefetch_related('crystal', 'protein__family','deletions','structure__state','insertions__insert_type')
    
    inserts = {}
    inserts['fusions'] = []
    inserts['other'] = {}
    for ins in ConstructInsertionType.objects.all().order_by('name','subtype'):
       # print(ins.name,ins.subtype,ins.sequence)
        if ins.name == 'fusion':
            inserts['fusions'].append(ins.subtype)
        else:
            if ins.name not in inserts['other']:
                inserts['other'][ins.name] = []
            if ins.subtype not in inserts['other'][ins.name]:
                inserts['other'][ins.name].append(ins.subtype)
        # fusion, f_results = c.fusion()
        # if fusion:
        #     f_protein = f_results[0][2]
        #     if f_protein not in inserts['fusions']:
        #         inserts['fusions'].append(f_protein)
        # else:
        #     for ins in c.insertions.all():
        #         print(ins)
    print(inserts)
    context['residues'] = residues
    context['residues_gn'] = residues_gn
    context['residues_pos'] = residues_pos
    context['class'] = c_level
    context['active_xtals'] = active_xtals
    context['inserts'] = inserts
    context['form'] = FileUploadForm
    #print(residues)

    return render(request,'tool.html',context)

@cache_page(60 * 60 * 24)
def json_fusion(request, slug, **response_kwargs):

    level = Protein.objects.filter(entry_name=slug).values_list('family__slug', flat = True).get()
    #proteins = Construct.objects.all().values_list('protein', flat = True)
    cons = Construct.objects.all().prefetch_related('crystal', 'protein__family','deletions')

    jsondata = "glyco"
    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)

@cache_page(60 * 60 * 24)
def json_palmi(request, slug, **response_kwargs):

    start_time = time.time()
    seq = Protein.objects.filter(entry_name=slug).values_list('sequence', flat = True).get()
    rs = Residue.objects.filter(protein_conformation__protein__entry_name=slug,protein_segment__slug__in=['H8','C-term']).order_by('sequence_number').prefetch_related('protein_segment')
    residues = {}
    seq = ''
    end_h8 = 0
    start_h8 = 0
    for r in rs:
        if not start_h8 and r.protein_segment.slug == 'H8':
            start_h8 = r.sequence_number
        if not end_h8 and r.protein_segment.slug == 'C-term':
            end_h8 = r.sequence_number-1 #end_h8 was prev residue
        elif end_h8 and r.sequence_number-10>end_h8:
            continue
        seq += r.amino_acid
        residues[r.sequence_number] = r.protein_segment.slug

    #No proline!
    p = re.compile("C")
    #print('all')
    mutations_all = []
    for m in p.finditer(seq):
        mutations_all.append([m.start()+start_h8,"A",'','',m.group(),residues[m.start()+start_h8]])


    palmi = OrderedDict()
    palmi['']= mutations_all

    jsondata = palmi
    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    end_time = time.time()
    diff = round(end_time - start_time,1)
    print("palmi",diff)
    return HttpResponse(jsondata, **response_kwargs)


@cache_page(60 * 60 * 24)
def json_glyco(request, slug, **response_kwargs):
    start_time = time.time()

    seq = Protein.objects.filter(entry_name=slug).values_list('sequence', flat = True).get()
    rs = Residue.objects.filter(protein_conformation__protein__entry_name=slug).prefetch_related('protein_segment')
    residues = {}
    for r in rs:
        residues[r.sequence_number] = r.protein_segment.slug

    #No proline!
    p = re.compile("N[^P][TS]")
    #print('all')
    mutations_all = []
    matches = re.finditer(r'(?=([N][^P][TS]))',seq)
    matches_seq = re.findall(r'(?=([N][^P][TS]))',seq)
    #{"all": [[39, "Q", "", "", "NTS", "N-term"], [203, "Q", "", "", "NNT", "ECL2"]], "mammalian": [[205, "V", 206, "V", "TTCVLNDPN", "ECL2"]]}
    for i,m in enumerate(matches):
        #print(matches_seq[i],m.start())
        #print(m.start(), m.group())
        if residues[m.start()+1] in ['N-term','ECL1','ECL2','ECL3']:
            mutations_all.append([m.start()+1,"Q",'','',matches_seq[i],residues[m.start()+1]])

    #print('mamalian')
    #p = re.compile("[TS]{2}[A-Z]{1,11}[N]", overlapped=True)

    matches = re.finditer(r'(?=([TS]{2}[A-Z]{1,10}[N]))',seq)
    matches_seq = re.findall(r'(?=([TS]{2}[A-Z]{1,10}[N]))',seq)
    #matches = re.findall(r'(?=(\w\w))', seq)
    #print(matches)
    mutations_mammalian = []
    for i,m in enumerate(matches):
        #print(matches_seq[i],m.start())
        if matches_seq[i][0]=="T":
            pos0 = "V"
        if matches_seq[i][1]=="T":
            pos1 = "V"
        if matches_seq[i][0]=="S":
            pos0 = "A"
        if matches_seq[i][1]=="S":
            pos1 = "A"

        if residues[m.start()+1] in ['N-term','ECL1','ECL2','ECL3']:
            mutations_mammalian.append([m.start()+1,pos0,m.start()+2,pos1,matches_seq[i],residues[m.start()+1]])

    glyco = OrderedDict()
    glyco['n-linked']= mutations_all
    glyco['o-linked'] = mutations_mammalian

    jsondata = glyco
    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    end_time = time.time()
    diff = round(end_time - start_time,1)
    print("glyco",diff)
    return HttpResponse(jsondata, **response_kwargs)

@cache_page(60 * 60 * 24)
def json_icl3(request, slug, **response_kwargs):
    start_time = time.time()
    level = Protein.objects.filter(entry_name=slug).values_list('family__slug', flat = True).get()

    ##PREPARE TM1 LOOKUP DATA
    proteins = Construct.objects.all().values_list('protein', flat = True)
    tm5_start = {}
    tm5_end = {}
    tm6_start = {}
    tm6_end = {}
    tm5_50 = {}
    tm6_50 = {}
    pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).prefetch_related('protein').filter(residue__protein_segment__slug='TM5').annotate(start=Min('residue__sequence_number'), end=Max('residue__sequence_number'))
    for pc in pconfs:
        tm5_start[pc.protein.entry_name] = pc.start
        tm5_end[pc.protein.entry_name] = pc.end
    pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).prefetch_related('protein').filter(residue__protein_segment__slug='TM6').annotate(start=Min('residue__sequence_number'), end=Max('residue__sequence_number'))
    for pc in pconfs:
        tm6_start[pc.protein.entry_name] = pc.start
        tm6_end[pc.protein.entry_name] = pc.end


    pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).prefetch_related('protein').filter(residue__generic_number__label__in=['5x50','6x50']).annotate(start=Min('residue__sequence_number'), end=Max('residue__sequence_number'))
    for pc in pconfs:
        tm5_50[pc.protein.entry_name] = pc.start
        tm6_50[pc.protein.entry_name] = pc.end


    cons = Construct.objects.all().prefetch_related('crystal', 'protein__family','deletions','structure__state','insertions__insert_type')
    
    deletions = OrderedDict()
    deletions['Receptor'] = {}
    deletions['Receptor Family'] = {}
    deletions['Ligand Type'] = {}
    deletions['Class'] = {}
    deletions['Different Class'] = {}
    states = {}
    for c in cons:
        p = c.protein
        entry_name = p.entry_name
        p_level = p.family.slug
        d_level, d_level_name = compare_family_slug(level,p_level)
        pdb = c.crystal.pdb_code
        state = c.structure.state.slug
        if pdb not in states:
            states[pdb] = state
        fusion, f_results = c.fusion()
        if fusion:
            f_protein = f_results[0][2]
        else:
            f_protein = ""
        for deletion in c.deletions.all():
            #print(pdb,deletion.start,deletion.end)
            if deletion.start > tm5_start[entry_name] and deletion.start < tm6_end[entry_name]:
                if p.entry_name not in deletions[d_level_name]:
                    deletions[d_level_name][entry_name] = {}
                #deletions[entry_name][pdb] = [tm5_end[entry_name],tm6_start[entry_name],deletion.start,deletion.end,deletion.start-tm5_end[entry_name],tm6_start[entry_name]-deletion.end]
                deletions[d_level_name][entry_name][pdb] = [deletion.start-tm5_50[entry_name]-1,tm6_50[entry_name]-deletion.end-1,state,str(fusion),f_protein]
                # if (str(fusion)=='icl3'):
                #     print(entry_name,pdb,50+deletion.start-tm5_50[entry_name],50-(tm6_50[entry_name]-deletion.end-1),str(fusion),f_protein)

    # for pdb,state in sorted(states.items()):
    #     print(pdb,"\t",state)

    jsondata = deletions
    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    end_time = time.time()
    diff = round(end_time - start_time,1)
    print("icl3",diff)
    return HttpResponse(jsondata, **response_kwargs)

@cache_page(60 * 60 * 24)
def json_icl2(request, slug, **response_kwargs):
    start_time = time.time()
    level = Protein.objects.filter(entry_name=slug).values_list('family__slug', flat = True).get()

    ##PREPARE TM1 LOOKUP DATA
    proteins = Construct.objects.all().values_list('protein', flat = True)
    tm3_start = {}
    tm3_end = {}
    tm4_start = {}
    tm4_end = {}
    tm3_50 = {}
    tm4_50 = {}
    pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).prefetch_related('protein').filter(residue__protein_segment__slug='TM3').annotate(start=Min('residue__sequence_number'), end=Max('residue__sequence_number'))
    for pc in pconfs:
        tm3_start[pc.protein.entry_name] = pc.start
        tm3_end[pc.protein.entry_name] = pc.end
    pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).prefetch_related('protein').filter(residue__protein_segment__slug='TM4').annotate(start=Min('residue__sequence_number'), end=Max('residue__sequence_number'))
    for pc in pconfs:
        tm4_start[pc.protein.entry_name] = pc.start
        tm4_end[pc.protein.entry_name] = pc.end


    pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).prefetch_related('protein').filter(residue__generic_number__label__in=['3x50','4x50']).annotate(start=Min('residue__sequence_number'), end=Max('residue__sequence_number'))
    for pc in pconfs:
        tm3_50[pc.protein.entry_name] = pc.start
        tm4_50[pc.protein.entry_name] = pc.end


    cons = Construct.objects.all().prefetch_related('crystal', 'protein__family','deletions','structure__state','insertions__insert_type')
    
    deletions = OrderedDict()
    deletions['Receptor'] = {}
    deletions['Receptor Family'] = {}
    deletions['Ligand Type'] = {}
    deletions['Class'] = {}
    deletions['Different Class'] = {}
    states = {}
    for c in cons:
        p = c.protein
        entry_name = p.entry_name
        p_level = p.family.slug
        d_level, d_level_name = compare_family_slug(level,p_level)
        pdb = c.crystal.pdb_code
        state = c.structure.state.slug
        if pdb not in states:
            states[pdb] = state
        fusion, f_results = c.fusion()
        if fusion:
            f_protein = f_results[0][2]
        else:
            f_protein = ""
        for deletion in c.deletions.all():
            #print(pdb,deletion.start,deletion.end)
            if deletion.start > tm3_start[entry_name] and deletion.start < tm4_end[entry_name]:
                if p.entry_name not in deletions[d_level_name]:
                    deletions[d_level_name][entry_name] = {}
                #deletions[entry_name][pdb] = [tm5_end[entry_name],tm6_start[entry_name],deletion.start,deletion.end,deletion.start-tm5_end[entry_name],tm6_start[entry_name]-deletion.end]
                deletions[d_level_name][entry_name][pdb] = [deletion.start-tm3_50[entry_name]-1,tm4_50[entry_name]-deletion.end-1,state,str(fusion),f_protein]

    # for pdb,state in sorted(states.items()):
    #     print(pdb,"\t",state)

    jsondata = deletions
    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    end_time = time.time()
    diff = round(end_time - start_time,1)
    print("icl2",diff)
    return HttpResponse(jsondata, **response_kwargs)

@cache_page(60 * 60 * 24)
def json_nterm(request, slug, **response_kwargs):
    start_time = time.time()

    level = Protein.objects.filter(entry_name=slug).values_list('family__slug', flat = True).get()

    ##PREPARE TM1 LOOKUP DATA
    proteins = Construct.objects.all().values_list('protein', flat = True)
    pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).prefetch_related('protein').filter(residue__protein_segment__slug='TM1').annotate(start=Min('residue__sequence_number'))
    #pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).filter(residue__generic_number__label__in=['1x50']).values_list('protein__entry_name','residue__sequence_number','residue__generic_number__label')
    tm1_start = {}
    for pc in pconfs:
        tm1_start[pc.protein.entry_name] = pc.start

    cons = Construct.objects.all().prefetch_related('crystal', 'protein__family','deletions','structure__state','insertions__insert_type')
    deletions = OrderedDict()
    deletions['Receptor'] = {}
    deletions['Receptor Family'] = {}
    deletions['Ligand Type'] = {}
    deletions['Class'] = {}
    deletions['Different Class'] = {}
    for c in cons:
        p = c.protein
        entry_name = p.entry_name
        p_level = p.family.slug
        d_level, d_level_name = compare_family_slug(level,p_level)
        pdb = c.crystal.pdb_code
        state = c.structure.state.slug
        fusion, f_results = c.fusion()
        if fusion:
            f_protein = f_results[0][2]
        else:
            f_protein = ""
        for deletion in c.deletions.all():
            if deletion.start < tm1_start[entry_name]:
                if p.entry_name not in deletions[d_level_name]:
                    deletions[d_level_name][entry_name] = {}
                deletions[d_level_name][entry_name][pdb] = [deletion.start,deletion.end-1, tm1_start[entry_name]-deletion.end-1,state,str(fusion),f_protein]

    jsondata = deletions
    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    end_time = time.time()
    diff = round(end_time - start_time,1)
    print("nterm",diff)
    return HttpResponse(jsondata, **response_kwargs)

@cache_page(60 * 60 * 24)
def json_cterm(request, slug, **response_kwargs):

    start_time = time.time()
    level = Protein.objects.filter(entry_name=slug).values_list('family__slug', flat = True).get()

    ##PREPARE TM1 LOOKUP DATA
    proteins = Construct.objects.all().values_list('protein', flat = True)
    # pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).filter(residue__protein_segment__slug='C-term').annotate(start=Min('residue__sequence_number'))
    # cterm_start = {}
    # for pc in pconfs:
    #     cterm_start[pc.protein.entry_name] = pc.start
    # pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).filter(residue__protein_segment__slug='C-term').annotate(start=Min('residue__sequence_number')).values_list('protein__entry_name','start','residue__generic_number__label')
    #pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).filter(residue__generic_number__label__in=['8x50']).values_list('protein__entry_name','residue__sequence_number','residue__generic_number__label')
    pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).prefetch_related('protein').filter(residue__protein_segment__slug='C-term').annotate(start=Min('residue__sequence_number'))
    cterm_start = {}
    for pc in pconfs:
        cterm_start[pc.protein.entry_name] = pc.start


    cons = Construct.objects.all().prefetch_related('crystal', 'protein__family','deletions','structure__state','insertions__insert_type')
    
    deletions = OrderedDict()
    deletions['Receptor'] = {}
    deletions['Receptor Family'] = {}
    deletions['Ligand Type'] = {}
    deletions['Class'] = {}
    deletions['Different Class'] = {}
    for c in cons:
        p = c.protein
        entry_name = p.entry_name
        p_level = p.family.slug
        d_level, d_level_name = compare_family_slug(level,p_level)
        pdb = c.crystal.pdb_code
        state = c.structure.state.slug
        fusion, f_results = c.fusion()
        if fusion:
            f_protein = f_results[0][2]
        else:
            f_protein = ""
        for deletion in c.deletions.all():
            if deletion.start >= cterm_start[entry_name]:
                if p.entry_name not in deletions[d_level_name]:
                    deletions[d_level_name][entry_name] = {}
                deletions[d_level_name][entry_name][pdb] = [deletion.start,deletion.end, deletion.start-cterm_start[entry_name],state,str(fusion),f_protein]

    jsondata = deletions
    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    end_time = time.time()
    diff = round(end_time - start_time,1)
    print("cterm",diff)
    return HttpResponse(jsondata, **response_kwargs)

@cache_page(60 * 60 * 24)
def thermostabilising(request, slug, **response_kwargs):

    start_time = time.time()

    rs = Residue.objects.filter(protein_conformation__protein__entry_name=slug).prefetch_related('protein_segment','display_generic_number','generic_number')

    wt_lookup = {}
    wt_lookup_pos = {}
    for r in rs:
        if r.generic_number:
            gn = r.generic_number.label
            wt_lookup[gn] = [r.amino_acid, r.sequence_number]
        pos = r.sequence_number
        wt_lookup_pos[pos] = [r.amino_acid]

    level = Protein.objects.filter(entry_name=slug).values_list('family__slug', flat = True).get()
    if level.split("_")[0]=='001':
        c_level = 'A'
    elif level.split("_")[0]=='002':
        c_level = 'B'
    elif level.split("_")[0]=='003':
        c_level = 'B'
    elif level.split("_")[0]=='004':
        c_level = 'C'
    elif level.split("_")[0]=='005':
        c_level = 'F'
    else:
        c_level = ''

    path = os.sep.join([settings.DATA_DIR, 'structure_data', 'construct_data', 'termo.xlsx'])
    d = parse_excel(path)
    if c_level in d:
        termo = d[c_level]
    else:
        termo = []

    results = OrderedDict()
    results['1'] = {}
    results['2'] = {} #fixed mut
    results['3'] = {} #fixed wt

    for mut in termo:
        gn = mut['GN']
        mut_aa = mut['MUT']
        wt_aa = mut['WT']
        entry_name = mut['UniProt']
        pos = int(mut['POS'])
        pdb = mut['PDB']
        if mut['Effect'] != 'Thermostabilising':
            continue #only thermo!
        if gn is "":
            continue
        if (entry_name == slug) or (entry_name.split('_')[0] == slug.split('_')[0] and wt_aa == wt_lookup[gn][0]):
            if gn not in results['1']:
                results['1'][gn] = {}
            if mut_aa not in results['1'][gn]:
                results['1'][gn][mut_aa] = {'pdbs':[], 'hits':0, 'wt':wt_lookup[gn]}
            if mut['PDB'] not in results['1'][gn][mut_aa]['pdbs']:
                results['1'][gn][mut_aa]['pdbs'].append(pdb)
            results['1'][gn][mut_aa]['hits'] += 1

        if gn:
            if gn in wt_lookup:
                if gn not in results['2']:
                    results['2'][gn] = {}
                if mut_aa not in results['2'][gn]:
                    results['2'][gn][mut_aa] = {'pdbs':[], 'proteins':[], 'hits':0, 'wt':wt_lookup[gn]}
                if entry_name not in results['2'][gn][mut_aa]['proteins']:
                    results['2'][gn][mut_aa]['proteins'].append(entry_name)
                    results['2'][gn][mut_aa]['hits'] += 1
                if wt_lookup[gn][0] == wt_aa:
                    if gn not in results['3']:
                        results['3'][gn] = {}
                    if wt_aa not in results['3'][gn]:
                        results['3'][gn][wt_aa] = {'pdbs':[], 'proteins':[], 'hits':0, 'wt':wt_lookup[gn], 'muts':[]}
                    if entry_name not in results['3'][gn][wt_aa]['proteins']:
                        results['3'][gn][wt_aa]['proteins'].append(entry_name)
                        results['3'][gn][wt_aa]['hits'] += 1
                        if mut_aa not in results['3'][gn][wt_aa]['muts']:
                            results['3'][gn][wt_aa]['muts'].append(mut_aa)
      
    temp = {}
    for gn, vals1 in results['2'].items():
        for mut_aa, vals2 in vals1.items():
            if vals2['hits']>1:
                if gn not in temp:
                    temp[gn] = {}
                if mut_aa not in temp[gn]:
                    temp[gn][mut_aa] = vals2
                #results['2'][gn].pop(mut_aa, None)
    results['2'] = temp

    temp_single = {}                  
    temp = {}
    for gn, vals1 in results['3'].items():
        for mut_aa, vals2 in vals1.items():
            if vals2['hits']>1:
                if gn not in temp:
                    temp[gn] = {}
                if mut_aa not in temp[gn]:
                    temp[gn][mut_aa] = vals2
                #results['2'][gn].pop(mut_aa, None)
            elif vals2['hits']==1:
                if gn not in temp_single:
                    temp_single[gn] = {}
                if mut_aa not in temp_single[gn]:
                    temp_single[gn][mut_aa] = vals2
    results['3'] = temp
    results['4'] = temp_single


    jsondata = results
    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    end_time = time.time()
    diff = round(end_time - start_time,1)
    print("termo",diff)
    return HttpResponse(jsondata, **response_kwargs)


@cache_page(60 * 60 * 24)
def structure_rules(request, slug, **response_kwargs):
    start_time = time.time()

    rs = Residue.objects.filter(protein_conformation__protein__entry_name=slug).prefetch_related('protein_segment','display_generic_number','generic_number')

    wt_lookup = {}
    wt_lookup_pos = {}
    for r in rs:
        if r.generic_number:
            gn = r.generic_number.label
            wt_lookup[gn] = [r.amino_acid, r.sequence_number]
        pos = r.sequence_number
        wt_lookup_pos[pos] = [r.amino_acid]

    level = Protein.objects.filter(entry_name=slug).values_list('family__slug', flat = True).get()
    if level.split("_")[0]=='001':
        c_level = 'A'
    elif level.split("_")[0]=='002':
        c_level = 'B'
    elif level.split("_")[0]=='003':
        c_level = 'B'
    elif level.split("_")[0]=='004':
        c_level = 'C'
    elif level.split("_")[0]=='005':
        c_level = 'F'
    else:
        c_level = ''

    # path = os.sep.join([settings.DATA_DIR, 'structure_data', 'construct_data', 'structure_rules.xlsx'])
    # d = parse_excel(path)
    # d_clean = {}
    # regex = r"(\d+)x(\d+)"
    # for rule_class, values in d.items():
    #     d_clean[rule_class] = []
    #     for rule in values:
    #         if rule['Type']!="Structure-based":
    #             # Only use structure based ones in this function
    #             continue
    #         if re.search(regex, rule['Definition']):
    #             match = re.search(regex, rule['Definition'])
    #             gn = match.group(1) + "x" + match.group(2)
    #             print(rule['Definition'],gn)
    #         else:
    #             continue
    #         regex = r"(\d+)x(\d+)"
    #         if re.search(regex, rule['Definition']):
    #             match = re.search(regex, rule['Definition'])
    #             rule['Generic Position'] = match.group(1) + "x" + match.group(2)
    #         else:
    #             continue
    #         d_clean[rule_class].append(rule)

    # # print(d)
    # print(json.dumps(d_clean,sort_keys=True, indent=4))
    d = STRUCTURAL_RULES
    # print(d)
    if c_level in d:
        rules = d[c_level]
    else:
        rules = []

    results = OrderedDict()
    results['active'] = {}
    results['inactive'] = {} #fixed mut

    for rule in rules:
        # if rule['Type']!="Structure-based":
        #     # Only use structure based ones in this function
        #     continue
        # regex = r"(\d+)x(\d+)"
        # if re.search(regex, rule['Definition']):
        #     match = re.search(regex, rule['Definition'])
        #     gn = match.group(1) + "x" + match.group(2)
        #     print(rule['Definition'],gn)
        # else:
        #     continue
        gn = rule['Generic Position']
        mut_aa = rule['Mut AA']
        wt_aas = rule['Wt AA'].split("/")
        definition = rule['Design Principle']+" "+rule['Addition / Removal']
        state = rule['State'].lower()
        valid = False
        if gn in wt_lookup:
            for wt_aa in wt_aas:
                if wt_aa=='X' and wt_lookup[gn][0]!=mut_aa: #if universal but not mut aa
                    valid = True
                elif wt_lookup[gn][0]==wt_aa:
                    valid = True
            if valid:
                mut = {'wt':wt_lookup[gn][0], 'gn': gn, 'pos':wt_lookup[gn][1], 'mut':mut_aa, 'definition':definition}
                if state=='all':
                    if gn not in results['active']: 
                        results['active'][gn] = []
                    if gn not in results['inactive']: 
                        results['inactive'][gn] = []
                    results['active'][gn].append(mut)
                    results['inactive'][gn].append(mut)
                else:
                    if gn not in results[state]: 
                        results[state][gn] = []
                    results[state][gn].append(mut)


    #     entry_name = mut['UniProt']
    #     pos = int(mut['POS'])
    #     pdb = mut['PDB']
    #     if mut['Effect'] != 'Thermostabilising':
    #         continue #only thermo!
    #     if entry_name == slug:
    #         if gn not in results['1']:
    #             results['1'][gn] = {}
    #         if mut_aa not in results['1'][gn]:
    #             results['1'][gn][mut_aa] = {'pdbs':[], 'hits':0, 'wt':wt_lookup[gn]}
    #         if mut['PDB'] not in results['1'][gn][mut_aa]['pdbs']:
    #             results['1'][gn][mut_aa]['pdbs'].append(pdb)
    #         results['1'][gn][mut_aa]['hits'] += 1

    #     if gn:
    #         if gn in wt_lookup:
    #             if gn not in results['2']:
    #                 results['2'][gn] = {}
    #             if mut_aa not in results['2'][gn]:
    #                 results['2'][gn][mut_aa] = {'pdbs':[], 'proteins':[], 'hits':0, 'wt':wt_lookup[gn]}
    #             if entry_name not in results['2'][gn][mut_aa]['proteins']:
    #                 results['2'][gn][mut_aa]['proteins'].append(entry_name)
    #                 results['2'][gn][mut_aa]['hits'] += 1
    #             if wt_lookup[gn][0] == wt_aa:
    #                 if gn not in results['3']:
    #                     results['3'][gn] = {}
    #                 if wt_aa not in results['3'][gn]:
    #                     results['3'][gn][wt_aa] = {'pdbs':[], 'proteins':[], 'hits':0, 'wt':wt_lookup[gn], 'muts':[]}
    #                 if entry_name not in results['3'][gn][wt_aa]['proteins']:
    #                     results['3'][gn][wt_aa]['proteins'].append(entry_name)
    #                     results['3'][gn][wt_aa]['hits'] += 1
    #                     if mut_aa not in results['3'][gn][wt_aa]['muts']:
    #                         results['3'][gn][wt_aa]['muts'].append(mut_aa)

    # temp = {}
    # for gn, vals1 in results['2'].items():
    #     for mut_aa, vals2 in vals1.items():
    #         if vals2['hits']>1:
    #             if gn not in temp:
    #                 temp[gn] = {}
    #             if mut_aa not in temp[gn]:
    #                 temp[gn][mut_aa] = vals2
    #             #results['2'][gn].pop(mut_aa, None)
    # results['2'] = temp

    # temp = {}
    # for gn, vals1 in results['3'].items():
    #     for mut_aa, vals2 in vals1.items():
    #         if vals2['hits']>1:
    #             if gn not in temp:
    #                 temp[gn] = {}
    #             if mut_aa not in temp[gn]:
    #                 temp[gn][mut_aa] = vals2
    #             #results['2'][gn].pop(mut_aa, None)
    # results['3'] = temp


    jsondata = results
    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    end_time = time.time()
    diff = round(end_time - start_time,1)
    print("rules",diff)
    return HttpResponse(jsondata, **response_kwargs)


# @cache_page(60 * 60 * 24)
def mutations(request, slug, **response_kwargs):
    from django.db import connection
    start_time = time.time()

    protein = Protein.objects.get(entry_name=slug)
    protein_class_slug = protein.family.slug.split("_")[0]
    protein_rf_name = protein.family.parent.name
    protein_rf_slug = protein.family.parent.slug
    protein_rf_count = ProteinFamily.objects.filter(parent__slug=protein_rf_slug).count()

    # Grab thermostabilising mutations
    key = "CD_all_thermo_mutations_class_%s" % protein_class_slug
    mutations = cache.get(key)
    if not mutations:
        mutations = []
        mutations_thermo = ConstructMutation.objects.filter(effects__slug='thermostabilising', construct__protein__family__parent__parent__parent__slug=protein_class_slug).all()\
                    .prefetch_related(
                        # "construct__structure__state",
                        "residue__generic_number",
                        # "residue__protein_segment",
                        "construct__protein__family__parent__parent__parent",
                        "construct__crystal"
                        )
        for mutant in mutations_thermo:
            if not mutant.residue.generic_number:
                continue
            prot = mutant.construct.protein
            p_receptor = prot.family.parent.name
            real_receptor = prot.entry_name
            pdb = mutant.construct.crystal.pdb_code
            gn = mutant.residue.generic_number.label
            mutations.append(([mutant.sequence_number,mutant.wild_type_amino_acid,mutant.mutated_amino_acid],real_receptor,pdb, p_receptor,gn))

        cache.set(key,mutations,60*60*24)

    # Build current target residue GN mapping
    rs = Residue.objects.filter(protein_conformation__protein__entry_name=slug, generic_number__isnull=False).prefetch_related('generic_number', 'protein_segment')
    
    # Build a dictionary to know how far a residue is from segment end/start
    # Used for propensity removals
    start_end_segments = {}
    for r in rs:
        if r.protein_segment.slug not in start_end_segments:
            start_end_segments[r.protein_segment.slug] = {'start':r.sequence_number}
        start_end_segments[r.protein_segment.slug]['end'] = r.sequence_number

    wt_lookup = {}
    GP_residues_in_target = []
    for r in rs:
        gn = r.generic_number.label
        from_start = r.sequence_number-start_end_segments[r.protein_segment.slug]['start']
        from_end = start_end_segments[r.protein_segment.slug]['end'] - r.sequence_number
        wt_lookup[gn] = [r.amino_acid, r.sequence_number,r.protein_segment.slug]
        if r.amino_acid in ["G","P"] and from_start>=4 and from_end>=4:
            # build a list of potential GP removals (ignore those close to helix borders)
            GP_residues_in_target.append(gn)

    # Go through all mutations and find groupings (common)
    mutation_list = OrderedDict()
    for mutation in mutations:
        pos = mutation[0][0]
        mut_wt = mutation[0][1]
        mut_mut = mutation[0][2]
        entry_name = mutation[1]
        pdb = mutation[2]
        family = mutation[3]
        gn = mutation[4]

        # First do the ones with same WT
        full_mutation = "%s_%s_%s" % (gn,mut_wt,"X")
        if gn in wt_lookup and wt_lookup[gn][0]==mut_wt:
            # Only use those that have the same WT residue at GN
            if full_mutation not in mutation_list:
                mutation_list[full_mutation] = {'proteins':[], 'hits':0, 'mutation':[[],[]], 'wt':'', 'pdbs':[], 'protein_families': []}

            entry_name = mutation[1].split("_")[0]

            if entry_name not in mutation_list[full_mutation]['proteins']:
                mutation_list[full_mutation]['proteins'].append(entry_name)
                mutation_list[full_mutation]['hits'] += 1  
                mutation_list[full_mutation]['mutation'][0].append(mut_wt)
                mutation_list[full_mutation]['mutation'][1].append(mut_mut)
                if gn in wt_lookup:
                    mutation_list[full_mutation]['wt'] = wt_lookup[gn]
                if family not in mutation_list[full_mutation]['protein_families']:
                    mutation_list[full_mutation]['protein_families'].append(family)

            if pdb not in mutation_list[full_mutation]['pdbs']:
                mutation_list[full_mutation]['pdbs'].append(pdb)

        # Second, check those with same mutated AA
        full_mutation = "%s_%s_%s" % (gn,"X",mut_mut)
        if gn in wt_lookup and wt_lookup[gn][0]!=mut_mut:
            if full_mutation not in mutation_list:
                mutation_list[full_mutation] = {'proteins':[], 'hits':0, 'mutation':[[],[]], 'wt':'', 'pdbs':[], 'protein_families': []}

            entry_name = mutation[1].split("_")[0]

            if entry_name not in mutation_list[full_mutation]['proteins']:
                mutation_list[full_mutation]['proteins'].append(entry_name)
                mutation_list[full_mutation]['hits'] += 1  
                mutation_list[full_mutation]['mutation'][0].append(mut_wt)
                mutation_list[full_mutation]['mutation'][1].append(mut_mut)
                if gn in wt_lookup:
                    mutation_list[full_mutation]['wt'] = wt_lookup[gn]
                if family not in mutation_list[full_mutation]['protein_families']:
                    mutation_list[full_mutation]['protein_families'].append(family)

            if pdb not in mutation_list[full_mutation]['pdbs']:
                mutation_list[full_mutation]['pdbs'].append(pdb)

    # Go through the previous list and filter with rules and add rule matches
    simple_list = OrderedDict()
    mutation_list = OrderedDict(sorted(mutation_list.items(), key=lambda x: x[1]['hits'],reverse=True))
    for gn, vals in mutation_list.items():
        definition_matches = []
        if gn.split("_")[1] == "X":
            # Below rules only apply the mutations that share the same mutation AA
            if slug.split("_")[0] in vals['proteins']:
                # Check if same receptor
                definition_matches.append([1,'same_receptor'])
            elif protein_rf_name in vals['protein_families']:
                # Check if same receptor receptor family
                definition_matches.append([2,'same_receptor_family'])
            elif len(vals['protein_families'])<2:
                # If not same receptor or receptor family and not in two receptor families,
                # it is just a single match on position used in B-F class
                if protein_class_slug!='001':
                    # If class A require two distinct receptor families
                    definition_matches.append([4,'same_pos'])

            if len(vals['protein_families'])>=2:
                # If mutation is seen in >=2 receptor families
                # Put this one outside the above logic, to allow multi definitions
                definition_matches.append([2,'common_mutation'])

            # Check for membrane binding
            if 'K' in vals['mutation'][1] or 'R' in vals['mutation'][1]:
                if vals['wt'][0] not in ['R','K']:
                    # Only if not R,K already
                    definition_matches.append([2,'membrane_binding'])
                elif vals['wt'][0] in ['K'] and 'R' in vals['mutation'][1]:
                    # If K
                    definition_matches.append([3,'membrane_binding_weak'])

        else:
            # Below rules is for the common WT (But different mut AA)
            if len(vals['protein_families'])>=2:
                definition_matches.append([2,'common_wt'])
            elif protein_rf_name not in vals['protein_families']:
                # if receptor family not the one, then check if it's a same wt match for B-F
                if protein_class_slug!='001':
                    # If class A require two distinct receptor families
                    definition_matches.append([3,'same_wt'])
        if definition_matches:
            min_priority = min(x[0] for x in definition_matches)
            pos = vals['wt'][1]
            wt_aa = vals['wt'][0]
            segment = vals['wt'][2]
            origin = {'pdbs': vals['pdbs'], 'protein_families': vals['protein_families'], 'proteins': vals['proteins'], 'hits':vals['hits']} 
            gpcrdb = gn.split("_")[0]
            for mut_aa in set(vals['mutation'][1]):
                if mut_aa!=wt_aa:
                    mut = {'wt_aa': wt_aa, 'segment': segment, 'pos': pos, 'gpcrdb':gpcrdb, 'mut_aa':mut_aa, 'definitions' : definition_matches, 'priority': min_priority, 'origin': [origin]}
                    key = '%s%s%s' % (wt_aa,pos,mut_aa)
                    # print(key,mut)
                    if key not in simple_list:
                        simple_list[key] = mut
                    else:
                        simple_list[key]['definitions'] += definition_matches
                        min_priority = min(x[0] for x in simple_list[key]['definitions'])
                        simple_list[key]['priority'] = min_priority
                        simple_list[key]['origin'].append(origin)


    # TODO : overlay with other types of mutations, e.g. surfacing expressing 

    # Conservation rules and Helix propensity rules

    rf_conservation = calculate_conservation(slug=protein_rf_slug)
    rf_cutoff = 5
    rf_conservation_priority = 3
    definition_matches = [rf_conservation_priority,'conservation_rf']
    for cons_gn, aa in rf_conservation.items():
        if int(aa[1])>=rf_cutoff and cons_gn in wt_lookup and wt_lookup[cons_gn][0]!=aa[0] and aa[0]!="+":
            # If cons_gn exist in target but AA is not the same
            mut = {'wt_aa': wt_lookup[cons_gn][0], 'segment': wt_lookup[cons_gn][2], 'pos': wt_lookup[cons_gn][1], 'gpcrdb':cons_gn, 'mut_aa':aa[0], 'definitions' : [definition_matches], 'priority': rf_conservation_priority}
            key = '%s%s%s' % (wt_lookup[cons_gn][0],wt_lookup[cons_gn][1],aa[0])
            if key not in simple_list:
                simple_list[key] = mut
            else:
                simple_list[key]['definitions'] += [definition_matches]
                min_priority = min(x[0] for x in simple_list[key]['definitions'])
                simple_list[key]['priority'] = min_priority

        # Apply helix propensity rules
        if cons_gn in GP_residues_in_target:
            if not (wt_lookup[cons_gn][0]==aa[0] and int(aa[1])>5):
                rule = [2,"remove_unconserved_%s" % wt_lookup[cons_gn][0]]
                mut = {'wt_aa': wt_lookup[cons_gn][0], 'segment': wt_lookup[cons_gn][2], 'pos': wt_lookup[cons_gn][1], 'gpcrdb':cons_gn, 'mut_aa':'A', 'definitions' : [rule], 'priority': 2}
                key = '%s%s%s' % (wt_lookup[cons_gn][0],wt_lookup[cons_gn][1],'A')
                if key not in simple_list:
                    simple_list[key] = mut
                else:
                    simple_list[key]['definitions'] += [rule]
                    min_priority = min(x[0] for x in simple_list[key]['definitions'])
                    simple_list[key]['priority'] = min_priority


    class_conservation = calculate_conservation(slug=protein_class_slug)
    class_cutoff = 7
    class_conservation_priority = 4
    definition_matches = [class_conservation_priority,'conservation_class']
    for cons_gn, aa in class_conservation.items():
        if int(aa[1])>=class_cutoff and cons_gn in wt_lookup and wt_lookup[cons_gn][0]!=aa[0] and aa[0]!="+":
            # If cons_gn exist in target but AA is not the same
            mut = {'wt_aa': wt_lookup[cons_gn][0], 'segment': wt_lookup[cons_gn][2], 'pos': wt_lookup[cons_gn][1], 'gpcrdb':cons_gn, 'mut_aa':aa[0], 'definitions' : [definition_matches], 'priority': class_conservation_priority}
            key = '%s%s%s' % (wt_lookup[cons_gn][0],wt_lookup[cons_gn][1],aa[0])
            if key not in simple_list:
                simple_list[key] = mut
            else:
                simple_list[key]['definitions'] += [definition_matches]
                min_priority = min(x[0] for x in simple_list[key]['definitions'])
                simple_list[key]['priority'] = min_priority

        # Apply helix propensity rules from class when receptor family only has one member or non-classA
        if (protein_rf_count==1 or protein_class_slug!='001') and cons_gn in GP_residues_in_target:
            if not (wt_lookup[cons_gn][0]==aa[0] and int(aa[1])>5):
                rule = [2,"remove_unconserved_%s" % wt_lookup[cons_gn][0]]
                mut = {'wt_aa': wt_lookup[cons_gn][0], 'segment': wt_lookup[cons_gn][2], 'pos': wt_lookup[cons_gn][1], 'gpcrdb':cons_gn, 'mut_aa':'A', 'definitions' : [rule], 'priority': 2}
                key = '%s%s%s' % (wt_lookup[cons_gn][0],wt_lookup[cons_gn][1],'A')
                if key not in simple_list:
                    simple_list[key] = mut
                else:
                    if rule not in simple_list[key]['definitions']:
                        # Do not add this rule if it is already there (From RF check)
                        simple_list[key]['definitions'] += [rule]
                        min_priority = min(x[0] for x in simple_list[key]['definitions'])
                        simple_list[key]['priority'] = min_priority


    xtals_conservation = cache.get("CD_xtal_cons_"+protein_class_slug)
    if not xtals_conservation:
        c_proteins = Construct.objects.filter(protein__family__slug__startswith = protein_class_slug).all().values_list('protein__pk', flat = True).distinct()
        xtal_proteins = Protein.objects.filter(pk__in=c_proteins)
        print(xtal_proteins)
        xtals_conservation = calculate_conservation(proteins=xtal_proteins)
        cache.set("CD_xtal_cons_"+protein_class_slug,xtals_conservation,60*60*24)

    xtals_cutoff = 5
    xtals_conservation_priority = 4
    definition_matches = [xtals_conservation_priority,'conservation_xtals']
    for cons_gn, aa in class_conservation.items():
        if int(aa[1])>=xtals_cutoff and cons_gn in wt_lookup and wt_lookup[cons_gn][0]!=aa[0] and aa[0]!="+":
            # If cons_gn exist in target but AA is not the same
            mut = {'wt_aa': wt_lookup[cons_gn][0], 'segment': wt_lookup[cons_gn][2], 'pos': wt_lookup[cons_gn][1], 'gpcrdb':cons_gn, 'mut_aa':aa[0], 'definitions' : [definition_matches], 'priority': xtals_conservation_priority}
            key = '%s%s%s' % (wt_lookup[cons_gn][0],wt_lookup[cons_gn][1],aa[0])
            if key not in simple_list:
                simple_list[key] = mut
            else:
                simple_list[key]['definitions'] += [definition_matches]
                min_priority = min(x[0] for x in simple_list[key]['definitions'])
                simple_list[key]['priority'] = min_priority

    simple_list = OrderedDict(sorted(simple_list.items(), key=lambda x: (x[1]['priority'],x[1]['pos']) ))
    for key, val in simple_list.items():
        val['definitions'] = [x[1] for x in val['definitions']]

    jsondata = simple_list
    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    diff = round(time.time() - start_time,1)
    print("muts",diff)
    return HttpResponse(jsondata, **response_kwargs)

@cache_page(60 * 60 * 24)
def cons_strucs(request, slug, **response_kwargs):
    start_time = time.time()

    level = Protein.objects.filter(entry_name=slug).values_list('family__slug', flat = True).get()
    ##PREPARE TM1 LOOKUP DATA
    c_proteins = Construct.objects.filter(protein__family__slug__startswith = level.split("_")[0]).all().values_list('protein__pk', flat = True).distinct()
    xtal_proteins = Protein.objects.filter(pk__in=c_proteins)
    align_segments = ProteinSegment.objects.all().filter(slug__in = list(settings.REFERENCE_POSITIONS.keys())).prefetch_related()

    amino_acids_stats = {}
    amino_acids_groups_stats = {}
        
    potentials = cache.get("CD_xtal_"+level.split("_")[0])

    if potentials==None:

        a = Alignment()

        a.load_proteins(xtal_proteins)

        a.load_segments(align_segments) #get all segments to make correct diagrams

        # build the alignment data matrix
        a.build_alignment()

        # calculate consensus sequence + amino acid and feature frequency
        a.calculate_statistics()

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
                    amino_acids_stats[n] = temp
                    amino_acids_groups_stats[n] = temp2
                a_id += 1
            s_id += 1

        potentials = {}
        for seg, aa_list in a.consensus.items():
            for gn, aa in aa_list.items():
                if int(aa[1])>5: #if conservations is >50%
                    potentials[gn] = [aa[0],aa[1]]
        cache.set("CD_xtal_"+level.split("_")[0],potentials,60*60*24)


    rs = Residue.objects.filter(protein_conformation__protein__entry_name=slug, generic_number__label__in=list(potentials.keys())).prefetch_related('protein_segment','display_generic_number','generic_number')

    results = {}
    for r in rs:
        gn = r.generic_number.label
        if r.amino_acid!=potentials[gn][0]:
            results[gn] = [r.amino_acid, r.sequence_number,potentials[gn][0],potentials[gn][1]]
    jsondata = json.dumps(results)
    response_kwargs['content_type'] = 'application/json'
    end_time = time.time()
    diff = round(end_time - start_time,1)
    print("cons_strucs",diff)
    return HttpResponse(jsondata, **response_kwargs)

@cache_page(60 * 60 * 24)
def cons_rf(request, slug, **response_kwargs):
    start_time = time.time()

    level = Protein.objects.filter(entry_name=slug).values_list('family__slug', flat = True).get()
    ##PREPARE TM1 LOOKUP DATA
    #c_proteins = Construct.objects.filter(protein__family__slug__startswith = level.split("_")[0]).all().values_list('protein__pk', flat = True).distinct()
    rf_proteins = Protein.objects.filter(family__slug__startswith="_".join(level.split("_")[0:3]), source__name='SWISSPROT',species__common_name='Human')
    align_segments = ProteinSegment.objects.all().filter(slug__in = list(settings.REFERENCE_POSITIONS.keys())).prefetch_related()

    amino_acids_stats = {}
    amino_acids_groups_stats = {}
        

    print(len(rf_proteins))

    try:

        # Load alignment 
        a = pickle.loads(AlignmentConsensus.objects.get(slug="_".join(level.split("_")[0:3])).alignment)
    except:
        print('failed!')
        a = Alignment()

        a.load_proteins(rf_proteins)

        a.load_segments(align_segments) #get all segments to make correct diagrams

        # build the alignment data matrix
        a.build_alignment()

        # calculate consensus sequence + amino acid and feature frequency
        a.calculate_statistics()

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
                amino_acids_stats[n] = temp
                amino_acids_groups_stats[n] = temp2
            a_id += 1
        s_id += 1

    potentials = {}
    for seg, aa_list in a.consensus.items():
        for gn, aa in aa_list.items():
            if int(aa[1])>5: #if conservations is >50%
                potentials[gn] = [aa[0],aa[1]]


    rs = Residue.objects.filter(protein_conformation__protein__entry_name=slug, generic_number__label__in=list(potentials.keys())).prefetch_related('protein_segment','display_generic_number','generic_number')

    results = {}
    for r in rs:
        gn = r.generic_number.label
        if r.amino_acid!=potentials[gn][0]:
            results[gn] = [r.amino_acid, r.sequence_number,potentials[gn][0],potentials[gn][1]]
    jsondata = json.dumps(results)
    response_kwargs['content_type'] = 'application/json'
    end_time = time.time()
    diff = round(end_time - start_time,1)
    print("cons_rf",diff)
    return HttpResponse(jsondata, **response_kwargs)

@cache_page(60 * 60 * 24)
def cons_rf_and_class(request, slug, **response_kwargs):
    start_time = time.time()

    level = Protein.objects.filter(entry_name=slug).values_list('family__slug', flat = True).get()
    ##PREPARE TM1 LOOKUP DATA
    #c_proteins = Construct.objects.filter(protein__family__slug__startswith = level.split("_")[0]).all().values_list('protein__pk', flat = True).distinct()
    rf_proteins = Protein.objects.filter(family__slug__startswith="_".join(level.split("_")[0:3]), source__name='SWISSPROT',species__common_name='Human')
    align_segments = ProteinSegment.objects.all().filter(slug__in = list(settings.REFERENCE_POSITIONS.keys())).prefetch_related()

    amino_acids_stats = {}
    amino_acids_groups_stats = {}
        
    try:
        # Load alignment 
        a = pickle.loads(AlignmentConsensus.objects.get(slug="_".join(level.split("_")[0:3])).alignment)
    except:
        print('failed!')
        a = Alignment()

        a.load_proteins(rf_proteins)

        a.load_segments(align_segments) #get all segments to make correct diagrams

        # build the alignment data matrix
        a.build_alignment()

        # calculate consensus sequence + amino acid and feature frequency
        a.calculate_statistics()

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
                amino_acids_stats[n] = temp
                amino_acids_groups_stats[n] = temp2
            a_id += 1
        s_id += 1

    potentials = {}
    for seg, aa_list in a.consensus.items():
        for gn, aa in aa_list.items():
            if int(aa[1])>5: #if conservations is >50%
                potentials[gn] = [aa[0],aa[1]]

    potentials2 = cache.get("CD_rfc_"+"_".join(level.split("_")[0:1]))

    if potentials2==None:
        class_proteins = Protein.objects.filter(family__slug__startswith="_".join(level.split("_")[0:1]), source__name='SWISSPROT',species__common_name='Human')
        align_segments = ProteinSegment.objects.all().filter(slug__in = list(settings.REFERENCE_POSITIONS.keys())).prefetch_related()

        amino_acids_stats = {}
        amino_acids_groups_stats = {}
            
        try:
            # Load alignment 
            a = pickle.loads(AlignmentConsensus.objects.get(slug="_".join(level.split("_")[0:1])).alignment)
        except:
            print('failed!')

            a = Alignment()

            a.load_proteins(class_proteins)

            a.load_segments(align_segments) #get all segments to make correct diagrams

            # build the alignment data matrix
            a.build_alignment()

            # calculate consensus sequence + amino acid and feature frequency
            a.calculate_statistics()

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
                    amino_acids_stats[n] = temp
                    amino_acids_groups_stats[n] = temp2
                a_id += 1
            s_id += 1

        potentials2 = {}
        for seg, aa_list in a.consensus.items():
            for gn, aa in aa_list.items():
                if int(aa[1])>5: #if conservations is >50%
                    potentials2[gn] = [aa[0],aa[1]]
        cache.set("CD_rfc_"+"_".join(level.split("_")[0:1]),potentials2,60*60*24)


    rs = Residue.objects.filter(protein_conformation__protein__entry_name=slug, generic_number__label__in=list(potentials.keys())).prefetch_related('protein_segment','display_generic_number','generic_number')

    results = {}
    for r in rs:
        gn = r.generic_number.label
        if r.amino_acid!=potentials[gn][0]:
            if gn in potentials2:
                results[gn] = [r.amino_acid, r.sequence_number,potentials[gn][0],potentials[gn][1]]
    jsondata = json.dumps(results)
    response_kwargs['content_type'] = 'application/json'
    end_time = time.time()
    diff = round(end_time - start_time,1)
    print("cons_rf_and_class",diff)
    return HttpResponse(jsondata, **response_kwargs)

@cache_page(60 * 60 * 24)
def cons_rm_GP(request, slug, **response_kwargs):
    start_time = time.time()
    level = Protein.objects.filter(entry_name=slug).values_list('family__slug', flat = True).get()
    ##PREPARE TM1 LOOKUP DATA
    #c_proteins = Construct.objects.filter(protein__family__slug__startswith = level.split("_")[0]).all().values_list('protein__pk', flat = True).distinct()
    rf_proteins = Protein.objects.filter(family__slug__startswith="_".join(level.split("_")[0:3]), source__name='SWISSPROT',species__common_name='Human')
    align_segments = ProteinSegment.objects.all().filter(slug__in = list(settings.REFERENCE_POSITIONS.keys())).prefetch_related()

    amino_acids_stats = {}
    amino_acids_groups_stats = {}
        
    a = Alignment()
    a.load_proteins(rf_proteins)
    a.load_segments(align_segments) #get all segments to make correct diagrams
    # build the alignment data matrix
    a.build_alignment()
    # calculate consensus sequence + amino acid and feature frequency
    a.calculate_statistics()

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
                amino_acids_stats[n] = temp
                amino_acids_groups_stats[n] = temp2
            a_id += 1
        s_id += 1

    potentials = {}
    for seg, aa_list in a.consensus.items():
        for gn, aa in aa_list.items():
            if int(aa[1])>5: #if conservations is >50%
                potentials[gn] = [aa[0],aa[1]]


    rs = Residue.objects.filter(protein_conformation__protein__entry_name=slug, generic_number__label__in=list(potentials.keys())).prefetch_related('protein_segment','display_generic_number','generic_number')

    results = {}
    results2 = {}
    for r in rs:
        gn = r.generic_number.label
        if r.amino_acid in ['G','P']:
            if r.amino_acid!=potentials[gn][0]:
                results[gn] = [r.amino_acid, r.sequence_number,potentials[gn][0],potentials[gn][1]]
            if r.amino_acid=='G' and potentials[gn][0]=='G':
                results2[gn] = [r.amino_acid, r.sequence_number,'A',potentials[gn][1]]
    jsondata = json.dumps({'non-conserved':results, 'conserved':results2})
    response_kwargs['content_type'] = 'application/json'
    end_time = time.time()
    diff = round(end_time - start_time,1)
    print("cons_rm_GP",diff)
    return HttpResponse(jsondata, **response_kwargs)

def calculate_conservation(proteins = None, slug = None):
    # Return a a dictionary of each generic number and the conserved residue and its frequency
    # Can either be used on a list of proteins or on a slug. If slug then use the cached alignment object.

    amino_acids_stats = {}
    amino_acids_groups_stats = {}

    if slug:
        try: 
            # Load alignment 
            alignment_consensus = AlignmentConsensus.objects.get(slug=slug)
            if alignment_consensus.gn_consensus:
                return pickle.loads(alignment_consensus.gn_consensus)
            a = pickle.loads(alignment_consensus.alignment)
        except: 
            proteins = Protein.objects.filter(family__slug__startswith=slug, source__name='SWISSPROT',species__common_name='Human')
            align_segments = ProteinSegment.objects.all().filter(slug__in = list(settings.REFERENCE_POSITIONS.keys())).prefetch_related()
            a = Alignment()
            a.load_proteins(proteins)
            a.load_segments(align_segments) 
            a.build_alignment()
            # calculate consensus sequence + amino acid and feature frequency
            a.calculate_statistics()
    elif proteins:
        align_segments = ProteinSegment.objects.all().filter(slug__in = list(settings.REFERENCE_POSITIONS.keys())).prefetch_related()
        a = Alignment()
        a.load_proteins(proteins)
        a.load_segments(align_segments) 
        a.build_alignment()
        # calculate consensus sequence + amino acid and feature frequency
        a.calculate_statistics()


    consensus = {}
    for seg, aa_list in a.consensus.items():
        for gn, aa in aa_list.items():
            if 'x' in gn: # only takes those GN positions that are actual 1x50 etc
                consensus[gn] = [aa[0],aa[1]]
    if slug and alignment_consensus:
        alignment_consensus.gn_consensus = pickle.dumps(consensus)
        alignment_consensus.save()

    return consensus
