from django.shortcuts import render
from django.http import HttpResponse
from django.db.models import Min, Count, Max
from django.conf import settings

from construct.models import *
from protein.models import ProteinConformation, Protein

import json
from collections import OrderedDict
import re
import xlrd
import yaml
import os

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

def tool(request):

    simple_selection = request.session.get('selection', False)
    proteins = []
    for target in simple_selection.targets:
        if target.type == 'protein':
            proteins.append(target.item)
    print(proteins)
    context = {}

    context['target'] = proteins[0]


    rs = Residue.objects.filter(protein_conformation__protein=proteins[0]).prefetch_related('protein_segment','display_generic_number','generic_number')

    residues = {}
    residues_gn = {}
    for r in rs:
        segment = r.protein_segment.slug
        segment = segment.replace("-","")
        if segment not in residues:
            residues[segment] = []
        residues[segment].append(r)
        if r.generic_number:
            residues_gn[r.generic_number.label] = r

    context['residues'] = residues
    context['residues_gn'] = residues_gn
    #print(residues)

    return render(request,'tool.html',context)

def json_fusion(request, slug, **response_kwargs):

    level = Protein.objects.filter(entry_name=slug).values_list('family__slug', flat = True).get()
    #proteins = Construct.objects.all().values_list('protein', flat = True)
    cons = Construct.objects.all().prefetch_related('crystal', 'protein__family','deletions')

    jsondata = "glyco"
    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)


def json_palmi(request, slug, **response_kwargs):


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
    palmi['all']= mutations_all

    jsondata = palmi
    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)



def json_glyco(request, slug, **response_kwargs):


    seq = Protein.objects.filter(entry_name=slug).values_list('sequence', flat = True).get()
    rs = Residue.objects.filter(protein_conformation__protein__entry_name=slug).prefetch_related('protein_segment')
    residues = {}
    for r in rs:
        residues[r.sequence_number] = r.protein_segment.slug

    #No proline!
    p = re.compile("N[^P][TS]")
    #print('all')
    mutations_all = []
    for m in p.finditer(seq):
        #print(m.start(), m.group())
        mutations_all.append([m.start()+1,"Q",'','',m.group(),residues[m.start()+1]])

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

        mutations_mammalian.append([m.start()+1,pos0,m.start()+2,pos1,matches_seq[i],residues[m.start()+1]])

    glyco = OrderedDict()
    glyco['all']= mutations_all
    glyco['mammalian'] = mutations_mammalian

    jsondata = glyco
    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)


def json_icl3(request, slug, **response_kwargs):
    level = Protein.objects.filter(entry_name=slug).values_list('family__slug', flat = True).get()

    ##PREPARE TM1 LOOKUP DATA
    proteins = Construct.objects.all().values_list('protein', flat = True)
    tm5_start = {}
    tm5_end = {}
    tm6_start = {}
    tm6_end = {}
    tm5_50 = {}
    tm6_50 = {}
    pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).filter(residue__protein_segment__slug='TM5').annotate(start=Min('residue__sequence_number'), end=Max('residue__sequence_number'))
    for pc in pconfs:
        tm5_start[pc.protein.entry_name] = pc.start
        tm5_end[pc.protein.entry_name] = pc.end
    pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).filter(residue__protein_segment__slug='TM6').annotate(start=Min('residue__sequence_number'), end=Max('residue__sequence_number'))
    for pc in pconfs:
        tm6_start[pc.protein.entry_name] = pc.start
        tm6_end[pc.protein.entry_name] = pc.end


    pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).filter(residue__generic_number__label__in=['5x50','6x50']).annotate(start=Min('residue__sequence_number'), end=Max('residue__sequence_number'))
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
        for deletion in c.deletions.all():
            #print(pdb,deletion.start,deletion.end)
            if deletion.start > tm5_start[entry_name] and deletion.start < tm6_end[entry_name]:
                if p.entry_name not in deletions[d_level_name]:
                    deletions[d_level_name][entry_name] = {}
                #deletions[entry_name][pdb] = [tm5_end[entry_name],tm6_start[entry_name],deletion.start,deletion.end,deletion.start-tm5_end[entry_name],tm6_start[entry_name]-deletion.end]
                deletions[d_level_name][entry_name][pdb] = [deletion.start-tm5_50[entry_name],tm6_50[entry_name]-deletion.end-1,state,str(fusion)]

    # for pdb,state in sorted(states.items()):
    #     print(pdb,"\t",state)

    jsondata = deletions
    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)

def json_nterm(request, slug, **response_kwargs):

    level = Protein.objects.filter(entry_name=slug).values_list('family__slug', flat = True).get()

    ##PREPARE TM1 LOOKUP DATA
    proteins = Construct.objects.all().values_list('protein', flat = True)
    #pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).filter(residue__protein_segment__slug='TM1').annotate(start=Min('residue__sequence_number'))
    pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).filter(residue__generic_number__label__in=['1x50']).values_list('protein__entry_name','residue__sequence_number','residue__generic_number__label')
    tm1_start = {}
    for pc in pconfs:
        tm1_start[pc[0]] = pc[1]


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
        for deletion in c.deletions.all():
            if deletion.start < tm1_start[entry_name]:
                if p.entry_name not in deletions[d_level_name]:
                    deletions[d_level_name][entry_name] = {}
                deletions[d_level_name][entry_name][pdb] = [deletion.start,deletion.end, tm1_start[entry_name]-deletion.end-1,state,str(fusion)]

    jsondata = deletions
    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)

def json_cterm(request, slug, **response_kwargs):
    level = Protein.objects.filter(entry_name=slug).values_list('family__slug', flat = True).get()

    ##PREPARE TM1 LOOKUP DATA
    proteins = Construct.objects.all().values_list('protein', flat = True)
    # pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).filter(residue__protein_segment__slug='C-term').annotate(start=Min('residue__sequence_number'))
    # cterm_start = {}
    # for pc in pconfs:
    #     cterm_start[pc.protein.entry_name] = pc.start
    pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).filter(residue__generic_number__label__in=['8x50']).values_list('protein__entry_name','residue__sequence_number','residue__generic_number__label')
    cterm_start = {}
    for pc in pconfs:
        cterm_start[pc[0]] = pc[1]


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
        for deletion in c.deletions.all():
            if deletion.start >= cterm_start[entry_name]:
                if p.entry_name not in deletions[d_level_name]:
                    deletions[d_level_name][entry_name] = {}
                deletions[d_level_name][entry_name][pdb] = [deletion.start,deletion.end, cterm_start[entry_name]-deletion.start,state,str(fusion)]

    jsondata = deletions
    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)

def thermostabilising(request, slug, **response_kwargs):


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
        if entry_name == slug:
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

    temp = {}
    for gn, vals1 in results['3'].items():
        for mut_aa, vals2 in vals1.items():
            if vals2['hits']>1:
                if gn not in temp:
                    temp[gn] = {}
                if mut_aa not in temp[gn]:
                    temp[gn][mut_aa] = vals2
                #results['2'][gn].pop(mut_aa, None)
    results['3'] = temp


    jsondata = results
    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)



def mutations(request, slug, **response_kwargs):

    ##PREPARE TM1 LOOKUP DATA
    # proteins = Construct.objects.all().values_list('protein', flat = True)
    # pconfs = ProteinConformation.objects.filter(protein_id__in=proteins).filter(residue__protein_segment__slug='C-term').annotate(start=Min('residue__sequence_number'))
    # cterm_start = {}
    # for pc in pconfs:
    #     cterm_start[pc.protein.entry_name] = pc.start

    cons = Construct.objects.all().prefetch_related('crystal', 'protein','mutations')
    mutations = []
    positions = []
    proteins = []
    for c in cons:
        p = c.protein
        entry_name = p.entry_name
        pdb = c.crystal.pdb_code
        for mutation in c.mutations.all():
            if p.entry_name not in proteins:
                proteins.append(entry_name)
            mutations.append((mutation,entry_name,pdb))
            if mutation.sequence_number not in positions:
                positions.append(mutation.sequence_number)

    #print(positions)
    #print(proteins)

    rs = Residue.objects.filter(protein_conformation__protein__entry_name__in=proteins, sequence_number__in=positions).prefetch_related('generic_number','protein_conformation__protein')

    rs_lookup = {}
    gns = []
    for r in rs:
        if not r.generic_number:
            continue #skip non GN
        entry_name = r.protein_conformation.protein.entry_name
        gn = r.generic_number.label
        pos = r.sequence_number
        if entry_name not in rs_lookup:
            rs_lookup[entry_name] = {}
        if pos not in rs_lookup[entry_name]:
            rs_lookup[entry_name][pos] = gn
        if gn not in gns:
            gns.append(gn)
    #print(rs_lookup)

    rs = Residue.objects.filter(protein_conformation__protein__entry_name=slug, generic_number__label__in=gns).prefetch_related('protein_segment','display_generic_number','generic_number')

    wt_lookup = {}
    for r in rs:
        gn = r.generic_number.label
        wt_lookup[gn] = [r.amino_acid, r.sequence_number]
  
    mutation_list = OrderedDict()
    for mutation in mutations:
        pos = mutation[0].sequence_number
        entry_name = mutation[1]
        if entry_name not in rs_lookup:
            continue
        if pos not in rs_lookup[entry_name]:
            continue
        gn = rs_lookup[entry_name][pos]

        if gn not in mutation_list:
            mutation_list[gn] = {'proteins':[], 'hits':0, 'mutation':[], 'wt':''}

        if entry_name not in mutation_list[gn]['proteins']:
            mutation_list[gn]['proteins'].append(entry_name)
            mutation_list[gn]['hits'] += 1  
            mutation_list[gn]['mutation'].append((mutation[0].wild_type_amino_acid,mutation[0].mutated_amino_acid))
            if gn in wt_lookup:
                mutation_list[gn]['wt'] = wt_lookup[gn]


    for gn, vals in mutation_list.items():
        if vals['hits']<2:
            mutation_list.pop(gn, None)

    mutation_list = OrderedDict(sorted(mutation_list.items(), key=lambda x: x[1]['hits'],reverse=True))
    #print(mutation_list)


    #print(json.dumps(mutation_list, indent=4))

    jsondata = mutation_list
    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)