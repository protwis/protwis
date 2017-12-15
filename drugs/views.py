from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.conf import settings
from django.db.models import Count, Max
from django.core.cache import cache
from django.views.decorators.cache import cache_page

from drugs.models import Drugs
from protein.models import Protein, ProteinFamily
from mutational_landscape.models import NHSPrescribings

import re
import json
import numpy as np
from collections import OrderedDict
from copy import deepcopy

def get_spaced_colors(n):
    max_value = 16581375 #255**3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]

    return ['#'+i for i in colors] # HEX colors
    # return [(int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors] # RGB colors

def striphtml(data):
    p = re.compile(r'<.*?>')
    return p.sub('', data)

@cache_page(60*60*24*30)
def drugstatistics(request):

    # ===== drugtargets =====
    drugtargets_raw_approved = Protein.objects.filter(drugs__status='approved').values('entry_name').annotate(value=Count('drugs__name', distinct = True)).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugtargets_raw_approved))
    drugtargets_approved = []
    for i, drugtarget in enumerate(drugtargets_raw_approved):
        drugtarget['label'] = drugtarget['entry_name'].replace("_human","").upper()
        # drugtarget['color'] = str(list_of_hec_colors[i])
        del drugtarget['entry_name']
        drugtargets_approved.append(drugtarget)

    drugtargets_raw_trials = Protein.objects.filter(drugs__status__in=['in trial'], drugs__clinicalstatus__in=['completed','not open yet','ongoing','recruiting','suspended']).values('entry_name').annotate(value=Count('drugs__name', distinct = True)).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugtargets_raw_trials))
    drugtargets_trials = []
    for i, drugtarget in enumerate(drugtargets_raw_trials):
        drugtarget['label'] = drugtarget['entry_name'].replace("_human","").upper()
        # drugtarget['color'] = str(list_of_hec_colors[i])
        del drugtarget['entry_name']
        drugtargets_trials.append(drugtarget)


    all_human_GPCRs = Protein.objects.filter(species_id=1, sequence_type_id=1, family__slug__startswith='00').distinct()

    in_trial = Protein.objects.filter(drugs__status__in=['in trial'] ).exclude(drugs__status='approved').distinct() #drugs__clinicalstatus__in=['completed','not open yet','ongoing','recruiting','suspended']

    not_targeted = len(all_human_GPCRs) - len(drugtargets_approved) - len(in_trial)

    # ===== drugfamilies =====
    drugfamilies_raw_approved = Protein.objects.filter(drugs__status='approved').values('family_id__parent__name').annotate(value=Count('drugs__name', distinct = True))

    list_of_hec_colors = get_spaced_colors(len(drugfamilies_raw_approved))
    drugfamilies_approved = []
    for i, drugfamily in enumerate(drugfamilies_raw_approved):
        drugfamily['label'] = striphtml(drugfamily['family_id__parent__name']).replace(" receptors","")
        drugfamily['color'] = str(list_of_hec_colors[i])
        del drugfamily['family_id__parent__name']
        drugfamilies_approved.append(drugfamily)

    drugfamilies_raw_trials = Protein.objects.exclude(drugs__status='approved').values('family_id__parent__name').annotate(value=Count('drugs__name', distinct = True))

    list_of_hec_colors = get_spaced_colors(len(drugfamilies_raw_trials))
    drugfamilies_trials = []
    for i, drugfamily in enumerate(drugfamilies_raw_trials):
        drugfamily['label'] = striphtml(drugfamily['family_id__parent__name']).replace(" receptors","")
        drugfamily['color'] = str(list_of_hec_colors[i])
        del drugfamily['family_id__parent__name']
        drugfamilies_trials.append(drugfamily)

    # ===== drugclas =====
    drugclasses_raw_approved = Protein.objects.filter(drugs__status='approved').values('family_id__parent__parent__parent__name').annotate(value=Count('drugs__name', distinct = True)).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugclasses_raw_approved)+1)
    drugClasses_approved = []
    for i, drugclas in enumerate(drugclasses_raw_approved):
        drugclas['label'] = drugclas['family_id__parent__parent__parent__name']
        drugclas['color'] = str(list_of_hec_colors[i+1])
        del drugclas['family_id__parent__parent__parent__name']
        drugClasses_approved.append(drugclas)

    drugclasses_raw_trials = Protein.objects.filter(drugs__status__in=['in trial'], drugs__clinicalstatus__in=['completed','not open yet','ongoing','recruiting','suspended']).values('family_id__parent__parent__parent__name').annotate(value=Count('drugs__name', distinct = True)).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugclasses_raw_trials)+1)
    drugClasses_trials = []
    for i, drugclas in enumerate(drugclasses_raw_trials):
        drugclas['label'] = drugclas['family_id__parent__parent__parent__name']
        drugclas['color'] = str(list_of_hec_colors[i+1])
        del drugclas['family_id__parent__parent__parent__name']
        drugClasses_trials.append(drugclas)

    # ===== drugtypes =====
    drugtypes_raw_approved = Drugs.objects.values('drugtype').filter(status='approved').annotate(value=Count('name', distinct = True)).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugtypes_raw_approved)+20)
    drugtypes_approved = []
    for i, drugtype in enumerate(drugtypes_raw_approved):
        drugtype['label'] = drugtype['drugtype']
        drugtype['color'] = str(list_of_hec_colors[i])
        del drugtype['drugtype']
        drugtypes_approved.append(drugtype)

    drugtypes_raw_trials = Drugs.objects.values('drugtype').filter(status='in trial', clinicalstatus__in=['completed','not open yet','ongoing','recruiting','suspended']).annotate(value=Count('name', distinct = True)).order_by('-value')

    # list_of_hec_colors = get_spaced_colors(len(drugtypes_raw_trials)+5)
    drugtypes_trials = []
    for i, drugtype in enumerate(drugtypes_raw_trials):
        drugtype['label'] = drugtype['drugtype']
        drugtype['color'] = str(list_of_hec_colors[i])
        del drugtype['drugtype']
        drugtypes_trials.append(drugtype)

    drugtypes_raw_not_estab = Drugs.objects.values('drugtype').filter(novelty='not established').annotate(value=Count('name', distinct = True)).order_by('-value')

    # list_of_hec_colors = get_spaced_colors(len(drugtypes_raw_not_estab)+5)
    drugtypes_not_estab = []
    for i, drugtype in enumerate(drugtypes_raw_not_estab):
        drugtype['label'] = drugtype['drugtype']
        drugtype['color'] = str(list_of_hec_colors[i])
        del drugtype['drugtype']
        drugtypes_not_estab.append(drugtype)

    drugtypes_raw_estab = Drugs.objects.values('drugtype').filter(novelty='established').annotate(value=Count('name', distinct = True)).order_by('-value')

    # list_of_hec_colors = get_spaced_colors(len(drugtypes_raw_estab)+5)
    drugtypes_estab = []
    for i, drugtype in enumerate(drugtypes_raw_estab):
        drugtype['label'] = drugtype['drugtype']
        drugtype['color'] = str(list_of_hec_colors[i])
        del drugtype['drugtype']
        drugtypes_estab.append(drugtype)

    # ===== modes of action =====
    moas_raw_approved = Drugs.objects.values('moa').filter(status='approved').annotate(value=Count('name', distinct = True)).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(moas_raw_approved)+5)
    moas_approved = []
    for i, moa in enumerate(moas_raw_approved):
        moa['label'] = moa['moa']
        moa['color'] = str(list_of_hec_colors[i])
        del moa['moa']
        moas_approved.append(moa)

    moa_raw_trials = Drugs.objects.values('moa').filter(status='in trial', clinicalstatus__in=['completed','not open yet','ongoing','recruiting','suspended']).annotate(value=Count('name', distinct = True)).order_by('-value')

    # list_of_hec_colors = get_spaced_colors(len(moa_raw_trials)+5)
    moas_trials = []
    for i, moa in enumerate(moa_raw_trials):
        moa['label'] = moa['moa']
        moa['color'] = str(list_of_hec_colors[i])
        del moa['moa']
        moas_trials.append(moa)

    # ===== Phase distributions =====
    # Distinguish between different Clinical Status
    phases_raw_active = Drugs.objects.values('phase').filter(status='in trial', clinicalstatus__in=['completed','not open yet','ongoing','recruiting','suspended']).annotate(value=Count('name', distinct = True)).order_by('-value')

    phase_trials = []
    list_of_hec_colors = ["#88df8c", "#43A047", "#b0f2b2"]
    for i, phase in enumerate(phases_raw_active):
        phase['label'] = 'Phase ' + phase['phase']
        phase['color'] = str(list_of_hec_colors[i])
        del phase['phase']
        phase_trials.append(phase)

    phases_raw_inactive = Drugs.objects.values('phase').filter(status='in trial', clinicalstatus__in=['terminated','discontinued','unknown','withdrawn']).annotate(value=Count('name', distinct = True)).order_by('-value')

    phase_trials_inactive = []
    list_of_hec_colors = ["#88df8c", "#43A047", "#b0f2b2"]
    for i, phase in enumerate(phases_raw_inactive):
        phase['label'] = 'Phase ' + phase['phase']
        phase['color'] = str(list_of_hec_colors[i])
        del phase['phase']
        phase_trials_inactive.append(phase)

    # ===== drugindications =====
    drugindications_raw_approved = Drugs.objects.values('indication').filter(status='approved').annotate(value=Count('name', distinct = True)).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugindications_raw_approved)+10)
    drugindications_approved = []
    for i, drugindication in enumerate(drugindications_raw_approved):
        drugindication['label'] = drugindication['indication']
        drugindication['color'] = str(list_of_hec_colors[i])
        del drugindication['indication']
        drugindications_approved.append(drugindication)

    drugindications_raw_trials = Drugs.objects.values('indication').filter(status='in trial', clinicalstatus__in=['completed','not open yet','ongoing','recruiting','suspended']).annotate(value=Count('name', distinct = True)).order_by('-value')

    # list_of_hec_colors = get_spaced_colors(len(drugindications_raw_trials))
    drugindications_trials = []
    for i, drugindication in enumerate(drugindications_raw_trials):
        drugindication['label'] = drugindication['indication']
        drugindication['color'] = str(list_of_hec_colors[i])
        del drugindication['indication']
        drugindications_trials.append(drugindication)

    # ===== drugtimes =====
    drugtime_raw = Drugs.objects.values('approval').filter(status='approved').annotate(y=Count('name', distinct = True)).order_by('approval')

    drugtimes = []
    running_total = 0

    for i, time in enumerate(range(1942,2017,1)):
        if str(time) in [i['approval'] for i in drugtime_raw]:
            y = [i['y'] for i in drugtime_raw if i['approval']==str(time)][0] + running_total
            x = time
            running_total = y
        else:
            x = time
            y = running_total
        if time % 2 == 0:
            drugtimes.append({'x':x,'y':y})

    drugs_over_time = [{"values": drugtimes, "yAxis": "1", "key": "GPCRs"}, {'values': [{'y': 2, 'x': '1942'}, {'x': '1944', 'y': 2}, {'y': 6, 'x': '1946'}, {'y': 9, 'x': '1948'}, {'y': 18, 'x': '1950'}, {'y': 30, 'x': '1952'}, {'y': 55, 'x': '1954'}, {'y': 72, 'x': '1956'}, {'y': 98, 'x': '1958'}, {'y': 131, 'x': '1960'}, {'y': 153, 'x': '1962'}, {'y': 171, 'x': '1964'}, {'y': 188, 'x': '1966'}, {'y': 205, 'x': '1968'}, {'y': 224, 'x': '1970'}, {'y': 242, 'x': '1972'}, {'y': 265, 'x': '1974'}, {'y': 300, 'x': '1976'}, {'y': 340, 'x': '1978'}, {'y': 361, 'x': '1980'}, {'y': 410, 'x': '1982'}, {'y': 442, 'x': '1984'}, {'y': 499, 'x': '1986'}, {'y': 542, 'x': '1988'}, {'y': 583, 'x': '1990'}, {'y': 639, 'x': '1992'}, {'y': 686, 'x': '1994'}, {'y': 779, 'x': '1996'}, {'y': 847, 'x': '1998'}, {'y': 909, 'x': '2000'}, {'y': 948, 'x': '2002'}, {'y': 1003, 'x': '2004'}, {'y': 1041, 'x': '2006'}, {'y': 1078, 'x': '2008'}, {'y': 1115, 'x': '2010'}, {'y': 1177, 'x': '2012'}, {'y': 1239, 'x': '2014'}, {'y': 1286, 'x': '2016'}], 'key': 'All FDA drugs', 'yAxis': '1'}]

    # ===== drugtimes =====


    return render(request, 'drugstatistics.html', {'drugtypes_approved':drugtypes_approved, 'drugtypes_trials':drugtypes_trials,  'drugtypes_estab':drugtypes_estab,  'drugtypes_not_estab':drugtypes_not_estab, 'drugindications_approved':drugindications_approved, 'drugindications_trials':drugindications_trials, 'drugtargets_approved':drugtargets_approved, 'drugtargets_trials':drugtargets_trials, 'phase_trials':phase_trials, 'phase_trials_inactive': phase_trials_inactive, 'moas_trials':moas_trials, 'moas_approved':moas_approved, 'drugfamilies_approved':drugfamilies_approved, 'drugfamilies_trials':drugfamilies_trials, 'drugClasses_approved':drugClasses_approved, 'drugClasses_trials':drugClasses_trials, 'drugs_over_time':drugs_over_time, 'in_trial':len(in_trial), 'not_targeted':not_targeted})

@cache_page(60*60*24*30)
def drugbrowser(request):
    # Get drugdata from here somehow

    name_of_cache = 'drug_browse3'

    context = cache.get(name_of_cache)

    if context == None:
        context = list()

        drugs = Drugs.objects.all().prefetch_related('target__family__parent__parent__parent')

        drugs_NHS_names = NHSPrescribings.objects.values_list('drugname__name', flat=True).distinct()

        for drug in drugs:
            drugname = drug.name
            drugtype = drug.drugtype
            status = drug.status
            approval = drug.approval
            targetlevel = drug.targetlevel
            phase = drug.phase
            moa = drug.moa
            indication = drug.indication
            novelty = drug.novelty
            clinicalstatus = drug.clinicalstatus
            references = [i for i in drug.references.split('|')]

            if drugname in drugs_NHS_names:
                NHS = 'yes'
            else:
                NHS = 'no'

            target_list = drug.target.all()
            for protein in target_list:

                clas = str(protein.family.parent.parent.parent.name)
                family = str(protein.family.parent.name)

                jsondata = {'name': drugname, 'target': str(protein), 'phase': phase, 'approval': approval, 'class': clas, 'family': family, 'indication': indication, 'status': status, 'drugtype': drugtype, 'moa': moa, 'novelty': novelty, 'targetlevel': targetlevel, 'clinicalstatus': clinicalstatus, 'references': references, 'NHS': NHS}
                context.append(jsondata)

            cache.set(name_of_cache, context, 60*60*24*25)   # two days timeout on cache

    return render(request, 'drugbrowser.html', {'drugdata': context})

@cache_page(60*60*24*30)
def drugmapping(request):
    context = dict()

    families = ProteinFamily.objects.all()
    lookup = {}
    for f in families:
        lookup[f.slug] = f.name.replace("receptors","").replace(" receptor","").replace(" hormone","").replace("/neuropeptide","/").replace(" (G protein-coupled)","").replace(" factor","").replace(" (LPA)","").replace(" (S1P)","").replace("GPR18, GPR55 and GPR119","GPR18/55/119").replace("-releasing","").replace(" peptide","").replace(" and oxytocin","/Oxytocin").replace("Adhesion class orphans","Adhesion orphans").replace("muscarinic","musc.").replace("-concentrating","-conc.")

    class_proteins = Protein.objects.filter(family__slug__startswith="00",source__name='SWISSPROT', species_id=1).prefetch_related('family').order_by('family__slug')

    temp = OrderedDict([
                    ('name',''),
                    ('trials', 0),
                    ('maxphase', 0),
                    ('approved', 0),
                    ('family_sum_approved', 0),
                    ('family_sum_trials' , 0),
                    ('establishment', 2),
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

    # # POULATE WITH DATA
    total_approved = 0
    drugtargets_approved_class = Protein.objects.filter(drugs__status='approved').values('family_id__parent__parent__parent__slug').annotate(value=Count('drugs__name', distinct = True))
    for i in drugtargets_approved_class:
        fid = i['family_id__parent__parent__parent__slug'].split("_")
        coverage[fid[0]]['family_sum_approved'] += i['value']
        total_approved += i['value']

    drugtargets_approved_type = Protein.objects.filter(drugs__status='approved').values('family_id__parent__parent__slug').annotate(value=Count('drugs__name', distinct = True))
    for i in drugtargets_approved_type:
        fid = i['family_id__parent__parent__slug'].split("_")
        coverage[fid[0]]['children'][fid[1]]['family_sum_approved'] += i['value']

    drugtargets_approved_family = Protein.objects.filter(drugs__status='approved').values('family_id__parent__slug').annotate(value=Count('drugs__name', distinct = True))
    for i in drugtargets_approved_family:
        fid = i['family_id__parent__slug'].split("_")
        coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['family_sum_approved'] += i['value']
    drugtargets_approved_target = Protein.objects.filter(drugs__status='approved').values('family_id__slug').annotate(value=Count('drugs__name', distinct = True)).annotate(maxphase=Max('drugs__phase'))
    for i in drugtargets_approved_target:
        fid = i['family_id__slug'].split("_")
        coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['approved'] += i['value']
        if int(i['maxphase']) > coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['maxphase']:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['maxphase'] = int(i['maxphase'])
        if i['value'] > 0:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['establishment'] = 4

    total_trials = 0
    drugtargets_trials_class = Protein.objects.filter(drugs__status__in=['in trial']).values('family_id__parent__parent__parent__slug').annotate(value=Count('drugs__name', distinct = True))
    for i in drugtargets_trials_class:
        fid = i['family_id__parent__parent__parent__slug'].split("_")
        coverage[fid[0]]['family_sum_trials'] += i['value']
        total_trials += i['value']

    drugtargets_trials_type = Protein.objects.filter(drugs__status__in=['in trial']).values('family_id__parent__parent__slug').annotate(value=Count('drugs__name', distinct = True))
    for i in drugtargets_trials_type:
        fid = i['family_id__parent__parent__slug'].split("_")
        coverage[fid[0]]['children'][fid[1]]['family_sum_trials'] += i['value']

    drugtargets_trials_family = Protein.objects.filter(drugs__status__in=['in trial']).values('family_id__parent__slug').annotate(value=Count('drugs__name', distinct = True))
    for i in drugtargets_trials_family:
        fid = i['family_id__parent__slug'].split("_")
        coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['family_sum_trials'] += i['value']

    drugtargets_trials_target = Protein.objects.filter(drugs__status__in=['in trial']).values('family_id__slug').annotate(value=Count('drugs__name', distinct = True)).annotate(maxphase=Max('drugs__phase'))
    # add highest reached trial here
    for i in drugtargets_trials_target:
        fid = i['family_id__slug'].split("_")
        coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['trials'] += i['value']
        if int(i['maxphase']) > coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['maxphase']:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['maxphase'] = int(i['maxphase'])
        if i['value'] > 0 and coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['establishment'] == 2:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['establishment'] = 7

    # MAKE THE TREE
    tree = OrderedDict({'name':'GPCRome', 'family_sum_approved': total_approved, 'family_sum_trials': total_trials,'children':[]})
    i = 0
    n = 0
    for c,c_v in coverage.items():
        c_v['name'] = c_v['name'].split("(")[0]
        if c_v['name'].strip() == 'Other GPCRs':
            # i += 1
            continue
            # pass
        children = []
        for lt,lt_v in c_v['children'].items():
            children_rf = []
            for rf,rf_v in lt_v['children'].items():
                rf_v['name'] = rf_v['name'].split("<")[0]
                # if rf_v['name'].strip() == 'Taste 2':
                    # continue
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

    jsontree = json.dumps(tree)

    context["drugdata"] = jsontree

    return render(request, 'drugmapping.html', {'drugdata':context})

@cache_page(60*60*24*15)
def nhs_drug(request, slug):

    nhs_data = NHSPrescribings.objects.filter(drugname__name=slug.lower()).order_by('date')

    data_dic = {}
    sections = []
    query_translate = {}
    for i in nhs_data:
        prescription_name = i.op_name +' (' + i.drugCode + ')'
        queryname = i.drugname.name

        if not prescription_name in data_dic:
            data_dic[prescription_name] = []
            sections.append(i.bnf_section)
        dic = {}
        dic['x'] = str(i.date)
        dic['y'] = int(i.actual_cost)
        data_dic[prescription_name].append(dic)

        if not prescription_name in query_translate:
            query_translate[prescription_name] = queryname

    data = []
    for nhs_name in data_dic.keys():
        data.append({'values': data_dic[nhs_name], 'query_key':str(query_translate[nhs_name]), 'key':nhs_name})

    return render(request, 'nhs.html', {'data':data, 'drug':slug, 'section':list(set(sections))})

@cache_page(60*60*24*15)
def nhs_section(request, slug):

    nhs_data = NHSPrescribings.objects.filter(bnf_section=slug).order_by('date')

    data_dic = {}
    sections = []
    query_translate = {}
    for i in nhs_data:
        prescription_name = i.op_name +' (' + i.drugCode + ')'
        queryname = i.drugname.name

        if not prescription_name in data_dic:
            data_dic[prescription_name] = []
            sections.append(i.bnf_section)

        dic = {}
        dic['x'] = str(i.date)
        dic['y'] = int(i.actual_cost)
        data_dic[prescription_name].append(dic)

        if not prescription_name in query_translate:
            query_translate[prescription_name] = queryname

    data = []
    for nhs_name in data_dic.keys():
        data.append({'values': data_dic[nhs_name], 'query_key':str(query_translate[nhs_name]), 'key':nhs_name})

    return render(request, 'nhs.html', {'data':data, 'drug':slug, 'section':list(set(sections))})
