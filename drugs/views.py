from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.conf import settings
from django.db.models import Count, Max, Q, F, Value, CharField, Case, When, IntegerField
from django.db.models import Count, Max
from django.core.cache import cache
from django.views.decorators.cache import cache_page

from drugs.models import Drugs,Drugs2024
from protein.models import Protein, ProteinFamily, TissueExpression
from structure.models import Structure
from drugs.models import Drugs, Drugs2024, Indication
from protein.models import Protein, ProteinFamily, Tissues, TissueExpression
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

@cache_page(60 * 60 * 24 * 28)
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

@cache_page(60 * 60 * 24 * 28)
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

            cache.set(name_of_cache, context, 60*60*24*28)

    return render(request, 'drugbrowser.html', {'drugdata': context})

class NewDrugsBrowser(TemplateView):
    # Template using this class #
    template_name = 'NewDrugsBrowser.html'
    # Get context for hmtl usage #
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        
        # Get data - server side - Queries #
        Drug_browser_data = Drugs2024.objects.all().prefetch_related('target','target__family__parent__parent','target__family__parent__parent__parent','indication','indication__code','ligand','ligand__ligand_type','moa')
        # Possible addons: 'target__family__parent','target__family__parent__parent','target__family__parent__parent__parent','indication','ligand','moa'
        # initialize context list for pushing data to html #
        context_data_browser = list()

        #proteins = list(TissueExpression.objects.all().values_list('protein__entry_name').distinct())
        #drugs = Drugs.objects.all().prefetch_related('target', 'target__family__parent__parent__parent', 'publication', 'publication__web_link', 'publication__web_link__web_resource', 'publication__journal')
        #drugs_NHS_names = list(NHSPrescribings.objects.values_list('drugname__name', flat=True).distinct())
        for entry in Drug_browser_data:
            ## For the drug browser table ##
            # Protein id, uniprot, and receptor name
            Protein_id = str(entry.target.id)
            Protein_uniprot = str(entry.target.accession)
            Protein_name = str(entry.target.entry_name).split("_")[0].upper()
            Protein_receptor_name = str(entry.target.name)
            Protein_family = str(entry.target.family.parent)
            Protein_class = str(entry.target.family.parent.parent.parent)
            Drug_name = str(entry.ligand.name)
            Indication = str(entry.indication.name)
            Indication_ID = str(entry.indication.code.index)
            Clinical_drug_status = str(entry.drug_status)
            Clinical_max_phase = str(entry.indication_max_phase)
            Clinical_approval_year = str(entry.approval_year)
            Clinical_drug_type = str(entry.ligand.ligand_type.name)
            #Clinical_drug_type = "bob"
            Clinical_moa = str(entry.moa.name)
            ### Columns that needs to be included ###
            
            #One row with target indication pair
            #Novelty score
            #IDG
            jsondata_browser = {
                'Index_number': Protein_id,
                'Drug': Drug_name,
                'Indication': Indication,
                'Indication_ID': Indication_ID,
                'Protein_uniprot': Protein_uniprot,
                'Protein_name': Protein_name,
                'Protein_receptor': Protein_receptor_name,
                'Protein_class': Protein_class,
                'Protein_family': Protein_family,
                'Drug_status': Clinical_drug_status,
                'Indication_max_phase': Clinical_max_phase,
                'Approval_year': Clinical_approval_year,
                'Drug_type': Clinical_drug_type,
                'Moa': Clinical_moa
            }
            context_data_browser.append(jsondata_browser)
        context['drug_data'] = context_data_browser
        return context

class TargetSelectionTool(TemplateView):
    # Template using this class #
    template_name = 'TargetSelectionTool.html'
    # Get context for hmtl usage #
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        # Get data - server side - Queries #
        TissueExp = TissueExpression.objects.all().prefetch_related('protein','tissue')
        Target_drug_data = Drugs2024.objects.all().prefetch_related('target','ligand','indication','indication__code')
        target_ids = [drug.target.id for drug in Target_drug_data if drug.target]
        #Structures_data = Structure.objects.all().prefetch_related('state', 'protein_conformation__protein')
        Structure_data = Structure.objects.values('protein_conformation_id', 'state_id').annotate(state_count=Count('state_id')).order_by('protein_conformation_id')
        filtered_structures = Structure.objects.annotate(target_id=Case(When(protein_conformation_id__protein_id__accession__isnull=False,
                                                                        then=F('protein_conformation_id__protein_id__id')),
                                                                        When(protein_conformation_id__protein_id__parent__accession__isnull=False,
                                                                        then=F('protein_conformation_id__protein_id__parent__id')),
                                                                        default=Value(None, output_field=IntegerField()),
                                                                        output_field=IntegerField())).filter(
                                                                            target_id__in=target_ids).values(
                                                                                'target_id').annotate(
                                                                                    active_count=Count(Case(When(state_id=1, then=1), output_field=IntegerField())),
                                                                                    inactive_count=Count(Case(When(state_id=2, then=1), output_field=IntegerField())),
                                                                                    intermediate_count=Count(Case(When(state_id=3, then=1), output_field=IntegerField())),
                                                                                    other_count=Count(Case(When(~Q(state_id__in=[1, 2, 3]), then=1), output_field=IntegerField())))
        #'target__family__parent__parent','target__family__parent__parent__parent','indication','indication__code','ligand','moa'
        # Context lists for pushing data #
        context_target_selection = list()
        context_data_tissue = list()
        # Dicts to modulate the data from the server side #
        structure_dict = {}
        target_selection_dict = {}
        target_indication_dict = {}
        Tissue_expression_dict = {}
        index_dict = {}
        # Create dict for structure total and state (active, inactive and intermediate)
        for entry in filtered_structures:
            id = entry['protein_conformation_id']
            
            if id not in structure_dict:
                structure_dict[id] = {}
                structure_dict[id]['Active']
                structure_dict[id]['Active']
                structure_dict[id]['Active']
            
        # Create Target selection browser #
        for entry in Target_drug_data:
            # Ids and keys
            Protein_id = str(entry.target.id)
            Indication_id = str(entry.indication.code.index)
            Target_indication_pair = "{}___{}".format(Protein_id,Indication_id)
            # Static values #
            Protein_uniprot = str(entry.target.accession)
            Protein_name = str(entry.target.entry_name).split("_")[0].upper()
            Indication = str(entry.indication.name)
            Drug_name  = str(entry.ligand.name)
            Drug_status = str(entry.drug_status)
            Novelty_score = float(entry.novelty_score)
            # Dicts #
            if Target_indication_pair not in target_indication_dict:
                target_indication_dict[Target_indication_pair] = {}
                target_indication_dict[Target_indication_pair]['information'] = [Protein_id,Indication_id,Indication]
                target_indication_dict[Target_indication_pair]['Novelty_score'] = Novelty_score
                #target_indication_dict[Target_indication_pair]['Drugs_total'] = 1
                target_indication_dict[Target_indication_pair]['Drugs'] = []
                target_indication_dict[Target_indication_pair]['Drugs'].append(Drug_name)
                # Handle drugs
                target_indication_dict[Target_indication_pair]['Drug__status'] = {}
                target_indication_dict[Target_indication_pair]['Drug__status']['Active'] = 0
                target_indication_dict[Target_indication_pair]['Drug__status']['Approved'] = 0
                target_indication_dict[Target_indication_pair]['Drug__status']['Discontinued'] = 0
                target_indication_dict[Target_indication_pair]['Drug__status'][Drug_status] += 1
                # Handle structures #
                target_indication_dict[Target_indication_pair]['Structures'] = {}
                target_indication_dict[Target_indication_pair]['Structures']['Total'] = structure_dict[Protein_id]['Total']
                target_indication_dict[Target_indication_pair]['Structures']['Active'] = structure_dict[Protein_id]['Active']
                target_indication_dict[Target_indication_pair]['Structures']['Inactive'] = structure_dict[Protein_id]['Inactive']
                target_indication_dict[Target_indication_pair]['Structures']['Intermediate'] = structure_dict[Protein_id]['Intermediate']
            else:
                # Drugs #
                target_indication_dict[Target_indication_pair]['Drugs'].append(Drug_name)
                target_indication_dict[Target_indication_pair]['Drug__status'][Drug_status] += 1

            if Protein_id not in target_selection_dict:
                target_selection_dict[Protein_id] = [Protein_uniprot,Protein_name]
        for key in target_indication_dict:
            key_id = str(target_indication_dict[key]['information'][0])
            jsondata_TargetSelectionTool = {
                    'Index_number': key_id,
                    'Target_name': target_selection_dict[key_id][1],
                    'Target_uniprot': target_selection_dict[key_id][0],
                    'Indication_name': target_indication_dict[key]['information'][2],
                    'Indication_id': target_indication_dict[key]['information'][1],
                    'Novelty_score': target_indication_dict[key]['Novelty_score'],
                    'IDG': "Coming soon",
                    'Drugs_approved_names': target_indication_dict[key]['Drugs'],
                    'Drugs_total': int(len(target_indication_dict[key]['Drugs'])),
                    'Drugs_approved': int(target_indication_dict[key]['Drug__status']['Approved']),
                    'Drugs_in_trial': int(target_indication_dict[key]['Drug__status']['Active']),
                    'Drugs_discontinued': int(target_indication_dict[key]['Drug__status']['Discontinued']),
                    'Structure_total': int(structure_dict[key_id]['Total']),
                    'Structure_active': int(target_indication_dict[key]['Structures']['Active']),
                    'Structure_inactive': int(target_indication_dict[key]['Structures']['Inactive']),
                    'Structure_intermediate': int(target_indication_dict[key]['Structures']['Intermediate'])
            }
            context_target_selection.append(jsondata_TargetSelectionTool)
        context['Target_data'] = context_target_selection
        # Go through server side data and modulate into a dict #
        for entry in TissueExp:
            # string values for Tissue expression table #
            protein_id = entry.protein.entry_name
            value = entry.value
            Tissue_id = entry.tissue.name
            # Index key #
            index_key = entry.protein.id
            if protein_id not in index_dict:
                index_dict[str(protein_id)] = index_key
            # Expression value linked to protein / target #
            if protein_id not in Tissue_expression_dict:
                Tissue_expression_dict[str(protein_id)] = {}
                Tissue_expression_dict[str(protein_id)][str(Tissue_id)] = float(value)
            else:
                Tissue_expression_dict[str(protein_id)][str(Tissue_id)] = float(value)
        # Run through dict and assign the correct values into the context data #
        for key in Tissue_expression_dict:
            jsondata_tissue = {
                    'Index_number': index_dict[key],
                    'ProteinID': key,
                    'adipose_tissue': Tissue_expression_dict[key]['adipose_tissue'],
                    'adrenal_gland': Tissue_expression_dict[key]['adrenal_gland'],
                    'amygdala': Tissue_expression_dict[key]['amygdala'],
                    'appendix': Tissue_expression_dict[key]['appendix'],
                    'basal_ganglia': Tissue_expression_dict[key]['basal_ganglia'],
                    'bone_marrow': Tissue_expression_dict[key]['bone_marrow'],
                    'breast': Tissue_expression_dict[key]['breast'],
                    'cerebellum': Tissue_expression_dict[key]['cerebellum'],
                    'cerebral_cortex': Tissue_expression_dict[key]['cerebral_cortex'],
                    'cervix': Tissue_expression_dict[key]['cervix'],
                    'choroid_plexus': Tissue_expression_dict[key]['choroid_plexus'],
                    'colon': Tissue_expression_dict[key]['colon'],
                    'duodenum': Tissue_expression_dict[key]['duodenum'],
                    'endometrium': Tissue_expression_dict[key]['endometrium_1'], #Should be updated
                    'epididymis': Tissue_expression_dict[key]['epididymis'],
                    'esophagus': Tissue_expression_dict[key]['esophagus'],
                    'fallopian_tube': Tissue_expression_dict[key]['fallopian_tube'],
                    'gallbladder': Tissue_expression_dict[key]['gallbladder'],
                    'heart_muscle': Tissue_expression_dict[key]['heart_muscle'],
                    'hippocampal_formation': Tissue_expression_dict[key]['hippocampal_formation'],
                    'hypothalamus': Tissue_expression_dict[key]['hypothalamus'],
                    'kidney': Tissue_expression_dict[key]['kidney'],
                    'liver': Tissue_expression_dict[key]['liver'],
                    'lung': Tissue_expression_dict[key]['lung'],
                    'lymph_node': Tissue_expression_dict[key]['lymph_node'],
                    'midbrain': Tissue_expression_dict[key]['midbrain'],
                    'ovary': Tissue_expression_dict[key]['ovary'],
                    'pancreas': Tissue_expression_dict[key]['pancreas'],
                    'parathyroid_gland': Tissue_expression_dict[key]['parathyroid_gland'],
                    'pituitary_gland': Tissue_expression_dict[key]['pituitary_gland'],
                    'placenta': Tissue_expression_dict[key]['placenta'],
                    'prostate': Tissue_expression_dict[key]['prostate'],
                    'rectum': Tissue_expression_dict[key]['rectum'],
                    'retina': Tissue_expression_dict[key]['retina'],
                    'salivary_gland': Tissue_expression_dict[key]['salivary_gland'],
                    'seminal_vesicle': Tissue_expression_dict[key]['seminal_vesicle'],
                    'skeletal_muscle': Tissue_expression_dict[key]['skeletal_muscle'],
                    'skin': Tissue_expression_dict[key]['skin_1'],
                    'small_intestine': Tissue_expression_dict[key]['small_intestine'],
                    'smooth_muscle': Tissue_expression_dict[key]['smooth_muscle'],
                    'spinal_cord': Tissue_expression_dict[key]['spinal_cord'],
                    'spleen': Tissue_expression_dict[key]['spleen'],
                    'stomach': Tissue_expression_dict[key]['stomach_1'],
                    'testis': Tissue_expression_dict[key]['testis'],
                    'thymus': Tissue_expression_dict[key]['thymus'],
                    'thyroid_gland': Tissue_expression_dict[key]['thyroid_gland'],
                    'tongue': Tissue_expression_dict[key]['tongue'],
                    'tonsil': Tissue_expression_dict[key]['tonsil'],
                    'urinary_bladder': Tissue_expression_dict[key]['urinary_bladder'],
                    'vagina': Tissue_expression_dict[key]['vagina']
                }
            # Append context data into list #
            context_data_tissue.append(jsondata_tissue)
        # Create context data for tissue expression data # 
        context['Tissue_data'] = context_data_tissue
        # Lastly return context for html usage #
        return context


@cache_page(60 * 60 * 24 * 28)
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

@cache_page(60 * 60 * 24 * 28)
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

# @cache_page(60 * 60 * 24 * 28)
def indication_detail(request, code):

    code = code.upper()
    context = dict()
    #code = 'EFO_0003843'
    indication_data = Drugs2024.objects.filter(indication__code__index=code).prefetch_related('ligand',
                                                                                              'target',
                                                                                              'indication',
                                                                                              'indication__code')

    indication_name = Indication.objects.filter(code__index=code).values_list('name', flat=True).distinct()[0]

    sankey = {"nodes": [],
              "links": []}
    caches = {'indication':[],
              'ligands': [],
              'targets': [],
              'entries': []}

    node_counter = 0
    for record in indication_data:
        #assess the values for indication/ligand/protein
        indication_code = record.indication.name.capitalize()
        ligand_name = record.ligand.name.capitalize()
        ligand_id = record.ligand.id
        protein_name = record.target.name
        target_name = record.target.entry_name
        #check for each value if it exists and retrieve the source node value
        if indication_code not in caches['indication']:
            sankey['nodes'].append({"node": node_counter, "name": indication_code, "url":'https://www.ebi.ac.uk/ols4/ontologies/efo/classes?short_form='+code})
            node_counter += 1
            caches['indication'].append(indication_code)
        indi_node = next((item['node'] for item in sankey['nodes'] if item['name'] == indication_code), None)
        if [ligand_name, ligand_id] not in caches['ligands']:
            sankey['nodes'].append({"node": node_counter, "name": ligand_name, "url":'/ligand/'+str(ligand_id)+'/info'})
            node_counter += 1
            caches['ligands'].append([ligand_name, ligand_id])
        lig_node = next((item['node'] for item in sankey['nodes'] if item['name'] == ligand_name), None)
        if protein_name not in caches['targets']:
            sankey['nodes'].append({"node": node_counter, "name": protein_name, "url":'/protein/'+str(target_name)})
            node_counter += 1
            caches['targets'].append(protein_name)
            caches['entries'].append(target_name)
        prot_node = next((item['node'] for item in sankey['nodes'] if item['name'] == protein_name), None)
        #append connection between indication and ligand
        sankey['links'].append({"source":indi_node, "target":lig_node, "value":1, "ligtrace": ligand_name, "prottrace": None})
        #append connection between ligand and target
        sankey['links'].append({"source":lig_node, "target":prot_node, "value":1, "ligtrace": ligand_name, "prottrace": protein_name})

    #Fixing redundancy in sankey['links']
    unique_combinations = {}

    for d in sankey['links']:
        # Create a key based on source and target for identifying unique combinations
        key = (d['source'], d['target'])

        if key in unique_combinations:
            # If the combination exists, add the value to the existing entry
            unique_combinations[key]['value'] += d['value']
        else:
            # If it's a new combination, add it to the dictionary
            unique_combinations[key] = d

    # Convert the unique_combinations back to a list of dictionaries
    sankey['links'] = list(unique_combinations.values())
    total_points = len(caches['targets']) + len(caches['targets']) + 1;
    if len(caches['ligands']) > len(caches['targets']):
        context['nodes_nr'] = len(caches['ligands'])
    else:
        context['nodes_nr'] = len(caches['targets'])
    context['indication_code'] = code
    context['indication'] = indication_name.capitalize()
    context['sankey'] = json.dumps(sankey)
    context['points'] = total_points
    context['targets'] = list(caches['entries'])
    context['ligands'] = list(caches['ligands'])
    return render(request, 'indication_detail.html', context)

@cache_page(60 * 60 * 24 * 28)
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
