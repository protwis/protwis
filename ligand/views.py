from django.db.models import Count, Avg, Min, Max
from collections import defaultdict
from django.shortcuts import render
from django.http import HttpResponse
from django.views.generic import TemplateView, View, ListView, DetailView

from common.models import ReleaseNotes
from common.phylogenetic_tree import PhylogeneticTreeGenerator
from common.selection import Selection, SelectionItem
from ligand.models import Ligand, BiasedExperiment, AssayExperiment, LigandProperities, LigandVendorLink, AnalyzedExperiment, AnnotationAssay
from protein.models import Protein, Species, ProteinFamily

from copy import deepcopy
import itertools
import operator
import json
import math

from itertools import chain
from operator import attrgetter
from django.db.models import Count, Min, Q


class LigandBrowser(TemplateView):
    """
    Per target summary of ligands.
    """
    template_name = 'ligand_browser.html'

    def get_context_data(self, **kwargs):

        context = super(LigandBrowser, self).get_context_data(**kwargs)

        ligands = AssayExperiment.objects.values(
            'protein__entry_name',
            'protein__species__common_name',
            'protein__family__name',
            'protein__family__parent__name',
            'protein__family__parent__parent__name',
            'protein__family__parent__parent__parent__name',
            'protein__species__common_name'
        ).annotate(num_ligands=Count('ligand', distinct=True))

        context['ligands'] = ligands

        return context


def LigandDetails(request, ligand_id):
    """
    The details of a ligand record. Lists all the assay experiments for a given ligand.
    """
    ligand_records = AssayExperiment.objects.filter(
        ligand__properities__web_links__index=ligand_id
    ).order_by('protein__entry_name')

    record_count = ligand_records.values(
        'protein',
    ).annotate(num_records=Count('protein__entry_name')
               ).order_by('protein__entry_name')

    ligand_data = []

    for record in record_count:
        per_target_data = ligand_records.filter(protein=record['protein'])
        protein_details = Protein.objects.get(pk=record['protein'])

        """
        A dictionary of dictionaries with a list of values.
        Assay_type
        |
        ->  Standard_type [list of values]
        """
        tmp = defaultdict(lambda: defaultdict(list))
        tmp_count = 0
        for data_line in per_target_data:
            tmp[data_line.assay_type][data_line.standard_type].append(
                data_line.standard_value)
            tmp_count += 1

        # Flattened list of lists of dict values
        values = list(itertools.chain(
            *[itertools.chain(*tmp[x].values()) for x in tmp.keys()]))

        ligand_data.append({
            'protein_name': protein_details.entry_name,
            'receptor_family': protein_details.family.parent.name,
            'ligand_type': protein_details.get_protein_family(),
            'class': protein_details.get_protein_class(),
            'record_count': tmp_count,
            'assay_type': ', '.join(tmp.keys()),
            # Flattened list of lists of dict keys:
            'value_types': ', '.join(itertools.chain(*(list(tmp[x]) for x in tmp.keys()))),
            'low_value': min(values),
            'average_value': sum(values) / len(values),
            'standard_units': ', '.join(list(set([x.standard_units for x in per_target_data])))
        })

    context = {'ligand_data': ligand_data, 'ligand': ligand_id}

    return render(request, 'ligand_details.html', context)


def TargetDetailsCompact(request, **kwargs):
    if 'slug' in kwargs:
        slug = kwargs['slug']
        if slug.count('_') == 0:
            ps = AssayExperiment.objects.filter(
                protein__family__parent__parent__parent__slug=slug, ligand__properities__web_links__web_resource__slug='chembl_ligand')
        elif slug.count('_') == 1 and len(slug) == 7:
            ps = AssayExperiment.objects.filter(
                protein__family__parent__parent__slug=slug, ligand__properities__web_links__web_resource__slug='chembl_ligand')
        elif slug.count('_') == 2:
            ps = AssayExperiment.objects.filter(
                protein__family__parent__slug=slug, ligand__properities__web_links__web_resource__slug='chembl_ligand')
        # elif slug.count('_') == 3:
        elif slug.count('_') == 1 and len(slug) != 7:
            ps = AssayExperiment.objects.filter(
                protein__entry_name=slug, ligand__properities__web_links__web_resource__slug='chembl_ligand')

        if slug.count('_') == 1 and len(slug) == 7:
            f = ProteinFamily.objects.get(slug=slug)
        else:
            f = slug

        context = {
            'target': f
        }
    else:
        simple_selection = request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)
        if selection.targets != []:
            prot_ids = [x.item.id for x in selection.targets]
            ps = AssayExperiment.objects.filter(
                protein__in=prot_ids, ligand__properities__web_links__web_resource__slug='chembl_ligand')
            context = {
                'target': ', '.join([x.item.entry_name for x in selection.targets])
            }

    ps = ps.prefetch_related(
        'protein', 'ligand__properities__web_links__web_resource', 'ligand__properities__vendors__vendor')
    d = {}
    for p in ps:
        if p.ligand not in d:
            d[p.ligand] = {}
        if p.protein not in d[p.ligand]:
            d[p.ligand][p.protein] = []
        d[p.ligand][p.protein].append(p)
    ligand_data = []
    for lig, records in d.items():

        links = lig.properities.web_links.all()
        chembl_id = [x for x in links if x.web_resource.slug ==
                     'chembl_ligand'][0].index

        vendors = lig.properities.vendors.all()
        purchasability = 'No'
        for v in vendors:
            if v.vendor.name not in ['ZINC', 'ChEMBL', 'BindingDB', 'SureChEMBL', 'eMolecules', 'MolPort', 'PubChem']:
                purchasability = 'Yes'

        for record, vals in records.items():
            per_target_data = vals
            protein_details = record

            """
            A dictionary of dictionaries with a list of values.
            Assay_type
            |
            ->  Standard_type [list of values]
            """
            tmp = defaultdict(list)
            tmp_count = 0
            for data_line in per_target_data:
                tmp["Bind" if data_line.assay_type == 'b' else "Funct"].append(
                    data_line.pchembl_value)
                tmp_count += 1
            values = list(itertools.chain(*tmp.values()))
            ligand_data.append({
                'ligand_id': chembl_id,
                'protein_name': protein_details.entry_name,
                'species': protein_details.species.common_name,
                'record_count': tmp_count,
                'assay_type': ', '.join(tmp.keys()),
                'purchasability': purchasability,
                # Flattened list of lists of dict keys:
                'low_value': min(values),
                'average_value': sum(values) / len(values),
                'high_value': max(values),
                'standard_units': ', '.join(list(set([x.standard_units for x in per_target_data]))),
                'smiles': lig.properities.smiles,
                'mw': lig.properities.mw,
                'rotatable_bonds': lig.properities.rotatable_bonds,
                'hdon': lig.properities.hdon,
                'hacc': lig.properities.hacc,
                'logp': lig.properities.logp,
            })
    context['ligand_data'] = ligand_data

    return render(request, 'target_details_compact.html', context)


def TargetDetails(request, **kwargs):

    if 'slug' in kwargs:
        slug = kwargs['slug']
        if slug.count('_') == 0:
            ps = AssayExperiment.objects.filter(
                protein__family__parent__parent__parent__slug=slug, ligand__properities__web_links__web_resource__slug='chembl_ligand')
        elif slug.count('_') == 1 and len(slug) == 7:
            ps = AssayExperiment.objects.filter(
                protein__family__parent__parent__slug=slug, ligand__properities__web_links__web_resource__slug='chembl_ligand')
        elif slug.count('_') == 2:
            ps = AssayExperiment.objects.filter(
                protein__family__parent__slug=slug, ligand__properities__web_links__web_resource__slug='chembl_ligand')
        # elif slug.count('_') == 3:
        elif slug.count('_') == 1 and len(slug) != 7:
            ps = AssayExperiment.objects.filter(
                protein__entry_name=slug, ligand__properities__web_links__web_resource__slug='chembl_ligand')

        if slug.count('_') == 1 and len(slug) == 7:
            f = ProteinFamily.objects.get(slug=slug)
        else:
            f = slug

        context = {
            'target': f
        }
    else:
        simple_selection = request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)
        if selection.targets != []:
            prot_ids = [x.item.id for x in selection.targets]
            ps = AssayExperiment.objects.filter(
                protein__in=prot_ids, ligand__properities__web_links__web_resource__slug='chembl_ligand')
            context = {
                'target': ', '.join([x.item.entry_name for x in selection.targets])
            }
    ps = ps.values('standard_type',
                   'standard_relation',
                   'standard_value',
                   'assay_description',
                   'assay_type',
                   # 'standard_units',
                   'pchembl_value',
                   'ligand__id',
                   'ligand__properities_id',
                   'ligand__properities__web_links__index',
                   # 'ligand__properities__vendors__vendor__name',
                   'protein__species__common_name',
                   'protein__entry_name',
                   'ligand__properities__mw',
                   'ligand__properities__logp',
                   'ligand__properities__rotatable_bonds',
                   'ligand__properities__smiles',
                   'ligand__properities__hdon',
                   'ligand__properities__hacc', 'protein'
                   ).annotate(num_targets=Count('protein__id', distinct=True))
    for record in ps:
        record['purchasability'] = 'Yes' if len(LigandVendorLink.objects.filter(lp=record['ligand__properities_id']).exclude(
            vendor__name__in=['ZINC', 'ChEMBL', 'BindingDB', 'SureChEMBL', 'eMolecules', 'MolPort', 'PubChem'])) > 0 else 'No'

    context['proteins'] = ps

    return render(request, 'target_details.html', context)


def TargetPurchasabilityDetails(request, **kwargs):

    simple_selection = request.session.get('selection', False)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)
    if selection.targets != []:
        prot_ids = [x.item.id for x in selection.targets]
        ps = AssayExperiment.objects.filter(
            protein__in=prot_ids, ligand__properities__web_links__web_resource__slug='chembl_ligand')
        context = {
            'target': ', '.join([x.item.entry_name for x in selection.targets])
        }

    ps = ps.values('standard_type',
                   'standard_relation',
                   'standard_value',
                   'assay_description',
                   'assay_type',
                   'standard_units',
                   'pchembl_value',
                   'ligand__id',
                   'ligand__properities_id',
                   'ligand__properities__web_links__index',
                   'ligand__properities__vendors__vendor__id',
                   'ligand__properities__vendors__vendor__name',
                   'protein__species__common_name',
                   'protein__entry_name',
                   'ligand__properities__mw',
                   'ligand__properities__logp',
                   'ligand__properities__rotatable_bonds',
                   'ligand__properities__smiles',
                   'ligand__properities__hdon',
                   'ligand__properities__hacc', 'protein'
                   ).annotate(num_targets=Count('protein__id', distinct=True))
    purchasable = []
    for record in ps:
        try:
            if record['ligand__properities__vendors__vendor__name'] in ['ZINC', 'ChEMBL', 'BindingDB', 'SureChEMBL', 'eMolecules', 'MolPort', 'PubChem', 'IUPHAR/BPS Guide to PHARMACOLOGY']:
                continue
            tmp = LigandVendorLink.objects.filter(
                vendor=record['ligand__properities__vendors__vendor__id'], lp=record['ligand__properities_id'])[0]
            record['vendor_id'] = tmp.vendor_external_id
            record['vendor_link'] = tmp.url
            purchasable.append(record)
        except:
            continue

    context['proteins'] = purchasable

    return render(request, 'target_purchasability_details.html', context)


class LigandStatistics(TemplateView):
    """
    Per class statistics of known ligands.
    """

    template_name = 'ligand_statistics.html'

    def get_context_data(self, **kwargs):

        context = super().get_context_data(**kwargs)
        assays = AssayExperiment.objects.all().prefetch_related(
            'protein__family__parent__parent__parent', 'protein__family')
        classes = ProteinFamily.objects.filter(
            slug__in=['001', '002', '003', '004', '005', '006'])  # ugly but fast
        proteins = Protein.objects.all().prefetch_related(
            'family__parent__parent__parent')
        ligands = []

        for fam in classes:
            lig_count = len(assays.filter(
                protein__family__parent__parent__parent=fam).distinct('ligand'))
            prot_count = len(proteins.filter(
                family__parent__parent__parent=fam).distinct('family'))
            target_count = len(assays.filter(
                protein__family__parent__parent__parent=fam).distinct('protein__family'))
            ligands.append({
                'name': fam.name,
                'num_ligands': lig_count,
                'avg_num_ligands': lig_count / prot_count,
                'target_percentage': target_count / prot_count * 100,
                'target_count': target_count
            })
        lig_count_total = sum([x['num_ligands'] for x in ligands])
        prot_count_total = len(proteins.distinct('family'))
        target_count_total = sum([x['target_count'] for x in ligands])
        lig_total = {
            'num_ligands': lig_count_total,
            'avg_num_ligands': lig_count_total / prot_count_total,
            'target_percentage': target_count_total / prot_count_total * 100,
            'target_count': target_count_total
        }
        # Elegant solution but kinda slow (6s querries):
        """
        ligands = AssayExperiment.objects.values(
            'protein__family__parent__parent__parent__name',
            'protein__family__parent__parent__parent',
            ).annotate(num_ligands=Count('ligand', distinct=True))
        for prot_class in ligands:
            class_subset = AssayExperiment.objects.filter(
                id=prot_class['protein__family__parent__parent__parent']).values(
                    'protein').annotate(
                        avg_num_ligands=Avg('ligand', distinct=True),
                        p_count=Count('protein')
                        )
            prot_class['avg_num_ligands']=class_subset[0]['avg_num_ligands']
            prot_class['p_count']=class_subset[0]['p_count']

        """
        context['ligands_total'] = lig_total
        context['ligands_by_class'] = ligands

        context['release_notes'] = ReleaseNotes.objects.all()[0]

        tree = PhylogeneticTreeGenerator()
        class_a_data = tree.get_tree_data(
            ProteinFamily.objects.get(name='Class A (Rhodopsin)'))
        context['class_a_options'] = deepcopy(tree.d3_options)
        context['class_a_options']['anchor'] = 'class_a'
        context['class_a_options']['leaf_offset'] = 50
        context['class_a_options']['label_free'] = []
        context['class_a'] = json.dumps(class_a_data.get_nodes_dict('ligands'))
        class_b1_data = tree.get_tree_data(
            ProteinFamily.objects.get(name__startswith='Class B1 (Secretin)'))
        context['class_b1_options'] = deepcopy(tree.d3_options)
        context['class_b1_options']['anchor'] = 'class_b1'
        context['class_b1_options']['branch_trunc'] = 60
        context['class_b1_options']['label_free'] = [1, ]
        context['class_b1'] = json.dumps(
            class_b1_data.get_nodes_dict('ligands'))
        class_b2_data = tree.get_tree_data(
            ProteinFamily.objects.get(name__startswith='Class B2 (Adhesion)'))
        context['class_b2_options'] = deepcopy(tree.d3_options)
        context['class_b2_options']['anchor'] = 'class_b2'
        context['class_b2_options']['label_free'] = [1, ]
        context['class_b2'] = json.dumps(
            class_b2_data.get_nodes_dict('ligands'))
        class_c_data = tree.get_tree_data(
            ProteinFamily.objects.get(name__startswith='Class C (Glutamate)'))
        context['class_c_options'] = deepcopy(tree.d3_options)
        context['class_c_options']['anchor'] = 'class_c'
        context['class_c_options']['branch_trunc'] = 50
        context['class_c_options']['label_free'] = [1, ]
        context['class_c'] = json.dumps(class_c_data.get_nodes_dict('ligands'))
        class_f_data = tree.get_tree_data(
            ProteinFamily.objects.get(name__startswith='Class F (Frizzled)'))
        context['class_f_options'] = deepcopy(tree.d3_options)
        context['class_f_options']['anchor'] = 'class_f'
        context['class_f_options']['label_free'] = [1, ]
        context['class_f'] = json.dumps(class_f_data.get_nodes_dict('ligands'))
        class_t2_data = tree.get_tree_data(
            ProteinFamily.objects.get(name='Taste 2'))
        context['class_t2_options'] = deepcopy(tree.d3_options)
        context['class_t2_options']['anchor'] = 'class_t2'
        context['class_t2_options']['label_free'] = [1, ]
        context['class_t2'] = json.dumps(
            class_t2_data.get_nodes_dict('ligands'))

        return context


class ExperimentEntryView(DetailView):
    context_object_name = 'experiment'
    model = BiasedExperiment
    template_name = 'biased_experiment_data.html'


def output_bias(request):
    couplings = BiasedExperiment.objects.prefetch_related('experiment_data', 'ligand', 'publication',
                                                          'receptor', 'mutation__protein', 'residue'
                                                          )

    context = {
        'data': couplings,

    }
    return render(request, 'biased_experiment.html', context)

'''
Bias browser between families
access data from db, fill empty fields with empty parse_children
'''

def bias_browser(request):
    context = dict()
    content = AnalyzedExperiment.objects.all().prefetch_related(
        'analyzed_data', 'ligand', 'receptor','receptor__family', 'receptor__species','publication', 'publication__web_link', 'publication__web_link__web_resource', 'publication__journal', 'analyzed_data__emax_ligand_reference')

    prepare_data = process_data1(content)

    keys = [k for k, v in prepare_data.items() if len(v['biasdata']) < 1]
    for x in keys:
        del prepare_data[x]

    multply_assay(prepare_data)
    context.update({'data': prepare_data})

    return render(request, 'bias_browser.html', context)


def process_data1(content):

    '''
    Merge BiasedExperiment with its children
    and pass it back to loop through dublicates
    '''
    rd = dict()
    increment = 0
    for instance in content:
        fin_obj = {}
        fin_obj['main'] = instance
        temp = dict()
        doubles = []

        temp['publication'] = instance.publication
        temp['ligand'] = instance.ligand
        temp['receptor'] = instance.receptor
        temp['biasdata'] = list()

        for entry in instance.analyzed_data.all():


            temp['reference_ligand'] = entry.emax_ligand_reference
            temp_dict = dict()
            temp_dict['family'] = entry.family
            temp_dict['signalling_protein'] = entry.signalling_protein
            temp_dict['cell_line'] = entry.cell_line
            temp_dict['assay_type'] = entry.assay_type
            temp_dict['assay_measure'] = entry.assay_measure
            temp_dict['assay_time_resolved'] = entry.assay_time_resolved
            temp_dict['ligand_function'] = entry.ligand_function
            temp_dict['quantitive_measure_type'] = entry.quantitive_measure_type
            temp_dict['quantitive_activity'] = entry.quantitive_activity
            temp_dict['quantitive_unit'] = entry.quantitive_unit
            temp_dict['qualitative_activity'] = entry.qualitative_activity
            temp_dict['quantitive_efficacy'] = entry.quantitive_efficacy
            temp_dict['efficacy_measure_type'] = entry.efficacy_measure_type
            temp_dict['efficacy_unit'] = entry.efficacy_unit
            temp_dict['potency'] =  entry.potency
            temp_dict['order_no'] =  int(entry.order_no)
            temp_dict['t_coefficient'] = entry.t_coefficient
            temp_dict['t_value'] = entry.t_value
            temp_dict['t_factor'] = entry.t_factor
            temp_dict['log_bias_factor'] = entry.log_bias_factor
            temp_dict['emax_ligand_reference'] = entry.emax_ligand_reference

            temp['biasdata'].append(temp_dict)

            doubles.append(temp_dict)

        #print('---data---', temp, '\n')

        rd[increment] = temp

        increment+=1
    return rd

def multply_assay(data):
    '''
    Used to fill empty template (5 in total, 5 families) spaces
    to hide columns im template
    '''
    for i in data.items():
        lenght = len(i[1]['biasdata'])
        for key in range(lenght,5):
            temp_dict = dict()
            temp_dict['pathway'] = ''
            temp_dict['bias'] = ''
            temp_dict['cell_line'] = ''
            temp_dict['assay_type'] = ''
            temp_dict['log_bias_factor'] = ''
            temp_dict['t_factor'] = ''
            temp_dict['ligand_function'] = ''
            temp_dict['order_no'] = lenght
            i[1]['biasdata'].append(temp_dict)
            lenght+=1
        test = sorted(i[1]['biasdata'], key=lambda x: x['order_no'],
                      reverse=False)
        i[1]['biasdata'] = test


'''
End  of Bias Browser
'''
