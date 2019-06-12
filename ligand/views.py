from django.db.models import Count, Avg, Min, Max
from collections import defaultdict
from django.shortcuts import render
from django.http import HttpResponse
from django.views.generic import TemplateView, View, ListView, DetailView

from common.models import ReleaseNotes
from common.phylogenetic_tree import PhylogeneticTreeGenerator
from common.selection import Selection, SelectionItem
from ligand.models import Ligand, BiasedExperiment, AssayExperiment, LigandProperities, LigandVendorLink
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


def process_data(content):
    '''
    Merge BiasedExperiment with its children
    and pass it back to loop through dublicates
    '''
    rd = []
    for instance in enumerate(content):
        temp_obj = []
        fin_obj = {}
        fin_obj['main'] = (instance[1])
        for entry in instance[1].experiment_data.all():
            temp_obj.append(entry)
        fin_obj['children'] = temp_obj
        rd.append(fin_obj)

    return rd


def change(rd):
    test = list()
    send = dict()
    increment = 0

    for j in rd:
        temp_dict = dict()
        temp = dict()
        doubles = []
        temp['publication'] = j['main'].publication
        temp['ligand'] = j['main'].ligand
        temp['receptor'] = j['main'].receptor
        if (j['children'][0].signalling_protein == 'β-arrestin' or
            j['children'][0].signalling_protein == 'β-arrestin-1 (non-visual arrestin-2)' or
                j['children'][0].signalling_protein == 'β-arrestin-2 (non-visual arrestin-3)'):
            temp_dict['pathway'] = 'β-arrestin'

        elif (j['children'][0].signalling_protein == 'gi/o-family' or
              j['children'][0].signalling_protein == 'gαi1' or
              j['children'][0].signalling_protein == 'gαi2' or
              j['children'][0].signalling_protein == 'gαi3' or
              j['children'][0].signalling_protein == 'gαo' or
              j['children'][0].signalling_protein == 'gαoA' or
              j['children'][0].signalling_protein == 'gαoB'):
            temp_dict['pathway'] = 'Gi/o-family'

        elif (j['children'][0].signalling_protein == 'gq-family' or
                j['children'][0].signalling_protein == 'gαq' or
                j['children'][0].signalling_protein == 'gαq11' or
                j['children'][0].signalling_protein == 'gαq14' or
                j['children'][0].signalling_protein == 'gαq14' or
                j['children'][0].signalling_protein == 'gαq16' or
                j['children'][0].signalling_protein == 'gαq14 (gαq16)'):
            temp_dict['pathway'] = 'Gq-family'

        elif (j['children'][0].signalling_protein == 'g12/13-family' or
              j['children'][0].signalling_protein == 'gα12' or
              j['children'][0].signalling_protein == 'gα13'):
            temp_dict['pathway'] = 'G12/13-family'

        elif (j['children'][0].signalling_protein == 'gs-family' or
              j['children'][0].signalling_protein == 'gαs' or
              j['children'][0].signalling_protein == 'gαolf'):
            temp_dict['pathway'] = 'Gs-family'
        else:
            temp_dict['pathway'] = 'No data'

        temp_dict['signalling_protein'] = j['children'][0].signalling_protein
        temp_dict['cell_line'] = j['children'][0].cell_line
        temp_dict['assay_type'] = j['children'][0].assay_type
        temp_dict['assay_measure_method'] = j['children'][0].assay_measure
        temp_dict['quantitive_activity'] = j['children'][0].quantitive_activity
        temp_dict['qualitative_activity'] = j['children'][0].qualitative_activity
        temp_dict['quantitive_efficacy'] = j['children'][0].quantitive_efficacy
        temp_dict['bias'] = 'Not yet'
        temp_dict['bias_reference'] = j['children'][0].bias_reference
        temp_dict['bias_value'] = j['children'][0].bias_value
        temp_dict['bias_value_initial'] = j['children'][0].bias_value_initial
        temp_dict['reference_ligand'] = j['children'][0].bias_ligand_reference
        temp_dict['ligand_function'] = j['children'][0].ligand_function

        doubles.append(temp_dict)
        temp['assay'] = doubles

        send[increment] = temp

        increment = increment + 1

    #    print('---temp---',temp,'\n')
    #    print('---send---',send[increment],'\n')

    return send


def process_references(dictionary):
    references = list()
    for item in dictionary.items():
        if item[1]['assay'][0]['bias_reference'].lower() == 'yes':
            references.append(item[1])
            #del item

    for item in dictionary.items():
        for i in references:
            if (i['publication'] == item[1]['publication'] and
                    i['receptor'] == item[1]['receptor']):
                if (item[1]['assay'][0]['signalling_protein'] == i['assay'][0]['signalling_protein'] and
                    item[1]['assay'][0]['assay_type'] == i['assay'][0]['assay_type'] and
                    item[1]['assay'][0]['cell_line'] == i['assay'][0]['cell_line'] and
                        item[1]['assay'][0]['assay_measure_method'] == i['assay'][0]['assay_measure_method']):
                    item[1]['reference'] = i['assay']

                else:
                    continue


def process_dublicates(dictionary):
    '''
    Recieve data from "process_data"
    search for objects with same publication ligand receptor
    Create new object of biased_experiment and all children
    '''
    doubles = []
    result = []
    context = {}
    send = {}
    for j in dictionary.items():
        name = str(j[1]['publication']) + \
            '/' + str(j[1]['ligand']) + '/' + str(j[1]['receptor'])
        temp_obj = list()
        reference_obj = list()
        if(name in context):
            temp_obj = context[name]['sub']
            reference_obj = context[name]['reference']

            for i in j[1]['assay']:
                temp_obj = temp_obj + j[1]['assay']
            for i in j[1]['reference']:
                reference_obj = reference_obj + j[1]['reference']
        else:
            for entry in j[1]['assay']:
                temp_obj.append(entry)

            for entry in j[1]['reference']:
                reference_obj.append(entry)

        context[name] = j[1]
        context[name].pop('assay')
        context[name]['sub'] = temp_obj
        context[name]['reference'] = reference_obj
        context[name]['experiment'] = parse_children(context[name])

        send[name] = parse_children(context[name])
        send[name]['reference'] = reference_obj
        doubles = send[name]['assay']
        doubles = [dict(t) for t in {tuple(d.items()) for d in doubles}]
        # potency sorting
        doubles = sorted(
            doubles, key=lambda k: k['quantitive_activity'] if k['quantitive_activity'] else 0,  reverse=True)

        send[name]['assay'] = doubles

    # assay_five(send)
    group_family(send)
    for i in send.items():
        test = dict()
        i[1].pop('assay')
        # family sorting
        test = sorted(i[1]['family'].items(), key=lambda x: x[1]['quantitive_activity']
                      if x[1]['quantitive_activity'] else 0,  reverse=False)
        i[1]['biasdata'] = test
        #verify_calc(i[1]['biasdata'], i[1]['reference'])
        calc_potency(i[1]['biasdata'])
        calc_bias_factor(i[1]['biasdata'], i[1]['reference'])
        print('---send---', i,'\n')
        for j in i[1]['family'].items():
            if(j[1]['quantitive_activity'] is not None and
            j[1]['quantitive_efficacy'] is not None):
                j[1]['quantitive_activity'] = round(math.log10(j[1]['quantitive_activity']),1)
                j[1]['quantitive_efficacy'] = int(j[1]['quantitive_efficacy'])
    assay_five(send)
    return send

def verify_calc(biasdata, reference):
    for item in enumerate(biasdata):
        print('---biasdata---', item,'\n')
    for item in enumerate(reference):
        print('---reference---', item,'\n')

def assay_five(send):
    '''Used to fill empty tempalte spaces '''
    for i in send.items():
        temp_dict = dict()
        bias_dict = 0

        if len(i[1]['biasdata']) < 5:
            temp_dict['pathway'] = ''
            temp_dict['bias'] = ''
            temp_dict['cell_line'] = ''
            temp_dict['assay_type'] = ''
            temp_dict['log_bias_factor'] = ''
            temp_dict['t_factor'] = ''
            temp_dict['ligand_function'] = ''

            for j in range(len(i[1]['biasdata']), 5):
                test = ('No_data', temp_dict)

                i[1]['biasdata'].append(test)


def parse_children(context):

    temp = dict()
    test = list()
    assay_list = list()

    for k in context.items():

        if k[0] == 'ligand':
            temp['ligand'] = k[1]
        elif k[0] == 'receptor':
            temp['receptor'] = k[1]
        elif k[0] == 'publication':
            temp['publication'] = k[1]

        elif k[0] == 'sub':
            counter = 0
            for i in k[1]:
                temp_dict = dict()
                # add logic of  families

                if (i['signalling_protein'] == 'β-arrestin' or
                    i['signalling_protein'] == 'β-arrestin-1 (non-visual arrestin-2)' or
                        i['signalling_protein'] == 'β-arrestin-2 (non-visual arrestin-3)'):
                    temp_dict['pathway'] = 'β-arrestin'

                elif (i['signalling_protein'] == 'gi/o-family' or
                      i['signalling_protein'] == 'gαi1' or
                      i['signalling_protein'] == 'gαi2' or
                      i['signalling_protein'] == 'gαi3' or
                      i['signalling_protein'] == 'gαo' or
                      i['signalling_protein'] == 'gαoA' or
                      i['signalling_protein'] == 'gαoB'):
                    temp_dict['pathway'] = 'Gi/o-family'

                elif (i['signalling_protein'] == 'gq-family' or
                        i['signalling_protein'] == 'gαq' or
                        i['signalling_protein'] == 'gαq11' or
                        i['signalling_protein'] == 'gαq14' or
                        i['signalling_protein'] == 'gαq14' or
                        i['signalling_protein'] == 'gαq16' or
                        i['signalling_protein'] == 'gαq14 (gαq16)'):
                    temp_dict['pathway'] = 'Gq-family'

                elif (i['signalling_protein'] == 'g12/13-family' or
                      i['signalling_protein'] == 'gα12' or
                      i['signalling_protein'] == 'gα13'):
                    temp_dict['pathway'] = 'G12/13-family'

                elif (i['signalling_protein'] == 'gs-family' or
                      i['signalling_protein'] == 'gαs' or
                      i['signalling_protein'] == 'gαolf'):
                    temp_dict['pathway'] = 'Gs-family'
                else:
                    temp_dict['pathway'] = 'No data'

                temp_dict['signalling_protein'] = i['signalling_protein']
                temp_dict['cell_line'] = i['cell_line']
                temp_dict['assay_type'] = i['assay_type']
                temp_dict['quantitive_activity'] = i['quantitive_activity']
                temp_dict['qualitative_activity'] = i['qualitative_activity']
                temp_dict['quantitive_efficacy'] = i['quantitive_efficacy']
                temp_dict['bias'] = 'Not yet'
                temp_dict['t_value'] = i['bias_value']
                temp_dict['t_coefficient'] = i['bias_value_initial']
                temp_dict['ligand_function'] = i['ligand_function']
                assay_list.append(temp_dict)
            temp['assay'] = assay_list

    context = temp
    return context
    for i in enumerate(biasdata):
        print('---biasdata---', biasdata, '\n')


def group_family(send):

    for item in send.items():

        context = dict()
        for i in item[1]['assay']:

            if i['pathway'] == 'β-arrestin':
                context['br'] = i
                #item[1]['family']= context
            elif i['pathway'] == 'Gi/o-family':
                context['gi'] = i
                #item[1]['gi']= i
            elif i['pathway'] == 'Gq-family':
                context['gq'] = i
                #item[1]['gq']= i
            elif i['pathway'] == 'G12/13-family':
                context['g12'] = i
                #item[1]['g12']= i
            elif i['pathway'] == 'Gs-family':
                context['gs'] = i
                #item[1]['gs']= i
            else:
                continue
            item[1]['family'] = context


def calc_potency(biasdata):
    most_potent = dict()
    for i in enumerate(biasdata):
        if i[0] == 0:
            most_potent = i[1][1]
        else:
            if i[1][1]['quantitive_activity'] is not None and most_potent['quantitive_activity'] is not None:
                i[1][1]['bias'] = round(
                    i[1][1]['quantitive_activity'] / most_potent['quantitive_activity'], 1)
            if i[1][1]['t_value'] is not None and most_potent['t_value'] is not None:
                i[1][1]['t_factor'] = round(
                    i[1][1]['t_value'] - most_potent['t_value'], 1)

        #print(ii,'---pre biasdata ---', data,'\n')


def calc_bias_factor(biasdata, reference):
    most_reference = dict()
    most_potent = dict()
    for i in enumerate(biasdata):
        temp_reference = dict()
        for j in reference:
            if i[0] == 0:
                most_potent = i[1][1]
                if (i[1][1]['assay_type'] == j['assay_type'] and
                    i[1][1]['cell_line'] == j['cell_line'] and
                        i[1][1]['signalling_protein'] == j['signalling_protein']):
                    most_reference = j

            else:
                if (i[1][1]['assay_type'] == j['assay_type'] and
                    i[1][1]['cell_line'] == j['cell_line'] and
                        i[1][1]['signalling_protein'] == j['signalling_protein']):
                    temp_reference = j
                if (i[1][1]['quantitive_activity'] is not None and most_potent['quantitive_activity'] is not None and
                    i[1][1]['quantitive_efficacy'] is not None and most_potent['quantitive_efficacy'] is not None and
                    j['quantitive_efficacy'] is not None and j['quantitive_efficacy'] is not None):
                    i[1][1]['log_bias_factor'] = round(((most_potent['quantitive_efficacy'] / most_potent['quantitive_activity']) - (most_reference['quantitive_efficacy'] / most_reference['quantitive_activity'])) -
                                                       ((i[1][1]['quantitive_efficacy'] / i[1][1]['quantitive_activity']) - (j['quantitive_efficacy'] / j['quantitive_activity'])), 1)


def bias_list(request):
    context = {}
    content = BiasedExperiment.objects.all().prefetch_related(
        'experiment_data', 'ligand', 'receptor', 'publication', 'experiment_data__bias_ligand_reference')
    # merge children
    pre_data = process_data(content)
    # transform to combined dictionary
    combined = change(pre_data)
    # combine references
    process_references(combined)

    # prepare dataset
    # process_dublicates(combined)

    #print("pre_data", pre_data)
    context.update({'data': process_dublicates(combined)})
    #print("context -----", context, '\n')
    return render(request, 'bias_list.html', context)
