from django.db.models import Count, Avg, Min, Max
from collections import defaultdict
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.views.generic import TemplateView, View, DetailView, ListView

from common.models import ReleaseNotes
from common.phylogenetic_tree import PhylogeneticTreeGenerator
from common.selection import Selection, SelectionItem
from ligand.models import *
from protein.models import Protein, Species, ProteinFamily

from django.views.decorators.csrf import csrf_exempt

from copy import deepcopy
import itertools
import json

class LigandBrowser(TemplateView):
    """
    Per target summary of ligands.
    """
    template_name = 'ligand_browser.html'

    def get_context_data (self, **kwargs):

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
        ).annotate(num_records = Count('protein__entry_name')
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
            tmp[data_line.assay_type][data_line.standard_type].append(data_line.standard_value)
            tmp_count += 1

        #Flattened list of lists of dict values
        values = list(itertools.chain(*[itertools.chain(*tmp[x].values()) for x in tmp.keys()]))
        # TEMPORARY workaround for handling string values
        values = [float(item) for item in values if float(item) ]

        if len(values) > 0:
            ligand_data.append({
                'protein_name': protein_details.entry_name,
                'receptor_family': protein_details.family.parent.name,
                'ligand_type': protein_details.get_protein_family(),
                'class': protein_details.get_protein_class(),
                'record_count': tmp_count,
                'assay_type': ', '.join(tmp.keys()),
                #Flattened list of lists of dict keys:
                'value_types': ', '.join(itertools.chain(*(list(tmp[x]) for x in tmp.keys()))),
                'low_value': min(values),
                'average_value': sum(values)/len(values),
                'standard_units': ', '.join(list(set([x.standard_units for x in per_target_data])))
                })

    context = {'ligand_data': ligand_data, 'ligand':ligand_id}

    return render(request, 'ligand_details.html', context)


def TargetDetailsCompact(request, **kwargs):
    if 'slug' in kwargs:
        slug = kwargs['slug']
        if slug.count('_') == 0 :
            ps = AssayExperiment.objects.filter(protein__family__parent__parent__parent__slug=slug, ligand__properities__web_links__web_resource__slug = 'chembl_ligand')
        elif slug.count('_') == 1 and len(slug) == 7:
            ps = AssayExperiment.objects.filter(protein__family__parent__parent__slug=slug, ligand__properities__web_links__web_resource__slug = 'chembl_ligand')
        elif slug.count('_') == 2:
            ps = AssayExperiment.objects.filter(protein__family__parent__slug=slug, ligand__properities__web_links__web_resource__slug = 'chembl_ligand')
        #elif slug.count('_') == 3:
        elif slug.count('_') == 1 and len(slug) != 7:
            ps = AssayExperiment.objects.filter(protein__entry_name = slug, ligand__properities__web_links__web_resource__slug = 'chembl_ligand')

        if slug.count('_') == 1 and len(slug) == 7:
            f = ProteinFamily.objects.get(slug=slug)
        else:
            f = slug

        context = {
            'target':f
            }
    else:
        simple_selection = request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)
        if selection.targets != []:
            prot_ids = [x.item.id for x in selection.targets]
            ps = AssayExperiment.objects.filter(protein__in=prot_ids, ligand__properities__web_links__web_resource__slug = 'chembl_ligand')
            context = {
                'target': ', '.join([x.item.entry_name for x in selection.targets])
                }

    ps = ps.prefetch_related('protein','ligand__properities__web_links__web_resource','ligand__properities__vendors__vendor')
    d = {}
    for p in ps:
        if p.ligand not in d:
            d[p.ligand] = {}
        if p.protein not in d[p.ligand]:
             d[p.ligand][p.protein] = []
        d[p.ligand][p.protein].append(p)
    ligand_data = []
    for lig, records  in d.items():

        links = lig.properities.web_links.all()
        chembl_id = [x for x in links if x.web_resource.slug=='chembl_ligand'][0].index

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
                tmp["Bind" if data_line.assay_type == 'b' else "Funct"].append(data_line.pchembl_value)
                tmp_count += 1

            # TEMPORARY workaround for handling string values
            values = [float(item) for item in itertools.chain(*tmp.values()) if float(item) ]

            if len(values)>0:
                ligand_data.append({
                    'ligand_id': chembl_id,
                    'protein_name': protein_details.entry_name,
                    'species': protein_details.species.common_name,
                    'record_count': tmp_count,
                    'assay_type': ', '.join(tmp.keys()),
                    'purchasability': purchasability,
                    #Flattened list of lists of dict keys:
                    'low_value': min(values),
                    'average_value': sum(values)/len(values),
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
        if slug.count('_') == 0 :
            ps = AssayExperiment.objects.filter(protein__family__parent__parent__parent__slug=slug, ligand__properities__web_links__web_resource__slug = 'chembl_ligand')
        elif slug.count('_') == 1 and len(slug) == 7:
            ps = AssayExperiment.objects.filter(protein__family__parent__parent__slug=slug, ligand__properities__web_links__web_resource__slug = 'chembl_ligand')
        elif slug.count('_') == 2:
            ps = AssayExperiment.objects.filter(protein__family__parent__slug=slug, ligand__properities__web_links__web_resource__slug = 'chembl_ligand')
        #elif slug.count('_') == 3:
        elif slug.count('_') == 1 and len(slug) != 7:
            ps = AssayExperiment.objects.filter(protein__entry_name = slug, ligand__properities__web_links__web_resource__slug = 'chembl_ligand')

        if slug.count('_') == 1 and len(slug) == 7:
            f = ProteinFamily.objects.get(slug=slug)
        else:
            f = slug

        context = {
            'target':f
            }
    else:
        simple_selection = request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)
        if selection.targets != []:
            prot_ids = [x.item.id for x in selection.targets]
            ps = AssayExperiment.objects.filter(protein__in=prot_ids, ligand__properities__web_links__web_resource__slug = 'chembl_ligand')
            context = {
                'target': ', '.join([x.item.entry_name for x in selection.targets])
                }
    ps = ps.values('standard_type',
                'standard_relation',
                'standard_value',
                'assay_description',
                'assay_type',
                #'standard_units',
                'pchembl_value',
                'ligand__id',
                'ligand__properities_id',
                'ligand__properities__web_links__index',
                #'ligand__properities__vendors__vendor__name',
                'protein__species__common_name',
                'protein__entry_name',
                'ligand__properities__mw',
                'ligand__properities__logp',
                'ligand__properities__rotatable_bonds',
                'ligand__properities__smiles',
                'ligand__properities__hdon',
                'ligand__properities__hacc','protein'
                ).annotate(num_targets = Count('protein__id', distinct=True))
    for record in ps:
        record['purchasability'] = 'Yes' if len(LigandVendorLink.objects.filter(lp=record['ligand__properities_id']).exclude(vendor__name__in=['ZINC', 'ChEMBL', 'BindingDB', 'SureChEMBL', 'eMolecules', 'MolPort', 'PubChem'])) > 0 else 'No'

    context['proteins'] = ps

    return render(request, 'target_details.html', context)


def TargetPurchasabilityDetails(request, **kwargs):

    simple_selection = request.session.get('selection', False)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)
    if selection.targets != []:
        prot_ids = [x.item.id for x in selection.targets]
        ps = AssayExperiment.objects.filter(protein__in=prot_ids, ligand__properities__web_links__web_resource__slug = 'chembl_ligand')
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
                'ligand__properities__hacc','protein'
                ).annotate(num_targets = Count('protein__id', distinct=True))
    purchasable = []
    for record in ps:
        try:
            if record['ligand__properities__vendors__vendor__name'] in ['ZINC', 'ChEMBL', 'BindingDB', 'SureChEMBL', 'eMolecules', 'MolPort', 'PubChem', 'IUPHAR/BPS Guide to PHARMACOLOGY']:
                continue
            tmp = LigandVendorLink.objects.filter(vendor=record['ligand__properities__vendors__vendor__id'], lp=record['ligand__properities_id'])[0]
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

    def get_context_data (self, **kwargs):

        context = super().get_context_data(**kwargs)
        assays = AssayExperiment.objects.all().prefetch_related('protein__family__parent__parent__parent', 'protein__family')

        lig_count_dict = {}
        assays_lig = list(AssayExperiment.objects.all().values('protein__family__parent__parent__parent__name').annotate(c=Count('ligand',distinct=True)))
        for a in assays_lig:
            lig_count_dict[a['protein__family__parent__parent__parent__name']] = a['c']
        target_count_dict = {}
        assays_target = list(AssayExperiment.objects.all().values('protein__family__parent__parent__parent__name').annotate(c=Count('protein__family',distinct=True)))
        for a in assays_target:
            target_count_dict[a['protein__family__parent__parent__parent__name']] = a['c']

        prot_count_dict = {}
        proteins_count = list(Protein.objects.all().values('family__parent__parent__parent__name').annotate(c=Count('family',distinct=True)))
        for pf in proteins_count:
            prot_count_dict[pf['family__parent__parent__parent__name']] = pf['c']

        classes = ProteinFamily.objects.filter(slug__in=['001', '002', '003', '004', '005', '006', '007']) #ugly but fast
        proteins = Protein.objects.all().prefetch_related('family__parent__parent__parent')
        ligands = []

        for fam in classes:
            if fam.name in lig_count_dict:
                lig_count = lig_count_dict[fam.name]
                target_count = target_count_dict[fam.name]
            else:
                lig_count = 0
                target_count = 0
            prot_count = prot_count_dict[fam.name]
            ligands.append({
                'name': fam.name,
                'num_ligands': lig_count,
                'avg_num_ligands': lig_count/prot_count,
                'target_percentage': target_count/prot_count*100,
                'target_count': target_count
                })
        lig_count_total = sum([x['num_ligands'] for x in ligands])
        prot_count_total = Protein.objects.filter(family__slug__startswith='00').all().distinct('family').count()
        target_count_total = sum([x['target_count'] for x in ligands])
        lig_total = {
            'num_ligands': lig_count_total,
            'avg_num_ligands': lig_count_total/prot_count_total,
            'target_percentage': target_count_total/prot_count_total*100,
            'target_count': target_count_total
            }
        #Elegant solution but kinda slow (6s querries):
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
        class_a_data = tree.get_tree_data(ProteinFamily.objects.get(name='Class A (Rhodopsin)'))
        context['class_a_options'] = deepcopy(tree.d3_options)
        context['class_a_options']['anchor'] = 'class_a'
        context['class_a_options']['leaf_offset'] = 50
        context['class_a_options']['label_free'] = []
        context['class_a'] = json.dumps(class_a_data.get_nodes_dict('ligands'))
        class_b1_data = tree.get_tree_data(ProteinFamily.objects.get(name__startswith='Class B1 (Secretin)'))
        context['class_b1_options'] = deepcopy(tree.d3_options)
        context['class_b1_options']['anchor'] = 'class_b1'
        context['class_b1_options']['branch_trunc'] = 60
        context['class_b1_options']['label_free'] = [1,]
        context['class_b1'] = json.dumps(class_b1_data.get_nodes_dict('ligands'))
        class_b2_data = tree.get_tree_data(ProteinFamily.objects.get(name__startswith='Class B2 (Adhesion)'))
        context['class_b2_options'] = deepcopy(tree.d3_options)
        context['class_b2_options']['anchor'] = 'class_b2'
        context['class_b2_options']['label_free'] = [1,]
        context['class_b2'] = json.dumps(class_b2_data.get_nodes_dict('ligands'))
        class_c_data = tree.get_tree_data(ProteinFamily.objects.get(name__startswith='Class C (Glutamate)'))
        context['class_c_options'] = deepcopy(tree.d3_options)
        context['class_c_options']['anchor'] = 'class_c'
        context['class_c_options']['branch_trunc'] = 50
        context['class_c_options']['label_free'] = [1,]
        context['class_c'] = json.dumps(class_c_data.get_nodes_dict('ligands'))
        class_f_data = tree.get_tree_data(ProteinFamily.objects.get(name__startswith='Class F (Frizzled)'))
        context['class_f_options'] = deepcopy(tree.d3_options)
        context['class_f_options']['anchor'] = 'class_f'
        context['class_f_options']['label_free'] = [1,]
        context['class_f'] = json.dumps(class_f_data.get_nodes_dict('ligands'))
        class_t2_data = tree.get_tree_data(ProteinFamily.objects.get(name__startswith='Class T (Taste 2)'))
        context['class_t2_options'] = deepcopy(tree.d3_options)
        context['class_t2_options']['anchor'] = 'class_t2'
        context['class_t2_options']['label_free'] = [1,]
        context['class_t2'] = json.dumps(class_t2_data.get_nodes_dict('ligands'))

        return context

#Biased Ligands part

class ExperimentEntryView(DetailView):
    context_object_name = 'experiment'
    model = AnalyzedExperiment
    template_name = 'biased_experiment_data.html'

#Biased pathways part

class PathwayExperimentEntryView(DetailView):
    context_object_name = 'experiment'
    model = BiasedPathways
    template_name = 'biased_pathways_data.html'


@csrf_exempt
def test_link(request):
    request.session['ids'] = ''
    # try:
    request.session['ids']
    if request.POST.get('action') == 'post':
        request.session.modified = True
        data = request.POST.get('ids')
        data = filter(lambda char: char not in " \"?.!/;:[]", data)
        datum = "".join(data)
        request.session['ids'] = datum
        request.session.set_expiry(15)
        # print('datum',datum )

    return HttpResponse(request)
    # except OSError as exc:
    #     raise



class BiasVendorBrowser(TemplateView):

    template_name = 'biased_ligand_vendor.html'
    #@cache_page(50000)
    def get_context_data(self, **kwargs):
        # try:
        context = dict()
        datum = self.request.session.get('ids')

        self.request.session.modified = True
        rd = list()
        for i in datum.split(','):
            ligand = Ligand.objects.filter(id=int(i))
            ligand = ligand.get()
            links = LigandVendorLink.objects.filter(lp=ligand.properities_id).prefetch_related('lp','vendor')
            for x in links:
                if x.vendor.name not in ['ZINC', 'ChEMBL', 'BindingDB', 'SureChEMBL', 'eMolecules', 'MolPort', 'PubChem']:
                    temp = dict()
                    vendor = LigandVendors.objects.filter(id=x.vendor_id)
                    vendor = vendor.get()
                    temp['ligand'] = ligand
                    temp['url'] = x.url
                    temp['vendor_id'] = x.vendor_external_id
                    temp['vendor'] = vendor

                    rd.append(temp)
            context['data'] = rd
        del self.request.session['ids']
        return context
        # except:
        #     raise

'''
Bias browser between families
access data from db, fill empty fields with empty parse_children
'''
class BiasBrowser(TemplateView):

    template_name = 'bias_browser.html'
    #@cache_page(50000)
    def get_context_data(self, *args, **kwargs  ):

        content = AnalyzedExperiment.objects.filter(source='different_family').prefetch_related(
        'analyzed_data', 'ligand','ligand__reference_ligand','reference_ligand',
        'endogenous_ligand' ,'ligand__properities','receptor','receptor','receptor__family',
        'receptor__family__parent','receptor__family__parent__parent__parent',
        'receptor__family__parent__parent','receptor__family', 'receptor__species',
        'publication', 'publication__web_link', 'publication__web_link__web_resource',
        'publication__journal', 'ligand__ref_ligand_bias_analyzed',
        'analyzed_data__emax_ligand_reference')
        context = dict()
        prepare_data = self.process_data(content)

        keys = [k for k, v in prepare_data.items() if len(v['biasdata']) < 2]
        for x in keys:
            del prepare_data[x]

        self.multply_assay(prepare_data)
        context.update({'data': prepare_data})

        return context

    def process_data(self, content):
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
            temp['experiment_id'] = instance.id
            temp['publication'] = instance.publication
            temp['ligand'] = instance.ligand
            temp['source'] = instance.source
            temp['chembl'] = instance.chembl
            temp['endogenous_ligand'] = instance.endogenous_ligand
            temp['vendor_quantity'] = instance.vendor_quantity
            temp['publication_quantity'] = instance.article_quantity
            temp['lab_quantity'] = instance.labs_quantity
            temp['reference_ligand'] = instance.reference_ligand
            temp['primary'] = instance.primary.replace(' family,','')
            temp['secondary'] = instance.secondary.replace(' family,','')

            if instance.receptor:
                temp['class'] = instance.receptor.family.parent.parent.parent.name.replace('Class','').strip()
                temp['receptor'] = instance.receptor
                temp['uniprot'] = instance.receptor.entry_short
                temp['IUPHAR'] = instance.receptor.name.split(' ', 1)[0].split('-adrenoceptor', 1)[0].strip()
            else:
                temp['receptor'] = 'Error appeared'
            temp['biasdata'] = list()
            increment_assay = 0
            for entry in instance.analyzed_data.all():
                if entry.order_no < 5:

                    temp_dict = dict()
                    temp_dict['emax_reference_ligand'] = entry.emax_ligand_reference
                    temp_dict['family'] = entry.family
                    temp_dict['show_family'] = entry.signalling_protein
                    temp_dict['signalling_protein'] = entry.signalling_protein
                    temp_dict['cell_line'] = entry.cell_line
                    temp_dict['assay_type'] = entry.assay_type
                    temp_dict['assay_measure'] = entry.assay_measure
                    temp_dict['assay_time_resolved'] = entry.assay_time_resolved
                    temp_dict['ligand_function'] = entry.ligand_function
                    temp_dict['quantitive_measure_type'] = entry.quantitive_measure_type
                    temp_dict['quantitive_activity'] = entry.quantitive_activity
                    temp_dict['quantitive_activity_initial'] = entry.quantitive_activity_initial
                    temp_dict['quantitive_unit'] = entry.quantitive_unit
                    temp_dict['qualitative_activity'] = entry.qualitative_activity
                    temp_dict['quantitive_efficacy'] = entry.quantitive_efficacy
                    temp_dict['efficacy_measure_type'] = entry.efficacy_measure_type
                    temp_dict['efficacy_unit'] = entry.efficacy_unit
                    temp_dict['order_no'] =  int(entry.order_no)
                    temp_dict['t_coefficient'] = entry.t_coefficient
                    if entry.t_value != None and entry.t_value !='None':
                        temp_dict['t_value'] = entry.t_value
                    else:
                        temp_dict['t_value'] = ''

                    if entry.t_factor != None and entry.t_factor !='None':
                        temp_dict['t_factor'] = entry.t_factor
                    else:
                        temp_dict['t_factor'] = ''

                    if entry.potency != None and entry.potency !='None':
                        temp_dict['potency'] =  entry.potency
                    else:
                        temp_dict['potency'] = ''

                    if entry.log_bias_factor != None and entry.log_bias_factor !='None':
                        temp_dict['log_bias_factor'] = entry.log_bias_factor
                    else:
                        temp_dict['log_bias_factor'] = ''

                    temp_dict['emax_ligand_reference'] = entry.emax_ligand_reference

                    temp['biasdata'].append(temp_dict)

                    doubles.append(temp_dict)
                    increment_assay+=1
                else:
                    continue

            rd[increment] = temp
            increment+=1


        return rd

    def multply_assay(self, data):

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
                          reverse=True)
            i[1]['biasdata'] = test

    '''
    End  of Bias Browser
    '''

class BiasBrowserGSubbtype(TemplateView):
    template_name = 'bias_browser_g.html'
    #@cache_page(50000)
    def get_context_data(self, *args, **kwargs  ):

        content = AnalyzedExperiment.objects.filter(source='same_family').prefetch_related(
    'analyzed_data', 'ligand','ligand__reference_ligand','reference_ligand',
    'endogenous_ligand' ,'ligand__properities','receptor','receptor__family__parent','receptor__family__parent__parent__parent',
    'receptor__family__parent__parent','receptor__species',
    'publication', 'publication__web_link', 'publication__web_link__web_resource',
    'publication__journal', 'ligand__ref_ligand_bias_analyzed',
    'analyzed_data__emax_ligand_reference')
        context = dict()
        prepare_data = self.process_data(content)

        keys = [k for k, v in prepare_data.items() if len(v['biasdata']) < 2]
        for x in keys:
            del prepare_data[x]

        self.multply_assay(prepare_data)
        context.update({'data': prepare_data})
        return context

    def process_data(self, content):
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
            # TODO: mutation residue
            temp['experiment_id'] = instance.id
            temp['publication'] = instance.publication
            temp['ligand'] = instance.ligand
            temp['source'] = instance.source
            temp['chembl'] = instance.chembl
            temp['endogenous_ligand'] = instance.endogenous_ligand
            temp['vendor_quantity'] = instance.vendor_quantity
            temp['publication_quantity'] = instance.article_quantity
            temp['lab_quantity'] = instance.labs_quantity
            temp['reference_ligand'] = instance.reference_ligand
            temp['primary'] =   instance.primary
            temp['secondary'] = instance.secondary
            if instance.receptor:
                temp['class'] = instance.receptor.family.parent.parent.parent.name.replace('Class','').strip()
                temp['receptor'] = instance.receptor
                temp['uniprot'] = instance.receptor.entry_short
                temp['IUPHAR'] = instance.receptor.name.split(' ', 1)[0].strip()
            else:
                temp['receptor'] = 'Error appeared'
            temp['biasdata'] = list()
            increment_assay = 0
            for entry in instance.analyzed_data.all():
                if entry.order_no < 5:

                    temp_dict = dict()
                    temp_dict['emax_reference_ligand'] = entry.emax_ligand_reference
                    temp_dict['family'] = entry.family
                    temp_dict['show_family'] = entry.signalling_protein
                    temp_dict['signalling_protein'] = entry.signalling_protein
                    temp_dict['cell_line'] = entry.cell_line
                    temp_dict['assay_type'] = entry.assay_type
                    temp_dict['assay_measure'] = entry.assay_measure
                    temp_dict['assay_time_resolved'] = entry.assay_time_resolved
                    temp_dict['ligand_function'] = entry.ligand_function
                    temp_dict['quantitive_measure_type'] = entry.quantitive_measure_type
                    temp_dict['quantitive_activity'] = entry.quantitive_activity
                    temp_dict['quantitive_activity_initial'] = entry.quantitive_activity_initial
                    temp_dict['quantitive_unit'] = entry.quantitive_unit
                    temp_dict['qualitative_activity'] = entry.qualitative_activity
                    temp_dict['quantitive_efficacy'] = entry.quantitive_efficacy
                    temp_dict['efficacy_measure_type'] = entry.efficacy_measure_type
                    temp_dict['efficacy_unit'] = entry.efficacy_unit
                    temp_dict['order_no'] =  int(entry.order_no)
                    temp_dict['t_coefficient'] = entry.t_coefficient
                    if entry.t_value != None and entry.t_value !='None':
                        temp_dict['t_value'] = entry.t_value
                    else:
                        temp_dict['t_value'] = ''

                    if entry.t_factor != None and entry.t_factor !='None':
                        temp_dict['t_factor'] = entry.t_factor
                    else:
                        temp_dict['t_factor'] = ''

                    if entry.potency != None and entry.potency !='None':
                        temp_dict['potency'] =  entry.potency
                    else:
                        temp_dict['potency'] = ''

                    if entry.log_bias_factor != None and entry.log_bias_factor !='None':
                        temp_dict['log_bias_factor'] = entry.log_bias_factor
                    else:
                        temp_dict['log_bias_factor'] = ''

                    temp_dict['emax_ligand_reference'] = entry.emax_ligand_reference

                    temp['biasdata'].append(temp_dict)

                    doubles.append(temp_dict)
                    increment_assay+=1
                else:
                    continue

            rd[increment] = temp
            increment+=1


        return rd

    def multply_assay(self, data):

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
                          reverse=True)
            i[1]['biasdata'] = test



'''
Bias browser between families
access data from db, fill empty fields with empty parse_children
'''
class BiasBrowserChembl(TemplateView):
    template_name = 'bias_browser_chembl.html'
    #@cache_page(50000)
    def get_context_data(self, *args, **kwargs  ):
        content = AnalyzedExperiment.objects.filter(source='chembl_data').prefetch_related(
            'analyzed_data', 'ligand','ligand__reference_ligand','reference_ligand',
            'endogenous_ligand' ,'ligand__properities','receptor','receptor__family',
            'receptor__family__parent','receptor__family__parent__parent__parent',
            'receptor__family__parent__parent','receptor__species',
            'publication', 'publication__web_link', 'publication__web_link__web_resource',
            'publication__journal', 'ligand__ref_ligand_bias_analyzed',
            'analyzed_data__emax_ligand_reference')
        context = dict()
        prepare_data = self.process_data(content)
        keys = [k for k, v in prepare_data.items() if len(v['biasdata']) < 2]
        for x in keys:
            del prepare_data[x]

        self.multply_assay(prepare_data)
        context.update({'data': prepare_data})
        return context

    def process_data(self, content):
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
            # TODO: mutation residue
            temp['experiment_id'] = instance.id
            temp['ligand'] = instance.ligand
            temp['source'] = instance.source
            temp['chembl'] = instance.chembl
            temp['endogenous_ligand'] = instance.endogenous_ligand
            temp['vendor_quantity'] = instance.vendor_quantity
            temp['primary'] =   instance.primary
            temp['secondary'] = instance.secondary

            if instance.receptor:
                temp['receptor'] = instance.receptor
            else:
                temp['receptor'] = 'Error appeared'
            temp['biasdata'] = list()
            increment_assay = 0
            for entry in instance.analyzed_data.all():
                if entry.order_no < 5:
                    temp_dict = dict()
                    temp_dict['family'] = entry.family
                    temp_dict['assay'] = entry.assay_type
                    temp_dict['assay_description'] = entry.assay_description
                    temp_dict['show_family'] = entry.signalling_protein
                    temp_dict['signalling_protein'] = entry.signalling_protein
                    temp_dict['quantitive_measure_type'] = entry.quantitive_measure_type
                    temp_dict['quantitive_activity'] = entry.quantitive_activity
                    temp_dict['quantitive_activity_initial'] = entry.quantitive_activity_initial
                    temp_dict['quantitive_unit'] = entry.quantitive_unit
                    temp_dict['qualitative_activity'] = entry.qualitative_activity
                    temp_dict['order_no'] =  int(entry.order_no)

                    if entry.potency != None and entry.potency !='None':
                        temp_dict['potency'] =  entry.potency
                    else:
                        temp_dict['potency'] = ''

                    temp['biasdata'].append(temp_dict)
                    doubles.append(temp_dict)
                    increment_assay+=1
                else:
                    continue
            rd[increment] = temp
            increment+=1
        return rd

    def multply_assay(self, data):

        for i in data.items():

            lenght = len(i[1]['biasdata'])
            for key in range(lenght,5):
                temp_dict = dict()
                temp_dict['pathway'] = ''
                temp_dict['order_no'] = lenght
                i[1]['biasdata'].append(temp_dict)
                lenght+=1
            test = sorted(i[1]['biasdata'], key=lambda x: x['order_no'],
                          reverse=True)
            i[1]['biasdata'] = test




    '''
    End  of Bias Browser
    '''

class BiasPathways(TemplateView):

    template_name = 'bias_browser_pathways.html'
    #@cache_page(50000)
    def get_context_data(self, *args, **kwargs  ):
        content = BiasedPathways.objects.all().prefetch_related(
        'biased_pathway', 'ligand','receptor','receptor','receptor__family',
        'receptor__family__parent','receptor__family__parent__parent__parent',
        'receptor__family__parent__parent','receptor__species',
        'publication', 'publication__web_link', 'publication__web_link__web_resource',
        'publication__journal')
        context = dict()
        prepare_data = self.process_data(content)
        context.update({'data': prepare_data})

        return context

    def process_data(self, content):
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
            # TODO: mutation residue
            temp['experiment_id'] = instance.id
            temp['publication'] = instance.publication
            temp['ligand'] = instance.ligand
            temp['rece'] = instance.chembl
            temp['chembl'] = instance.chembl
            temp['relevance'] = instance.relevance
            temp['signalling_protein'] = instance.signalling_protein

            if instance.receptor:
                temp['receptor'] = instance.receptor
                temp['uniprot'] = instance.receptor.entry_short
                temp['IUPHAR'] = instance.receptor.name.split(' ', 1)[0].strip()
            else:
                temp['receptor'] = 'Error appeared'
            # at the moment, there is only 1 pathways for every biased_pathway
            # change if more pathways added (logic changed)
            for entry in instance.biased_pathway.all():
                temp['pathway_outcome_high'] = entry.pathway_outcome_high
                temp['pathway_outcome_summary'] = entry.pathway_outcome_summary
                temp['pathway_outcome_detail'] = entry.pathway_outcome_detail
                temp['experiment_pathway_distinction'] = entry.experiment_pathway_distinction
                temp['experiment_system'] = entry.experiment_system
                temp['experiment_outcome_method'] = entry.experiment_outcome_method

            rd[increment] = temp
            increment+=1


        return rd




    '''
    End  of Bias Browser
    '''
