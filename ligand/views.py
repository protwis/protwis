import hashlib
import itertools
import json
import re

from random import SystemRandom
from copy import deepcopy
from collections import defaultdict, OrderedDict

from django.shortcuts import render, redirect
from django.http import HttpResponse
from django.views.generic import TemplateView, DetailView, ListView

from django.db.models import Count, Subquery, OuterRef
from django.views.decorators.csrf import csrf_exempt

from django.core.cache import cache

from common.views import AbsTargetSelectionTable, Alignment, AbsReferenceSelectionTable, getReferenceTable
from common.models import ReleaseNotes
from common.phylogenetic_tree import PhylogeneticTreeGenerator
from common.selection import Selection
from ligand.models import Ligand, LigandVendorLink,LigandVendors, AnalyzedExperiment, AnalyzedAssay, BiasedPathways, AssayExperiment
from protein.models import Protein, ProteinFamily, ProteinCouplings
from interaction.models import StructureLigandInteraction
from mutation.models import MutationExperiment

class LigandTargetSelection(AbsTargetSelectionTable):
    step = 1
    number_of_steps = 1
    filter_tableselect = False
    docs = 'sequences.html#structure-based-alignments'
    title = "SELECT RECEPTORS with ligand assays"
    description = 'Select receptors in the table (below) or browse the classification tree (right). You can select entire' \
        + ' families or individual receptors.\n\nOnce you have selected all your receptors, click the green button.'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Next',
            'onclick': "submitSelection('/ligand/browser');",
            'color': 'success',
        },
    }

class LigandBrowser(TemplateView):
    """
    Per target summary of ligands.
    """
    template_name = 'ligand_browser.html'

    def get_context_data(self, **kwargs):
        protein_list = list()
        try:
            simple_selection = self.request.session.get('selection', False)
            a = Alignment()
            # load data from selection into the alignment
            a.load_proteins_from_selection(simple_selection)
            for items in a.proteins:
                protein_list.append(items.protein)
        except:
            protein_list.append(1)
        context = super(LigandBrowser, self).get_context_data(**kwargs)
        assay_b = LigandReceptorStatistics.objects.filter(
            type='b',
            protein=OuterRef('protein_id'),
        )
        assay_f = LigandReceptorStatistics.objects.filter(
            type='f',
            protein=OuterRef('protein_id'),
        )
        ligands = AssayExperiment.objects.filter(protein__in=protein_list,).values(
            'protein',
            'protein__entry_name',
            'protein__species__common_name',
            'protein__family__name',
            'protein__family__parent__name',
            'protein__family__parent__parent__name',
            'protein__family__parent__parent__parent__name',
            'protein__species__common_name'
        ).annotate(num_ligands=Count('ligand', distinct=True),
            assay_f_count = Subquery(assay_f.values('value')),
            assay_b_count = Subquery(assay_b.values('value')),
            primary = Subquery(assay_f.values('primary')),
            primary1 = Subquery(assay_b.values('primary')),
            secondary = Subquery(assay_f.values('secondary')),
            secondary1 = Subquery(assay_b.values('secondary'))
        ).prefetch_related('protein', 'LigandReceptorStatistics')
        # import pdb; pdb.set_trace()
        context['ligands'] = ligands

        return context

    def fetch_receptor_transducers(self, receptor):
        primary = set()
        temp = str()
        temp1 = str()
        secondary = set()
        try:
            gprotein = ProteinCouplings.objects.filter(protein=receptor)
            for x in gprotein:
                if x.transduction and x.transduction == 'primary':
                    primary.add(x.g_protein.name)
                elif x.transduction and x.transduction == 'secondary':
                    secondary.add(x.g_protein.name)
            for i in primary:
                temp += str(i.replace(' family', '')) + str(', ')

            for i in secondary:
                temp1 += str(i.replace('family', '')) + str(', ')
            return temp, temp1
        except:
            self.logger.info('receptor not found error')
            return None, None

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
        Assay_type|->  Standard_type [list of values]
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
        # TEMPORARY workaround for handling string values
        values = [float(item) for item in values if float(item)]

        if len(values) > 0:
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
        elif slug.count('_') == 3:
            ps = AssayExperiment.objects.filter(
                protein__family__slug=slug, ligand__properities__web_links__web_resource__slug='chembl_ligand')
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
        if simple_selection == False or not simple_selection.targets :
            return redirect("ligand_browser")
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
    # if queryset is empty redirect to ligand browser
    if not ps:
        return redirect("ligand_browser")


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

            # TEMPORARY workaround for handling string values
            values = [float(item) for item in itertools.chain(
                *tmp.values()) if float(item)]

            if len(values) > 0:
                ligand_data.append({
                    'lig_id': lig.id,
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
        elif slug.count('_') == 3:
            ps = AssayExperiment.objects.filter(
                protein__family__slug=slug, ligand__properities__web_links__web_resource__slug='chembl_ligand')
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
        if simple_selection == False or not simple_selection.targets :
            return redirect("ligand_browser")
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
    # if queryset is empty redirect to ligand browser
    if not ps:
        return redirect("ligand_browser")

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

#Biased Effector Family Browser
class BiasTargetSelection(AbsReferenceSelectionTable):
    step = 1
    number_of_steps = 1
    filters = False
    filter_tableselect = False
    family_tree = False
    docs = 'sequences.html#structure-based-alignments'
    title = "SELECT RECEPTORS with ligands biased for a G protein or arrestin family (relative to an endogenous reference ligand)"
    description = 'Select a receptor in the table (below) or browse the classification tree (center).\nThen click a green button (right)'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Physiology-biased ligands<br>(endogenous agonist reference)',
            'onclick': "submitSelection('/ligand/biased');",
            'color': 'success',
        },
        "pathway": {
            "label": "Pathway-biased ligands<br>(balanced reference ligand)",
            "modal": "path_bias",
            "color": 'default active',
            "sameSize": True,
        },
        "biased": {
            "label": "Biased ligands<br>(any reference ligand)",
            "modal": "biased",
            "color": 'default active',
            "sameSize": True,
        },
    }

    table_data = getReferenceTable("different_family", "tested_assays")

#Biased Pathway Preferred Browser
class BiasPredictionTargetSelection(AbsReferenceSelectionTable):
    step = 1
    number_of_steps = 1
    filters = False
    filter_tableselect = False
    family_tree = False
    docs = 'sequences.html#structure-based-alignments'
    title = "SELECT RECEPTORS to retrieve ligands with a preferred G protein or arrestin pathway (ΔLog(Emax/EC50  values across pathways for one ligand (no reference ligand)))"
    description = 'Select a receptor in the table (below) or browse the classification tree (center).\nThen click the green button (right)'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Physiology-biased ligands<br>(no reference ligand)',
            'onclick': "submitSelection('/ligand/biasedpredicted');",
            'color': 'success',
        },
    }

    table_data = getReferenceTable('predicted_family', 'predicted_tested_assays')

#Biased Subtype Browser
class BiasGTargetSelection(AbsReferenceSelectionTable):
    step = 1
    number_of_steps = 1
    filters = False
    filter_tableselect = False
    family_tree = False
    docs = 'sequences.html#structure-based-alignments'
    title = "SELECT RECEPTORS with ligands biased for a G protein or arrestin subtype"
    description = 'Select a receptor in the table (below) or browse the classification tree (center).\nThen click a green button (right)'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        "continue": {
            "label": "Physiology-biased ligands<br>(endogenous agonist reference)",
            "onclick": "submitSelection('/ligand/biasedsubtypes');",
            "color": 'success',
            "sameSize": True,
        },
        "pathway": {
            "label": "Pathway-biased ligands<br>(balanced reference ligand)",
            "modal": "path_bias",
            "color": 'default active',
            "sameSize": True,
        },
        "biased": {
            "label": "Biased ligands<br>(any reference ligand)",
            "modal": "biased",
            "color": 'default active',
            "sameSize": True,
        },
    }

    table_data = getReferenceTable("sub_different_family", "sub_tested_assays")

class RankOrderSelection(AbsReferenceSelectionTable):
    step = 1
    number_of_steps = 1
    filters = False
    filter_tableselect = False
    family_tree = False
    title = "SELECT RECEPTOR for Ligand Bias Rank Order ΔΔLog(Emax/EC50)"
    description = 'Select a receptor in the table (below).' \
        + ' \n\nand click the green button (upper right).'
    selection_boxes = OrderedDict([
        ('reference', True),
        ('targets', False),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Physiology-biased ligands<br>(endogenous agonist reference)',
            'onclick': "submitSelection('/ligand/emax_rankorder');",
            'color': 'success',
        },
        "pathway": {
            "label": "Pathway-biased ligands<br>(balanced reference ligand)",
            "modal": "path_bias",
            "color": 'default active',
            "sameSize": True,
        },
        "biased": {
            "label": "Biased ligands<br>(any reference ligand)",
            "modal": "biased",
            "color": 'default active',
            "sameSize": True,
        },
    }

    table_data = getReferenceTable("different_family", "tested_assays")

class TauRankOrderSelection(AbsReferenceSelectionTable):
    step = 1
    number_of_steps = 1
    filters = False
    filter_tableselect = False
    family_tree = False
    title = "SELECT RECEPTOR for Ligand Bias Rank Order ΔΔLog(Tau/KA)"
    description = 'Select a receptor in the table (below).' \
        + ' \n\nand click the green button (upper right).'
    selection_boxes = OrderedDict([
        ('reference', True),
        ('targets', False),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Physiology-biased ligands<br>(endogenous agonist reference)',
            'onclick': "submitSelection('/ligand/tau_rankorder');",
            'color': 'success',
        },
        "pathway": {
            "label": "Pathway-biased ligands<br>(balanced reference ligand)",
            "modal": "path_bias",
            "color": 'default active',
            "sameSize": True,
        },
        "biased": {
            "label": "Biased ligands<br>(any reference ligand)",
            "modal": "biased",
            "color": 'default active',
            "sameSize": True,
        },
    }

    table_data = getReferenceTable("different_family", "tested_assays")

class TauSubtypeRankOrderSelection(AbsReferenceSelectionTable):
    step = 1
    number_of_steps = 1
    filters = False
    filter_tableselect = False
    family_tree = False
    title = "SELECT RECEPTOR for Ligand Bias (subtype) Rank Order ΔΔLog(Tau/KA)"
    description = 'Select a receptor in the table (below).' \
        + ' \n\nand click the green button (upper right).'
    selection_boxes = OrderedDict([
        ('reference', True),
        ('targets', False),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Physiology-biased ligands<br>(endogenous agonist reference)',
            'onclick': "submitSelection('/ligand/subtype_tau_rankorder');",
            'color': 'success',
        },
        "pathway": {
            "label": "Pathway-biased ligands<br>(balanced reference ligand)",
            "modal": "path_bias",
            "color": 'default active',
            "sameSize": True,
        },
        "biased": {
            "label": "Biased ligands<br>(any reference ligand)",
            "modal": "biased",
            "color": 'default active',
            "sameSize": True,
        },
    }

    table_data = getReferenceTable("sub_different_family", "sub_tested_assays")

class EmaxPathProfileSelection(AbsReferenceSelectionTable):
    step = 1
    number_of_steps = 1
    filters = False
    filter_tableselect = False
    family_tree = False
    title = "SELECT RECEPTOR for Ligand Pathway Profiles ΔLog(Emax/EC50)"
    description = 'Select a receptor in the table (below).' \
        + ' \n\nand click the green button (upper right).'
    selection_boxes = OrderedDict([
        ('reference', True),
        ('targets', False),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Physiology-biased ligands<br>(endogenous agonist reference)',
            'onclick': "submitSelection('/ligand/emax_path_profiles');",
            'color': 'success',
        },
        "pathway": {
            "label": "Pathway-biased ligands<br>(balanced reference ligand)",
            "modal": "path_bias",
            "color": 'default active',
            "sameSize": True,
        },
        "biased": {
            "label": "Biased ligands<br>(any reference ligand)",
            "modal": "biased",
            "color": 'default active',
            "sameSize": True,
        },
    }

    table_data = getReferenceTable("different_family", "tested_assays")

class TauPathProfileSelection(AbsReferenceSelectionTable):
    step = 1
    number_of_steps = 1
    filters = False
    filter_tableselect = False
    family_tree = False
    title = "SELECT RECEPTOR for Ligand Pathway Profiles ΔΔLog(Tau/KA)"
    description = 'Select a receptor in the table (below).' \
        + ' \n\nand click the green button (upper right).'
    selection_boxes = OrderedDict([
        ('reference', True),
        ('targets', False),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Physiology-biased ligands<br>(endogenous agonist reference)',
            'onclick': "submitSelection('/ligand/tau_path_profiles');",
            'color': 'success',
        },
        "pathway": {
            "label": "Pathway-biased ligands<br>(balanced reference ligand)",
            "modal": "path_bias",
            "color": 'default active',
            "sameSize": True,
        },
        "biased": {
            "label": "Biased ligands<br>(any reference ligand)",
            "modal": "biased",
            "color": 'default active',
            "sameSize": True,
        },
    }

    table_data = getReferenceTable("different_family", "tested_assays")

class TauSubtypePathProfileSelection(AbsReferenceSelectionTable):
    step = 1
    number_of_steps = 1
    filters = False
    filter_tableselect = False
    family_tree = False
    title = "SELECT RECEPTOR for Ligand (subtype) Pathway Profiles ΔΔLog(Tau/KA)"
    description = 'Select a receptor in the table (below).' \
        + ' \n\nand click the green button (upper right).'
    selection_boxes = OrderedDict([
        ('reference', True),
        ('targets', False),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Physiology-biased ligands<br>(endogenous agonist reference)',
            'onclick': "submitSelection('/ligand/subtype_tau_path_profiles');",
            'color': 'success',
        },
        "pathway": {
            "label": "Pathway-biased ligands<br>(balanced reference ligand)",
            "modal": "path_bias",
            "color": 'default active',
            "sameSize": True,
        },
        "biased": {
            "label": "Biased ligands<br>(any reference ligand)",
            "modal": "biased",
            "color": 'default active',
            "sameSize": True,
        },
    }

    table_data = getReferenceTable("sub_different_family", "sub_tested_assays")

class EmaxPathPrefRankOrderSelection(AbsReferenceSelectionTable):
    step = 1
    number_of_steps = 1
    filters = False
    filter_tableselect = False
    family_tree = False
    title = "SELECT RECEPTOR for Ligand Pathway Preference rank orders (ΔLog(Emax/EC50))"
    description = 'Select a receptor in the table (below).' \
        + ' \n\nand click the green button (upper right).'
    selection_boxes = OrderedDict([
        ('reference', True),
        ('targets', False),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Physiology-biased ligands<br>(no reference ligand)',
            'onclick': "submitSelection('/ligand/path_preference_emax_rankorder');",
            'color': 'success',
        },
    }

    table_data = getReferenceTable('predicted_family', 'predicted_tested_assays')


class EmaxPathPrefPathProfilesSelection(AbsReferenceSelectionTable):
    step = 1
    number_of_steps = 1
    filters = False
    filter_tableselect = False
    family_tree = False
    title = "SELECT RECEPTOR for Ligand Pathway Profiles (log(Emax/EC50))"
    description = 'Select a receptor in the table (below).' \
        + ' \n\nand click the green button (upper right).'
    selection_boxes = OrderedDict([
        ('reference', True),
        ('targets', False),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Physiology-biased ligands<br>(no reference ligand)',
            'onclick': "submitSelection('/ligand/path_preference_emax_path_profiles');",
            'color': 'success',
        },
    }

    table_data = getReferenceTable('predicted_family', 'predicted_tested_assays')


class EmaxSubtypeRankOrderSelection(AbsReferenceSelectionTable):
    step = 1
    number_of_steps = 1
    filters = False
    filter_tableselect = False
    family_tree = False
    title = "SELECT RECEPTOR for Ligand biased (subtype) rank orders (ΔΔLog(Emax/EC50))"
    description = 'Select a receptor in the table (below).' \
        + ' \n\nand click the green button (upper right).'
    selection_boxes = OrderedDict([
        ('reference', True),
        ('targets', False),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Physiology-biased ligands<br>(endogenous agonist reference)',
            'onclick': "submitSelection('/ligand/subtype_emax_rankorder');",
            'color': 'success',
        },
        "pathway": {
            "label": "Pathway-biased ligands<br>(balanced reference ligand)",
            "modal": "path_bias",
            "color": 'default active',
            "sameSize": True,
        },
        "biased": {
            "label": "Biased ligands<br>(any reference ligand)",
            "modal": "biased",
            "color": 'default active',
            "sameSize": True,
        },
    }

    table_data = getReferenceTable('sub_different_family', 'sub_tested_assays')

class EmaxSubtypePathProfilesSelection(AbsReferenceSelectionTable):
    step = 1
    number_of_steps = 1
    filters = False
    filter_tableselect = False
    family_tree = False
    title = "SELECT RECEPTOR for Ligand biased (subtype) Pathway Profiles (Δlog(Emax/EC50))"
    description = 'Select a receptor in the table (below).' \
        + ' \n\nand click the green button (upper right).'
    selection_boxes = OrderedDict([
        ('reference', True),
        ('targets', False),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Physiology-biased ligands<br>(endogenous agonist reference)',
            'onclick': "submitSelection('/ligand/subtype_emax_path_profiles');",
            'color': 'success',
        },
        "pathway": {
            "label": "Pathway-biased ligands<br>(balanced reference ligand)",
            "modal": "path_bias",
            "color": 'default active',
            "sameSize": True,
        },
        "biased": {
            "label": "Biased ligands<br>(any reference ligand)",
            "modal": "biased",
            "color": 'default active',
            "sameSize": True,
        },
    }

    table_data = getReferenceTable('sub_different_family', 'sub_tested_assays')


class BiasedRankOrder(TemplateView):
    #set a global variable for different pages
    page = "rankorder"
    label = "emax"
    template_name = "biased_rank_orders.html"
    source = "different_family"
    assay = "tested_assays"

    @staticmethod
    def jitter_tooltip(page, assay, ligand, value, headers, prefix='', small_data=None, large_data=None, small_ref=None):
        #small and large data has to structured
        #small --> pathway/value/value
        #large --> pathway1/delta/value/value pathway2/delta/value/value
        ref_small = ''
        small = ''
        large = ''
        head = "<b>Compound Name:</b> " + str(ligand) + \
               "<br><b>Plotted Value:</b> " + str(value) + \
               "<hr class='solid'>"
        if small_data:
            #small table to show reference ligand or single datapoint
            small =  "<table>" + \
                     "      <tr>" + \
                     "        <th>" + str(small_data[3]) + "</th>" + \
                     "        <th>" + headers[0] + "</th>" + \
                     "        <th>" + headers[1] + "</th>" + \
                     "      </tr>" + \
                     "      <tr>" + \
                     "        <td>" + str(small_data[0]) + "</td>" + \
                     "        <td>" + str(small_data[1]) + "</td>" + \
                     "        <td>" + str(small_data[2]) + "</td>" + \
                     "      </tr>" + \
                     "</table>" + \
                     "<hr class='solid'>"
        if small_ref:
            #small table to show reference ligand or single datapoint
            ref_small =  "<table>" + \
                         "      <tr>" + \
                         "        <th>" + str(small_ref[3]) + "</th>" + \
                         "        <th>" + headers[0] + "</th>" + \
                         "        <th>" + headers[1] + "</th>" + \
                         "      </tr>" + \
                         "      <tr>" + \
                         "        <td>" + str(small_ref[0]) + "</td>" + \
                         "        <td>" + str(small_ref[1]) + "</td>" + \
                         "        <td>" + str(small_ref[2]) + "</td>" + \
                         "      </tr>" + \
                         "</table>"
        if large_data:
            #large table showing also Δ data
            large =  "<table>" + \
                     "      <tr>" + \
                     "        <th>" + str(ligand) + "</th>" + \
                     "        <th>" + prefix + "Log(" + headers[0] + "/" + headers[1] + ") </th>" + \
                     "        <th>" + headers[0] + "</th>" + \
                     "        <th>" + headers[1] + "</th>" + \
                     "      </tr>" + \
                     "      <tr>" + \
                     "        <td>" + str(large_data[0]) + "</td>" + \
                     "        <td>" + str(large_data[1]) + "</td>" + \
                     "        <td>" + str(large_data[2]) + "</td>" + \
                     "        <td>" + str(large_data[3]) + "</td>" + \
                     "      </tr>" + \
                     "      <tr>" + \
                     "        <td>" + str(large_data[4]) + "</td>" + \
                     "        <td>" + str(large_data[5]) + "</td>" + \
                     "        <td>" + str(large_data[6]) + "</td>" + \
                     "        <td>" + str(large_data[7]) + "</td>" + \
                     "      </tr>" + \
                     "</table>" + \
                     "<hr class='solid'>"

        if assay == "predicted_tested_assays":
            #no reference values
            if page == "rankorder":
                tip = head + large
                #dot plots without reference values
            else:
                tip = head + small
                #line charts without reference values
        else:
            #reference values required
            if page == "rankorder":
                tip = head + large + small
                #dot plots with reference values
            else:
                tip = head + small + ref_small
                #line chart wtih reference values
        return tip
        #line charts with reference values

    @staticmethod
    def create_rgb_color(): #, name, power): # pseudo-randomization function
        # h = hash( name + str(power) ) # hash string and int together
        # if h < 0: # ensure positive number
        #     h = h * -1
        # random.seed(h) # set the seed to use for randomization
        cryptogen = SystemRandom()          # added cryptography random number creation
        output = cryptogen.randrange(0,225) # 225 instead of 255 to avoid all very pale colors
        return output

    def get_context_data(self, **kwargs):

        context = super().get_context_data(**kwargs)
        simple_selection = self.request.session.get('selection', False)
        #I know it's a for cycle, but it should be just one element
        #since it's coming from a reference
        for item in simple_selection.reference:
            receptor = item.item

        rec_name = list(Protein.objects.filter(id=receptor)
                                        .values_list('entry_name',
                                                     'name'))

        tooltip_dict =  {'G12/13': 'G<sub>12/13</sub>',
                         'Gi/o': 'G<sub>i/o</sub>',
                         'Gq/11': 'G<sub>q/11</sub>',
                         'Gs': 'G<sub>s</sub>'}

        sign_prot_conversion = {'-' : 'None',
                                'arrestin-2 (b-arrestin-1)': 'Arrestin 2 (&beta;<sub>1</sub>)',
                                'arrestin-3 (b-arrestin-2)': 'Arrestin 3 (&beta;<sub>2</sub>)',
                                'g11': 'G<sub>11</sub>',
                                'g12': 'G<sub>12</sub>',
                                'g13': 'G<sub>13</sub>',
                                'g14': 'G<sub>14</sub>',
                                'g15': 'G<sub>15</sub>',
                                'g16': 'G<sub>16</sub>',
                                'gaqi5': 'G&alpha;qi5',
                                'gaqδ6i4myr': 'G&alpha;qδ6i4myr',
                                'gi1': 'G<sub>i1</sub>',
                                'gi2': 'G<sub>i2</sub>',
                                'gi3': 'G<sub>i3</sub>',
                                'goa': 'G<sub>o&alpha;</sub>',
                                'gob': 'G<sub>o&beta;</sub>',
                                'golf': 'G<sub>olfactory</sub>',
                                'gq': 'G<sub>q</sub>',
                                'gs': 'G<sub>s</sub>',
                                'gz': 'G<sub>z</sub>',
                                'minigi': 'Mini-G'}

        if self.assay != "predicted_tested_assays":
            prefix = 'Δ'
        else:
            prefix = ''

        upgrade_value = ["High activity", "High activity (Potency and Emax)", "Full agonism"]
        downgrade_value = ["Low activity", "No activity", "Inverse agonism/antagonism"]
        exclude_list = ["Agonism", "Partial agonism", "Medium activity"] #full agonism should be removed
        publications = list(AnalyzedAssay.objects.filter(
                        family__isnull=False,
                        experiment__receptor=receptor,
                        assay_description=self.assay,
                        experiment__source=self.source).exclude(
                        qualitative_activity__in=exclude_list,
                        ).values_list(
                        "family", #pathway                                              0
                        "experiment__ligand__name", #name                               1
                        "experiment__publication__web_link__index", # DOI               2
                        "experiment__publication__year", #year                          3
                        "experiment__publication__journal__name", #journal              4
                        "experiment__publication__authors",  #authors                   5
                        "experiment__ligand",    #ligand_id for hash                    6
                        "experiment__endogenous_ligand__name", #endogenous              7
                        "qualitative_activity",  #activity values                       8
                        "log_bias_factor_a", #ΔLog(Emax/EC50)                           9
                        "log_bias_factor",  #ΔΔLog(Emax/EC50)                           10
                        "order_no",         #ranking                                    11
                        "t_coefficient",    #ΔLog(TAU/Ka)                               12
                        "t_factor",         #ΔΔLog(TAU/Ka)                              13
                        "quantitive_activity",   #EC 50                                 14
                        "quantitive_efficacy",   #Emax                                  15
                        "reference_assay_initial_id__quantitive_activity", #            16 EC 50 compared drug
                        "reference_assay_initial_id__quantitive_efficacy", #            17 Emax compared drug
                        "reference_assay_initial_id__family", # Pathway compared drug   18
                        "reference_assay_initial_id__biased_experiment__ligand__name", #19 Name compared drug
                        "signalling_protein",   #subtype                                20
                        "reference_assay_initial_id__signalling_protein", #             21 Compared signalling protein
                        ).distinct()) #check

        list_of_ligands = []
        list_of_publications = []
        full_data = {}
        full_ligands = {}
        SpiderOptions = {}
        jitterDict = {}
        jitterPlot = {}
        jitterLegend = {}
        Colors = {}
        pathway_nr = {}
        labels_dict = {}

        for result in publications:
            if self.assay == 'sub_tested_assays':
                reference_path = sign_prot_conversion[result[21]]
            else:
                try:
                    reference_path = tooltip_dict[result[18]]
                except KeyError:
                    reference_path = result[18]

            #checking the value to plot
            #based on the landing page
            if self.label == 'emax':
                single_delta = result[9]
                double_delta = result[10]
                emax_tau = result[15]
                EC50_ka = result[14]
                components = ['Emax', 'EC50']
                reference_emax_tau = result[17]
                reference_EC50_ka = result[16]
            else:
                single_delta = result[12]
                double_delta = result[13]
                emax_tau = "NA" #need to be updated IF datacolumn for TAU will be added
                EC50_ka = "NA"  #need to be updated IF datacolumn for KA will be added
                components = ['Tau', 'KA']
                reference_emax_tau = "NA" #need to be updated IF datacolumn for TAU will be added
                reference_EC50_ka = "NA"  #need to be updated IF datacolumn for KA will be added

            #fixing ligand name (hash hash baby)
            lig_name = result[1]
            if result[1][0].isdigit():
                lig_name = "Ligand-"+result[1]

            hashed_lig_name = 'L' + hashlib.md5((str(result[6])).encode('utf-8')).hexdigest()
            # replace second white space with closing and opening tspan for svg
            journal_name = result[4]
            if result[4]:
                if ' ' in result[4]:
                    journal_name  = re.sub(r'(\s\S*?)\s', r'\1 closeTS openTS ', result[4])
            else:
                journal_name = "Not listed"
            if result[5] == None:
                authors = "Authors not listed, (" + str(result[3]) + ')'
                shortAuthors = "Authors not listed, (" + str(result[3]) + ')'
                jitterAuthors = 'openTS ' + shortAuthors + ' closeTS openTS ' + journal_name + ' closeTS openTS (' + str(result[3]) + ') closeTS'
                labels_dict[jitterAuthors] = shortAuthors
            else:
                authors = result[5].split(',')[0] + ', et al., ' + str(result[4]) + ', (' + str(result[3]) + ')'
                shortAuthors = result[5].split(',')[0] + ', et al.'
                jitterAuthors = 'openTS ' + shortAuthors + ' closeTS openTS ' + journal_name + ' closeTS openTS (' + str(result[3]) + ') closeTS'
                labels_dict[jitterAuthors] = shortAuthors

            list_of_ligands.append(tuple((hashed_lig_name, lig_name)))
            list_of_publications.append(authors)
            #start parsing the data to create the big dictionary
            if single_delta == None:               #∆Emax / ∆Tau flex
                value = 0
            else:
                value = float(single_delta)

            if result[0] not in jitterPlot.keys():
                jitterLegend[result[0]] = []
                jitterPlot[result[0]] = []

            if jitterAuthors not in jitterDict.keys():
                jitterDict[jitterAuthors] =  {}

            if lig_name not in jitterDict[jitterAuthors].keys():
                jitterDict[jitterAuthors][lig_name] = {}

            try:
                DD = [float(double_delta), "REAL"]
            except (ValueError, TypeError):
                DD = [0, double_delta]

            if result[11] == 1:
                try:
                    jitterDict[jitterAuthors][lig_name]['2nd_Pathway'] = tooltip_dict[result[0]]
                    jitterDict[jitterAuthors][lig_name]['2nd_Pathway_delta'] = value
                    jitterDict[jitterAuthors][lig_name]['2nd_Pathway_emax_tau'] = emax_tau
                    jitterDict[jitterAuthors][lig_name]['2nd_Pathway_EC50_KA'] = EC50_ka
                except KeyError:
                    jitterDict[jitterAuthors][lig_name]['2nd_Pathway'] = result[0]
                    jitterDict[jitterAuthors][lig_name]['2nd_Pathway_delta'] = value
                    jitterDict[jitterAuthors][lig_name]['2nd_Pathway_emax_tau'] = emax_tau
                    jitterDict[jitterAuthors][lig_name]['2nd_Pathway_EC50_KA'] = EC50_ka
                jitterDict[jitterAuthors][lig_name]['deltadelta'] = DD
                jitterDict[jitterAuthors][lig_name]['signalling_prot'] = result[20]

            if result[11] == 0:
                jitterDict[jitterAuthors][lig_name]["Pathway"] = result[0]
                jitterDict[jitterAuthors][lig_name]["delta"] = value
                jitterDict[jitterAuthors][lig_name]["Emax_Tau"] = emax_tau
                jitterDict[jitterAuthors][lig_name]["EC50_KA"] = EC50_ka

            tooltip_info = BiasedRankOrder.jitter_tooltip(self.page, self.assay, lig_name, value, components,
                                                          small_data=[result[0], emax_tau, EC50_ka, lig_name],
                                                          small_ref=[result[18], reference_emax_tau, reference_EC50_ka, result[19]])

            # initialization of the dictionary for new publication
            if result[2] not in full_data.keys():
                full_ligands[result[2]] = []
                full_data[result[2]] = {"Authors": authors,
                                        "Journal": result[4],
                                        "Year": result[3],
                                        "Endogenous": result[7],
                                        "Ligand": [lig_name],
                                        "Qualitative": [result[8]],
                                        "Data": [{"name": hashed_lig_name,
                                                 "PathwaysData":
                                                    [{"Pathway": result[0],
                                                      "value": [value, "REAL"],
                                                      "tooltip": tooltip_info}]}]}
                SpiderOptions[authors] = {}
                SpiderOptions[authors][hashed_lig_name] = {"Data":[[{'axis':result[0],
                                                                     'value':value}]],
                                                             "Options": {'levels': 4,
                                                                         'maxValue': 4,
                                                                         'roundStrokes': False,
                                                                         'title': lig_name}}
            #new ligand pushed into existing publication
            if hashed_lig_name not in [d['name'] for d in full_data[result[2]]["Data"]]:
                full_data[result[2]]["Ligand"].append(lig_name)
                full_data[result[2]]["Qualitative"].append(result[8])
                full_data[result[2]]["Data"].append(
                                            {"name": hashed_lig_name,
                                             "PathwaysData":
                                                [{"Pathway": result[0],
                                                  "value": [value,"REAL"],
                                                  "tooltip": tooltip_info}]})

                SpiderOptions[authors][hashed_lig_name] = {"Data":[[{'axis':result[0],
                                                                     'value':value}]],
                                                             "Options": {'levels': 4,
                                                                         'maxValue': 4,
                                                                         'roundStrokes': False,
                                                                         'title': lig_name}}
            else:
                ID = [d['name'] for d in full_data[result[2]]["Data"]].index(hashed_lig_name)
                if result[0] not in [d["Pathway"] for d in full_data[result[2]]["Data"][ID]["PathwaysData"]]:
                    full_data[result[2]]["Data"][ID]["PathwaysData"].append(
                                            {"Pathway": result[0],
                                             "value": [value, "REAL"],
                                             "tooltip": tooltip_info})
                    SpiderOptions[authors][hashed_lig_name]["Data"][0].append({'axis':result[0],
                                                                               'value':value})

        for item in full_data.keys():
            vals = []
            for name in full_data[item]["Data"]:
                to_fix = [subdict['value'][0] for subdict in name['PathwaysData']]
                vals = vals + to_fix
            MIN = min(vals) - 1
            MAX = max(vals) + 1
            for name in full_data[item]["Data"]:
                to_fix = [subdict['value'][0] for subdict in name['PathwaysData']]
                indices = [i for i, x in enumerate(to_fix) if x == 0]
                quality = full_data[item]["Qualitative"][full_data[item]["Data"].index(name)]
                if quality in upgrade_value:
                    try:
                        for i in indices:
                            name["PathwaysData"][i]["value"] = [MAX,"ARTIFICIAL"]
                            name["PathwaysData"][i]["tooltip"] = BiasedRankOrder.jitter_tooltip(self.page, self.assay, lig_name, value, components,
                                                                                                small_data=[result[0], emax_tau, 'High', lig_name],
                                                                                                small_ref=[result[18], reference_emax_tau, reference_EC50_ka, result[19]])
                    except ValueError:
                        continue
                if quality in downgrade_value:
                    try:
                        for i in indices:
                            name["PathwaysData"][i]["value"] = [MIN,"ARTIFICIAL"]
                            name["PathwaysData"][i]["tooltip"] = BiasedRankOrder.jitter_tooltip(self.page, self.assay, lig_name, value, components,
                                                                                                small_data=[result[0], emax_tau, 'Low', lig_name],
                                                                                                small_ref=[result[18], reference_emax_tau, reference_EC50_ka, result[19]])
                    except ValueError:
                        continue

        #The big dictionary is created, now it needs to be ordered
        #the publication with more pathway (or tied) should be first

        #From this, reorder the original one
        for key in full_data.keys():
            for item in full_data[key]["Data"]:
                if len(item['PathwaysData']) not in pathway_nr.keys():
                    pathway_nr[len(item['PathwaysData'])] = []
                pathway_nr[len(item['PathwaysData'])].append(tuple((key, len(full_data[key]["Data"]))))
        #get the number of pathways available
        keys = sorted(list(pathway_nr.keys()))[::-1]
        #starting from the highest value
        sorted_full_data = {}
        for value in keys:
            pairs = sorted(list(set(pathway_nr[value])), key=lambda x: x[1], reverse=True)
            for item in pairs:
                sorted_full_data[item[0]] = deepcopy(full_data[item[0]])
                # test_data[item[0]] = deepcopy(full_data[item[0]])
                max_values = []
                for name in full_data[item[0]]['Data']:
                    max_values.append(tuple((full_data[item[0]]['Data'].index(name), max([subdict['value'] for subdict in name['PathwaysData']]), name['name'])))
                max_values = sorted(max_values, key=lambda x: x[1], reverse=True)
                sorted_full_data[item[0]]["Data"] = [] #need to be sorted_full_data
                for couple in max_values:
                    sorted_full_data[item[0]]["Data"].append(full_data[item[0]]["Data"][couple[0]]) #need to be sorted_full_data
                    if tuple((couple[2], full_data[item[0]]["Ligand"][couple[0]])) not in full_ligands[item[0]]:
                        full_ligands[item[0]].append(tuple((couple[2], full_data[item[0]]["Ligand"][couple[0]])))

        #now the sorted dict is done, we can clear cache the og one
        del full_data

        for pub in jitterDict.keys():
            for ligand in jitterDict[pub]:
                try:
                    if ligand not in Colors.keys():
                        color = '#%02x%02x%02x' % (BiasedRankOrder.create_rgb_color(), BiasedRankOrder.create_rgb_color(), BiasedRankOrder.create_rgb_color())
                        Colors[ligand] = color
                    little = [reference_path, reference_EC50_ka, reference_emax_tau, result[19]]
                    if self.assay == "sub_tested_assays":
                        big = [sign_prot_conversion[jitterDict[pub][ligand]["signalling_prot"]], jitterDict[pub][ligand]['delta'], jitterDict[pub][ligand]['Emax_Tau'], jitterDict[pub][ligand]['EC50_KA'],
                               sign_prot_conversion[jitterDict[pub][ligand]["signalling_prot"]], jitterDict[pub][ligand]['2nd_Pathway_delta'], jitterDict[pub][ligand]['2nd_Pathway_emax_tau'], jitterDict[pub][ligand]['2nd_Pathway_EC50_KA']]
                    else:
                        big = [jitterDict[pub][ligand]["Pathway"], jitterDict[pub][ligand]['delta'], jitterDict[pub][ligand]['Emax_Tau'], jitterDict[pub][ligand]['EC50_KA'],
                               jitterDict[pub][ligand]['2nd_Pathway'], jitterDict[pub][ligand]['2nd_Pathway_delta'], jitterDict[pub][ligand]['2nd_Pathway_emax_tau'], jitterDict[pub][ligand]['2nd_Pathway_EC50_KA']]
                    if (jitterDict[pub][ligand]['deltadelta'][1] == 'High Bias') or (jitterDict[pub][ligand]['deltadelta'][1] == 'Full Bias'):
                        tooltip = BiasedRankOrder.jitter_tooltip(self.page, self.assay, ligand, jitterDict[pub][ligand]['deltadelta'][1], components, prefix, small_data=little, large_data=big)
                    else:
                        tooltip = BiasedRankOrder.jitter_tooltip(self.page, self.assay, ligand, jitterDict[pub][ligand]['deltadelta'][0], components, prefix, small_data=little, large_data=big)
                    jitterPlot[jitterDict[pub][ligand]["Pathway"]].append([pub, jitterDict[pub][ligand]['deltadelta'][0], Colors[ligand], ligand, jitterDict[pub][ligand]['deltadelta'][1], tooltip])
                    jitterLegend[jitterDict[pub][ligand]["Pathway"]].append(tuple((ligand, jitterDict[pub][ligand]['deltadelta'][0])))
                except KeyError:
                    continue
                jitterLegend[jitterDict[pub][ligand]["Pathway"]] = sorted(list(set(jitterLegend[jitterDict[pub][ligand]["Pathway"]])), key=lambda x: x[1], reverse=True)

        #addressing qualitative points in the ∆∆ data (full bias/high bias/none)
        for pathway in jitterPlot.keys():
            highest = 0
            for datapoint in jitterPlot[pathway]:
                if datapoint[1] > highest:
                    highest = datapoint[1]
            change = {'None' : 0, "High Bias": highest + 1, "Full Bias": highest + 1.5}
            for datapoint in jitterPlot[pathway]:
                if datapoint[4] in change:
                    datapoint[1] = change[datapoint[4]]

        jitterPlot = {k: v for k, v in jitterPlot.items() if len(v) > 0}

        for key in jitterLegend.keys():
            jitterLegend[key] = list(dict.fromkeys([name[0] for name in jitterLegend[key]]))[:20]

        context['column_dict'] = json.dumps(labels_dict)
        context['source'] = self.source
        context['label'] = self.label
        context['page'] = self.page
        context['scatter_legend'] = json.dumps(jitterLegend)
        context['colors'] = json.dumps(Colors)
        context['scatter'] = json.dumps(jitterPlot)
        context['spider'] = json.dumps(SpiderOptions)
        context['all_ligands'] = list(set(list_of_ligands))
        context['all_publications'] = list(set(list_of_publications))
        context['query'] = rec_name[0][0]
        context['IUPHAR'] = rec_name[0][1]
        context['full_data'] = json.dumps(sorted_full_data)
        context['full_ligands'] = json.dumps(full_ligands)
        return context

class LigandStatistics(TemplateView):
    """
    Per class statistics of known ligands.
    """

    template_name = 'ligand_statistics.html'
    page = 'ligands' #either ligand_bias or pathway_pref OR subtype (NEW PAGE)

    def get_context_data(self, **kwargs):

        context = super().get_context_data(**kwargs)
        # assays = AssayExperiment.objects.all().prefetch_related('protein__family__parent__parent__parent', 'protein__family')
        lig_count_dict = {}
        # create a dict for conversion of input data,
        # if new families will be added to db just add values to this dict
        # values are: experiment_source, assay_description, render for html
        conversion = {
        'ligand_bias': ['different_family', 'tested_assays', 'bias'],
        'pathway_pref': ['predicted_family', 'predicted_tested_assays', 'pathway'],
        'subtype': ['sub_different_family', 'sub_tested_assays', 'subtype']
        }
        if self.page == 'ligands':
            assays_lig = list(AssayExperiment.objects.all().values(
                'protein__family__parent__parent__parent__name').annotate(c=Count('ligand', distinct=True)))
            for a in assays_lig:
                lig_count_dict[a['protein__family__parent__parent__parent__name']] = a['c']
        else:
            assays_lig = list(AnalyzedAssay.objects
                .filter(log_bias_factor__gte=1,
                        experiment__source=conversion[self.page][0])
                .values('experiment__receptor__family__parent__parent__parent__name')
                .annotate(c=Count('experiment__ligand_id', distinct=True)))
            for a in assays_lig:
                lig_count_dict[a['experiment__receptor__family__parent__parent__parent__name']] = a['c']

        target_count_dict = {}
        assays_target = list(AssayExperiment.objects.all().values(
            'protein__family__parent__parent__parent__name').annotate(c=Count('protein__family', distinct=True)))
        for a in assays_target:
            target_count_dict[a['protein__family__parent__parent__parent__name']] = a['c']

        prot_count_dict = {}
        proteins_count = list(Protein.objects.all().values(
            'family__parent__parent__parent__name').annotate(c=Count('family', distinct=True)))
        for pf in proteins_count:
            prot_count_dict[pf['family__parent__parent__parent__name']] = pf['c']

        if self.page == 'ligands':
            classes = ProteinFamily.objects.filter(
                slug__in=['001', '002', '003', '004', '005', '006', '007'])  # ugly but fast
        else:
            classes = ProteinFamily.objects.filter(
                slug__in=['001', '002', '003', '004', '006', '007'])  # ugly but fast

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
                'name': fam.name.replace('Class', ''),
                'num_ligands': lig_count,
                'avg_num_ligands': lig_count / prot_count,
                'target_percentage': target_count / prot_count * 100,
                'target_count': target_count
            })
        lig_count_total = sum([x['num_ligands'] for x in ligands])
        prot_count_total = Protein.objects.filter(
            family__slug__startswith='00').all().distinct('family').count()
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
        class_a_data = tree.get_tree_data(ProteinFamily.objects.get(name='Class A (Rhodopsin)'))
        context['class_a_options'] = deepcopy(tree.d3_options)
        context['class_a_options']['anchor'] = 'class_a'
        context['class_a_options']['leaf_offset'] = 50
        context['class_a_options']['label_free'] = []
        # section to remove Orphan from Class A tree and apply to a different tree
        whole_class_a = class_a_data.get_nodes_dict(self.page)
        for item in whole_class_a['children']:
            if item['name'] == 'Orphan':
                orphan_data = OrderedDict(
                    [('name', ''), ('value', 3000), ('color', ''), ('children', [item])])
                whole_class_a['children'].remove(item)
                break
        context['class_a'] = json.dumps(whole_class_a)
        class_b1_data = tree.get_tree_data(
            ProteinFamily.objects.get(name__startswith='Class B1 (Secretin)'))
        context['class_b1_options'] = deepcopy(tree.d3_options)
        context['class_b1_options']['anchor'] = 'class_b1'
        context['class_b1_options']['branch_trunc'] = 60
        context['class_b1_options']['label_free'] = [1, ]
        context['class_b1'] = json.dumps(
            class_b1_data.get_nodes_dict(self.page))
        class_b2_data = tree.get_tree_data(
            ProteinFamily.objects.get(name__startswith='Class B2 (Adhesion)'))
        context['class_b2_options'] = deepcopy(tree.d3_options)
        context['class_b2_options']['anchor'] = 'class_b2'
        context['class_b2_options']['label_free'] = [1, ]
        context['class_b2'] = json.dumps(
            class_b2_data.get_nodes_dict(self.page))
        class_c_data = tree.get_tree_data(
            ProteinFamily.objects.get(name__startswith='Class C (Glutamate)'))
        context['class_c_options'] = deepcopy(tree.d3_options)
        context['class_c_options']['anchor'] = 'class_c'
        context['class_c_options']['branch_trunc'] = 50
        context['class_c_options']['label_free'] = [1, ]
        context['class_c'] = json.dumps(class_c_data.get_nodes_dict(self.page))
        class_f_data = tree.get_tree_data(
            ProteinFamily.objects.get(name__startswith='Class F (Frizzled)'))
        context['class_f_options'] = deepcopy(tree.d3_options)
        context['class_f_options']['anchor'] = 'class_f'
        context['class_f_options']['label_free'] = [1, ]
        context['class_f'] = json.dumps(class_f_data.get_nodes_dict(self.page))
        class_t2_data = tree.get_tree_data(
            ProteinFamily.objects.get(name__startswith='Class T (Taste 2)'))
        context['class_t2_options'] = deepcopy(tree.d3_options)
        context['class_t2_options']['anchor'] = 'class_t2'
        context['class_t2_options']['label_free'] = [1, ]
        context['class_t2'] = json.dumps(
            class_t2_data.get_nodes_dict(self.page))
        # definition of the class a orphan tree
        context['orphan_options'] = deepcopy(tree.d3_options)
        context['orphan_options']['anchor'] = 'orphan'
        context['orphan_options']['label_free'] = [1, ]
        context['orphan'] = json.dumps(orphan_data)

        whole_receptors = Protein.objects.prefetch_related(
            "family", "family__parent__parent__parent")
        whole_rec_dict = {}
        for rec in whole_receptors:
            rec_uniprot = rec.entry_short()
            rec_iuphar = rec.family.name.replace("receptor", '').replace(
                "<i>", "").replace("</i>", "").strip()
            if (rec_iuphar[0].isupper()) or (rec_iuphar[0].isdigit()):
                whole_rec_dict[rec_uniprot] = [rec_iuphar]
            else:
                whole_rec_dict[rec_uniprot] = [rec_iuphar.capitalize()]

        context["whole_receptors"] = json.dumps(whole_rec_dict)
        if self.page == 'ligands':
            context["render"] = "not_bias"

        else:
            assay_qs = AnalyzedAssay.objects.filter(
                assay_description=conversion[self.page][1]).values_list(
                "family", "experiment__receptor__entry_name").order_by(
                "family", "experiment__receptor__entry_name").distinct(
                "family", "experiment__receptor__entry_name")

            ligand_qs = AnalyzedAssay.objects.filter(
                order_no=0,
                assay_description=conversion[self.page][1]).values_list(
                "family", "experiment__receptor__entry_name", "experiment__ligand").order_by(
                "family", "experiment__receptor__entry_name", "experiment__ligand").distinct(
                "family", "experiment__receptor__entry_name", "experiment__ligand")

            circles = {}
            for data in assay_qs:
                if data[1].split('_')[1] == 'human':
                    key = data[1].split('_')[0].upper()
                    if key not in circles.keys():
                        circles[key] = {}
                        circles[key][data[0]] = 0
                    else:
                        circles[key][data[0]] = 0

            for data in ligand_qs:
                if data[1].split('_')[1] == 'human':
                    key = data[1].split('_')[0].upper()
                    circles[key][data[0]] += 1

            context["circles_data"] = json.dumps(circles)
            context["render"] = conversion[self.page][2]

        return context

class ExperimentEntryView(DetailView):
    context_object_name = 'experiment'
    model = AnalyzedExperiment
    template_name = 'biased_experiment_data.html'

# Biased pathways part


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

    return HttpResponse(request)


class BiasVendorBrowser(TemplateView):

    template_name = 'biased_ligand_vendor.html'

    def get_context_data(self, **kwargs):
        # try:
        context = dict()
        datum = self.request.session.get('ids')

        self.request.session.modified = True
        rd = list()
        for i in datum.split(','):
            ligand = Ligand.objects.filter(id=i)
            ligand = ligand.get()
            links = LigandVendorLink.objects.filter(
                lp=ligand.properities_id).prefetch_related('lp', 'vendor')
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

class LigandInformationView(TemplateView):
    template_name = 'ligand_info.html'

    def get_context_data(self, *args, **kwargs):
        context = super(LigandInformationView, self).get_context_data(**kwargs)
        ligand_id = self.kwargs['pk']
        ligand_data = Ligand.objects.get(id=ligand_id)
        assay_data = list(AssayExperiment.objects.filter(ligand=ligand_id).prefetch_related(
            'ligand', 'ligand__properities', 'protein', 'protein__family',
            'protein__family__parent', 'protein__family__parent__parent__parent',
            'protein__family__parent__parent', 'protein__family', 'protein__species',
            'publication', 'publication__web_link', 'publication__web_link__web_resource',
            'publication__journal'))
        context = dict()
        structures = LigandInformationView.get_structure(ligand_data)
        ligand_data = LigandInformationView.process_ligand(ligand_data)
        assay_data = LigandInformationView.process_assay(assay_data)
        assay_data = LigandInformationView.process_values(assay_data)
        mutations = LigandInformationView.get_mutations(ligand_data)
        context.update({'structure': structures})
        context.update({'ligand': ligand_data})
        context.update({'assay': assay_data})
        context.update({'mutations': mutations})
        return context

    @staticmethod
    def get_structure(ligand):
        return_list = list()
        structures = list(
            StructureLigandInteraction.objects.filter(ligand=ligand))
        for i in structures:
            structure_dict = dict()
            structure_dict['structure_pdb'] = i.structure.pdb_code.index
            return_list.append(structure_dict)
        return return_list

    @staticmethod
    def get_mutations(ligand):
        return_set = set()
        return_list = list()
        mutations = list(
            MutationExperiment.objects.filter(ligand=ligand['ligand_id']).only('protein').order_by("protein__name"))
        for i in mutations:
            if i.protein.family_id in return_set:
                pass
            else:
                return_list.append({"id":i.protein.family_id, "name": i.protein.name.split(' ', 1)[0].split('-adrenoceptor', 1)[0].strip()})
                return_set.add(i.protein.family_id)

        return return_list

    @staticmethod
    def get_min_max_values(value):
        value = list(map(float, value))
        maximum = max(value)
        minimum = min(value)
        avg = sum(value) / len(value)
        return minimum, avg, maximum

    @staticmethod
    def process_assay(assays):
        return_dict = dict()
        for i in assays:
            name = str(i.protein)
            if name in return_dict:
                if i.standard_type == 'EC50' or i.standard_type == 'potency':
                    return_dict[name]['potency_values'].append(
                        i.standard_value)
                elif i.standard_type == 'IC50':
                    return_dict[name]['affinity_values'].append(
                        i.standard_value)
                # TODO: recalculate min max avg_num_ligands
            else:
                return_dict[name] = dict()
                return_dict[name]['potency_values'] = list()
                return_dict[name]['affinity_values'] = list()
                return_dict[name]['receptor_gtp'] = i.protein.name.split(
                    ' ', 1)[0].split('-adrenoceptor', 1)[0].strip()
                return_dict[name]['receptor_uniprot'] = i.protein.entry_name
                return_dict[name]['receptor_species'] = i.protein.species.common_name
                return_dict[name]['receptor_family'] = i.protein.family.parent.name
                return_dict[name]['receptor_class'] = i.protein.family.parent.parent.parent.name
                if i.standard_type == 'EC50' or i.standard_type == 'potency':
                    return_dict[name]['potency_values'].append(
                        i.standard_value)
                elif i.standard_type == 'IC50':
                    return_dict[name]['affinity_values'].append(
                        i.standard_value)
        return return_dict

    @staticmethod
    def process_values(return_dict):
        return_list = list()
        for item in return_dict.items():
            temp_dict = dict()
            temp_dict = item[1]

            if item[1]['potency_values']:
                temp_dict['potency_min'],temp_dict['potency_avg'],temp_dict['potency_max'] = LigandInformationView.get_min_max_values(
                    item[1]['potency_values'])
                return_list.append(temp_dict)
            if item[1]['affinity_values']:
                temp_dict['affinity_min'],temp_dict['affinity_avg'],temp_dict['affinity_max'] = LigandInformationView.get_min_max_values(
                    item[1]['affinity_values'])
                return_list.append(temp_dict)
        return return_list

    @staticmethod
    def process_ligand(ligand_data):
        ld = dict()
        ld['ligand_id'] = ligand_data.id
        ld['ligand_name'] = ligand_data.name
        ld['ligand_smiles'] = ligand_data.properities.smiles
        ld['ligand_inchikey'] = ligand_data.properities.inchikey
        try:
            ld['type'] = ligand_data.properities.ligand_type.name
        except:
            ld['type'] = None
        ld['rotatable'] = ligand_data.properities.rotatable_bonds
        ld['sequence'] = ligand_data.properities.sequence
        ld['hacc'] = ligand_data.properities.hacc
        ld['hdon'] = ligand_data.properities.hdon
        ld['logp'] = ligand_data.properities.logp
        ld['mw'] = ligand_data.properities.mw
        ld['wl'] = list()
        ld['picture'] = None
        for i in ligand_data.properities.web_links.all():
            ld['wl'].append({'name': i.web_resource.name, "link": str(i)})
            if i.web_resource.slug == 'chembl_ligand':
                ld['picture'] = i.index
        return ld


class BiasPathways(TemplateView):
    template_name = 'bias_browser_pathways.html'
    # @cache_page(50000)
    def get_context_data(self, *args, **kwargs):
        content = BiasedPathways.objects.all().prefetch_related(
            'biased_pathway', 'ligand', 'receptor', 'receptor', 'receptor__family',
            'receptor__family__parent', 'receptor__family__parent__parent__parent',
            'receptor__family__parent__parent', 'receptor__species',
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
        self.logger.info('receptor not found error')
        for instance in content:
            fin_obj = {}
            fin_obj['main'] = instance
            temp = dict()
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
                temp['IUPHAR'] = instance.receptor.name.split(' ', 1)[
                    0].strip()
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
            increment += 1
        return rd

    '''
    End  of Bias Browser
    '''

def CachedBiasBrowser(request):
    return CachedBiasBrowsers("biasbrowser", request)

def CachedBiasGBrowser(request):
    return CachedBiasBrowsers("biasgbrowser", request)

def CachedBiasPredictBrowser(request):
    return CachedBiasBrowsers("biasprecictedbrowser", request)

def CachedBiasBrowsers(browser_type, request):
    protein_ids = []
    try:
        simple_selection = request.session.get('selection', False)

        families = []
        for target in simple_selection.targets:
            if target.type == 'protein':
                protein_ids.append(target.item.entry_name)
            elif target.type == 'family':
                families.append(target.item.slug)

        if len(families) > 0:
            # species filter
            species_list = []
            for species in simple_selection.species:
                species_list.append(species.item)

            # annotation filter
            protein_source_list = []
            for protein_source in simple_selection.annotation:
                protein_source_list.append(protein_source.item)

            family_proteins = Protein.objects.filter(family__slug__in=families, source__in=(protein_source_list))
            if len(species_list) > 0:
                family_proteins.filter(species__in=(species_list))

            for prot_name in family_proteins.values_list("entry_name", flat=True):
                protein_ids.append(prot_name)

    except:
        protein_ids = ["NOSELECTION"]

    protein_ids.sort()
    cache_key = "BIASBROWSER_" + browser_type + "_" + hashlib.md5("_".join(protein_ids).encode('utf-8')).hexdigest()
    return_html = cache.get(cache_key)
    if return_html == None:
        if browser_type == "biasbrowser":
            return_html = BiasBrowser.as_view()(request).render()
        elif browser_type == "biasgbrowser":
            return_html = BiasGBrowser.as_view()(request).render()
        else:
            return_html = BiasPredictionBrowser.as_view()(request).render()
        cache.set(cache_key, return_html, 60*60*24*7)
    return return_html


'''
Bias browser between families
access data from db, fill empty fields with empty parse_children
'''
class BiasBrowser(ListView):
    # serializer_class = AnalyzedExperimentSerializer
    template_name = 'bias_browser.html'
    context_object_name = 'data_test'
    # @cache_page(50000)
    def get_queryset(self):
        protein_list = list()

        try:
            simple_selection = self.request.session.get('selection', False)
            a = Alignment()
            # load data from selection into the alignment
            a.load_proteins_from_selection(simple_selection)
            for items in a.proteins:
                protein_list.append(items.protein)
        except:
            protein_list.append(1)
        assay_qs = AnalyzedAssay.objects.filter(
            order_no__lte=5,
            assay_description='tested_assays',

            experiment=OuterRef('pk'),
        ).order_by('order_no')

        ref_assay_qs = AnalyzedAssay.objects.filter(
            order_no__lte=5,
            assay_description='endogenous',
            experiment=OuterRef('pk'),
        ).order_by('order_no')

        queryset = AnalyzedExperiment.objects.filter(
            source='different_family',
            receptor__in=protein_list,
        ).prefetch_related(
            'analyzed_data', 'ligand', 'ligand__reference_ligand', 'reference_ligand',
            'endogenous_ligand', 'ligand__properities', 'receptor', 'receptor', 'receptor__family',
            'receptor__family__parent', 'receptor__family__parent__parent__parent',
            'receptor__family__parent__parent', 'receptor__family', 'receptor__species',
            'publication', 'publication__web_link', 'publication__web_link__web_resource',
            'publication__journal', 'ligand__ref_ligand_bias_analyzed',
            'analyzed_data__emax_ligand_reference'
        ).annotate(
            # pathways
            pathways_p1=Subquery(assay_qs.values('family')[:1]),
            pathways_p2=Subquery(assay_qs.values('family')[1:2]),
            pathways_p3=Subquery(assay_qs.values('family')[2:3]),
            pathways_p4=Subquery(assay_qs.values('family')[3:4]),
            pathways_p5=Subquery(assay_qs.values('family')[4:5]),
            # t_factor
            opmodel_p2_p1=Subquery(assay_qs.values('t_factor')[1:2]),
            opmodel_p3_p1=Subquery(assay_qs.values('t_factor')[2:3]),
            opmodel_p4_p1=Subquery(assay_qs.values('t_factor')[3:4]),
            opmodel_p5_p1=Subquery(assay_qs.values('t_factor')[4:5]),
            # log bias factor
            lbf_p2_p1=Subquery(assay_qs.values('log_bias_factor')[1:2]),
            lbf_p3_p1=Subquery(assay_qs.values('log_bias_factor')[2:3]),
            lbf_p4_p1=Subquery(assay_qs.values('log_bias_factor')[3:4]),
            lbf_p5_p1=Subquery(assay_qs.values('log_bias_factor')[4:5]),
            # Potency ratio
            potency_p2_p1=Subquery(assay_qs.values('potency')[1:2]),
            potency_p3_p1=Subquery(assay_qs.values('potency')[2:3]),
            potency_p4_p1=Subquery(assay_qs.values('potency')[3:4]),
            potency_p5_p1=Subquery(assay_qs.values('potency')[4:5]),
            # Potency
            activity_p1=Subquery(assay_qs.values(
                'quantitive_activity_initial')[:1]),
            activity_p2=Subquery(assay_qs.values(
                'quantitive_activity_initial')[1:2]),
            activity_p3=Subquery(assay_qs.values(
                'quantitive_activity_initial')[2:3]),
            activity_p4=Subquery(assay_qs.values(
                'quantitive_activity_initial')[3:4]),
            activity_p5=Subquery(assay_qs.values(
                'quantitive_activity_initial')[4:5]),
            quality_activity_p1=Subquery(
                assay_qs.values('qualitative_activity')[:1]),
            quality_activity_p2=Subquery(
                assay_qs.values('qualitative_activity')[1:2]),
            quality_activity_p3=Subquery(
                assay_qs.values('qualitative_activity')[2:3]),
            quality_activity_p4=Subquery(
                assay_qs.values('qualitative_activity')[3:4]),
            quality_activity_p5=Subquery(
                assay_qs.values('qualitative_activity')[4:5]),
            # quality_activity
            standard_type_p1=Subquery(
                assay_qs.values('quantitive_measure_type')[:1]),
            standard_type_p2=Subquery(assay_qs.values(
                'quantitive_measure_type')[1:2]),
            standard_type_p3=Subquery(assay_qs.values(
                'quantitive_measure_type')[2:3]),
            standard_type_p4=Subquery(assay_qs.values(
                'quantitive_measure_type')[3:4]),
            standard_type_p5=Subquery(assay_qs.values(
                'quantitive_measure_type')[4:5]),
            # E Max
            emax_p1=Subquery(assay_qs.values('quantitive_efficacy')[:1]),
            emax_p2=Subquery(assay_qs.values('quantitive_efficacy')[1:2]),
            emax_p3=Subquery(assay_qs.values('quantitive_efficacy')[2:3]),
            emax_p4=Subquery(assay_qs.values('quantitive_efficacy')[3:4]),
            emax_p5=Subquery(assay_qs.values('quantitive_efficacy')[4:5]),

            lbf_part_p1=Subquery(assay_qs.values('log_bias_factor_a')[:1]),
            lbf_part_p2=Subquery(assay_qs.values('log_bias_factor_a')[1:2]),
            lbf_part_p3=Subquery(assay_qs.values('log_bias_factor_a')[2:3]),
            lbf_part_p4=Subquery(assay_qs.values('log_bias_factor_a')[3:4]),
            lbf_part_p5=Subquery(assay_qs.values('log_bias_factor_a')[4:5]),
            # reference assay
            reference_ligand_p1=Subquery(assay_qs.values('reference_ligand_id')[:1]),
            reference_ligand_p2=Subquery(assay_qs.values('reference_ligand_id')[1:2]),
            reference_ligand_p3=Subquery(assay_qs.values('reference_ligand_id')[2:3]),
            reference_ligand_p4=Subquery(assay_qs.values('reference_ligand_id')[3:4]),
            reference_ligand_p5=Subquery(assay_qs.values('reference_ligand_id')[4:5]),

            # T factor
            tfactor_p1=Subquery(assay_qs.values('t_value')[:1]),
            tfactor_p2=Subquery(assay_qs.values('t_value')[1:2]),
            tfactor_p3=Subquery(assay_qs.values('t_value')[2:3]),
            tfactor_p4=Subquery(assay_qs.values('t_value')[3:4]),
            tfactor_p5=Subquery(assay_qs.values('t_value')[4:5]),

            molecule1_p1=Subquery(assay_qs.values('molecule_1')[:1]),
            molecule1_p2=Subquery(assay_qs.values('molecule_1')[1:2]),
            molecule1_p3=Subquery(assay_qs.values('molecule_1')[2:3]),
            molecule1_p4=Subquery(assay_qs.values('molecule_1')[3:4]),
            molecule1_p5=Subquery(assay_qs.values('molecule_1')[4:5]),
            #molecule
            molecule2_p1=Subquery(assay_qs.values('molecule_2')[:1]),
            molecule2_p2=Subquery(assay_qs.values('molecule_2')[1:2]),
            molecule2_p3=Subquery(assay_qs.values('molecule_2')[2:3]),
            molecule2_p4=Subquery(assay_qs.values('molecule_2')[3:4]),
            molecule2_p5=Subquery(assay_qs.values('molecule_2')[4:5]),
            # Assay
            assay_p1=Subquery(assay_qs.values('signalling_protein')[:1]),
            assay_p2=Subquery(assay_qs.values('signalling_protein')[1:2]),
            assay_p3=Subquery(assay_qs.values('signalling_protein')[2:3]),
            assay_p4=Subquery(assay_qs.values('signalling_protein')[3:4]),
            assay_p5=Subquery(assay_qs.values('signalling_protein')[4:5]),

            # Cell Line
            cell_p1=Subquery(assay_qs.values('cell_line')[:1]),
            cell_p2=Subquery(assay_qs.values('cell_line')[1:2]),
            cell_p3=Subquery(assay_qs.values('cell_line')[2:3]),
            cell_p4=Subquery(assay_qs.values('cell_line')[3:4]),
            cell_p5=Subquery(assay_qs.values('cell_line')[4:5]),

            # Time
            time_p1=Subquery(assay_qs.values('assay_time_resolved')[:1]),
            time_p2=Subquery(assay_qs.values('assay_time_resolved')[1:2]),
            time_p3=Subquery(assay_qs.values('assay_time_resolved')[2:3]),
            time_p4=Subquery(assay_qs.values('assay_time_resolved')[3:4]),
            time_p5=Subquery(assay_qs.values('assay_time_resolved')[4:5]),

            measured_biological_process_p1=Subquery(assay_qs.values('measured_biological_process')[:1]),
            measured_biological_process_p2=Subquery(assay_qs.values('measured_biological_process')[1:2]),
            measured_biological_process_p3=Subquery(assay_qs.values('measured_biological_process')[2:3]),
            measured_biological_process_p4=Subquery(assay_qs.values('measured_biological_process')[3:4]),
            measured_biological_process_p5=Subquery(assay_qs.values('measured_biological_process')[4:5]),


            reference_quantitive_activity_initial_p1=Subquery(ref_assay_qs.values('quantitive_activity_initial')[:1]),
            reference_quantitive_activity_initial_p2=Subquery(ref_assay_qs.values('quantitive_activity_initial')[1:2]),
            reference_quantitive_activity_initial_p3=Subquery(ref_assay_qs.values('quantitive_activity_initial')[2:3]),
            reference_quantitive_activity_initial_p4=Subquery(ref_assay_qs.values('quantitive_activity_initial')[3:4]),
            reference_quantitive_activity_initial_p5=Subquery(ref_assay_qs.values('quantitive_activity_initial')[4:5]),

            reference_qualitative_activity_p1=Subquery(ref_assay_qs.values('qualitative_activity')[:1]),
            reference_qualitative_activity_p2=Subquery(ref_assay_qs.values('qualitative_activity')[1:2]),
            reference_qualitative_activity_p3=Subquery(ref_assay_qs.values('qualitative_activity')[2:3]),
            reference_qualitative_activity_p4=Subquery(ref_assay_qs.values('qualitative_activity')[3:4]),
            reference_qualitative_activity_p5=Subquery(ref_assay_qs.values('qualitative_activity')[4:5]),

            reference_quantitive_efficacy_p1=Subquery(ref_assay_qs.values('quantitive_efficacy')[:1]),
            reference_quantitive_efficacy_p2=Subquery(ref_assay_qs.values('quantitive_efficacy')[1:2]),
            reference_quantitive_efficacy_p3=Subquery(ref_assay_qs.values('quantitive_efficacy')[2:3]),
            reference_quantitive_efficacy_p4=Subquery(ref_assay_qs.values('quantitive_efficacy')[3:4]),
            reference_quantitive_efficacy_p5=Subquery(ref_assay_qs.values('quantitive_efficacy')[4:5]),

            reference_assay_type_p1=Subquery(ref_assay_qs.values('assay_type')[:1]),
            reference_assay_type_p2=Subquery(ref_assay_qs.values('assay_type')[1:2]),
            reference_assay_type_p3=Subquery(ref_assay_qs.values('assay_type')[2:3]),
            reference_assay_type_p4=Subquery(ref_assay_qs.values('assay_type')[3:4]),
            reference_assay_type_p5=Subquery(ref_assay_qs.values('assay_type')[4:5]),

            reference_a_p1=Subquery(assay_qs.values('log_bias_factor_a')[:1]),
            reference_a_p2=Subquery(assay_qs.values('log_bias_factor_a')[1:2]),
            reference_a_p3=Subquery(assay_qs.values('log_bias_factor_a')[2:3]),
            reference_a_p4=Subquery(assay_qs.values('log_bias_factor_a')[3:4]),
            reference_a_p5=Subquery(assay_qs.values('log_bias_factor_a')[4:5]),
        )
        return queryset


class BiasGBrowser(ListView):
    # serializer_class = AnalyzedExperimentSerializer
    template_name = 'bias_browser_subtypes.html'
    context_object_name = 'data_test'
    # @cache_page(50000)
    def get_queryset(self):
        protein_list = list()

        try:
            simple_selection = self.request.session.get('selection', False)
            a = Alignment()
            # load data from selection into the alignment
            a.load_proteins_from_selection(simple_selection)
            for items in a.proteins:
                protein_list.append(items.protein)
        except:
            protein_list.append(1)
        assay_qs = AnalyzedAssay.objects.filter(
            order_no__lte=5,
            assay_description='sub_tested_assays',

            experiment=OuterRef('pk'),
        ).order_by('order_no')

        ref_assay_qs = AnalyzedAssay.objects.filter(
            order_no__lte=5,
            assay_description='sub_endogenous',
            experiment=OuterRef('pk'),
        ).order_by('order_no')

        queryset = AnalyzedExperiment.objects.filter(
            source='sub_different_family',
            receptor__in=protein_list,
        ).prefetch_related(
            'analyzed_data', 'ligand', 'ligand__reference_ligand', 'reference_ligand',
            'endogenous_ligand', 'ligand__properities', 'receptor', 'receptor', 'receptor__family',
            'receptor__family__parent', 'receptor__family__parent__parent__parent',
            'receptor__family__parent__parent', 'receptor__family', 'receptor__species',
            'publication', 'publication__web_link', 'publication__web_link__web_resource',
            'publication__journal', 'ligand__ref_ligand_bias_analyzed',
            'analyzed_data__emax_ligand_reference'
        ).annotate(
            # pathways
            pathways_p1=Subquery(assay_qs.values('family')[:1]),
            pathways_p2=Subquery(assay_qs.values('family')[1:2]),
            pathways_p3=Subquery(assay_qs.values('family')[2:3]),
            pathways_p4=Subquery(assay_qs.values('family')[3:4]),
            pathways_p5=Subquery(assay_qs.values('family')[4:5]),
            # t_factor
            opmodel_p2_p1=Subquery(assay_qs.values('t_factor')[1:2]),
            opmodel_p3_p1=Subquery(assay_qs.values('t_factor')[2:3]),
            opmodel_p4_p1=Subquery(assay_qs.values('t_factor')[3:4]),
            opmodel_p5_p1=Subquery(assay_qs.values('t_factor')[4:5]),
            # log bias factor
            lbf_p2_p1=Subquery(assay_qs.values('log_bias_factor')[1:2]),
            lbf_p3_p1=Subquery(assay_qs.values('log_bias_factor')[2:3]),
            lbf_p4_p1=Subquery(assay_qs.values('log_bias_factor')[3:4]),
            lbf_p5_p1=Subquery(assay_qs.values('log_bias_factor')[4:5]),
            # Potency ratio
            potency_p2_p1=Subquery(assay_qs.values('potency')[1:2]),
            potency_p3_p1=Subquery(assay_qs.values('potency')[2:3]),
            potency_p4_p1=Subquery(assay_qs.values('potency')[3:4]),
            potency_p5_p1=Subquery(assay_qs.values('potency')[4:5]),
            # Potency
            activity_p1=Subquery(assay_qs.values(
                'quantitive_activity_initial')[:1]),
            activity_p2=Subquery(assay_qs.values(
                'quantitive_activity_initial')[1:2]),
            activity_p3=Subquery(assay_qs.values(
                'quantitive_activity_initial')[2:3]),
            activity_p4=Subquery(assay_qs.values(
                'quantitive_activity_initial')[3:4]),
            activity_p5=Subquery(assay_qs.values(
                'quantitive_activity_initial')[4:5]),
            quality_activity_p1=Subquery(
                assay_qs.values('qualitative_activity')[:1]),
            quality_activity_p2=Subquery(
                assay_qs.values('qualitative_activity')[1:2]),
            quality_activity_p3=Subquery(
                assay_qs.values('qualitative_activity')[2:3]),
            quality_activity_p4=Subquery(
                assay_qs.values('qualitative_activity')[3:4]),
            quality_activity_p5=Subquery(
                assay_qs.values('qualitative_activity')[4:5]),
            # quality_activity
            standard_type_p1=Subquery(
                assay_qs.values('quantitive_measure_type')[:1]),
            standard_type_p2=Subquery(assay_qs.values(
                'quantitive_measure_type')[1:2]),
            standard_type_p3=Subquery(assay_qs.values(
                'quantitive_measure_type')[2:3]),
            standard_type_p4=Subquery(assay_qs.values(
                'quantitive_measure_type')[3:4]),
            standard_type_p5=Subquery(assay_qs.values(
                'quantitive_measure_type')[4:5]),
            # E Max
            emax_p1=Subquery(assay_qs.values('quantitive_efficacy')[:1]),
            emax_p2=Subquery(assay_qs.values('quantitive_efficacy')[1:2]),
            emax_p3=Subquery(assay_qs.values('quantitive_efficacy')[2:3]),
            emax_p4=Subquery(assay_qs.values('quantitive_efficacy')[3:4]),
            emax_p5=Subquery(assay_qs.values('quantitive_efficacy')[4:5]),

            lbf_part_p1=Subquery(assay_qs.values('log_bias_factor_a')[:1]),
            lbf_part_p2=Subquery(assay_qs.values('log_bias_factor_a')[1:2]),
            lbf_part_p3=Subquery(assay_qs.values('log_bias_factor_a')[2:3]),
            lbf_part_p4=Subquery(assay_qs.values('log_bias_factor_a')[3:4]),
            lbf_part_p5=Subquery(assay_qs.values('log_bias_factor_a')[4:5]),
            # reference assay
            reference_ligand_p1=Subquery(assay_qs.values('reference_ligand_id')[:1]),
            reference_ligand_p2=Subquery(assay_qs.values('reference_ligand_id')[1:2]),
            reference_ligand_p3=Subquery(assay_qs.values('reference_ligand_id')[2:3]),
            reference_ligand_p4=Subquery(assay_qs.values('reference_ligand_id')[3:4]),
            reference_ligand_p5=Subquery(assay_qs.values('reference_ligand_id')[4:5]),

            # T factor
            tfactor_p1=Subquery(assay_qs.values('t_value')[:1]),
            tfactor_p2=Subquery(assay_qs.values('t_value')[1:2]),
            tfactor_p3=Subquery(assay_qs.values('t_value')[2:3]),
            tfactor_p4=Subquery(assay_qs.values('t_value')[3:4]),
            tfactor_p5=Subquery(assay_qs.values('t_value')[4:5]),

            molecule1_p1=Subquery(assay_qs.values('molecule_1')[:1]),
            molecule1_p2=Subquery(assay_qs.values('molecule_1')[1:2]),
            molecule1_p3=Subquery(assay_qs.values('molecule_1')[2:3]),
            molecule1_p4=Subquery(assay_qs.values('molecule_1')[3:4]),
            molecule1_p5=Subquery(assay_qs.values('molecule_1')[4:5]),
            #molecule
            molecule2_p1=Subquery(assay_qs.values('molecule_2')[:1]),
            molecule2_p2=Subquery(assay_qs.values('molecule_2')[1:2]),
            molecule2_p3=Subquery(assay_qs.values('molecule_2')[2:3]),
            molecule2_p4=Subquery(assay_qs.values('molecule_2')[3:4]),
            molecule2_p5=Subquery(assay_qs.values('molecule_2')[4:5]),
            # signalling protein
            assay_p1=Subquery(assay_qs.values('signalling_protein')[:1]),
            assay_p2=Subquery(assay_qs.values('signalling_protein')[1:2]),
            assay_p3=Subquery(assay_qs.values('signalling_protein')[2:3]),
            assay_p4=Subquery(assay_qs.values('signalling_protein')[3:4]),
            assay_p5=Subquery(assay_qs.values('signalling_protein')[4:5]),

            # Cell Line
            cell_p1=Subquery(assay_qs.values('cell_line')[:1]),
            cell_p2=Subquery(assay_qs.values('cell_line')[1:2]),
            cell_p3=Subquery(assay_qs.values('cell_line')[2:3]),
            cell_p4=Subquery(assay_qs.values('cell_line')[3:4]),
            cell_p5=Subquery(assay_qs.values('cell_line')[4:5]),

            # Time
            time_p1=Subquery(assay_qs.values('assay_time_resolved')[:1]),
            time_p2=Subquery(assay_qs.values('assay_time_resolved')[1:2]),
            time_p3=Subquery(assay_qs.values('assay_time_resolved')[2:3]),
            time_p4=Subquery(assay_qs.values('assay_time_resolved')[3:4]),
            time_p5=Subquery(assay_qs.values('assay_time_resolved')[4:5]),

            measured_biological_process_p1=Subquery(assay_qs.values('measured_biological_process')[:1]),
            measured_biological_process_p2=Subquery(assay_qs.values('measured_biological_process')[1:2]),
            measured_biological_process_p3=Subquery(assay_qs.values('measured_biological_process')[2:3]),
            measured_biological_process_p4=Subquery(assay_qs.values('measured_biological_process')[3:4]),
            measured_biological_process_p5=Subquery(assay_qs.values('measured_biological_process')[4:5]),


            reference_quantitive_activity_initial_p1=Subquery(ref_assay_qs.values('quantitive_activity_initial')[:1]),
            reference_quantitive_activity_initial_p2=Subquery(ref_assay_qs.values('quantitive_activity_initial')[1:2]),
            reference_quantitive_activity_initial_p3=Subquery(ref_assay_qs.values('quantitive_activity_initial')[2:3]),
            reference_quantitive_activity_initial_p4=Subquery(ref_assay_qs.values('quantitive_activity_initial')[3:4]),
            reference_quantitive_activity_initial_p5=Subquery(ref_assay_qs.values('quantitive_activity_initial')[4:5]),

            reference_qualitative_activity_p1=Subquery(ref_assay_qs.values('qualitative_activity')[:1]),
            reference_qualitative_activity_p2=Subquery(ref_assay_qs.values('qualitative_activity')[1:2]),
            reference_qualitative_activity_p3=Subquery(ref_assay_qs.values('qualitative_activity')[2:3]),
            reference_qualitative_activity_p4=Subquery(ref_assay_qs.values('qualitative_activity')[3:4]),
            reference_qualitative_activity_p5=Subquery(ref_assay_qs.values('qualitative_activity')[4:5]),

            reference_quantitive_efficacy_p1=Subquery(ref_assay_qs.values('quantitive_efficacy')[:1]),
            reference_quantitive_efficacy_p2=Subquery(ref_assay_qs.values('quantitive_efficacy')[1:2]),
            reference_quantitive_efficacy_p3=Subquery(ref_assay_qs.values('quantitive_efficacy')[2:3]),
            reference_quantitive_efficacy_p4=Subquery(ref_assay_qs.values('quantitive_efficacy')[3:4]),
            reference_quantitive_efficacy_p5=Subquery(ref_assay_qs.values('quantitive_efficacy')[4:5]),

            reference_assay_type_p1=Subquery(ref_assay_qs.values('assay_type')[:1]),
            reference_assay_type_p2=Subquery(ref_assay_qs.values('assay_type')[1:2]),
            reference_assay_type_p3=Subquery(ref_assay_qs.values('assay_type')[2:3]),
            reference_assay_type_p4=Subquery(ref_assay_qs.values('assay_type')[3:4]),
            reference_assay_type_p5=Subquery(ref_assay_qs.values('assay_type')[4:5]),

            reference_a_p1=Subquery(assay_qs.values('log_bias_factor_a')[:1]),
            reference_a_p2=Subquery(assay_qs.values('log_bias_factor_a')[1:2]),
            reference_a_p3=Subquery(assay_qs.values('log_bias_factor_a')[2:3]),
            reference_a_p4=Subquery(assay_qs.values('log_bias_factor_a')[3:4]),
            reference_a_p5=Subquery(assay_qs.values('log_bias_factor_a')[4:5]),
        )
        return queryset


class BiasPredictionBrowser(ListView):
    # serializer_class = AnalyzedExperimentSerializer
    template_name = 'bias_browser_predict.html'
    context_object_name = 'data_test'
    # @cache_page(50000)
    def get_queryset(self):
        protein_list = list()

        try:
            simple_selection = self.request.session.get('selection', False)
            a = Alignment()
            # load data from selection into the alignment
            a.load_proteins_from_selection(simple_selection)
            for items in a.proteins:
                protein_list.append(items.protein)
        except:
            protein_list.append(1)
        assay_qs = AnalyzedAssay.objects.filter(
            order_no__lte=5,
            assay_description='predicted_tested_assays',

            experiment=OuterRef('pk'),
        ).order_by('order_no')

        queryset = AnalyzedExperiment.objects.filter(
            source='predicted_family',
            receptor__in=protein_list,
        ).prefetch_related(
            'analyzed_data', 'ligand', 'ligand__reference_ligand', 'reference_ligand',
            'endogenous_ligand', 'ligand__properities', 'receptor', 'receptor', 'receptor__family',
            'receptor__family__parent', 'receptor__family__parent__parent__parent',
            'receptor__family__parent__parent', 'receptor__family', 'receptor__species',
            'publication', 'publication__web_link', 'publication__web_link__web_resource',
            'publication__journal', 'ligand__ref_ligand_bias_analyzed',
            'analyzed_data__emax_ligand_reference'
        ).annotate(
            # pathways
            pathways_p1=Subquery(assay_qs.values('family')[:1]),
            pathways_p2=Subquery(assay_qs.values('family')[1:2]),
            pathways_p3=Subquery(assay_qs.values('family')[2:3]),
            pathways_p4=Subquery(assay_qs.values('family')[3:4]),
            pathways_p5=Subquery(assay_qs.values('family')[4:5]),
            # t_factor
            opmodel_p2_p1=Subquery(assay_qs.values('t_factor')[1:2]),
            opmodel_p3_p1=Subquery(assay_qs.values('t_factor')[2:3]),
            opmodel_p4_p1=Subquery(assay_qs.values('t_factor')[3:4]),
            opmodel_p5_p1=Subquery(assay_qs.values('t_factor')[4:5]),
            # log bias factor
            lbf_p2_p1=Subquery(assay_qs.values('log_bias_factor')[1:2]),
            lbf_p3_p1=Subquery(assay_qs.values('log_bias_factor')[2:3]),
            lbf_p4_p1=Subquery(assay_qs.values('log_bias_factor')[3:4]),
            lbf_p5_p1=Subquery(assay_qs.values('log_bias_factor')[4:5]),
            # Potency ratio
            potency_p2_p1=Subquery(assay_qs.values('potency')[1:2]),
            potency_p3_p1=Subquery(assay_qs.values('potency')[2:3]),
            potency_p4_p1=Subquery(assay_qs.values('potency')[3:4]),
            potency_p5_p1=Subquery(assay_qs.values('potency')[4:5]),
            # Potency
            activity_p1=Subquery(assay_qs.values(
                'quantitive_activity_initial')[:1]),
            activity_p2=Subquery(assay_qs.values(
                'quantitive_activity_initial')[1:2]),
            activity_p3=Subquery(assay_qs.values(
                'quantitive_activity_initial')[2:3]),
            activity_p4=Subquery(assay_qs.values(
                'quantitive_activity_initial')[3:4]),
            activity_p5=Subquery(assay_qs.values(
                'quantitive_activity_initial')[4:5]),
            quality_activity_p1=Subquery(
                assay_qs.values('qualitative_activity')[:1]),
            quality_activity_p2=Subquery(
                assay_qs.values('qualitative_activity')[1:2]),
            quality_activity_p3=Subquery(
                assay_qs.values('qualitative_activity')[2:3]),
            quality_activity_p4=Subquery(
                assay_qs.values('qualitative_activity')[3:4]),
            quality_activity_p5=Subquery(
                assay_qs.values('qualitative_activity')[4:5]),
            # quality_activity
            standard_type_p1=Subquery(
                assay_qs.values('quantitive_measure_type')[:1]),
            standard_type_p2=Subquery(assay_qs.values(
                'quantitive_measure_type')[1:2]),
            standard_type_p3=Subquery(assay_qs.values(
                'quantitive_measure_type')[2:3]),
            standard_type_p4=Subquery(assay_qs.values(
                'quantitive_measure_type')[3:4]),
            standard_type_p5=Subquery(assay_qs.values(
                'quantitive_measure_type')[4:5]),
            # E Max
            emax_p1=Subquery(assay_qs.values('quantitive_efficacy')[:1]),
            emax_p2=Subquery(assay_qs.values('quantitive_efficacy')[1:2]),
            emax_p3=Subquery(assay_qs.values('quantitive_efficacy')[2:3]),
            emax_p4=Subquery(assay_qs.values('quantitive_efficacy')[3:4]),
            emax_p5=Subquery(assay_qs.values('quantitive_efficacy')[4:5]),

            lbf_part_p1=Subquery(assay_qs.values('log_bias_factor_a')[:1]),
            lbf_part_p2=Subquery(assay_qs.values('log_bias_factor_a')[1:2]),
            lbf_part_p3=Subquery(assay_qs.values('log_bias_factor_a')[2:3]),
            lbf_part_p4=Subquery(assay_qs.values('log_bias_factor_a')[3:4]),
            lbf_part_p5=Subquery(assay_qs.values('log_bias_factor_a')[4:5]),
            # reference assay
            reference_ligand_p1=Subquery(assay_qs.values('reference_ligand_id')[:1]),
            reference_ligand_p2=Subquery(assay_qs.values('reference_ligand_id')[1:2]),
            reference_ligand_p3=Subquery(assay_qs.values('reference_ligand_id')[2:3]),
            reference_ligand_p4=Subquery(assay_qs.values('reference_ligand_id')[3:4]),
            reference_ligand_p5=Subquery(assay_qs.values('reference_ligand_id')[4:5]),

            # T factor
            tfactor_p1=Subquery(assay_qs.values('t_value')[:1]),
            tfactor_p2=Subquery(assay_qs.values('t_value')[1:2]),
            tfactor_p3=Subquery(assay_qs.values('t_value')[2:3]),
            tfactor_p4=Subquery(assay_qs.values('t_value')[3:4]),
            tfactor_p5=Subquery(assay_qs.values('t_value')[4:5]),

            molecule1_p1=Subquery(assay_qs.values('molecule_1')[:1]),
            molecule1_p2=Subquery(assay_qs.values('molecule_1')[1:2]),
            molecule1_p3=Subquery(assay_qs.values('molecule_1')[2:3]),
            molecule1_p4=Subquery(assay_qs.values('molecule_1')[3:4]),
            molecule1_p5=Subquery(assay_qs.values('molecule_1')[4:5]),
            #molecule
            molecule2_p1=Subquery(assay_qs.values('molecule_2')[:1]),
            molecule2_p2=Subquery(assay_qs.values('molecule_2')[1:2]),
            molecule2_p3=Subquery(assay_qs.values('molecule_2')[2:3]),
            molecule2_p4=Subquery(assay_qs.values('molecule_2')[3:4]),
            molecule2_p5=Subquery(assay_qs.values('molecule_2')[4:5]),
            # signalling protein
            assay_p1=Subquery(assay_qs.values('signalling_protein')[:1]),
            assay_p2=Subquery(assay_qs.values('signalling_protein')[1:2]),
            assay_p3=Subquery(assay_qs.values('signalling_protein')[2:3]),
            assay_p4=Subquery(assay_qs.values('signalling_protein')[3:4]),
            assay_p5=Subquery(assay_qs.values('signalling_protein')[4:5]),

            # Cell Line
            cell_p1=Subquery(assay_qs.values('cell_line')[:1]),
            cell_p2=Subquery(assay_qs.values('cell_line')[1:2]),
            cell_p3=Subquery(assay_qs.values('cell_line')[2:3]),
            cell_p4=Subquery(assay_qs.values('cell_line')[3:4]),
            cell_p5=Subquery(assay_qs.values('cell_line')[4:5]),

            # Time
            time_p1=Subquery(assay_qs.values('assay_time_resolved')[:1]),
            time_p2=Subquery(assay_qs.values('assay_time_resolved')[1:2]),
            time_p3=Subquery(assay_qs.values('assay_time_resolved')[2:3]),
            time_p4=Subquery(assay_qs.values('assay_time_resolved')[3:4]),
            time_p5=Subquery(assay_qs.values('assay_time_resolved')[4:5]),

            measured_biological_process_p1=Subquery(assay_qs.values('measured_biological_process')[:1]),
            measured_biological_process_p2=Subquery(assay_qs.values('measured_biological_process')[1:2]),
            measured_biological_process_p3=Subquery(assay_qs.values('measured_biological_process')[2:3]),
            measured_biological_process_p4=Subquery(assay_qs.values('measured_biological_process')[3:4]),
            measured_biological_process_p5=Subquery(assay_qs.values('measured_biological_process')[4:5]),

            reference_a_p1=Subquery(assay_qs.values('log_bias_factor_a')[:1]),
            reference_a_p2=Subquery(assay_qs.values('log_bias_factor_a')[1:2]),
            reference_a_p3=Subquery(assay_qs.values('log_bias_factor_a')[2:3]),
            reference_a_p4=Subquery(assay_qs.values('log_bias_factor_a')[3:4]),
            reference_a_p5=Subquery(assay_qs.values('log_bias_factor_a')[4:5]),
        )
        return queryset

class BiasGuidelines(TemplateView):

    template_name = 'bias_guidelines.html'

    def get_context_data(self, **kwargs):

        context = super().get_context_data(**kwargs)
        # assays = AssayExperiment.objects.all().prefetch_related('protein__family__parent__parent__parent', 'protein__family')
        return context
