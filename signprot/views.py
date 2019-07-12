from django.shortcuts import get_object_or_404, render
from django.http import HttpResponse
from django.http import JsonResponse
from django.core.cache import cache
from django.views.decorators.cache import cache_page

from protein.models import Protein, ProteinConformation, ProteinAlias, ProteinSegment, ProteinFamily, Gene, ProteinGProtein, ProteinGProteinPair
from residue.models import Residue, ResiduePositionSet, ResidueGenericNumberEquivalent

from structure.models import Structure
from mutation.models import MutationExperiment
from common.selection import Selection
from common.diagrams_gpcr import DrawSnakePlot
from common.diagrams_gprotein import DrawGproteinPlot
from common.diagrams_arrestin import DrawArrestinPlot
from common.definitions import AMINO_ACIDS, AMINO_ACID_GROUPS, AMINO_ACID_GROUP_NAMES, AMINO_ACID_GROUP_PROPERTIES

from seqsign.sequence_signature import SignatureMatch
from seqsign.sequence_signature import SequenceSignature
from signprot.models import SignprotStructure, SignprotBarcode, SignprotInteractions, SignprotComplex
from signprot.interactions import (
    get_entry_names,
    get_ignore_info,
    get_protein_segments,
    get_generic_numbers,
    get_signature_features,
    group_signature_features,
    get_signature_consensus,
    prepare_signature_match,
)

from common import definitions
from collections import OrderedDict
from collections import Counter
from common.views import AbsTargetSelection

import json
import re
import time
import pickle
from itertools import chain

from django.core.exceptions import ObjectDoesNotExist
from decimal import Decimal


class BrowseSelection(AbsTargetSelection):
    step = 1
    number_of_steps = 1
    psets = False
    filters = True
    filter_gprotein = True

    type_of_selection = 'browse_gprot'

    description = 'Select a G protein or family by searching or browsing in the right column.'
    description = 'Select a G protein (family) by searching or browsing in the middle. The selection is viewed to' \
        + ' the right.'
    docs = 'receptors.html'
    target_input=False

    selection_boxes = OrderedDict([
        ('reference', False), ('targets', True),
        ('segments', False),
    ])
    try:
        ppf_g = ProteinFamily.objects.get(slug="100_001")
        # ppf_a = ProteinFamily.objects.get(slug="200_000")
        # pfs = ProteinFamily.objects.filter(parent__in=[ppf_g.id,ppf_a.id])
        pfs = ProteinFamily.objects.filter(parent__in=[ppf_g.id])
        ps = Protein.objects.filter(family__in=[ppf_g]) # ,ppf_a
        tree_indent_level = []
        # action = 'expand'
        # remove the parent family (for all other families than the root of the tree, the parent should be shown)
        # del ppf_g
        # del ppf_a
    except Exception as e:
        pass

#@cache_page(60*60*24*2) # 2 days caching
def GProtein(request, dataset = "GuideToPharma"):

    name_of_cache = 'gprotein_statistics_{}'.format(dataset)

    context = cache.get(name_of_cache)
    if context==None:

        context = OrderedDict()
        i=0
        gproteins = ProteinGProtein.objects.all().prefetch_related('proteingproteinpair_set')
        slugs = ['001','002','004','005']
        slug_translate = {'001':"ClassA", '002':"ClassB1",'004':"ClassC", '005':"ClassF"}
        selectivitydata = {}
        for slug in slugs:
            jsondata = {}
            for gp in gproteins:
                # ps = gp.proteingproteinpair_set.all()
                ps = gp.proteingproteinpair_set.filter(protein__family__slug__startswith=slug, source=dataset).prefetch_related('protein')
                # print(ps,len(ps))
                if ps:
                    jsondata[str(gp)] = []
                    for p in ps:
                        if dataset=="Aska" and p.log_rai_mean<-1:
                            continue
                        if str(p.protein.entry_name).split('_')[0].upper() not in selectivitydata:
                            selectivitydata[str(p.protein.entry_name).split('_')[0].upper()] = []
                        selectivitydata[str(p.protein.entry_name).split('_')[0].upper()].append(str(gp))
                        # print(p.protein.family.parent.parent.parent)
                        jsondata[str(gp)].append(str(p.protein.entry_name)+'\n')

                    jsondata[str(gp)] = ''.join(jsondata[str(gp)])

            context[slug_translate[slug]] = jsondata

        context["selectivitydata"] = selectivitydata


    cache.set(name_of_cache, context, 60*60*24*2) #two days timeout on cache

    return render(request, 'signprot/gprotein.html', context)

#@cache_page(60*60*24*2) # 2 days caching
def Couplings(request):

    context = OrderedDict()

    threshold_primary = -0.1
    threshold_secondary = -1


    proteins = Protein.objects.filter(sequence_type__slug='wt',family__slug__startswith='00',species__common_name='Human').all().prefetch_related('family')
    data = {}
    class_names = {}
    for p in proteins:
        p_class = p.family.slug.split('_')[0]
        if p_class not in class_names:
            class_names[p_class] =  re.sub(r'\([^)]*\)', '', p.family.parent.parent.parent.name)
        p_class_name = class_names[p_class].strip()
        data[p.entry_short()] = {'class':p_class_name,'pretty':p.short()[:15],'GuideToPharma':{},'Aska':{}}

    distinct_g_families = []
    distinct_g_subunit_families = {}
    distinct_sources = ['GuideToPharma','Aska']

    couplings = ProteinGProteinPair.objects.all().prefetch_related('protein','g_protein_subunit','g_protein')
    for c in couplings:
        p = c.protein.entry_short()
        s = c.source
        t = c.transduction
        m = c.log_rai_mean
        gf = c.g_protein.name
        # print(gf)
        gf = gf.replace(" family","")

        if gf not in distinct_g_families:
            distinct_g_families.append(gf)
            distinct_g_subunit_families[gf] = []

        if c.g_protein_subunit:
            g = c.g_protein_subunit.entry_name
            g = g.replace("_human","")
            # print("g",g)
            if g not in distinct_g_subunit_families[gf]:
                distinct_g_subunit_families[gf].append(g)
                distinct_g_subunit_families[gf] = sorted(distinct_g_subunit_families[gf])

        if s not in data[p]:
            data[p][s] = {}

        if gf not in data[p][s]:
            data[p][s][gf] = {}

        # If transduction in GuideToPharma data
        if t:
            data[p][s][gf] = t
        else:
            if 'subunits' not in data[p][s][gf]:
                data[p][s][gf] = {'subunits':{},'best':-2.00}
            data[p][s][gf]['subunits'][g] = round(Decimal(m),2)
            if round(Decimal(m),2)== -0.00:
                data[p][s][gf]['subunits'][g] = 0.00
            # get the lowest number into 'best'
            if m>data[p][s][gf]['best']:
                data[p][s][gf]['best'] = round(Decimal(m),2)

    fd = {} #final data

    distinct_g_families = sorted(distinct_g_families)
    distinct_g_families = ['Gs','Gi/Go', 'Gq/G11', 'G12/G13', ]
    distinct_g_subunit_families = OrderedDict([('Gs',['gnas2','gnal']), ('Gi/Go',['gnai1', 'gnai3', 'gnao', 'gnaz']), ('Gq/G11',['gnaq', 'gna14', 'gna15']), ('G12/G13',['gna12', 'gna13'])])

    for p,v in data.items():
        fd[p] = [v['class'],p,v['pretty']]

        s = 'GuideToPharma'
        #Merge
        for gf in distinct_g_families:
            values = []
            if 'GuideToPharma' in v and gf in v['GuideToPharma']:
                values.append(v['GuideToPharma'][gf])
            if 'Aska' in v and gf in v['Aska']:
                best = v['Aska'][gf]['best']
                if best > threshold_primary:
                    values.append('primary')
                elif best > threshold_secondary:
                    values.append('secondary')
            if 'primary' in values:
                fd[p].append('primary')
            elif 'secondary' in values:
                fd[p].append('secondary')
            else:
                fd[p].append('')

        s = 'GuideToPharma'
        #First loop over GuideToPharma
        for gf in distinct_g_families:
            if gf in v[s]:
                fd[p].append(v[s][gf])
            else:
                fd[p].append("")

        s = 'Aska'
        for gf in distinct_g_families:
            if gf in v[s]:
                if v[s][gf]['best']>threshold_primary:
                    fd[p].append("primary")
                elif v[s][gf]['best']>threshold_secondary:
                    fd[p].append("secondary")
                else:
                    fd[p].append("No coupling")
            else:
                fd[p].append("")

        for gf,sfs in distinct_g_subunit_families.items():
            for sf in sfs:
                if gf in v[s]:
                    if sf in v[s][gf]['subunits']:
                        fd[p].append(v[s][gf]['subunits'][sf])
                    else:
                        fd[p].append("")
                else:
                    fd[p].append("")


    context['data'] = fd
    context['distinct_gf'] = distinct_g_families
    context['distinct_sf'] = distinct_g_subunit_families

    return render(request, 'signprot/browser.html', context)

@cache_page(60*60*24*2)
def familyDetail(request, slug):
    # get family
    pf = ProteinFamily.objects.get(slug=slug)

    # get family list
    ppf = pf
    families = [ppf.name]
    while ppf.parent.parent:
        families.append(ppf.parent.name)
        ppf = ppf.parent
    families.reverse()

    # number of proteins
    proteins = Protein.objects.filter(family__slug__startswith=pf.slug, sequence_type__slug='wt')
    no_of_proteins = proteins.count()
    no_of_human_proteins = Protein.objects.filter(family__slug__startswith=pf.slug, species__id=1,
        sequence_type__slug='wt').count()
    list_proteins = list(proteins.values_list('pk',flat=True))

    # get structures of this family
    structures = SignprotStructure.objects.filter(protein__family__slug__startswith=slug
        )

    mutations = MutationExperiment.objects.filter(protein__in=proteins).prefetch_related('residue__generic_number', 'exp_qual', 'ligand')

    mutations_list = {}
    for mutation in mutations:
        if not mutation.residue.generic_number: continue #cant map those without display numbers
        if mutation.residue.generic_number.label not in mutations_list: mutations_list[mutation.residue.generic_number.label] = []
        if mutation.ligand:
            ligand = mutation.ligand.name
        else:
            ligand = ''
        if mutation.exp_qual:
            qual = mutation.exp_qual.qual
        else:
            qual = ''
        mutations_list[mutation.residue.generic_number.label].append([mutation.foldchange,ligand.replace("'", "\\'"),qual])

    # Update to consensus sequence in protein confirmation!
    try:
        pc = ProteinConformation.objects.filter(protein__family__slug=slug, protein__sequence_type__slug='consensus')
    except ProteinConformation.DoesNotExist:
        pc = ProteinConformation.objects.get(protein__family__slug=slug, protein__species_id=1,
            protein__sequence_type__slug='wt')

    residues = Residue.objects.filter(protein_conformation=pc).order_by('sequence_number').prefetch_related(
        'protein_segment', 'generic_number', 'display_generic_number')

    jsondata = {}
    jsondata_interaction = {}
    for r in residues:
        if r.generic_number:
            if r.generic_number.label in mutations_list:
                jsondata[r.sequence_number] = [mutations_list[r.generic_number.label]]
            if r.generic_number.label in interaction_list:
                jsondata_interaction[r.sequence_number] = interaction_list[r.generic_number.label]

    # process residues and return them in chunks of 10
    # this is done for easier scaling on smaller screens
    chunk_size = 10
    r_chunks = []
    r_buffer = []
    last_segment = False
    border = False
    title_cell_skip = 0
    for i, r in enumerate(residues):
        # title of segment to be written out for the first residue in each segment
        segment_title = False

        # keep track of last residues segment (for marking borders)
        if r.protein_segment.slug != last_segment:
            last_segment = r.protein_segment.slug
            border = True

        # if on a border, is there room to write out the title? If not, write title in next chunk
        if i == 0 or (border and len(last_segment) <= (chunk_size - i % chunk_size)):
            segment_title = True
            border = False
            title_cell_skip += len(last_segment) # skip cells following title (which has colspan > 1)

        if i and i % chunk_size == 0:
            r_chunks.append(r_buffer)
            r_buffer = []

        r_buffer.append((r, segment_title, title_cell_skip))

        # update cell skip counter
        if title_cell_skip > 0:
            title_cell_skip -= 1
    if r_buffer:
        r_chunks.append(r_buffer)

    context = {'pf': pf, 'families': families, 'structures': structures, 'no_of_proteins': no_of_proteins,
        'no_of_human_proteins': no_of_human_proteins, 'mutations':mutations, 'r_chunks': r_chunks, 'chunk_size': chunk_size}

    return render(request, 'signprot/family_details.html', context)

class TargetSelection(AbsTargetSelection):
    step = 1
    number_of_steps = 1
    filters = False
    psets = False
    target_input = False
    redirect_on_select = True
    type_of_selection = 'ginterface'
    title = 'SELECT TARGET for Gs INTERFACE'
    description = 'Select a reference target by searching or browsing.' \
        + '\n\nThe Gs interface from adrb2 (PDB: 3SN6) will be superposed onto the selected target.' \
        + '\n\nAn interaction browser for the adrb2 Gs interface will be given for comparison"'

    # template_name = 'common/targetselection.html'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])

    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '#',
            'color': 'success',
        },
    }

@cache_page(60*60*24*2)
def Ginterface(request, protein = None):

    residuelist = Residue.objects.filter(protein_conformation__protein__entry_name=protein).prefetch_related('protein_segment','display_generic_number','generic_number')
    SnakePlot = DrawSnakePlot(
                residuelist, "Class A (Rhodopsin)", protein, nobuttons=1)

    # TEST
    gprotein_residues = Residue.objects.filter(protein_conformation__protein__entry_name='gnaz_human').prefetch_related('protein_segment','display_generic_number','generic_number')
    gproteinplot = DrawGproteinPlot(
                gprotein_residues, "Gprotein", protein)

    crystal = Structure.objects.get(pdb_code__index="3SN6")
    aa_names = definitions.AMINO_ACID_GROUP_NAMES_OLD
    names_aa = dict(zip(aa_names.values(),aa_names.keys()))
    names_aa['Polar (S/T)'] = 'pol_short'
    names_aa['Polar (N/Q/H)'] = 'pol_long'

    residues_browser = [{'pos': 135, 'aa': 'I', 'gprotseg': "H5",'segment': 'TM3', 'ligand': 'Gs', 'type': aa_names['hp'], 'gpcrdb': '3.54x54', 'gpnum': 'G.H5.16', 'gpaa': 'Q384', 'availability': 'interacting'},{'pos': 136, 'aa': 'T', 'gprotseg': "H5",'segment': 'TM3', 'ligand': 'Gs', 'type': 'Polar (S/T)', 'gpcrdb': '3.55x55', 'gpnum': 'G.H5.12', 'gpaa': 'R380', 'availability': 'interacting'},{'pos': 139, 'aa': 'F', 'gprotseg': "H5",'segment': 'ICL2', 'ligand': 'Gs', 'type': 'Aromatic', 'gpcrdb': '34.51x51', 'gpnum': 'G.H5.8', 'gpaa': 'F376', 'availability': 'interacting'},{'pos': 139, 'aa': 'F', 'gprotseg': "S1",'segment': 'ICL2', 'ligand': 'Gs', 'type': 'Aromatic', 'gpcrdb': '34.51x51', 'gpnum': 'G.S1.2', 'gpaa': 'H41', 'availability': 'interacting'},{'pos': 141, 'aa': 'Y', 'gprotseg': "H5",'segment': 'ICL2', 'ligand': 'Gs', 'type': 'Aromatic', 'gpcrdb': '34.53x53', 'gpnum': 'G.H5.19', 'gpaa': 'H387', 'availability': 'interacting'},{'pos': 225, 'aa': 'E', 'gprotseg': "H5",'segment': 'TM5', 'ligand': 'Gs', 'type': 'Negative charge', 'gpcrdb': '5.64x64', 'gpnum': 'G.H5.12', 'gpaa': 'R380', 'availability': 'interacting'},{'pos': 225, 'aa': 'E', 'gprotseg': "H5",'segment': 'TM5', 'ligand': 'Gs', 'type': 'Negative charge', 'gpcrdb': '5.64x64', 'gpnum': 'G.H5.16', 'gpaa': 'Q384', 'availability': 'interacting'},{'pos': 229, 'aa': 'Q', 'gprotseg': "H5",'segment': 'TM5', 'ligand': 'Gs', 'type': 'Polar (N/Q/H)', 'gpcrdb': '5.68x68', 'gpnum': 'G.H5.13', 'gpaa': 'D381', 'availability': 'interacting'},{'pos': 229, 'aa': 'Q', 'gprotseg': "H5",'segment': 'TM5', 'ligand': 'Gs', 'type': 'Polar (N/Q/H)', 'gpcrdb': '5.68x68', 'gpnum': 'G.H5.16', 'gpaa': 'Q384', 'availability': 'interacting'},{'pos': 229, 'aa': 'Q', 'gprotseg': "H5",'segment': 'TM5', 'ligand': 'Gs', 'type': 'Polar (N/Q/H)', 'gpcrdb': '5.68x68', 'gpnum': 'G.H5.17', 'gpaa': 'R385', 'availability': 'interacting'},{'pos': 274, 'aa': 'T', 'gprotseg': "H5",'segment': 'TM6', 'ligand': 'Gs', 'type': 'Polar (S/T)', 'gpcrdb': '6.36x36', 'gpnum': 'G.H5.24', 'gpaa': 'E392', 'availability': 'interacting'},{'pos': 328, 'aa': 'R', 'gprotseg': "H5",'segment': 'TM7', 'ligand': 'Gs', 'type': 'Positive charge', 'gpcrdb': '7.55x55', 'gpnum': 'G.H5.24', 'gpaa': 'E392', 'availability': 'interacting'}, {'pos': 232, 'aa': 'K', 'segment': 'TM5', 'ligand': 'Gs', 'type': 'Positive charge', 'gpcrdb': '5.71x71', 'gprotseg': "H5", 'gpnum': 'G.H5.13', 'gpaa': 'D381', 'availability': 'interacting'}]

    # accessible_gn = ['3.50x50', '3.53x53', '3.54x54', '3.55x55', '34.50x50', '34.51x51', '34.53x53', '34.54x54', '5.61x61', '5.64x64', '5.65x65', '5.67x67', '5.68x68', '5.71x71', '5.72x72', '5.74x74', '5.75x75', '6.29x29', '6.32x32', '6.33x33', '6.36x36', '6.37x37', '7.55x55', '8.48x48', '8.49x49']

    accessible_gn = ['3.50x50', '3.53x53', '3.54x54', '3.55x55', '3.56x56', '34.50x50', '34.51x51', '34.52x52', '34.53x53', '34.54x54', '34.55x55', '34.56x56', '34.57x57', '5.61x61', '5.64x64', '5.65x65', '5.66x66', '5.67x67', '5.68x68', '5.69x69', '5.71x71', '5.72x72', '5.74x74', '5.75x75', '6.25x25', '6.26x26', '6.28x28', '6.29x29', '6.32x32', '6.33x33', '6.36x36', '6.37x37', '6.40x40', '7.55x55', '7.56x56', '8.47x47', '8.48x48', '8.49x49', '8.51x51']

    exchange_table = OrderedDict([('hp', ('V','I', 'L', 'M')),
                                 ('ar', ('F', 'H', 'W', 'Y')),
                                 ('pol_short', ('S', 'T')), # Short/hydroxy
                                 ('pol_long', ('N', 'Q', 'H')), # Amino-like (both donor and acceptor
                                 ('neg', ('D', 'E')),
                                 ('pos', ('K', 'R'))])

    interacting_gn = []

    accessible_pos = list(residuelist.filter(display_generic_number__label__in=accessible_gn).values_list('sequence_number', flat=True))

    # Which of the Gs interacting_pos are conserved?
    GS_none_equivalent_interacting_pos = []
    GS_none_equivalent_interacting_gn = []

    for interaction in residues_browser:
        interacting_gn.append(interaction['gpcrdb'])
        gs_b2_interaction_type_long = (next((item['type'] for item in residues_browser if item['gpcrdb'] == interaction['gpcrdb']), None))

        interacting_aa = residuelist.filter(display_generic_number__label__in=[interaction['gpcrdb']]).values_list('amino_acid', flat=True)

        if interacting_aa:
            interaction['aa'] = interacting_aa[0]
            pos = residuelist.filter(display_generic_number__label__in=[interaction['gpcrdb']]).values_list('sequence_number', flat=True)[0]
            interaction['pos'] = pos

            feature = names_aa[gs_b2_interaction_type_long]

            if interacting_aa[0] not in exchange_table[feature]:
                GS_none_equivalent_interacting_pos.append(pos)
                GS_none_equivalent_interacting_gn.append(interaction['gpcrdb'])

    GS_equivalent_interacting_pos = list(residuelist.filter(display_generic_number__label__in=interacting_gn).values_list('sequence_number', flat=True))

    gProteinData = ProteinGProteinPair.objects.filter(protein__entry_name=protein)

    primary = []
    secondary = []

    for entry in gProteinData:
        if entry.transduction == 'primary':
            primary.append((entry.g_protein.name.replace("Gs","G<sub>s</sub>").replace("Gi","G<sub>i</sub>").replace("Go","G<sub>o</sub>").replace("G11","G<sub>11</sub>").replace("G12","G<sub>12</sub>").replace("G13","G<sub>13</sub>").replace("Gq","G<sub>q</sub>").replace("G","G&alpha;"),entry.g_protein.slug))
        elif entry.transduction == 'secondary':
            secondary.append((entry.g_protein.name.replace("Gs","G<sub>s</sub>").replace("Gi","G<sub>i</sub>").replace("Go","G<sub>o</sub>").replace("G11","G<sub>11</sub>").replace("G12","G<sub>12</sub>").replace("G13","G<sub>13</sub>").replace("Gq","G<sub>q</sub>").replace("G","G&alpha;"),entry.g_protein.slug))


    return render(request, 'signprot/ginterface.html', {'pdbname': '3SN6', 'snakeplot': SnakePlot, 'gproteinplot': gproteinplot, 'crystal': crystal, 'interacting_equivalent': GS_equivalent_interacting_pos, 'interacting_none_equivalent': GS_none_equivalent_interacting_pos, 'accessible': accessible_pos, 'residues': residues_browser, 'mapped_protein': protein, 'interacting_gn': GS_none_equivalent_interacting_gn, 'primary_Gprotein': set(primary), 'secondary_Gprotein': set(secondary)} )


def ajaxInterface(request, slug, **response_kwargs):

    name_of_cache = 'ajaxInterface_' + slug

    jsondata = cache.get(name_of_cache)

    if jsondata == None:

        p = Protein.objects.filter(entry_name=slug).get()

        if p.family.slug.startswith('200'):
            rsets = ResiduePositionSet.objects.get(name="Arrestin interface")
        else:
            rsets = ResiduePositionSet.objects.get(name="Gprotein Barcode")

        jsondata = {}
        for x, residue in enumerate(rsets.residue_position.all()):
            try:
                pos = str(list(Residue.objects.filter(protein_conformation__protein__entry_name=slug, display_generic_number__label=residue.label))[0])
            except:
                print("Protein has no residue position at", residue.label)
            a = pos[1:]

            jsondata[a] = [5, 'Receptor interface position', residue.label]

        jsondata = json.dumps(jsondata)

    cache.set(name_of_cache, jsondata, 60*60*24*2) #two days timeout on cache

    response_kwargs['content_type'] = 'application/json'

    return HttpResponse(jsondata, **response_kwargs)

def ajaxBarcode(request, slug, cutoff, **response_kwargs):

    name_of_cache = 'ajaxBarcode_'+slug+cutoff

    jsondata = cache.get(name_of_cache)

    if jsondata == None:
        jsondata = {}

        selectivity_pos = list(SignprotBarcode.objects.filter(protein__entry_name=slug, seq_identity__gte=cutoff).values_list('residue__display_generic_number__label', flat=True))

        conserved = list(SignprotBarcode.objects.filter(protein__entry_name=slug, paralog_score__gte=cutoff, seq_identity__gte=cutoff).prefetch_related('residue__display_generic_number').values_list('residue__display_generic_number__label', flat=True))

        na_data = list(SignprotBarcode.objects.filter(protein__entry_name=slug, seq_identity=0, paralog_score=0).values_list('residue__display_generic_number__label', flat=True))

        all_positions = Residue.objects.filter(protein_conformation__protein__entry_name=slug).prefetch_related('display_generic_number')

        for res in all_positions:
            cgn = str(res.generic_number)
            res = str(res.sequence_number)
            if cgn in conserved:
                jsondata[res] = [0, 'Conserved', cgn]
            elif cgn in selectivity_pos and cgn not in conserved:
                jsondata[res] = [1, 'Selectivity determining', cgn]
            elif cgn in na_data:
                jsondata[res] = [3, 'NA', cgn]
            else:
                jsondata[res] = [2, 'Evolutionary neutral', cgn]

        jsondata = json.dumps(jsondata)
        response_kwargs['content_type'] = 'application/json'

        cache.set(name_of_cache, jsondata, 60*60*24*2) #two days timeout on cache

    return HttpResponse(jsondata, **response_kwargs)

@cache_page(60*60*24*2)
def StructureInfo(request, pdbname):
    """
    Show structure details
    """
    protein = Protein.objects.get(signprotstructure__PDB_code=pdbname)

    crystal = SignprotStructure.objects.get(PDB_code=pdbname)

    return render(request,'signprot/structure_info.html',{'pdbname': pdbname, 'protein': protein, 'crystal': crystal})

# @cache_page(60*60*24*2)
def signprotdetail(request, slug):
    # get protein

    slug = slug.lower()
    p = Protein.objects.prefetch_related('web_links__web_resource').get(entry_name=slug, sequence_type__slug='wt')

    # get family list
    pf = p.family
    families = [pf.name]
    while pf.parent.parent:
        families.append(pf.parent.name)
        pf = pf.parent
    families.reverse()

    # get protein aliases
    aliases = ProteinAlias.objects.filter(protein=p).values_list('name', flat=True)

    # get genes
    genes = Gene.objects.filter(proteins=p).values_list('name', flat=True)
    gene = genes[0]
    alt_genes = genes[1:]

    # get structures of this signal protein
    structures = SignprotStructure.objects.filter(protein=p)
    complex_structures = SignprotComplex.objects.filter(protein=p)

    # mutations
    mutations = MutationExperiment.objects.filter(protein=p)


    # get residues
    pc = ProteinConformation.objects.get(protein=p)

    residues = Residue.objects.filter(protein_conformation=pc).order_by('sequence_number').prefetch_related(
        'protein_segment', 'generic_number', 'display_generic_number')

    # process residues and return them in chunks of 10
    # this is done for easier scaling on smaller screens
    chunk_size = 10
    r_chunks = []
    r_buffer = []
    last_segment = False
    border = False
    title_cell_skip = 0
    for i, r in enumerate(residues):
        # title of segment to be written out for the first residue in each segment
        segment_title = False

        # keep track of last residues segment (for marking borders)
        if r.protein_segment.slug != last_segment:
            last_segment = r.protein_segment.slug
            border = True

        # if on a border, is there room to write out the title? If not, write title in next chunk
        if i == 0 or (border and len(last_segment) <= (chunk_size - i % chunk_size)):
            segment_title = True
            border = False
            title_cell_skip += len(last_segment) # skip cells following title (which has colspan > 1)

        if i and i % chunk_size == 0:
            r_chunks.append(r_buffer)
            r_buffer = []

        r_buffer.append((r, segment_title, title_cell_skip))

        # update cell skip counter
        if title_cell_skip > 0:
            title_cell_skip -= 1
    if r_buffer:
        r_chunks.append(r_buffer)
    context = {'p': p, 'families': families, 'r_chunks': r_chunks, 'chunk_size': chunk_size, 'aliases': aliases,
        'gene': gene, 'alt_genes': alt_genes, 'structures': structures, 'complex_structures': complex_structures, 'mutations': mutations}

    return render(request, 'signprot/signprot_details.html', context)


def interface_dataset():
    dataset = {
        '4x1h' : [
        ['A','N',73,'2.40x40','C','C',347,'G.H5.23', ["water-mediated"]],
        ['A','A',233,'5.68x68','C','V',340,'G.H5.16', ["hydrophobic"]],
        ['A','V',138,'3.53x53','C','D',343,'G.H5.19', ["polar-backbone-sidechain", "hydrophobic"]],
        ['A','V',250,'6.33x33','C','L',349,'G.H5.25', ["hydrophobic"]],
        ['A','V',139,'3.54x54','C','V',340,'G.H5.16', ["hydrophobic", "water-mediated"]],
        ['A','V',250,'6.33x33','C','L',344,'G.H5.20', ["hydrophobic"]],
        ['A','V',139,'3.54x54','C','L',344,'G.H5.20', ["hydrophobic"]],
        ['A','M',253,'6.36x36','C','G',348,'G.H5.24', ["water-mediated"]],
        ['A','T',243,'6.26x26','C','L',341,'G.H5.17', ["hydrophobic"]],
        ['A','K',311,'8.48x48','C','F',350,'G.H5.26', ["polar-sidechain-sidechain", "polar-sidechain-backbone", "cation-pi", "van-der-waals"]],
        ['A','N',310,'8.47x47','C','C',347,'G.H5.23', ["polar-sidechain-backbone"]],
        ['A','A',246,'6.29x29','C','F',350,'G.H5.26', ["van-der-waals", "hydrophobic"]],
        ['A','L',226,'5.61x61','C','L',349,'G.H5.25', ["hydrophobic"]],
        ['A','K',245,'6.28x28','C','F',350,'G.H5.26', ["hydrophobic"]],
        ['A','N',310,'8.47x47','C','G',348,'G.H5.24', ["polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['A','L',226,'5.61x61','C','L',344,'G.H5.20', ["hydrophobic"]],
        ['A','T',242,'6.25x25','C','F',350,'G.H5.26', ["van-der-waals", "hydrophobic"]],
        ['A','A',246,'6.29x29','C','L',341,'G.H5.17', ["van-der-waals", "hydrophobic"]],
        ['A','M',309,'7.56x56','C','G',348,'G.H5.24', ["water-mediated"]],
        ['A','T',242,'6.25x25','C','L',341,'G.H5.17', ["hydrophobic"]],
        ['A','V',230,'5.65x65','C','L',344,'G.H5.20', ["hydrophobic"]],
        ['A','L',72,'2.39x39','C','S',346,'G.H5.22', ["van-der-waals", "hydrophobic"]],
        ['A','E',249,'6.32x32','C','F',350,'G.H5.26', ["hydrophobic"]],
        ['A','R',135,'3.50x50','C','L',349,'G.H5.25', ["van-der-waals", "hydrophobic"]],
        ['A','K',141,'3.56x56','C','D',343,'G.H5.19', ["ionic", "polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ['A','L',72,'2.39x39','C','C',347,'G.H5.23', ["water-mediated"]],
        ['A','T',229,'5.64x64','C','V',340,'G.H5.16', ["hydrophobic"]],
        ['A','V',139,'3.54x54','C','D',343,'G.H5.19', ["hydrophobic", "water-mediated"]],
        ['A','A',246,'6.29x29','C','L',344,'G.H5.20', ["hydrophobic"]],
        ['A','R',135,'3.50x50','C','C',347,'G.H5.23', ["polar-sidechain-backbone", "van-der-waals", "hydrophobic", "water-mediated"]],
        ],
        '3sn6' : [
        ['R','F',139,'34.51x51','A','R',380,'G.H5.12', ["van-der-waals", "hydrophobic"]],
        ['R','Q',142,'34.54x54','A','I',383,'G.H5.15', ["hydrophobic"]],
        ['R','R',228,'5.67x67','A','D',381,'G.H5.13', ["ionic", "polar-sidechain-sidechain"]],
        ['R','E',225,'5.64x64','A','R',380,'G.H5.12', ["ionic", "polar-sidechain-sidechain"]],
        ['R','F',139,'34.51x51','A','V',217,'G.S3.01', ["van-der-waals", "hydrophobic"]],
        ['R','I',135,'3.54x54','A','L',388,'G.H5.20', ["van-der-waals", "hydrophobic"]],
        ['R','K',140,'34.52x52','A','R',380,'G.H5.12', ["polar-sidechain-sidechain"]],
        ['R','I',233,'5.72x72','A','Y',358,'G.h4s6.20', ["van-der-waals", "hydrophobic"]],
        ['R','P',138,'34.50x50','A','R',380,'G.H5.12', ["hydrophobic"]],
        ['R','A',271,'6.33x33','A','L',393,'G.H5.25', ["van-der-waals", "hydrophobic"]],
        ['R','L',230,'5.69x69','A','L',394,'G.H5.26', ["van-der-waals", "hydrophobic"]],
        ['R','P',138,'34.50x50','A','Q',384,'G.H5.16', ["van-der-waals", "hydrophobic"]],
        ['R','T',274,'6.36x36','A','L',393,'G.H5.25', ["polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['R','E',225,'5.64x64','A','Q',384,'G.H5.16', ["polar-sidechain-sidechain", "van-der-waals"]],
        ['R','K',232,'5.71x71','A','D',381,'G.H5.13', ["ionic", "h-bond donor-acceptor", "polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','T',274,'6.36x36','A','E',392,'G.H5.24', ["polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['R','P',138,'34.50x50','A','I',383,'G.H5.15', ["van-der-waals", "hydrophobic"]],
        ['R','R',239,'-','A','R',347,'G.H4.17', ["polar-sidechain-backbone", "hydrophobic"]],
        ['R','R',239,'-','A','T',350,'G.h4s6.03', ["polar-sidechain-sidechain", "hydrophobic"]],
        ['R','P',138,'34.50x50','A','H',387,'G.H5.19', ["hydrophobic"]],
        ['R','F',139,'34.51x51','A','F',376,'G.H5.08', ["edge-to-face", "van-der-waals", "hydrophobic"]],
        ['R','R',239,'-','A','D',343,'G.H4.13', ["polar-sidechain-backbone"]],
        ['R','I',233,'5.72x72','A','L',394,'G.H5.26', ["hydrophobic"]],
        ['R','F',139,'34.51x51','A','H',41,'G.S1.02', ["edge-to-face", "hydrophobic"]],
        ['R','R',131,'3.50x50','A','Y',391,'G.H5.23', ["cation-pi", "van-der-waals", "hydrophobic"]],
        ['R','Q',142,'34.54x54','A','H',387,'G.H5.19', ["polar-sidechain-sidechain"]],
        ['R','T',136,'3.55x55','A','R',380,'G.H5.12', ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','Q',229,'5.68x68','A','D',381,'G.H5.13', ["polar-sidechain-sidechain", "polar-sidechain-backbone"]],
        ['R','Q',229,'5.68x68','A','L',388,'G.H5.20', ["hydrophobic"]],
        ['R','S',143,'34.55x55','A','A',39,'G.hns1.03', ["van-der-waals", "hydrophobic"]],
        ['R','Q',229,'5.68x68','A','R',385,'G.H5.17', ["h-bond acceptor-donor", "polar-sidechain-sidechain", "polar-backbone-sidechain", "hydrophobic"]],
        ['R','L',275,'6.37x37','A','L',393,'G.H5.25', ["van-der-waals", "hydrophobic"]],
        ['R','I',135,'3.54x54','A','Y',391,'G.H5.23', ["hydrophobic"]],
        ['R','I',135,'3.54x54','A','H',387,'G.H5.19', ["polar-backbone-sidechain", "hydrophobic"]],
        ['R','R',228,'5.67x67','A','Q',384,'G.H5.16', ["polar-sidechain-sidechain"]],
        ['R','V',222,'5.61x61','A','L',393,'G.H5.25', ["van-der-waals", "hydrophobic"]],
        ['R','I',135,'3.54x54','A','Q',384,'G.H5.16', ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','T',274,'6.36x36','A','Y',391,'G.H5.23', ["polar-sidechain-backbone"]],
        ['R','F',139,'34.51x51','A','I',383,'G.H5.15', ["hydrophobic"]],
        ['R','R',239,'-','A','L',346,'G.H4.16', ["hydrophobic"]],
        ['R','F',139,'34.51x51','A','C',379,'G.H5.11', ["hydrophobic"]],
        ['R','A',226,'5.65x65','A','L',388,'G.H5.20', ["van-der-waals", "hydrophobic"]],
        ['R','K',235,'5.74x74','A','D',323,'G.hgh4.13', ["polar-backbone-sidechain"]],
        ['R','D',130,'3.49x49','A','Y',391,'G.H5.23', ["polar-sidechain-sidechain"]],
        ['R','I',233,'5.72x72','A','R',385,'G.H5.17', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','Y',141,'34.53x53','A','H',387,'G.H5.19', ["pi-cation", "edge-to-face", "hydrophobic"]],
        ['R','A',134,'3.53x53','A','H',387,'G.H5.19', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','I',135,'3.54x54','A','L',393,'G.H5.25', ["van-der-waals", "hydrophobic"]],
        ['R','Q',229,'5.68x68','A','Q',384,'G.H5.16', ["polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ],
        '5g53' : [
        ['A','Y',112,'34.53x53','C','H',387,'G.H5.19', ["pi-cation", "edge-to-face", "hydrophobic"]],
        ['A','A',105,'3.53x53','C','H',387,'G.H5.19', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['A','I',106,'3.54x54','C','Y',391,'G.H5.23', ["hydrophobic"]],
        ['A','R',111,'34.52x52','C','R',380,'G.H5.12', ["polar-sidechain-sidechain"]],
        ['A','M',211,'5.72x72','C','Y',358,'G.h4s6.20', ["polar-sidechain-sidechain", "hydrophobic"]],
        ['A','P',109,'34.50x50','C','Q',384,'G.H5.16', ["hydrophobic"]],
        ['A','R',107,'3.55x55','C','R',380,'G.H5.12', ["polar-backbone-sidechain", "van-der-waals"]],
        ['A','I',106,'3.54x54','C','L',388,'G.H5.20', ["van-der-waals", "hydrophobic"]],
        ['A','R',293,'8.48x48','C','E',392,'G.H5.24', ["polar-backbone-sidechain", "hydrophobic"]],
        ['A','R',111,'34.52x52','C','V',217,'G.S3.01', ["van-der-waals", "hydrophobic"]],
        ['A','A',203,'5.64x64','C','L',388,'G.H5.20', ["hydrophobic"]],
        ['A','Q',207,'5.68x68','C','Q',384,'G.H5.16', ["h-bond donor-acceptor", "polar-sidechain-sidechain"]],
        ['A','Q',210,'5.71x71','C','D',381,'G.H5.13', ["h-bond donor-acceptor", "polar-sidechain-sidechain", "van-der-waals"]],
        ['A','A',203,'5.64x64','C','Q',384,'G.H5.16', ["polar-backbone-sidechain", "hydrophobic"]],
        ['A','L',110,'34.51x51','C','V',217,'G.S3.01', ["van-der-waals", "hydrophobic"]],
        ['A','Q',207,'5.68x68','C','D',381,'G.H5.13', ["polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['A','L',110,'34.51x51','C','R',380,'G.H5.12', ["hydrophobic"]],
        ['A','R',111,'34.52x52','C','D',215,'G.s2s3.01', ["polar-sidechain-backbone", "van-der-waals"]],
        ['A','L',208,'5.69x69','C','L',394,'G.H5.26', ["van-der-waals", "hydrophobic"]],
        ['A','R',291,'7.56x56','C','Y',391,'G.H5.23', ["polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['A','I',106,'3.54x54','C','Q',384,'G.H5.16', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['A','I',108,'3.56x56','C','R',380,'G.H5.12', ["polar-backbone-sidechain"]],
        ['A','R',102,'3.50x50','C','Y',391,'G.H5.23', ["polar-backbone-sidechain", "cation-pi", "van-der-waals", "hydrophobic"]],
        ['A','A',231,'6.33x33','C','L',393,'G.H5.25', ["hydrophobic"]],
        ['A','R',291,'7.56x56','C','E',392,'G.H5.24', ["polar-backbone-sidechain"]],
        ['A','L',110,'34.51x51','C','F',219,'G.S3.03', ["hydrophobic"]],
        ['A','I',200,'5.61x61','C','L',393,'G.H5.25', ["hydrophobic"]],
        ['A','Q',207,'5.68x68','C','L',388,'G.H5.20', ["hydrophobic"]],
        ['A','Q',207,'5.68x68','C','R',385,'G.H5.17', ["polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['A','P',109,'34.50x50','C','I',383,'G.H5.15', ["van-der-waals", "hydrophobic"]],
        ['A','L',110,'34.51x51','C','H',41,'H.HD.11', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['A','K',227,'6.29x29','C','L',394,'G.H5.26', ["polar-sidechain-sidechain", "hydrophobic"]],
        ['A','I',106,'3.54x54','C','H',387,'G.H5.19', ["hydrophobic"]],
        ['A','L',235,'6.37x37','C','L',393,'G.H5.25', ["van-der-waals", "hydrophobic"]],
        ['A','L',110,'34.51x51','C','C',379,'G.H5.11', ["hydrophobic"]],
        ['A','R',296,'8.51x51','C','E',392,'G.H5.24', ["ionic", "polar-sidechain-sidechain", "van-der-waals"]],
        ['A','I',200,'5.61x61','C','L',388,'G.H5.20', ["hydrophobic"]],
        ['A','Q',207,'5.68x68','C','Y',360,'G.S6.02', ["polar-sidechain-sidechain"]],
        ['A','L',110,'34.51x51','C','F',376,'G.H5.08', ["van-der-waals", "hydrophobic"]],
        ['A','L',110,'34.51x51','C','I',383,'G.H5.15', ["hydrophobic"]],
        ['A','A',204,'5.65x65','C','L',388,'G.H5.20', ["hydrophobic"]],
        ['A','P',109,'34.50x50','C','R',380,'G.H5.12', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ],
        '5uz7' : [
        ['R','L',247,'3.57x57','A','H',387,'G.H5.19', ["van-der-waals", "hydrophobic"]],
        ['R','K',326,'5.64x64','A','R',380,'G.H5.12', ["polar-sidechain-sidechain"]],
        ['R','L',244,'3.54x54','A','Y',391,'G.H5.23', ["hydrophobic"]],
        ['R','L',348,'6.45x45','A','E',392,'G.H5.24', ["hydrophobic"]],
        ['R','V',252,'-','A','I',383,'G.H5.15', ["van-der-waals", "hydrophobic"]],
        ['R','K',326,'5.64x64','A','Q',384,'G.H5.16', ["polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','M',327,'5.65x65','A','L',394,'G.H5.26', ["van-der-waals", "hydrophobic"]],
        ['R','N',396,'8.48x48','A','E',392,'G.H5.24', ["polar-backbone-sidechain"]],
        ['R','F',253,'-','A','H',41,'G.S1.02', ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','I',248,'3.58x58','A','L',388,'G.H5.20', ["hydrophobic"]],
        ['R','N',396,'8.48x48','A','R',356,'G.h4s6.12', ["polar-sidechain-sidechain"]],
        ['R','R',180,'2.46x46','A','Y',391,'G.H5.23', ["polar-sidechain-backbone", "cation-pi", "van-der-waals", "hydrophobic"]],
        ['R','L',348,'6.45x45','A','L',393,'G.H5.25', ["hydrophobic"]],
        ['R','V',252,'-','A','R',380,'G.H5.12', ["van-der-waals", "hydrophobic"]],
        ['R','C',394,'7.60x60','A','E',392,'G.H5.24', ["polar-backbone-sidechain"]],
        ['R','V',249,'3.59x59','A','Q',384,'G.H5.16', ["polar-backbone-sidechain"]],
        ['R','R',180,'2.46x46','A','Q',390,'G.H5.22', ["h-bond donor-acceptor", "polar-sidechain-sidechain", "polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['R','L',247,'3.57x57','A','Y',391,'G.H5.23', ["hydrophobic"]],
        ['R','I',248,'3.58x58','A','Q',384,'G.H5.16', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','T',330,'-','A','Y',358,'G.h4s6.20', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','Y',243,'3.53x53','A','Y',391,'G.H5.23', ["van-der-waals"]],
        ['R','H',331,'-','A','Y',358,'G.h4s6.20', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','L',323,'5.61x61','A','L',388,'G.H5.20', ["van-der-waals", "hydrophobic"]],
        ['R','F',253,'-','A','V',217,'G.S3.01', ["hydrophobic"]],
        ['R','V',252,'-','A','Q',384,'G.H5.16', ["van-der-waals", "hydrophobic"]],
        ['R','T',254,'-','A','H',387,'G.H5.19', ["polar-sidechain-sidechain"]],
        ['R','E',329,'-','A','R',385,'G.H5.17', ["polar-backbone-sidechain"]],
        ['R','I',248,'3.58x58','A','H',387,'G.H5.19', ["hydrophobic"]],
        ['R','T',345,'6.42x42','A','L',393,'G.H5.25', ["van-der-waals", "hydrophobic"]],
        ['R','K',326,'5.64x64','A','R',385,'G.H5.17', ["polar-backbone-sidechain"]],
        ['R','H',184,'2.50x50','A','Y',391,'G.H5.23', ["hydrophobic"]],
        ],
        '5vai' : [
        ['R','L',339,'-','A','Y',358,'G.h4s6.20', ["van-der-waals", "hydrophobic"]],
        ['R','S',352,'6.41x41','A','L',393,'G.H5.25', ["polar-sidechain-backbone", "van-der-waals"]],
        ['R','L',401,'7.56x56','A','E',392,'G.H5.24', ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','S',352,'6.41x41','A','E',392,'G.H5.24', ["polar-sidechain-backbone"]],
        ['R','E',262,'4.38x39','A','Q',35,'G.HN.52', ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','V',331,'5.61x61','A','L',388,'G.H5.20', ["van-der-waals", "hydrophobic"]],
        ['R','T',353,'6.42x42','A','L',393,'G.H5.25', ["hydrophobic"]],
        ['R','L',255,'3.58x58','A','R',380,'G.H5.12', ["polar-backbone-sidechain"]],
        ['R','L',359,'6.48x48','A','Y',391,'G.H5.23', ["van-der-waals", "hydrophobic"]],
        ['R','L',251,'3.54x54','A','Y',391,'G.H5.23', ["van-der-waals", "hydrophobic"]],
        ['R','V',405,'7.60x60','A','E',392,'G.H5.24', ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','R',176,'2.46x46','A','Q',390,'G.H5.22', ["polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','N',407,'8.48x48','A','E',392,'G.H5.24', ["polar-backbone-sidechain"]],
        ['R','L',255,'3.58x58','A','Q',384,'G.H5.16', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','A',256,'3.59x59','A','R',380,'G.H5.12', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','Q',263,'4.39x40','A','Q',35,'G.HN.52', ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','Q',263,'4.39x40','A','Q',31,'G.HN.48', ["polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','S',352,'6.41x41','A','L',394,'G.H5.26', ["polar-sidechain-backbone"]],
        ['R','H',180,'2.50x50','A','Y',391,'G.H5.23', ["polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','K',334,'5.64x64','A','R',385,'G.H5.17', ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','L',339,'-','A','L',394,'G.H5.26', ["hydrophobic"]],
        ['R','E',408,'8.49x49','A','Q',390,'G.H5.22', ["h-bond acceptor-donor", "polar-sidechain-sidechain", "polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['R','K',334,'5.64x64','A','D',381,'G.H5.13', ["ionic", "h-bond donor-acceptor", "polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','N',406,'8.47x47','A','E',392,'G.H5.24', ["polar-sidechain-sidechain", "polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','L',339,'-','A','R',385,'G.H5.17', ["polar-backbone-sidechain"]],
        ['R','N',338,'-','A','C',359,'G.S6.01', ["polar-sidechain-backbone"]],
        ['R','N',338,'-','A','Y',360,'G.S6.02', ["van-der-waals", "hydrophobic"]],
        ['R','L',356,'6.45x45','A','L',393,'G.H5.25', ["hydrophobic"]],
        ['R','N',406,'8.47x47','A','Q',390,'G.H5.22', ["polar-sidechain-backbone"]],
        ['R','L',254,'3.57x57','A','H',387,'G.H5.19', ["van-der-waals", "hydrophobic"]],
        ['R','Y',402,'7.57x57','A','E',392,'G.H5.24', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','A',256,'3.59x59','A','Q',384,'G.H5.16', ["polar-backbone-sidechain"]],
        ['R','S',261,'4.37x38','A','Q',35,'G.HN.52', ["polar-sidechain-sidechain", "polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','L',356,'6.45x45','A','Y',391,'G.H5.23', ["hydrophobic"]],
        ['R','K',342,'-','A','T',350,'G.h4s6.03', ["h-bond donor-acceptor", "polar-sidechain-sidechain", "polar-sidechain-backbone", "van-der-waals"]],
        ['R','R',176,'2.46x46','A','Y',391,'G.H5.23', ["hydrophobic"]],
        ['R','R',264,'4.40x41','A','Q',35,'G.HN.52', ["polar-backbone-sidechain"]],
        ],
        '6b3j' : [
        ['R','V',331,'5.61x61','A','L',393,'G.H5.25', ["hydrophobic"]],
        ['R','S',352,'6.41x41','A','L',393,'G.H5.25', ["polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['R','E',262,'4.38x39','A','K',34,'G.HN.51', ["ionic", "polar-sidechain-sidechain"]],
        ['R','V',331,'5.61x61','A','L',388,'G.H5.20', ["hydrophobic"]],
        ['R','L',251,'3.54x54','A','Y',391,'G.H5.23', ["van-der-waals", "hydrophobic"]],
        ['R','S',261,'-','A','Q',35,'G.HN.52', ["polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','R',176,'2.46x46','A','Q',390,'G.H5.22', ["polar-sidechain-sidechain", "polar-sidechain-backbone", "hydrophobic"]],
        ['R','N',407,'8.48x48','A','E',392,'G.H5.24', ["polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','L',255,'3.58x58','A','Q',384,'G.H5.16', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','S',352,'6.41x41','A','L',394,'G.H5.26', ["polar-sidechain-backbone"]],
        ['R','H',180,'2.50x50','A','Y',391,'G.H5.23', ["hydrophobic"]],
        ['R','S',258,'-','A','I',383,'G.H5.15', ["hydrophobic"]],
        ['R','R',348,'6.37x37','A','L',394,'G.H5.26', ["hydrophobic"]],
        ['R','K',334,'5.64x64','A','Q',384,'G.H5.16', ["polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','K',334,'5.64x64','A','D',381,'G.H5.13', ["ionic", "h-bond donor-acceptor", "polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','V',331,'5.61x61','A','L',394,'G.H5.26', ["hydrophobic"]],
        ['R','V',259,'-','A','V',217,'G.S3.01', ["hydrophobic"]],
        ['R','Y',250,'3.53x53','A','Y',391,'G.H5.23', ["van-der-waals"]],
        ['R','L',255,'3.58x58','A','H',387,'G.H5.19', ["hydrophobic"]],
        ['R','L',356,'6.45x45','A','L',393,'G.H5.25', ["hydrophobic"]],
        ['R','K',334,'5.64x64','A','R',385,'G.H5.17', ["polar-backbone-sidechain", "hydrophobic"]],
        ['R','V',327,'5.57x57','A','L',393,'G.H5.25', ["hydrophobic"]],
        ['R','L',254,'3.57x57','A','H',387,'G.H5.19', ["van-der-waals", "hydrophobic"]],
        ['R','E',262,'4.38x39','A','R',38,'G.hns1.02', ["ionic", "polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','L',255,'3.58x58','A','L',388,'G.H5.20', ["hydrophobic"]],
        ['R','R',176,'2.46x46','A','Y',391,'G.H5.23', ["cation-pi", "hydrophobic"]],
        ['R','K',334,'5.64x64','A','L',388,'G.H5.20', ["hydrophobic"]],
        ],
        '6cmo' : [
        ['R','V',139,'3.54x54','A','L',348,'G.H5.20', ["hydrophobic"]],
        ['R','K',311,'8.48x48','A','F',354,'G.H5.26', ["polar-sidechain-sidechain", "polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['R','Q',237,'5.72x72','A','D',341,'G.H5.13', ["hydrophobic"]],
        ['R','K',311,'8.48x48','A','G',352,'G.H5.24', ["polar-sidechain-backbone"]],
        ['R','M',309,'7.56x56','A','G',352,'G.H5.24', ["van-der-waals"]],
        ['R','M',253,'6.36x36','A','L',353,'G.H5.25', ["van-der-waals", "hydrophobic"]],
        ['R','Q',237,'5.72x72','A','Y',320,'G.h4s6.20', ["polar-sidechain-sidechain"]],
        ['R','A',241,'6.24x24','A','E',318,'G.h4s6.12', ["polar-backbone-sidechain"]],
        ['R','A',246,'6.29x29','A','L',348,'G.H5.20', ["hydrophobic"]],
        ['R','R',147,'34.55x55','A','A',31,'G.hns1.02', ["polar-sidechain-backbone"]],
        ['R','K',245,'6.28x28','A','F',354,'G.H5.26', ["hydrophobic"]],
        ['R','V',250,'6.33x33','A','L',353,'G.H5.25', ["hydrophobic"]],
        ['R','K',311,'8.48x48','A','K',349,'G.H5.21', ["polar-sidechain-backbone"]],
        ['R','T',242,'6.25x25','A','D',315,'G.h4s6.09', ["polar-sidechain-backbone"]],
        ['R','S',240,'-','A','K',345,'G.H5.17', ["h-bond acceptor-donor", "polar-sidechain-sidechain"]],
        ['R','E',239,'-','A','Y',320,'G.h4s6.20', ["van-der-waals", "hydrophobic"]],
        ['R','R',147,'34.55x55','A','R',32,'G.hns1.03', ["polar-sidechain-backbone", "hydrophobic"]],
        ['R','V',139,'3.54x54','A','N',347,'G.H5.19', ["hydrophobic"]],
        ['R','R',135,'3.50x50','A','L',353,'G.H5.25', ["hydrophobic"]],
        ['R','K',141,'3.56x56','A','D',193,'G.s2s3.02', ["ionic", "polar-sidechain-sidechain"]],
        ['R','T',243,'6.26x26','A','D',341,'G.H5.13', ["polar-sidechain-sidechain"]],
        ['R','E',239,'-','A','E',318,'G.h4s6.12', ["polar-backbone-sidechain"]],
        ['R','E',249,'6.32x32','A','F',354,'G.H5.26', ["hydrophobic"]],
        ['R','S',240,'-','A','E',318,'G.h4s6.12', ["h-bond donor-acceptor", "polar-sidechain-sidechain", "polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','K',311,'8.48x48','A','L',353,'G.H5.25', ["polar-sidechain-backbone"]],
        ['R','N',310,'8.47x47','A','G',352,'G.H5.24', ["polar-sidechain-backbone", "hydrophobic"]],
        ['R','T',242,'6.25x25','A','F',354,'G.H5.26', ["van-der-waals", "hydrophobic"]],
        ['R','E',249,'6.32x32','A','L',353,'G.H5.25', ["polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['R','A',246,'6.29x29','A','F',354,'G.H5.26', ["van-der-waals", "hydrophobic"]],
        ],
        '6d9h' : [
        ['R','V',203,'5.61x61','A','L',354,'G.H5.25', ["hydrophobic"]],
        ['R','K',224,'6.25x25','A','D',316,'G.h4s6.09', ["ionic", "h-bond donor-acceptor", "polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','I',292,'8.47x47','A','G',353,'G.H5.24', ["hydrophobic"]],
        ['R','K',213,'5.71x71','A','D',342,'G.H5.13', ["ionic", "polar-sidechain-sidechain"]],
        ['R','K',231,'6.32x32','A','F',355,'G.H5.26', ["polar-sidechain-backbone"]],
        ['R','L',113,'34.51x51','A','I',344,'G.H5.15', ["hydrophobic"]],
        ['R','R',108,'3.53x53','A','N',348,'G.H5.19', ["polar-sidechain-sidechain", "polar-backbone-sidechain", "hydrophobic"]],
        ['R','R',108,'3.53x53','A','D',351,'G.H5.22', ["ionic", "polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','I',232,'6.33x33','A','L',354,'G.H5.25', ["van-der-waals", "hydrophobic"]],
        ['R','L',113,'34.51x51','A','F',337,'G.H5.08', ["hydrophobic"]],
        ['R','R',105,'3.50x50','A','C',352,'G.H5.23', ["polar-sidechain-backbone", "hydrophobic"]],
        ['R','Q',210,'5.68x68','A','I',345,'G.H5.16', ["van-der-waals", "hydrophobic"]],
        ['R','V',109,'3.54x54','A','L',349,'G.H5.20', ["van-der-waals", "hydrophobic"]],
        ['R','Q',210,'5.68x68','A','D',342,'G.H5.13', ["polar-sidechain-sidechain", "van-der-waals"]],
        ['R','K',228,'6.29x29','A','D',316,'G.h4s6.09', ["ionic", "polar-sidechain-sidechain"]],
        ['R','L',211,'5.69x69','A','K',346,'G.H5.17', ["hydrophobic"]],
        ['R','K',294,'8.49x49','A','D',351,'G.H5.22', ["ionic", "h-bond donor-acceptor", "polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','K',228,'6.29x29','A','F',355,'G.H5.26', ["polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','L',236,'6.37x37','A','L',354,'G.H5.25', ["van-der-waals", "hydrophobic"]],
        ['R','I',232,'6.33x33','A','F',355,'G.H5.26', ["hydrophobic"]],
        ['R','D',42,'2.37x37','A','D',351,'G.H5.22', ["polar-sidechain-sidechain"]],
        ['R','R',105,'3.50x50','A','L',354,'G.H5.25', ["hydrophobic"]],
        ['R','I',207,'5.65x65','A','L',349,'G.H5.20', ["hydrophobic"]],
        ['R','F',45,'2.40x40','A','D',351,'G.H5.22', ["van-der-waals"]],
        ['R','K',224,'6.25x25','A','E',319,'G.h4s6.12', ["ionic", "polar-sidechain-sidechain"]],
        ['R','P',112,'34.50x50','A','N',348,'G.H5.19', ["polar-backbone-sidechain"]],
        ['R','R',108,'3.53x53','A','C',352,'G.H5.23', ["polar-sidechain-sidechain", "van-der-waals"]],
        ['R','P',112,'34.50x50','A','I',344,'G.H5.15', ["van-der-waals", "hydrophobic"]],
        ['R','L',113,'34.51x51','A','L',195,'G.S3.01', ["hydrophobic"]],
        ['R','I',292,'8.47x47','A','C',352,'G.H5.23', ["van-der-waals", "hydrophobic"]],
        ['R','L',113,'34.51x51','A','T',341,'G.H5.12', ["hydrophobic"]],
        ['R','I',207,'5.65x65','A','L',354,'G.H5.25', ["hydrophobic"]],
        ['R','R',291,'7.56x56','A','G',353,'G.H5.24', ["hydrophobic"]],
        ['R','P',112,'34.50x50','A','I',345,'G.H5.16', ["van-der-waals", "hydrophobic"]],
        ],
        '6dde' : [
        ['R','M',264,'-','A','D',341,'G.H5.13', ["polar-backbone-sidechain", "hydrophobic"]],
        ['R','S',268,'6.23x23','A','D',315,'G.h4s6.09', ["polar-backbone-sidechain"]],
        ['R','K',271,'6.26x26','A','K',314,'G.h4s6.08', ["polar-sidechain-backbone", "van-der-waals"]],
        ['R','T',103,'2.39x39','A','C',351,'G.H5.23', ["polar-sidechain-sidechain", "polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['R','R',179,'34.57x57','A','N',347,'G.H5.19', ["polar-sidechain-sidechain", "van-der-waals"]],
        ['R','D',177,'34.55x55','A','R',32,'G.hns1.03', ["ionic", "polar-sidechain-sidechain", "polar-backbone-sidechain"]],
        ['R','T',103,'2.39x39','A','D',350,'G.H5.22', ["polar-sidechain-backbone"]],
        ['R','L',259,'5.65x65','A','L',348,'G.H5.20', ["hydrophobic"]],
        ['R','L',176,'34.54x54','A','R',32,'G.hns1.03', ["polar-backbone-sidechain", "hydrophobic"]],
        ['R','I',278,'6.33x33','A','L',353,'G.H5.25', ["hydrophobic"]],
        ['R','E',270,'6.25x25','A','D',315,'G.h4s6.09', ["polar-sidechain-sidechain", "polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['R','R',165,'3.50x50','A','C',351,'G.H5.23', ["polar-sidechain-sidechain", "van-der-waals"]],
        ['R','A',168,'3.53x53','A','N',347,'G.H5.19', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','R',263,'-','A','I',319,'G.h4s6.13', ["polar-sidechain-backbone"]],
        ['R','M',255,'5.61x61','A','L',353,'G.H5.25', ["hydrophobic"]],
        ['R','R',258,'5.64x64','A','I',344,'G.H5.16', ["van-der-waals", "hydrophobic"]],
        ['R','I',278,'6.33x33','A','F',354,'G.H5.26', ["van-der-waals", "hydrophobic"]],
        ['R','M',264,'-','A','K',345,'G.H5.17', ["hydrophobic"]],
        ['R','L',259,'5.65x65','A','I',344,'G.H5.16', ["hydrophobic"]],
        ['R','M',281,'6.36x36','A','L',353,'G.H5.25', ["van-der-waals", "hydrophobic"]],
        ['R','V',173,'34.51x51','A','D',193,'G.s2s3.02', ["van-der-waals", "hydrophobic"]],
        ['R','V',173,'34.51x51','A','F',336,'G.H5.08', ["van-der-waals", "hydrophobic"]],
        ['R','R',277,'6.32x32','A','L',353,'G.H5.25', ["polar-sidechain-backbone"]],
        ['R','V',262,'5.68x68','A','I',344,'G.H5.16', ["hydrophobic"]],
        ['R','P',172,'34.50x50','A','I',343,'G.H5.15', ["van-der-waals", "hydrophobic"]],
        ['R','R',165,'3.50x50','A','L',353,'G.H5.25', ["hydrophobic"]],
        ['R','V',173,'34.51x51','A','L',194,'G.S3.01', ["hydrophobic"]],
        ['R','L',176,'34.54x54','A','I',343,'G.H5.15', ["van-der-waals", "hydrophobic"]],
        ['R','M',264,'-','A','T',316,'G.h4s6.10', ["polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['R','P',172,'34.50x50','A','T',340,'G.H5.12', ["hydrophobic"]],
        ['R','R',182,'4.40x40','A','R',24,'G.HN.48', ["polar-sidechain-sidechain", "hydrophobic"]],
        ['R','R',263,'-','A','Y',320,'G.h4s6.20', ["van-der-waals", "hydrophobic"]],
        ['R','P',172,'34.50x50','A','I',344,'G.H5.16', ["van-der-waals", "hydrophobic"]],
        ['R','D',340,'8.47x47','A','G',352,'G.H5.24', ["polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['R','V',262,'5.68x68','A','D',341,'G.H5.13', ["van-der-waals", "hydrophobic"]],
        ['R','L',176,'34.54x54','A','L',194,'G.S3.01', ["hydrophobic"]],
        ['R','K',271,'6.26x26','A','D',315,'G.h4s6.09', ["van-der-waals", "hydrophobic"]],
        ['R','V',169,'3.54x54','A','L',348,'G.H5.20', ["van-der-waals", "hydrophobic"]],
        ['R','R',263,'-','A','D',341,'G.H5.13', ["polar-backbone-sidechain"]],
        ['R','I',278,'6.33x33','A','L',348,'G.H5.20', ["hydrophobic"]],
        ],
        '6ddf' : [
        ['R','R',263,'-','A','I',319,'G.h4s6.13', ["polar-sidechain-backbone", "van-der-waals"]],
        ['R','M',264,'-','A','K',345,'G.H5.17', ["hydrophobic"]],
        ['R','R',258,'5.64x64','A','I',344,'G.H5.16', ["van-der-waals", "hydrophobic"]],
        ['R','R',179,'34.57x57','A','N',347,'G.H5.19', ["polar-sidechain-sidechain"]],
        ['R','L',259,'5.65x65','A','L',348,'G.H5.20', ["van-der-waals", "hydrophobic"]],
        ['R','I',278,'6.33x33','A','L',348,'G.H5.20', ["hydrophobic"]],
        ['R','I',278,'6.33x33','A','L',353,'G.H5.25', ["hydrophobic"]],
        ['R','E',270,'6.25x25','A','D',315,'G.h4s6.09', ["polar-sidechain-sidechain", "polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['R','R',165,'3.50x50','A','C',351,'G.H5.23', ["polar-sidechain-backbone", "hydrophobic"]],
        ['R','T',103,'2.39x39','A','C',351,'G.H5.23', ["polar-sidechain-sidechain", "polar-sidechain-backbone", "hydrophobic"]],
        ['R','I',278,'6.33x33','A','F',354,'G.H5.26', ["van-der-waals", "hydrophobic"]],
        ['R','R',179,'34.57x57','A','C',351,'G.H5.23', ["polar-sidechain-sidechain", "van-der-waals"]],
        ['R','A',168,'3.53x53','A','N',347,'G.H5.19', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','E',341,'8.48x48','A','F',354,'G.H5.26', ["hydrophobic"]],
        ['R','P',172,'34.50x50','A','I',343,'G.H5.15', ["van-der-waals", "hydrophobic"]],
        ['R','T',103,'2.39x39','A','D',350,'G.H5.22', ["polar-sidechain-backbone", "van-der-waals"]],
        ['R','V',173,'34.51x51','A','L',194,'G.S3.01', ["hydrophobic"]],
        ['R','L',176,'34.54x54','A','I',343,'G.H5.15', ["hydrophobic"]],
        ['R','M',264,'-','A','T',316,'G.h4s6.10', ["polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['R','P',172,'34.50x50','A','T',340,'G.H5.12', ["van-der-waals", "hydrophobic"]],
        ['R','R',263,'-','A','Y',320,'G.h4s6.20', ["polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','P',172,'34.50x50','A','I',344,'G.H5.16', ["van-der-waals", "hydrophobic"]],
        ['R','D',340,'8.47x47','A','G',352,'G.H5.24', ["polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['R','V',262,'5.68x68','A','D',341,'G.H5.13', ["van-der-waals", "hydrophobic"]],
        ['R','L',176,'34.54x54','A','L',194,'G.S3.01', ["hydrophobic"]],
        ['R','E',341,'8.48x48','A','G',352,'G.H5.24', ["polar-sidechain-backbone", "van-der-waals"]],
        ['R','M',255,'5.61x61','A','L',353,'G.H5.25', ["hydrophobic"]],
        ['R','S',268,'6.23x23','A','D',315,'G.h4s6.09', ["polar-sidechain-sidechain", "polar-backbone-sidechain", "van-der-waals"]],
        ['R','K',271,'6.26x26','A','K',314,'G.h4s6.08', ["polar-sidechain-backbone", "van-der-waals"]],
        ['R','L',259,'5.65x65','A','I',344,'G.H5.16', ["hydrophobic"]],
        ['R','D',177,'34.55x55','A','R',32,'G.hns1.03', ["ionic", "polar-sidechain-sidechain"]],
        ['R','D',340,'8.47x47','A','L',353,'G.H5.25', ["polar-sidechain-backbone"]],
        ['R','L',176,'34.54x54','A','R',32,'G.hns1.03', ["hydrophobic"]],
        ['R','M',281,'6.36x36','A','L',353,'G.H5.25', ["van-der-waals", "hydrophobic"]],
        ['R','R',277,'6.32x32','A','F',354,'G.H5.26', ["polar-sidechain-backbone", "van-der-waals"]],
        ['R','A',168,'3.53x53','A','C',351,'G.H5.23', ["van-der-waals"]],
        ['R','V',173,'34.51x51','A','D',193,'G.s2s3.02', ["hydrophobic"]],
        ['R','V',173,'34.51x51','A','F',336,'G.H5.08', ["van-der-waals", "hydrophobic"]],
        ['R','R',277,'6.32x32','A','L',353,'G.H5.25', ["polar-sidechain-backbone", "van-der-waals"]],
        ['R','R',165,'3.50x50','A','L',353,'G.H5.25', ["hydrophobic"]],
        ['R','R',263,'-','A','E',318,'G.h4s6.12', ["van-der-waals"]],
        ['R','E',341,'8.48x48','A','L',353,'G.H5.25', ["polar-sidechain-backbone"]],
        ['R','V',262,'5.68x68','A','I',344,'G.H5.16', ["hydrophobic"]],
        ['R','K',271,'6.26x26','A','K',317,'G.h4s6.11', ["polar-sidechain-backbone"]],
        ['R','R',263,'-','A','D',341,'G.H5.13', ["polar-backbone-sidechain", "hydrophobic"]],
        ['R','M',264,'-','A','D',341,'G.H5.13', ["polar-backbone-sidechain"]],
        ['R','K',174,'34.52x52','A','D',193,'G.s2s3.02', ["ionic", "polar-sidechain-sidechain"]],
        ['R','D',164,'3.49x49','A','C',351,'G.H5.23', ["polar-sidechain-sidechain"]],
        ['R','K',271,'6.26x26','A','D',315,'G.h4s6.09', ["van-der-waals", "hydrophobic"]],
        ['R','V',169,'3.54x54','A','L',348,'G.H5.20', ["van-der-waals", "hydrophobic"]],
        ['R','P',172,'34.50x50','A','N',347,'G.H5.19', ["polar-backbone-sidechain"]],
        ],
        '6e3y' : [
        ['R','I',241,'3.58x58','A','Y',391,'G.H5.23', ["hydrophobic"]],
        ['R','V',242,'3.59x59','A','R',380,'G.H5.12', ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','F',387,'7.60x60','A','E',392,'G.H5.24', ["polar-backbone-sidechain"]],
        ['R','I',241,'3.58x58','A','H',387,'G.H5.19', ["hydrophobic"]],
        ['R','F',246,'-','A','H',41,'G.S1.02', ["hydrophobic"]],
        ['R','L',316,'5.61x61','A','L',393,'G.H5.25', ["hydrophobic"]],
        ['R','R',173,'2.46x46','A','Q',390,'G.H5.22', ["polar-sidechain-sidechain", "polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['R','N',388,'8.47x47','A','E',392,'G.H5.24', ["polar-sidechain-sidechain", "van-der-waals"]],
        ['R','G',389,'8.48x48','A','E',392,'G.H5.24', ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','V',242,'3.59x59','A','Q',384,'G.H5.16', ["polar-backbone-sidechain"]],
        ['R','K',319,'5.64x64','A','L',388,'G.H5.20', ["hydrophobic"]],
        ['R','F',246,'-','A','F',376,'G.H5.08', ["van-der-waals", "hydrophobic"]],
        ['R','L',316,'5.61x61','A','L',388,'G.H5.20', ["van-der-waals", "hydrophobic"]],
        ['R','L',240,'3.57x57','A','Y',391,'G.H5.23', ["hydrophobic"]],
        ['R','R',336,'6.40x40','A','L',393,'G.H5.25', ["polar-sidechain-backbone"]],
        ['R','K',319,'5.64x64','A','R',385,'G.H5.17', ["polar-backbone-sidechain", "van-der-waals"]],
        ['R','F',246,'-','A','I',383,'G.H5.15', ["van-der-waals", "hydrophobic"]],
        ['R','L',237,'3.54x54','A','Y',391,'G.H5.23', ["van-der-waals", "hydrophobic"]],
        ['R','L',320,'5.65x65','A','L',394,'G.H5.26', ["van-der-waals", "hydrophobic"]],
        ['R','T',323,'-','A','Y',358,'G.h4s6.20', ["hydrophobic"]],
        ['R','R',173,'2.46x46','A','Y',391,'G.H5.23', ["polar-sidechain-backbone", "cation-pi", "van-der-waals", "hydrophobic"]],
        ['R','R',336,'6.40x40','A','L',394,'G.H5.26', ["polar-sidechain-sidechain"]],
        ['R','H',177,'2.50x50','A','Y',391,'G.H5.23', ["hydrophobic"]],
        ['R','E',248,'-','A','A',39,'G.hns1.03', ["hydrophobic"]],
        ['R','V',245,'-','A','H',387,'G.H5.19', ["polar-backbone-sidechain"]],
        ['R','F',246,'-','A','F',219,'G.S3.03', ["hydrophobic"]],
        ['R','V',245,'-','A','Q',384,'G.H5.16', ["van-der-waals", "hydrophobic"]],
        ['R','K',319,'5.64x64','A','Q',384,'G.H5.16', ["h-bond donor-acceptor", "polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','V',245,'-','A','I',383,'G.H5.15', ["van-der-waals", "hydrophobic"]],
        ['R','F',246,'-','A','R',380,'G.H5.12', ["van-der-waals", "hydrophobic"]],
        ['R','I',312,'5.57x57','A','L',393,'G.H5.25', ["hydrophobic"]],
        ['R','E',248,'-','A','R',38,'G.hns1.02', ["van-der-waals", "hydrophobic"]],
        ['R','F',246,'-','A','V',217,'G.S3.01', ["van-der-waals", "hydrophobic"]],
        ['R','V',243,'3.60x60','A','R',380,'G.H5.12', ["polar-backbone-sidechain"]],
        ['R','V',245,'-','A','R',380,'G.H5.12', ["van-der-waals", "hydrophobic"]],
        ['R','L',341,'6.45x45','A','L',393,'G.H5.25', ["van-der-waals", "hydrophobic"]],
        ['R','I',241,'3.58x58','A','L',388,'G.H5.20', ["hydrophobic"]],
        ['R','R',336,'6.40x40','A','E',392,'G.H5.24', ["ionic", "polar-sidechain-sidechain", "polar-sidechain-backbone", "van-der-waals"]],
        ['R','L',240,'3.57x57','A','H',387,'G.H5.19', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['R','L',316,'5.61x61','A','L',394,'G.H5.26', ["van-der-waals", "hydrophobic"]],
        ['R','F',246,'-','A','C',379,'G.H5.11', ["hydrophobic"]],
        ['R','I',241,'3.58x58','A','L',393,'G.H5.25', ["hydrophobic"]],
        ['R','K',333,'6.37x37','A','L',394,'G.H5.26', ["polar-sidechain-sidechain", "polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['R','K',319,'5.64x64','A','D',381,'G.H5.13', ["van-der-waals", "hydrophobic"]],
        ['R','I',241,'3.58x58','A','Q',384,'G.H5.16', ["polar-backbone-sidechain", "van-der-waals"]],
        ],
        '6g79' : [
        ['S','R',308,'6.29x29','A','Y',354,'G.H5.26', ["polar-sidechain-sidechain", "polar-sidechain-backbone", "cation-pi", "van-der-waals", "hydrophobic"]],
        ['S','N',373,'8.47x47','A','G',350,'G.H5.22', ["polar-sidechain-backbone"]],
        ['S','K',311,'6.32x32','A','G',352,'G.H5.24', ["polar-sidechain-backbone"]],
        ['S','I',231,'5.61x61','A','L',353,'G.H5.25', ["van-der-waals", "hydrophobic"]],
        ['S','V',155,'34.51x51','A','F',336,'G.H5.08', ["hydrophobic"]],
        ['S','I',151,'3.54x54','A','I',344,'G.H5.16', ["hydrophobic"]],
        ['S','R',147,'3.50x50','A','C',351,'G.H5.23', ["polar-sidechain-sidechain", "polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['S','R',161,'34.57x57','A','N',347,'G.H5.19', ["polar-sidechain-sidechain"]],
        ['S','I',239,'5.69x69','A','Y',354,'G.H5.26', ["hydrophobic"]],
        ['S','T',315,'6.36x36','A','L',353,'G.H5.25', ["polar-sidechain-backbone", "hydrophobic"]],
        ['S','S',372,'7.56x56','A','C',351,'G.H5.23', ["polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['S','I',151,'3.54x54','A','N',347,'G.H5.19', ["hydrophobic"]],
        ['S','A',235,'5.65x65','A','L',353,'G.H5.25', ["hydrophobic"]],
        ['S','R',238,'5.68x68','A','L',348,'G.H5.20', ["hydrophobic"]],
        ['S','T',315,'6.36x36','A','G',352,'G.H5.24', ["polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['S','I',151,'3.54x54','A','L',348,'G.H5.20', ["van-der-waals", "hydrophobic"]],
        ['S','A',312,'6.33x33','A','L',353,'G.H5.25', ["van-der-waals", "hydrophobic"]],
        ['S','K',311,'6.32x32','A','Y',354,'G.H5.26', ["polar-sidechain-sidechain", "polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['S','R',238,'5.68x68','A','D',341,'G.H5.13', ["ionic", "polar-sidechain-sidechain", "polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['S','S',372,'7.56x56','A','G',352,'G.H5.24', ["hydrophobic"]],
        ['S','L',316,'6.37x37','A','L',353,'G.H5.25', ["hydrophobic"]],
        ['S','R',238,'5.68x68','A','I',344,'G.H5.16', ["van-der-waals", "hydrophobic"]],
        ['S','A',154,'34.50x50','A','I',344,'G.H5.16', ["van-der-waals", "hydrophobic"]],
        ['S','A',235,'5.65x65','A','L',348,'G.H5.20', ["hydrophobic"]],
        ['S','A',154,'34.50x50','A','I',343,'G.H5.15', ["hydrophobic"]],
        ['S','A',150,'3.53x53','A','N',347,'G.H5.19', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['S','R',238,'5.68x68','A','Y',354,'G.H5.26', ["hydrophobic"]],
        ['S','V',155,'34.51x51','A','T',340,'G.H5.12', ["hydrophobic"]],
        ['S','R',308,'6.29x29','A','N',316,'G.h4s6.10', ["polar-sidechain-backbone"]],
        ],
        '6gdg' : [
        ['A','Q',207,'5.68x68','D','R',375,'G.H5.17', ["polar-sidechain-backbone"]],
        ['A','I',106,'3.54x54','D','Y',381,'G.H5.23', ["hydrophobic"]],
        ['A','R',111,'34.52x52','D','D',215,'G.s2s3.01', ["polar-sidechain-backbone"]],
        ['A','I',200,'5.61x61','D','L',383,'G.H5.25', ["van-der-waals", "hydrophobic"]],
        ['A','R',111,'34.52x52','D','K',216,'G.s2s3.02', ["polar-sidechain-backbone", "hydrophobic"]],
        ['A','N',113,'34.54x54','D','A',39,'H.HD.09', ["hydrophobic"]],
        ['A','R',293,'8.48x48','D','E',382,'G.H5.24', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['A','E',294,'8.49x49','D','Q',380,'G.H5.22', ["polar-backbone-sidechain"]],
        ['A','A',105,'3.53x53','D','H',377,'G.H5.19', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['A','L',110,'34.51x51','D','R',370,'G.H5.12', ["hydrophobic"]],
        ['A','R',291,'7.56x56','D','Y',381,'G.H5.23', ["polar-sidechain-backbone"]],
        ['A','Q',207,'5.68x68','D','L',384,'G.H5.26', ["hydrophobic"]],
        ['A','A',231,'6.33x33','D','L',383,'G.H5.25', ["van-der-waals", "hydrophobic"]],
        ['A','Q',207,'5.68x68','D','L',378,'G.H5.20', ["hydrophobic"]],
        ['A','Q',207,'5.68x68','D','D',371,'G.H5.13', ["h-bond donor-acceptor", "polar-sidechain-sidechain", "polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['A','R',293,'8.48x48','D','Q',380,'G.H5.22', ["polar-backbone-sidechain"]],
        ['A','H',230,'6.32x32','D','E',382,'G.H5.24', ["polar-sidechain-backbone"]],
        ['A','P',109,'34.50x50','D','Q',374,'G.H5.16', ["van-der-waals", "hydrophobic"]],
        ['A','L',208,'5.69x69','D','L',384,'G.H5.26', ["hydrophobic"]],
        ['A','G',114,'34.55x55','D','A',39,'H.HD.09', ["hydrophobic"]],
        ['A','S',234,'6.36x36','D','L',383,'G.H5.25', ["van-der-waals", "hydrophobic"]],
        ['A','L',110,'34.51x51','D','V',217,'G.S3.01', ["hydrophobic"]],
        ['A','R',291,'7.56x56','D','E',382,'G.H5.24', ["polar-sidechain-backbone", "polar-backbone-sidechain", "van-der-waals"]],
        ['A','R',111,'34.52x52','D','V',217,'G.S3.01', ["van-der-waals", "hydrophobic"]],
        ['A','Y',112,'34.53x53','D','H',377,'G.H5.19', ["edge-to-face", "hydrophobic"]],
        ['A','Q',210,'5.71x71','D','D',371,'G.H5.13', ["polar-sidechain-sidechain"]],
        ['A','L',110,'34.51x51','D','H',41,'H.HD.11', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['A','L',110,'34.51x51','D','I',373,'G.H5.15', ["hydrophobic"]],
        ['A','I',106,'3.54x54','D','H',377,'G.H5.19', ["van-der-waals", "hydrophobic"]],
        ['A','A',203,'5.64x64','D','L',378,'G.H5.20', ["hydrophobic"]],
        ['A','L',110,'34.51x51','D','C',369,'G.H5.11', ["hydrophobic"]],
        ['A','P',109,'34.50x50','D','R',370,'G.H5.12', ["hydrophobic"]],
        ['A','I',106,'3.54x54','D','L',378,'G.H5.20', ["van-der-waals", "hydrophobic"]],
        ['A','N',113,'34.54x54','D','R',38,'H.HD.08', ["polar-sidechain-backbone", "van-der-waals", "hydrophobic"]],
        ['A','L',110,'34.51x51','D','F',366,'G.H5.08', ["van-der-waals", "hydrophobic"]],
        ['A','M',211,'5.72x72','D','Y',348,'G.h4s6.20', ["hydrophobic"]],
        ['A','S',234,'6.36x36','D','E',382,'G.H5.24', ["polar-sidechain-backbone"]],
        ['A','I',106,'3.54x54','D','Q',374,'G.H5.16', ["polar-backbone-sidechain", "van-der-waals", "hydrophobic"]],
        ['A','Q',207,'5.68x68','D','Q',374,'G.H5.16', ["polar-sidechain-sidechain", "van-der-waals", "hydrophobic"]],
        ['A','P',109,'34.50x50','D','I',373,'G.H5.15', ["van-der-waals", "hydrophobic"]],
        ['A','R',102,'3.50x50','D','Y',381,'G.H5.23', ["cation-pi", "van-der-waals", "hydrophobic"]],
        ['A','L',235,'6.37x37','D','L',383,'G.H5.25', ["van-der-waals", "hydrophobic"]],
        ['A','A',204,'5.65x65','D','L',378,'G.H5.20', ["hydrophobic"]],
        ['A','R',107,'3.55x55','D','R',370,'G.H5.12', ["polar-sidechain-sidechain", "polar-backbone-sidechain", "van-der-waals"]],
        ],
    }

    return dataset


def InteractionMatrix(request):
    from django.db.models import F
    from django.db.models import Q
    from signprot.views import interface_dataset
    import requests

    dataset = interface_dataset()

    # generate complex info dataset
    filt = [e.upper() for e in list(dataset)]
    struc = Structure.objects.filter(pdb_code__index__in=filt).prefetch_related('protein_conformation__protein__parent')

    complex_info = []
    for s in struc:
        r = {}
        r['pdb_id'] = str.lower(s.pdb_code.index)
        r['name'] = s.protein_conformation.protein.parent.name
        r['entry_name'] = s.protein_conformation.protein.parent.entry_name
        r['class'] = s.protein_conformation.protein.get_protein_class()
        r['family'] = s.protein_conformation.protein.get_protein_family()
        r['conf_id'] = s.protein_conformation.id
        r['organism'] = s.protein_conformation.protein.species.common_name
        try:
            r['gprot'] = s.get_stab_agents_gproteins()
        except Exception:
            r['gprot'] = ''
        try:
            r['gprot_class'] = s.get_signprot_gprot_family()
        except Exception:
            r['gprot_class'] = ''
        complex_info.append(r)

    data = Protein.objects.filter(
        sequence_type__slug='wt',
        species__common_name="Human",
        family__slug__startswith='00',  # receptors, no gproteins
        ).prefetch_related(
            'parent__protein_conformation',
            'family__parent__parent__parent',
        )

    # proteins = []
    # for s in data:
    #     r = {}
    #     r['name'] = s.name
    #     r['entry_name'] = s.entry_name
    #     r['protein_family'] = s.family.parent.short()
    #     r['protein_class'] = s.family.parent.parent.parent.short()
    #     r['ligand'] = s.family.parent.parent.short()
    #     proteins.append(r)

    interactions_metadata = complex_info
    gprotein_order = ProteinSegment.objects.filter(proteinfamily='Alpha').values('id', 'slug')
    prot_conf_ids = [i['conf_id'] for i in complex_info]
    remaining_residues = Residue.objects.filter(
            protein_conformation_id__in=prot_conf_ids,
            ).values(
                rec_id = F('protein_conformation__protein__id'),
                name = F('protein_conformation__protein__parent__name'),
                entry_name = F('protein_conformation__protein__parent__entry_name'),
                rec_aa = F('amino_acid'),
                rec_gn = F('display_generic_number__label'),
            ).exclude(
                Q(rec_gn=None)
            )

    new_dataset = []
    for pdb_key in dataset:
        for residue_list in dataset[pdb_key]:
            curr_meta = None
            while curr_meta is None:
                for meta in interactions_metadata:
                    if meta['pdb_id'].upper() == pdb_key.upper():
                        curr_meta = meta
                break
            if curr_meta is not None:
                gprot = curr_meta['gprot']
                entry_name = curr_meta['entry_name']
                pdb_id = curr_meta['pdb_id']
                residue_list.extend([gprot, entry_name, pdb_id])
                new_dataset.append(residue_list)

    context = {
        'interactions': new_dataset,
        # 'interactions': json.dumps(dataset),
        'non_interactions': json.dumps(list(remaining_residues)),
        'interactions_metadata': json.dumps(interactions_metadata),
        # 'ps': json.dumps(list(proteins)),
        'gprot': json.dumps(list(gprotein_order)),
        }

    request.session['signature'] = None
    request.session.modified = True

    return render(request, 'signprot/matrix.html', context)


def IMSequenceSignature(request):
    '''Accept set of proteins + generic numbers and calculate the signature for those'''
    t1 = time.time()

    pos_set_in = get_entry_names(request)
    ignore_in_alignment = get_ignore_info(request)
    segments = get_protein_segments(request)

    # get pos objects
    pos_set = Protein.objects.filter(entry_name__in=pos_set_in).select_related('residue_numbering_scheme', 'species')

    # Calculate Sequence Signature
    signature = SequenceSignature()
    signature.setup_alignments_signprot(segments, pos_set, ignore_in_alignment=ignore_in_alignment)
    signature.calculate_signature_onesided()
    # preprocess data for return
    signature_data = signature.prepare_display_data_onesided()

    # FEATURES AND REGIONS
    feats = [feature for feature in signature_data['a_pos'].features_combo]

    # GET GENERIC NUMBERS
    generic_numbers = get_generic_numbers(signature_data)

    # FEATURE FREQUENCIES
    signature_features = get_signature_features(signature_data, generic_numbers, feats)
    grouped_features = group_signature_features(signature_features)

    # FEATURE CONSENSUS
    generic_numbers_flat = list(chain.from_iterable(generic_numbers))
    sigcons = get_signature_consensus(signature_data, generic_numbers_flat)

    rec_class = pos_set[0].get_protein_class()

    # dump = {
    #     'rec_class': rec_class,
    #     'signature': signature,
    #     'consensus': signature_data,
    #     }
    # with open('signprot/notebooks/interface_pickles/{}.p'.format(rec_class), 'wb+') as out_file:
    #     pickle.dump(dump, out_file)

    # pass back to front
    res = {
        'cons': sigcons,
        'feat': grouped_features,
    }

    request.session['signature'] = signature.prepare_session_data()
    request.session.modified = True

    t2 = time.time()
    print('Runtime: {}'.format((t2-t1)*1000.0))

    return JsonResponse(res, safe=False)


def IMSignatureMatch(request):
    '''Take the signature stored in the session and query the db'''
    signature_data = request.session.get('signature')
    ss_pos = request.POST.getlist('pos[]')
    cutoff = request.POST.get('cutoff')
    request.session['ss_pos'] = ss_pos
    request.session['cutoff'] = cutoff

    pos_set = Protein.objects.filter(entry_name__in=ss_pos).select_related('residue_numbering_scheme', 'species')
    pos_set = [protein for protein in pos_set]
    pfam = [protein.family.slug[:3] for protein in pos_set]

    signature_match = SignatureMatch(
        signature_data['common_positions'],
        signature_data['numbering_schemes'],
        signature_data['common_segments'],
        signature_data['diff_matrix'],
        pos_set,
        pos_set,
        cutoff = 0
    )

    maj_pfam = Counter(pfam).most_common()[0][0]
    signature_match.score_protein_class(maj_pfam)
    # request.session['signature_match'] = signature_match

    signature_match = {
        'scores': signature_match.protein_report,
        'scores_pos': signature_match.scores_pos,
        'scores_neg': signature_match.scores_neg,
        'protein_signatures': signature_match.protein_signatures,
        'signatures_pos': signature_match.signatures_pos,
        'signatures_neg': signature_match.signatures_neg,
        'signature_filtered': signature_match.signature_consensus,
        'relevant_gn': signature_match.relevant_gn,
        'relevant_segments': signature_match.relevant_segments,
        'numbering_schemes': signature_match.schemes,
    }

    signature_match = prepare_signature_match(signature_match)
    return JsonResponse(signature_match, safe=False)


def render_IMSigMat(request):

    # signature_match = request.session.get('signature_match')
    signature_data = request.session.get('signature')
    ss_pos = request.session.get('ss_pos')
    cutoff = request.session.get('cutoff')

    pos_set = Protein.objects.filter(entry_name__in=ss_pos).select_related('residue_numbering_scheme', 'species')
    pos_set = [protein for protein in pos_set]
    pfam = [protein.family.slug[:3] for protein in pos_set]

    signature_match = SignatureMatch(
        signature_data['common_positions'],
        signature_data['numbering_schemes'],
        signature_data['common_segments'],
        signature_data['diff_matrix'],
        pos_set,
        pos_set,
        cutoff = 0
    )

    maj_pfam = Counter(pfam).most_common()[0][0]
    signature_match.score_protein_class(maj_pfam)


    response = render(
        request,
        'signature_match.html',
        {'scores': signature_match}
        )
    return response
