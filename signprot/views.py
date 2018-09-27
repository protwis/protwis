from django.shortcuts import get_object_or_404, render
from django.http import HttpResponse
from django.core.cache import cache
from django.views.decorators.cache import cache_page

from protein.models import Protein, ProteinConformation, ProteinAlias, ProteinFamily, Gene, ProteinGProtein, ProteinGProteinPair
from residue.models import Residue, ResiduePositionSet

from structure.models import Structure
from mutation.models import MutationExperiment
from common.selection import Selection
from common.diagrams_gpcr import DrawSnakePlot
from common.diagrams_gprotein import DrawGproteinPlot
from common.diagrams_arrestin import DrawArrestinPlot

from signprot.models import SignprotStructure, SignprotBarcode, SignprotInteractions

from common import definitions
from collections import OrderedDict
from common.views import AbsTargetSelection

import json
# Create your views here.
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
        ppf_g = ProteinFamily.objects.get(slug="100_000")
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

@cache_page(60*60*24*2) # 2 days caching
def GProtein(request):

    name_of_cache = 'gprotein_statistics'

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
                ps = gp.proteingproteinpair_set.filter(protein__family__slug__startswith=slug)

                if ps:
                    jsondata[str(gp)] = []
                    for p in ps:
                        if str(p.protein.entry_name).split('_')[0].upper() not in selectivitydata:
                            selectivitydata[str(p.protein.entry_name).split('_')[0].upper()] = []
                        selectivitydata[str(p.protein.entry_name).split('_')[0].upper()].append(str(gp))
                        # print(p.protein.family.parent.parent.parent)
                        jsondata[str(gp)].append(str(p.protein.entry_name)+'\n')

                    jsondata[str(gp)] = ''.join(jsondata[str(gp)])

            context[slug_translate[slug]] = jsondata

        context["selectivitydata"] = selectivitydata



    return render(request, 'signprot/gprotein.html', context)

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
    structures = SignprotStructure.objects.filter(origin__family__slug__startswith=slug
        )

    mutations = MutationExperiment.objects.filter(protein__in=proteins).prefetch_related('residue__generic_number',
                                'exp_qual', 'ligand')

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

    name_of_cache = 'ajaxInterface_'+slug

    jsondata = cache.get(name_of_cache)

    if jsondata == None:

        if slug == "arrs_human":
            rsets = ResiduePositionSet.objects.get(name="Arrestin interface")
        else:
            rsets = ResiduePositionSet.objects.get(name="Gprotein Barcode")
        # residues = Residue.objects.filter(protein_conformation__protein__entry_name=slug, display_generic_number__label=residue.label)

        jsondata = {}
        positions = []
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
    structures = SignprotStructure.objects.filter(origin=p)

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
        'gene': gene, 'alt_genes': alt_genes, 'structures': structures, 'mutations': mutations}

    return render(request, 'signprot/signprot_details.html', context)

def InteractionMatrix(request):

    dataset = {
        '3sn6' : [
            ['R','139','34.51x51','A','376', ["edge-to-face", "face-to-edge", "hydrophobic", "van-der-waals"]],
            ['R','233','5.72x72','A','394', ["hydrophobic"]],
            ['R','233','5.72x72','A','385', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['R','135','3.54x54','A','388', ["hydrophobic", "van-der-waals"]],
            ['R','135','3.54x54','A','393', ["hydrophobic", "van-der-waals"]],
            ['R','62','12.48x48','B','312', ["polar-sidechain-sidechain"]],
            ['R','229','5.68x68','A','388', ["hydrophobic"]],
            ['R','63','12.49x49','B','312', ["polar-backbone-sidechain"]],
            ['R','142','34.54x54','A','387', ["polar-sidechain-sidechain"]],
            ['R','274','6.36x36','A','391', ["polar-sidechain-backbone"]],
            ['R','239','-','A','350', ["hydrophobic", "polar-sidechain-sidechain"]],
            ['R','232','5.71x71','A','381', ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
            ['R','135','3.54x54','A','384', ["polar-backbone-sidechain", "van-der-waals"]],
            ['R','139','34.51x51','A','383', ["hydrophobic"]],
            ['R','142','34.54x54','A','383', ["hydrophobic"]],
            ['R','235','5.74x74','A','323', ["polar-backbone-sidechain"]],
            ['R','271','6.33x33','A','393', ["hydrophobic", "van-der-waals"]],
            ['R','225','5.64x64','A','380', ["polar-sidechain-sidechain"]],
            ['R','138','34.50x50','A','380', ["hydrophobic"]],
            ['R','233','5.72x72','A','358', ["hydrophobic", "van-der-waals"]],
            ['R','139','34.51x51','A','217', ["hydrophobic", "van-der-waals"]],
            ['R','138','34.50x50','A','384', ["hydrophobic", "van-der-waals"]],
            ['R','228','5.67x67','A','381', ["polar-sidechain-sidechain"]],
            ['R','239','-','A','347', ["hydrophobic", "polar-sidechain-backbone"]],
            ['R','139','34.51x51','A','379', ["hydrophobic"]],
            ['R','135','3.54x54','A','387', ["hydrophobic", "polar-backbone-sidechain"]],
            ['R','141','34.53x53','A','387', ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
            ['R','140','34.52x52','A','380', ["polar-sidechain-sidechain"]],
            ['R','239','-','A','346', ["hydrophobic"]],
            ['R','136','3.55x55','A','380', ["polar-backbone-sidechain", "van-der-waals"]],
            ['R','274','6.36x36','A','392', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
            ['R','225','5.64x64','A','384', ["polar-sidechain-sidechain", "van-der-waals"]],
            ['R','230','5.69x69','A','394', ["hydrophobic", "van-der-waals"]],
            ['R','229','5.68x68','A','384', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','239','-','A','343', ["polar-sidechain-backbone"]],
            ['R','138','34.50x50','A','387', ["hydrophobic"]],
            ['R','135','3.54x54','A','391', ["hydrophobic"]],
            ['R','139','34.51x51','A','41', ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
            ['R','130','3.49x49','A','391', ["polar-sidechain-sidechain"]],
            ['R','138','34.50x50','A','383', ["hydrophobic", "van-der-waals"]],
            ['R','275','6.37x37','A','393', ["hydrophobic", "van-der-waals"]],
            ['R','228','5.67x67','A','384', ["polar-sidechain-sidechain"]],
            ['R','229','5.68x68','A','385', ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain", "h-bond"]],
            ['R','134','3.53x53','A','387', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['R','226','5.65x65','A','388', ["hydrophobic", "van-der-waals"]],
            ['R','274','6.36x36','A','393', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
            ['R','139','34.51x51','A','380', ["hydrophobic", "van-der-waals"]],
            ['R','143','34.55x55','A','39', ["hydrophobic", "van-der-waals"]],
            ['R','222','5.61x61','A','393', ["hydrophobic", "van-der-waals"]],
            ['R','229','5.68x68','A','381', ["polar-sidechain-backbone", "polar-sidechain-sidechain"]],
            ['R','131','3.50x50','A','391', ["cation-pi", "hydrophobic", "van-der-waals"]],
           ],
        '4x1h' : [
            ['A','135','3.50x50','C','347', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals", "water-mediated", "water-mediated"]],
            ['A','72','2.39x39','C','347', ["water-mediated", "water-mediated"]],
            ['A','310','8.47x47','C','346', ["water-mediated"]],
            ['A','246','6.29x29','C','341', ["hydrophobic", "van-der-waals"]],
            ['A','135','3.50x50','C','346', ["water-mediated"]],
            ['A','309','7.56x56','C','348', ["water-mediated", "water-mediated"]],
            ['A','305','7.52x52','C','348', ["water-mediated"]],
            ['A','138','3.53x53','C','343', ["hydrophobic", "polar-backbone-sidechain", "water-mediated"]],
            ['A','139','3.54x54','C','340', ["hydrophobic", "water-mediated"]],
            ['A','73','2.40x40','C','346', ["water-mediated", "water-mediated"]],
            ['A','306','7.53x53','C','348', ["water-mediated", "water-mediated"]],
            ['A','312','8.49x49','C','346', ["water-mediated"]],
            ['A','139','3.54x54','C','344', ["hydrophobic"]],
            ['A','250','6.33x33','C','344', ["hydrophobic"]],
            ['A','250','6.33x33','C','349', ["hydrophobic"]],
            ['A','135','3.50x50','C','348', ["water-mediated", "water-mediated"]],
            ['A','243','6.26x26','C','341', ["hydrophobic"]],
            ['A','141','3.56x56','C','343', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals", "water-mediated"]],
            ['A','312','8.49x49','C','348', ["water-mediated"]],
            ['A','141','3.56x56','C','340', ["water-mediated"]],
            ['A','229','5.64x64','C','340', ["hydrophobic"]],
            ['A','73','2.40x40','C','347', ["water-mediated", "water-mediated"]],
            ['A','71','2.38x38','C','346', ["water-mediated"]],
            ['A','242','6.25x25','C','350', ["hydrophobic", "van-der-waals"]],
            ['A','246','6.29x29','C','344', ["hydrophobic"]],
            ['A','311','8.48x48','C','350', ["polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
            ['A','226','5.61x61','C','344', ["hydrophobic"]],
            ['A','70','2.37x37','C','346', ["water-mediated"]],
            ['A','249','6.32x32','C','350', ["hydrophobic"]],
            ['A','139','3.54x54','C','343', ["hydrophobic", "water-mediated"]],
            ['A','226','5.61x61','C','349', ["hydrophobic"]],
            ['A','230','5.65x65','C','344', ["hydrophobic"]],
            ['A','310','8.47x47','C','348', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals", "water-mediated", "water-mediated", "water-mediated", "water-mediated"]],
            ['A','233','5.68x68','C','340', ["hydrophobic"]],
            ['A','246','6.29x29','C','350', ["hydrophobic", "van-der-waals"]],
            ['A','138','3.53x53','C','340', ["water-mediated"]],
            ['A','310','8.47x47','C','347', ["polar-sidechain-backbone", "water-mediated"]],
            ['A','135','3.50x50','C','349', ["hydrophobic", "van-der-waals"]],
            ['A','73','2.40x40','C','348', ["water-mediated"]],
            ['A','72','2.39x39','C','346', ["hydrophobic", "van-der-waals", "water-mediated", "water-mediated"]],
            ['A','245','6.28x28','C','350', ["hydrophobic"]],
            ['A','253','6.36x36','C','348', ["water-mediated"]],
            ['A','312','8.49x49','C','347', ["water-mediated"]],
            ['A','242','6.25x25','C','341', ["hydrophobic"]],
            ],
        '5g53' : [
            ['A','106','3.54x54','C','388', ["hydrophobic", "van-der-waals"]],
            ['A','111','34.52x52','C','215', ["polar-sidechain-backbone", "van-der-waals"]],
            ['A','207','5.68x68','C','384', ["polar-sidechain-sidechain", "h-bond"]],
            ['A','200','5.61x61','C','388', ["hydrophobic"]],
            ['A','110','34.51x51','C','41', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['A','109','34.50x50','C','384', ["hydrophobic"]],
            ['A','107','3.55x55','C','380', ["polar-backbone-sidechain", "van-der-waals"]],
            ['A','293','8.48x48','C','392', ["hydrophobic", "polar-backbone-sidechain"]],
            ['A','200','5.61x61','C','393', ["hydrophobic"]],
            ['A','204','5.65x65','C','388', ["hydrophobic"]],
            ['A','210','5.71x71','C','381', ["polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
            ['A','291','7.56x56','C','392', ["polar-backbone-sidechain"]],
            ['A','110','34.51x51','C','379', ["hydrophobic"]],
            ['A','110','34.51x51','C','376', ["hydrophobic", "van-der-waals"]],
            ['A','106','3.54x54','C','391', ["hydrophobic"]],
            ['A','207','5.68x68','C','385', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
            ['A','207','5.68x68','C','381', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
            ['A','105','3.53x53','C','387', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['A','110','34.51x51','C','217', ["hydrophobic", "van-der-waals"]],
            ['A','203','5.64x64','C','384', ["hydrophobic", "polar-backbone-sidechain"]],
            ['A','110','34.51x51','C','380', ["hydrophobic"]],
            ['A','106','3.54x54','C','384', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['A','111','34.52x52','C','217', ["hydrophobic", "van-der-waals"]],
            ['A','291','7.56x56','C','391', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
            ['A','106','3.54x54','C','387', ["hydrophobic"]],
            ['A','111','34.52x52','C','380', ["polar-sidechain-sidechain"]],
            ['A','110','34.51x51','C','219', ["hydrophobic"]],
            ['A','211','5.72x72','C','358', ["hydrophobic", "polar-sidechain-sidechain"]],
            ['A','110','34.51x51','C','383', ["hydrophobic"]],
            ['A','108','3.56x56','C','380', ["polar-backbone-sidechain"]],
            ['A','207','5.68x68','C','388', ["hydrophobic"]],
            ['A','231','6.33x33','C','393', ["hydrophobic"]],
            ['A','109','34.50x50','C','380', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['A','208','5.69x69','C','394', ["hydrophobic", "van-der-waals"]],
            ['A','207','5.68x68','C','360', ["polar-sidechain-sidechain"]],
            ['A','109','34.50x50','C','383', ["hydrophobic", "van-der-waals"]],
            ['A','112','34.53x53','C','387', ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
            ['A','235','6.37x37','C','393', ["hydrophobic", "van-der-waals"]],
            ['A','203','5.64x64','C','388', ["hydrophobic"]],
            ['A','102','3.50x50','C','391', ["cation-pi", "hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['A','296','8.51x51','C','392', ["polar-sidechain-sidechain", "van-der-waals"]],
            ['A','227','6.29x29','C','394', ["hydrophobic", "polar-sidechain-sidechain"]],
            ],
        '5g53' : [
            ['B','200','5.61x61','D','388', ["hydrophobic"]],
            ['B','203','5.64x64','D','388', ["hydrophobic"]],
            ['B','235','6.37x37','D','393', ["hydrophobic", "van-der-waals"]],
            ['B','114','34.55x55','D','39', ["hydrophobic", "van-der-waals"]],
            ['B','207','5.68x68','D','384', ["polar-sidechain-sidechain"]],
            ['B','109','34.50x50','D','384', ["hydrophobic"]],
            ['B','293','8.48x48','D','392', ["polar-backbone-sidechain", "van-der-waals"]],
            ['B','110','34.51x51','D','217', ["hydrophobic", "van-der-waals"]],
            ['B','203','5.64x64','D','384', ["hydrophobic", "polar-backbone-sidechain"]],
            ['B','207','5.68x68','D','381', ["polar-sidechain-backbone"]],
            ['B','106','3.54x54','D','384', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['B','113','34.54x54','D','39', ["hydrophobic", "van-der-waals"]],
            ['B','296','8.51x51','D','392', ["polar-sidechain-sidechain", "van-der-waals"]],
            ['B','207','5.68x68','D','385', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
            ['B','102','3.50x50','D','391', ["hydrophobic", "polar-backbone-sidechain"]],
            ['B','111','34.52x52','D','380', ["polar-sidechain-sidechain"]],
            ['B','110','34.51x51','D','41', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['B','231','6.33x33','D','393', ["hydrophobic", "van-der-waals"]],
            ['B','106','3.54x54','D','388', ["hydrophobic", "van-der-waals"]],
            ['B','110','34.51x51','D','379', ["hydrophobic"]],
            ['B','204','5.65x65','D','393', ["hydrophobic", "van-der-waals"]],
            ['B','204','5.65x65','D','388', ["hydrophobic"]],
            ['B','105','3.53x53','D','387', ["hydrophobic", "polar-backbone-sidechain"]],
            ['B','110','34.51x51','D','219', ["hydrophobic"]],
            ['B','291','7.56x56','D','392', ["polar-backbone-sidechain"]],
            ['B','110','34.51x51','D','376', ["hydrophobic", "van-der-waals"]],
            ['B','111','34.52x52','D','217', ["hydrophobic", "van-der-waals"]],
            ['B','106','3.54x54','D','387', ["hydrophobic"]],
            ['B','207','5.68x68','D','388', ["hydrophobic"]],
            ['B','110','34.51x51','D','380', ["hydrophobic"]],
            ['B','111','34.52x52','D','215', ["polar-sidechain-backbone", "van-der-waals"]],
            ['B','109','34.50x50','D','383', ["hydrophobic", "van-der-waals"]],
            ['B','200','5.61x61','D','393', ["hydrophobic", "van-der-waals"]],
            ['B','109','34.50x50','D','380', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['B','112','34.53x53','D','387', ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
            ['B','105','3.53x53','D','391', ["van-der-waals"]],
            ['B','106','3.54x54','D','391', ["hydrophobic"]],
            ['B','107','3.55x55','D','380', ["polar-backbone-sidechain", "van-der-waals"]],
            ['B','110','34.51x51','D','383', ["hydrophobic"]],
            ['B','108','3.56x56','D','380', ["polar-backbone-sidechain"]],
            ],
        '5uz7' : [
            ['R','326','5.64x64','A','384', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','415','8.67x67','B','307', ["polar-sidechain-backbone", "van-der-waals"]],
            ['R','408','8.60x60','B','311', ["polar-sidechain-backbone", "van-der-waals"]],
            ['R','331','-','A','358', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['R','253','-','A','41', ["polar-backbone-sidechain", "van-der-waals"]],
            ['R','348','6.45x45','A','393', ["hydrophobic"]],
            ['R','329','-','A','385', ["polar-backbone-sidechain"]],
            ['R','408','8.60x60','B','309', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
            ['R','243','3.53x53','A','391', ["van-der-waals"]],
            ['R','247','3.57x57','A','387', ["hydrophobic", "van-der-waals"]],
            ['R','249','3.59x59','A','384', ["polar-backbone-sidechain"]],
            ['R','411','8.63x63','B','307', ["hydrophobic"]],
            ['R','394','7.60x60','A','392', ["polar-backbone-sidechain"]],
            ['R','348','6.45x45','A','392', ["hydrophobic"]],
            ['R','323','5.61x61','A','388', ["hydrophobic", "van-der-waals"]],
            ['R','404','8.56x56','B','312', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','254','-','A','387', ["polar-sidechain-sidechain"]],
            ['R','180','2.46x46','A','391', ["cation-pi", "hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
            ['R','252','-','A','384', ["hydrophobic", "van-der-waals"]],
            ['R','415','8.67x67','B','44', ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
            ['R','396','8.48x48','A','356', ["polar-sidechain-sidechain"]],
            ['R','326','5.64x64','A','385', ["polar-backbone-sidechain"]],
            ['R','252','-','A','380', ["hydrophobic", "van-der-waals"]],
            ['R','252','-','A','383', ["hydrophobic", "van-der-waals"]],
            ['R','408','8.60x60','B','310', ["polar-sidechain-backbone"]],
            ['R','247','3.57x57','A','391', ["hydrophobic"]],
            ['R','253','-','A','217', ["hydrophobic"]],
            ['R','327','5.65x65','A','394', ["hydrophobic", "van-der-waals"]],
            ['R','345','6.42x42','A','393', ["hydrophobic", "van-der-waals"]],
            ['R','180','2.46x46','A','390', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
            ['R','248','3.58x58','A','387', ["hydrophobic"]],
            ['R','184','2.50x50','A','391', ["hydrophobic"]],
            ['R','244','3.54x54','A','391', ["hydrophobic"]],
            ['R','330','-','A','358', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['R','248','3.58x58','A','388', ["hydrophobic"]],
            ['R','396','8.48x48','A','392', ["polar-backbone-sidechain"]],
            ['R','326','5.64x64','A','380', ["polar-sidechain-sidechain"]],
            ['R','248','3.58x58','A','384', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ],
        '5vai' : [
            ['R','415','8.56x56','B','292', ["cation-pi", "hydrophobic", "van-der-waals"]],
            ['R','262','4.38x39','A','35', ["polar-backbone-sidechain", "van-der-waals"]],
            ['R','356','6.45x45','A','393', ["hydrophobic"]],
            ['R','264','4.40x41','A','35', ["polar-backbone-sidechain"]],
            ['R','176','2.46x46','A','390', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','261','4.37x38','A','35', ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','419','8.60x60','B','310', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
            ['R','263','4.39x40','A','31', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','338','-','A','359', ["polar-sidechain-backbone"]],
            ['R','405','7.60x60','A','392', ["polar-backbone-sidechain", "van-der-waals"]],
            ['R','256','3.59x59','A','384', ["polar-backbone-sidechain"]],
            ['R','339','-','A','385', ["polar-backbone-sidechain"]],
            ['R','353','6.42x42','A','393', ["hydrophobic"]],
            ['R','419','8.60x60','B','293', ["polar-sidechain-sidechain"]],
            ['R','339','-','A','358', ["hydrophobic", "van-der-waals"]],
            ['R','408','8.49x49','A','390', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
            ['R','419','8.60x60','B','312', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','255','3.58x58','A','384', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['R','171','12.49x49','B','312', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','352','6.41x41','A','393', ["polar-sidechain-backbone", "van-der-waals"]],
            ['R','406','8.47x47','A','390', ["polar-sidechain-backbone"]],
            ['R','256','3.59x59','A','380', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['R','339','-','A','394', ["hydrophobic"]],
            ['R','419','8.60x60','B','309', ["hydrophobic", "van-der-waals"]],
            ['R','331','5.61x61','A','388', ["hydrophobic", "van-der-waals"]],
            ['R','407','8.48x48','A','392', ["polar-backbone-sidechain"]],
            ['R','412','8.53x53','B','312', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','254','3.57x57','A','387', ["hydrophobic", "van-der-waals"]],
            ['R','352','6.41x41','A','392', ["polar-sidechain-backbone"]],
            ['R','263','4.39x40','A','35', ["polar-backbone-sidechain", "van-der-waals"]],
            ['R','342','-','A','350', ["polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
            ['R','176','2.46x46','A','391', ["hydrophobic"]],
            ['R','334','5.64x64','A','385', ["polar-backbone-sidechain", "van-der-waals"]],
            ['R','359','6.48x48','A','391', ["hydrophobic", "van-der-waals"]],
            ['R','255','3.58x58','A','380', ["polar-backbone-sidechain"]],
            ['R','334','5.64x64','A','381', ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
            ['R','352','6.41x41','A','394', ["polar-sidechain-backbone"]],
            ['R','180','2.50x50','A','391', ["cation-pi", "hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','356','6.45x45','A','391', ["hydrophobic"]],
            ['R','338','-','A','360', ["hydrophobic", "van-der-waals"]],
            ['R','415','8.56x56','B','291', ["polar-sidechain-backbone"]],
            ['R','419','8.60x60','B','311', ["polar-sidechain-backbone"]],
            ['R','406','8.47x47','A','392', ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','402','7.57x57','A','392', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['R','401','7.56x56','A','392', ["polar-backbone-sidechain", "van-der-waals"]],
            ['R','251','3.54x54','A','391', ["hydrophobic", "van-der-waals"]],
            ['R','419','8.60x60','B','292', ["polar-sidechain-backbone"]],
            ['R','170','12.48x48','B','52', ["polar-sidechain-sidechain"]],
            ],
        '6b3j' : [
            ['R','356','6.45x45','A','393', ["hydrophobic"]],
            ['R','262','4.38x39','A','38', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','250','3.53x53','A','391', ["van-der-waals"]],
            ['R','334','5.64x64','A','384', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','327','5.57x57','A','393', ["hydrophobic"]],
            ['R','258','-','A','383', ["hydrophobic"]],
            ['R','261','-','A','35', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','176','2.46x46','A','391', ["cation-pi", "hydrophobic"]],
            ['R','352','6.41x41','A','393', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
            ['R','334','5.64x64','A','385', ["hydrophobic", "polar-backbone-sidechain"]],
            ['R','262','4.38x39','A','34', ["polar-sidechain-sidechain"]],
            ['R','423','8.64x64','B','44', ["hydrophobic"]],
            ['R','259','-','A','217', ["hydrophobic"]],
            ['R','334','5.64x64','A','381', ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
            ['R','352','6.41x41','A','394', ["polar-sidechain-backbone"]],
            ['R','348','6.37x37','A','394', ["hydrophobic"]],
            ['R','254','3.57x57','A','387', ["hydrophobic", "van-der-waals"]],
            ['R','180','2.50x50','A','391', ["hydrophobic"]],
            ['R','331','5.61x61','A','393', ["hydrophobic"]],
            ['R','334','5.64x64','A','388', ["hydrophobic"]],
            ['R','415','8.56x56','B','312', ["polar-sidechain-sidechain", "h-bond"]],
            ['R','255','3.58x58','A','384', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['R','171','12.49x49','B','312', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','251','3.54x54','A','391', ["hydrophobic", "van-der-waals"]],
            ['R','331','5.61x61','A','388', ["hydrophobic"]],
            ['R','176','2.46x46','A','390', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain"]],
            ['R','331','5.61x61','A','394', ["hydrophobic"]],
            ['R','255','3.58x58','A','387', ["hydrophobic"]],
            ['R','419','8.60x60','B','309', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
            ['R','255','3.58x58','A','388', ["hydrophobic"]],
            ['R','407','8.48x48','A','392', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
            ],
        '6cmo' : [
            ['R','246','6.29x29','A','348', ["hydrophobic"]],
            ['R','141','3.56x56','A','193', ["polar-sidechain-sidechain"]],
            ['R','237','5.72x72','A','320', ["polar-sidechain-sidechain"]],
            ['R','241','6.24x24','A','318', ["polar-backbone-sidechain"]],
            ['R','245','6.28x28','A','354', ["hydrophobic"]],
            ['R','242','6.25x25','A','354', ["hydrophobic", "van-der-waals"]],
            ['R','249','6.32x32','A','354', ["hydrophobic"]],
            ['R','246','6.29x29','A','354', ["hydrophobic", "van-der-waals"]],
            ['R','240','-','A','345', ["polar-sidechain-sidechain", "h-bond"]],
            ['R','239','-','A','318', ["polar-backbone-sidechain"]],
            ['R','135','3.50x50','A','353', ["hydrophobic"]],
            ['R','66','12.48x48','B','312', ["polar-sidechain-sidechain"]],
            ['R','240','-','A','318', ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
            ['R','139','3.54x54','A','347', ["hydrophobic"]],
            ['R','237','5.72x72','A','341', ["hydrophobic"]],
            ['R','243','6.26x26','A','341', ["polar-sidechain-sidechain"]],
            ['R','311','8.48x48','A','353', ["polar-sidechain-backbone"]],
            ['R','311','8.48x48','A','352', ["polar-sidechain-backbone"]],
            ['R','242','6.25x25','A','315', ["polar-sidechain-backbone"]],
            ['R','253','6.36x36','A','353', ["hydrophobic", "van-der-waals"]],
            ['R','249','6.32x32','A','353', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
            ['R','147','34.55x55','A','31', ["polar-sidechain-backbone"]],
            ['R','147','34.55x55','A','32', ["hydrophobic", "polar-sidechain-backbone"]],
            ['R','309','7.56x56','A','352', ["van-der-waals"]],
            ['R','139','3.54x54','A','348', ["hydrophobic"]],
            ['R','311','8.48x48','A','349', ["polar-sidechain-backbone"]],
            ['R','250','6.33x33','A','353', ["hydrophobic"]],
            ['R','310','8.47x47','A','352', ["hydrophobic", "polar-sidechain-backbone"]],
            ['R','239','-','A','320', ["hydrophobic", "van-der-waals"]],
            ['R','311','8.48x48','A','354', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
            ],
        '6d9h' : [
            ['R','45','2.40x40','A','351', ["van-der-waals"]],
            ['R','38','12.48x48','B','335', ["hydrophobic", "van-der-waals"]],
            ['R','207','5.65x65','A','354', ["hydrophobic"]],
            ['R','232','6.33x33','A','354', ["hydrophobic", "van-der-waals"]],
            ['R','113','34.51x51','A','337', ["hydrophobic"]],
            ['R','228','6.29x29','A','316', ["polar-sidechain-sidechain"]],
            ['R','294','8.49x49','A','351', ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
            ['R','203','5.61x61','A','354', ["hydrophobic"]],
            ['R','37','1.60x60','B','312', ["polar-sidechain-sidechain"]],
            ['R','113','34.51x51','A','344', ["hydrophobic"]],
            ['R','228','6.29x29','A','355', ["cation-pi", "hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','210','5.68x68','A','345', ["hydrophobic", "van-der-waals"]],
            ['R','210','5.68x68','A','342', ["polar-sidechain-sidechain", "van-der-waals"]],
            ['R','207','5.65x65','A','349', ["hydrophobic"]],
            ['R','292','8.47x47','A','353', ["hydrophobic"]],
            ['R','113','34.51x51','A','341', ["hydrophobic"]],
            ['R','224','6.25x25','A','319', ["polar-sidechain-sidechain"]],
            ['R','236','6.37x37','A','354', ["hydrophobic", "van-der-waals"]],
            ['R','224','6.25x25','A','316', ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
            ['R','112','34.50x50','A','345', ["hydrophobic", "van-der-waals"]],
            ['R','105','3.50x50','A','354', ["hydrophobic"]],
            ['R','109','3.54x54','A','349', ["hydrophobic", "van-der-waals"]],
            ['R','113','34.51x51','A','195', ["hydrophobic"]],
            ['R','292','8.47x47','A','352', ["hydrophobic", "van-der-waals"]],
            ['R','213','5.71x71','A','342', ["polar-sidechain-sidechain"]],
            ['R','231','6.32x32','A','355', ["polar-sidechain-backbone"]],
            ['R','112','34.50x50','A','344', ["hydrophobic", "van-der-waals"]],
            ['R','291','7.56x56','A','353', ["hydrophobic"]],
            ['R','301','8.56x56','B','312', ["polar-sidechain-backbone"]],
            ['R','232','6.33x33','A','355', ["hydrophobic"]],
            ['R','42','2.37x37','A','351', ["polar-sidechain-sidechain"]],
            ['R','105','3.50x50','A','352', ["hydrophobic", "polar-sidechain-backbone"]],
            ['R','112','34.50x50','A','348', ["polar-backbone-sidechain"]],
            ['R','108','3.53x53','A','352', ["polar-sidechain-sidechain", "van-der-waals"]],
            ['R','108','3.53x53','A','351', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','108','3.53x53','A','348', ["hydrophobic", "polar-backbone-sidechain", "polar-sidechain-sidechain"]],
            ['R','211','5.69x69','A','346', ["hydrophobic"]],
            ['R','301','8.56x56','B','292', ["cation-pi", "hydrophobic", "van-der-waals"]],
            ],
        '6dde' : [
            ['R','263','-','A','341', ["polar-backbone-sidechain"]],
            ['R','176','34.54x54','A','343', ["hydrophobic", "van-der-waals"]],
            ['R','173','34.51x51','A','194', ["hydrophobic"]],
            ['R','264','-','A','341', ["hydrophobic", "polar-backbone-sidechain"]],
            ['R','258','5.64x64','A','344', ["hydrophobic", "van-der-waals"]],
            ['R','268','6.23x23','A','315', ["polar-backbone-sidechain"]],
            ['R','264','-','A','316', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
            ['R','281','6.36x36','A','353', ["hydrophobic", "van-der-waals"]],
            ['R','271','6.26x26','A','314', ["polar-sidechain-backbone", "van-der-waals"]],
            ['R','177','34.55x55','A','32', ["polar-backbone-sidechain", "polar-sidechain-sidechain"]],
            ['R','263','-','A','319', ["polar-sidechain-backbone"]],
            ['R','165','3.50x50','A','353', ["hydrophobic"]],
            ['R','278','6.33x33','A','354', ["hydrophobic", "van-der-waals"]],
            ['R','165','3.50x50','A','351', ["polar-sidechain-sidechain", "van-der-waals"]],
            ['R','277','6.32x32','A','353', ["polar-sidechain-backbone"]],
            ['R','259','5.65x65','A','348', ["hydrophobic"]],
            ['R','172','34.50x50','A','343', ["hydrophobic", "van-der-waals"]],
            ['R','255','5.61x61','A','353', ["hydrophobic"]],
            ['R','103','2.39x39','A','351', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','176','34.54x54','A','194', ["hydrophobic"]],
            ['R','169','3.54x54','A','348', ["hydrophobic", "van-der-waals"]],
            ['R','103','2.39x39','A','350', ["polar-sidechain-backbone"]],
            ['R','264','-','A','345', ["hydrophobic"]],
            ['R','262','5.68x68','A','344', ["hydrophobic"]],
            ['R','172','34.50x50','A','344', ["hydrophobic", "van-der-waals"]],
            ['R','278','6.33x33','A','348', ["hydrophobic"]],
            ['R','182','4.40x40','A','24', ["hydrophobic", "polar-sidechain-sidechain"]],
            ['R','176','34.54x54','A','32', ["hydrophobic", "polar-backbone-sidechain"]],
            ['R','278','6.33x33','A','353', ["hydrophobic"]],
            ['R','271','6.26x26','A','315', ["hydrophobic", "van-der-waals"]],
            ['R','259','5.65x65','A','344', ["hydrophobic"]],
            ['R','179','34.57x57','A','347', ["polar-sidechain-sidechain", "van-der-waals"]],
            ['R','340','8.47x47','A','352', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
            ['R','263','-','A','320', ["cation-pi", "hydrophobic", "van-der-waals"]],
            ['R','168','3.53x53','A','347', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['R','173','34.51x51','A','336', ["hydrophobic", "van-der-waals"]],
            ['R','270','6.25x25','A','315', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','172','34.50x50','A','340', ["hydrophobic"]],
            ['R','173','34.51x51','A','193', ["hydrophobic", "van-der-waals"]],
            ['R','262','5.68x68','A','341', ["hydrophobic", "van-der-waals"]],
            ],
        '6ddf' : [
            ['R','263','-','A','341', ["hydrophobic", "polar-backbone-sidechain"]],
            ['R','176','34.54x54','A','343', ["hydrophobic"]],
            ['R','271','6.26x26','A','317', ["polar-sidechain-backbone"]],
            ['R','263','-','A','319', ["polar-sidechain-backbone", "van-der-waals"]],
            ['R','259','5.65x65','A','348', ["hydrophobic", "van-der-waals"]],
            ['R','271','6.26x26','A','314', ["polar-sidechain-backbone", "van-der-waals"]],
            ['R','177','34.55x55','A','32', ["polar-sidechain-sidechain"]],
            ['R','341','8.48x48','A','354', ["hydrophobic"]],
            ['R','165','3.50x50','A','353', ["hydrophobic"]],
            ['R','165','3.50x50','A','351', ["hydrophobic", "polar-sidechain-backbone"]],
            ['R','271','6.26x26','A','315', ["hydrophobic", "van-der-waals"]],
            ['R','103','2.39x39','A','351', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain"]],
            ['R','263','-','A','318', ["van-der-waals"]],
            ['R','179','34.57x57','A','347', ["polar-sidechain-sidechain"]],
            ['R','103','2.39x39','A','350', ["polar-sidechain-backbone", "van-der-waals"]],
            ['R','164','3.49x49','A','351', ["polar-sidechain-sidechain"]],
            ['R','278','6.33x33','A','348', ["hydrophobic"]],
            ['R','340','8.47x47','A','352', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
            ['R','173','34.51x51','A','194', ["hydrophobic"]],
            ['R','278','6.33x33','A','353', ["hydrophobic"]],
            ['R','172','34.50x50','A','347', ["polar-backbone-sidechain"]],
            ['R','174','34.52x52','A','193', ["polar-sidechain-sidechain"]],
            ['R','179','34.57x57','A','351', ["polar-sidechain-sidechain", "van-der-waals"]],
            ['R','263','-','A','320', ["cation-pi", "hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','172','34.50x50','A','343', ["hydrophobic", "van-der-waals"]],
            ['R','172','34.50x50','A','340', ["hydrophobic", "van-der-waals"]],
            ['R','169','3.54x54','A','348', ["hydrophobic", "van-der-waals"]],
            ['R','264','-','A','341', ["polar-backbone-sidechain"]],
            ['R','258','5.64x64','A','344', ["hydrophobic", "van-der-waals"]],
            ['R','176','34.54x54','A','32', ["hydrophobic"]],
            ['R','264','-','A','316', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
            ['R','281','6.36x36','A','353', ["hydrophobic", "van-der-waals"]],
            ['R','277','6.32x32','A','353', ["polar-sidechain-backbone", "van-der-waals"]],
            ['R','173','34.51x51','A','336', ["hydrophobic", "van-der-waals"]],
            ['R','264','-','A','345', ["hydrophobic"]],
            ['R','277','6.32x32','A','354', ["polar-sidechain-backbone", "van-der-waals"]],
            ['R','255','5.61x61','A','353', ["hydrophobic"]],
            ['R','341','8.48x48','A','353', ["polar-sidechain-backbone"]],
            ['R','168','3.53x53','A','351', ["van-der-waals"]],
            ['R','176','34.54x54','A','194', ["hydrophobic"]],
            ['R','262','5.68x68','A','344', ["hydrophobic"]],
            ['R','172','34.50x50','A','344', ["hydrophobic", "van-der-waals"]],
            ['R','340','8.47x47','A','353', ["polar-sidechain-backbone"]],
            ['R','278','6.33x33','A','354', ["hydrophobic", "van-der-waals"]],
            ['R','341','8.48x48','A','352', ["polar-sidechain-backbone", "van-der-waals"]],
            ['R','173','34.51x51','A','193', ["hydrophobic"]],
            ['R','270','6.25x25','A','315', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
            ['R','168','3.53x53','A','347', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['R','259','5.65x65','A','344', ["hydrophobic"]],
            ['R','262','5.68x68','A','341', ["hydrophobic", "van-der-waals"]],
            ['R','268','6.23x23','A','315', ["polar-backbone-sidechain", "polar-sidechain-sidechain", "van-der-waals"]],
            ],
        '6g79' : [
            ['S','308','6.29x29','A','306', ["polar-sidechain-backbone"]],
            ['S','373','8.47x47','A','340', ["polar-sidechain-backbone"]],
            ['S','151','3.54x54','A','338', ["hydrophobic", "van-der-waals"]],
            ['S','311','6.32x32','A','342', ["polar-sidechain-backbone"]],
            ['S','312','6.33x33','A','343', ["hydrophobic", "van-der-waals"]],
            ['S','154','34.50x50','A','333', ["hydrophobic"]],
            ['S','372','7.56x56','A','342', ["hydrophobic"]],
            ['S','147','3.50x50','A','341', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
            ['S','239','5.69x69','A','344', ["hydrophobic"]],
            ['S','150','3.53x53','A','337', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['S','155','34.51x51','A','330', ["hydrophobic"]],
            ['S','372','7.56x56','A','341', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
            ['S','315','6.36x36','A','342', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
            ['S','161','34.57x57','A','337', ["polar-sidechain-sidechain"]],
            ['S','231','5.61x61','A','343', ["hydrophobic", "van-der-waals"]],
            ['S','151','3.54x54','A','334', ["hydrophobic"]],
            ['S','315','6.36x36','A','343', ["hydrophobic", "polar-sidechain-backbone"]],
            ['S','316','6.37x37','A','343', ["hydrophobic"]],
            ['S','238','5.68x68','A','338', ["hydrophobic"]],
            ['S','238','5.68x68','A','334', ["hydrophobic", "van-der-waals"]],
            ['S','238','5.68x68','A','344', ["hydrophobic"]],
            ['S','235','5.65x65','A','343', ["hydrophobic"]],
            ['S','311','6.32x32','A','344', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
            ['S','154','34.50x50','A','334', ["hydrophobic", "van-der-waals"]],
            ['S','235','5.65x65','A','338', ["hydrophobic"]],
            ['S','238','5.68x68','A','331', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
            ['S','155','34.51x51','A','326', ["hydrophobic"]],
            ['S','308','6.29x29','A','344', ["cation-pi", "hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "van-der-waals"]],
            ['S','151','3.54x54','A','337', ["hydrophobic"]],
            ],
        '6gdg' : [
            ['A','106','3.54x54','D','378', ["hydrophobic", "van-der-waals"]],
            ['A','112','34.53x53','D','377', ["edge-to-face", "face-to-edge", "pi-cation", "hydrophobic"]],
            ['A','207','5.68x68','D','374', ["hydrophobic", "polar-sidechain-sidechain", "van-der-waals"]],
            ['A','106','3.54x54','D','374', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['A','110','34.51x51','D','373', ["hydrophobic"]],
            ['A','203','5.64x64','D','378', ["hydrophobic"]],
            ['A','109','34.50x50','D','374', ["hydrophobic", "van-der-waals"]],
            ['A','102','3.50x50','D','381', ["cation-pi", "hydrophobic", "van-der-waals"]],
            ['A','207','5.68x68','D','384', ["hydrophobic"]],
            ['A','113','34.54x54','D','39', ["hydrophobic"]],
            ['A','105','3.53x53','D','377', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['A','204','5.65x65','D','378', ["hydrophobic"]],
            ['A','200','5.61x61','D','383', ["hydrophobic", "van-der-waals"]],
            ['A','110','34.51x51','D','369', ["hydrophobic"]],
            ['A','106','3.54x54','D','381', ["hydrophobic"]],
            ['A','109','34.50x50','D','370', ["hydrophobic"]],
            ['A','211','5.72x72','D','348', ["hydrophobic"]],
            ['A','293','8.48x48','D','380', ["polar-backbone-sidechain"]],
            ['A','207','5.68x68','D','375', ["polar-sidechain-backbone"]],
            ['A','291','7.56x56','D','381', ["polar-sidechain-backbone"]],
            ['A','110','34.51x51','D','370', ["hydrophobic"]],
            ['A','113','34.54x54','D','38', ["hydrophobic", "polar-sidechain-backbone", "van-der-waals"]],
            ['A','106','3.54x54','D','377', ["hydrophobic", "van-der-waals"]],
            ['A','294','8.49x49','D','380', ["polar-backbone-sidechain"]],
            ['A','110','34.51x51','D','217', ["hydrophobic"]],
            ['A','293','8.48x48','D','382', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['A','207','5.68x68','D','371', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
            ['A','111','34.52x52','D','217', ["hydrophobic", "van-der-waals"]],
            ['A','207','5.68x68','D','378', ["hydrophobic"]],
            ['A','34','1.60x60','B','312', ["polar-sidechain-sidechain"]],
            ['A','114','34.55x55','D','39', ["hydrophobic"]],
            ['A','230','6.32x32','D','382', ["polar-sidechain-backbone"]],
            ['A','107','3.55x55','D','370', ["polar-backbone-sidechain", "polar-sidechain-sidechain", "van-der-waals"]],
            ['A','234','6.36x36','D','383', ["hydrophobic", "van-der-waals"]],
            ['A','109','34.50x50','D','373', ["hydrophobic", "van-der-waals"]],
            ['A','110','34.51x51','D','366', ["hydrophobic", "van-der-waals"]],
            ['A','36','12.49x49','B','312', ["hydrophobic", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
            ['A','208','5.69x69','D','384', ["hydrophobic"]],
            ['A','35','12.48x48','B','333', ["hydrophobic", "polar-sidechain-backbone", "polar-sidechain-sidechain", "h-bond", "van-der-waals"]],
            ['A','35','12.48x48','B','335', ["hydrophobic", "van-der-waals"]],
            ['A','35','12.48x48','B','312', ["hydrophobic"]],
            ['A','291','7.56x56','D','382', ["polar-backbone-sidechain", "polar-sidechain-backbone", "van-der-waals"]],
            ['A','231','6.33x33','D','383', ["hydrophobic", "van-der-waals"]],
            ['A','111','34.52x52','D','215', ["polar-sidechain-backbone"]],
            ['A','111','34.52x52','D','216', ["hydrophobic", "polar-sidechain-backbone"]],
            ['A','110','34.51x51','D','41', ["hydrophobic", "polar-backbone-sidechain", "van-der-waals"]],
            ['A','235','6.37x37','D','383', ["hydrophobic", "van-der-waals"]],
            ['A','210','5.71x71','D','371', ["polar-sidechain-sidechain"]],
            ['A','234','6.36x36','D','382', ["polar-sidechain-backbone"]],
            ['A','38','12.51x51','B','52', ["polar-sidechain-sidechain", "van-der-waals"]],
            ]
}

    complex_info = [
        {
            'pdb_id': '3sn6',
            'receptor': 'beta2',
            'gprotein': 'GNAS2_BOVIN',
            'alternative_gprotein': 'GNAS2_HUMAN'
            },
        {
            'pdb_id': '4x1h',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '5g53',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '5g53',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '5uz7',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '5vai',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '6b3j',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '6cmo',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '6d9h',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '6dde',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '6ddf',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '6g79',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
        {
            'pdb_id': '6gdg',
            'receptor': '',
            'gprotein': '',
            'alternative_protein': ''
            },
    ]

    # rs = Residue.objects.filter(protein_conformation__protein=proteins[0]).prefetch_related('protein_segment','display_generic_number','generic_number')

    interactions = SignprotInteractions.objects.all().values_list(
        'gpcr_residue__sequence_number',
        'gpcr_residue__display_generic_number__label',
        'structure__pdb_code__index','interaction_type',
        'signprot_residue__sequence_number',
        'signprot_residue__display_generic_number__label'
        )

    interactions_metadata = complex_info
    context = {
        'interactions': json.dumps(list(interactions)),
        'interactions_metadata': interactions_metadata
        }

    return render(request, 'signprot/matrix.html', context)
