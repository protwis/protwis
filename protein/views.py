from django.shortcuts import get_object_or_404, render, redirect
from django.views import generic
from django.http import HttpResponse
from django.db.models import Q, F, Func, Value
from django.core.cache import cache
from django.views.decorators.cache import cache_page
from django.urls import reverse

from protein.models import Protein, ProteinConformation, ProteinAlias, ProteinFamily, Gene,ProteinGProteinPair
from residue.models import Residue
from structure.models import Structure, StructureModel
from mutation.models import MutationExperiment
from common.selection import Selection
from common.views import AbsBrowseSelection

import json
from collections import OrderedDict


class BrowseSelection(AbsBrowseSelection):
    title = 'SELECT A RECEPTOR (FAMILY)'
    description = 'Select a target or family by searching or browsing in the right column.'
    description = 'Select a receptor (family) by searching or browsing in the middle. The selection is viewed to' \
        + ' the right.'
    docs = 'receptors.html'
    target_input=False


@cache_page(60 * 60 * 24)
def detail(request, slug):
    # get protein
    slug = slug.lower()
    try:
        if Protein.objects.filter(entry_name=slug).exists():
            p = Protein.objects.prefetch_related('web_links__web_resource').get(entry_name=slug, sequence_type__slug='wt')
        else:
            p = Protein.objects.prefetch_related('web_links__web_resource').get(accession=slug.upper(), sequence_type__slug='wt')
    except:
        context = {'protein_no_found': slug}

        return render(request, 'protein/protein_detail.html', context)


    if p.family.slug.startswith('100') or p.family.slug.startswith('200'):
        # If this protein is a gprotein, redirect to that page.
        return redirect(reverse('signprotdetail', kwargs={'slug': slug}))

    # get family list
    pf = p.family
    families = [pf.name]
    while pf.parent.parent:
        families.append(pf.parent.name)
        pf = pf.parent
    families.reverse()

    # get default conformation
    pc = ProteinConformation.objects.get(protein=p)

    # get protein aliases
    aliases = ProteinAlias.objects.filter(protein=p).values_list('name', flat=True)

    # get genes
    genes = Gene.objects.filter(proteins=p).values_list('name', flat=True)
    gene = genes[0]
    alt_genes = genes[1:]

    # get structures of this protein
    structures = Structure.objects.filter(protein_conformation__protein__parent=p).order_by('-representative',
        'resolution')

    # get residues
    residues = Residue.objects.filter(protein_conformation=pc).order_by('sequence_number').prefetch_related(
        'protein_segment', 'generic_number', 'display_generic_number')

    mutations = MutationExperiment.objects.filter(protein=p)

    protein_links = p.web_links.all().distinct('web_resource__slug')

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

    homology_models = StructureModel.objects.filter(protein=p)

    context = {'p': p, 'families': families, 'r_chunks': r_chunks, 'chunk_size': chunk_size, 'aliases': aliases,
        'gene': gene, 'alt_genes': alt_genes, 'structures': structures, 'mutations': mutations, 'protein_links': protein_links,'homology_models': homology_models}

    return render(request, 'protein/protein_detail.html', context)

def SelectionAutocomplete(request):

    if request.is_ajax():
        q = request.GET.get('term')
        type_of_selection = request.GET.get('type_of_selection')
        selection_only_receptors = request.GET.get('selection_only_receptors')
        referer = request.META.get('HTTP_REFERER')

        if 'gproteinselection' in str(referer) or 'signprot' in str(referer) and not 'ginterface' in str(referer):
            exclusion_slug = '00'
        else:
            exclusion_slug = '100'

        results = []

        # session
        simple_selection = request.session.get('selection')
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)

        # species filter
        species_list = []
        for species in selection.species:
            species_list.append(species.item)

        # annotation filter
        protein_source_list = []
        for protein_source in selection.annotation:
            protein_source_list.append(protein_source.item)

        # find proteins
        if type_of_selection!='navbar':
            ps = Protein.objects.filter(Q(name__icontains=q) | Q(entry_name__icontains=q),
                species__in=(species_list),
                source__in=(protein_source_list)).exclude(family__slug__startswith=exclusion_slug)[:10]
        else:
            ps = Protein.objects.filter(Q(name__icontains=q) | Q(entry_name__icontains=q) | Q(family__name__icontains=q) | Q(accession=q),
                species__common_name='Human', source__name='SWISSPROT').exclude(family__slug__startswith=exclusion_slug)[:10]

        # Try matching protein name after stripping html tags
        if ps.count() == 0:
            ps = Protein.objects.annotate(filtered=Func(F('name'), Value('<[^>]+>'), Value(''), Value('gi'), function='regexp_replace')).filter(Q(filtered__icontains=q), species__common_name='Human', source__name='SWISSPROT')

        for p in ps:
            p_json = {}
            p_json['id'] = p.id
            p_json['label'] = p.name + " [" + p.species.common_name + "]"
            p_json['slug'] = p.entry_name
            p_json['type'] = 'protein'
            p_json['category'] = 'Targets'
            results.append(p_json)


        if type_of_selection!='navbar':
            # find protein aliases
            pas = ProteinAlias.objects.prefetch_related('protein').filter(name__icontains=q,
                protein__species__in=(species_list),
                protein__source__in=(protein_source_list)).exclude(protein__family__slug__startswith=exclusion_slug)[:10]
            for pa in pas:
                pa_json = {}
                pa_json['id'] = pa.protein.id
                pa_json['label'] = pa.protein.name  + " [" + pa.protein.species.common_name + "]"
                pa_json['slug'] = pa.protein.entry_name
                pa_json['type'] = 'protein'
                pa_json['category'] = 'Targets'
                if pa_json not in results:
                    results.append(pa_json)

            # protein families
            if (type_of_selection == 'targets' or type_of_selection == 'browse' or type_of_selection == 'gproteins') and selection_only_receptors!="True":
                # find protein families
                pfs = ProteinFamily.objects.filter(name__icontains=q).exclude(slug='000').exclude(slug__startswith=exclusion_slug)[:10]
                for pf in pfs:
                    pf_json = {}
                    pf_json['id'] = pf.id
                    pf_json['label'] = pf.name
                    pf_json['slug'] = pf.slug
                    pf_json['type'] = 'family'
                    pf_json['category'] = 'Target families'
                    results.append(pf_json)

        data = json.dumps(results)
    else:
        data = 'fail'
    mimetype = 'application/json'
    return HttpResponse(data, mimetype)


def g_proteins(request, **response_kwargs):
    ''' Example of g_proteins '''
    proteins = Protein.objects.filter(source__name='SWISSPROT').prefetch_related('proteingproteinpair_set')
    jsondata = {}
    for p in proteins:
        gps = p.proteingproteinpair_set.all()
        if gps:
            jsondata[str(p)] = []
            for gp in gps:
                jsondata[str(p)].append(str(gp))
    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)
