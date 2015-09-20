from django.shortcuts import get_object_or_404, render
from django.views import generic
from django.http import HttpResponse
from django.db.models import Q

from protein.models import Protein, ProteinConformation, ProteinAlias, ProteinFamily, Gene
from residue.models import Residue
from structure.models import Structure
from common.selection import Selection
from common.views import AbsBrowseSelection

import json
from collections import OrderedDict


class BrowseSelection(AbsBrowseSelection):
    title = 'SELECT A RECEPTOR (FAMILY)'
    description = 'Select a target or family by searching or browsing in the right column.'
    description = 'Select a receptor (family) by searching or browsing in the middle. The selection is viewed to' \
        + ' the right.'
    docs = '/docs/browse'
    buttons = {}
        

def detail(request, slug):
    # get protein
    p = Protein.objects.prefetch_related('web_links__web_resource').get(entry_name=slug, sequence_type__slug='wt')

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
    structures = Structure.objects.filter(protein_conformation__protein__parent=p)

    # get residues
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
        'gene': gene, 'alt_genes': alt_genes, 'structures': structures}

    return render(request, 'protein/protein_detail.html', context)

def SelectionAutocomplete(request):
    if request.is_ajax():
        q = request.GET.get('term')
        type_of_selection = request.GET.get('type_of_selection')
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

        if type_of_selection == 'targets' or type_of_selection == 'browse':
            # find protein families
            pfs = ProteinFamily.objects.filter(name__icontains=q).exclude(slug='000')[:10]
            for pf in pfs:
                pf_json = {}
                pf_json['id'] = pf.id
                pf_json['label'] = pf.name
                pf_json['slug'] = pf.slug
                pf_json['type'] = 'family'
                pf_json['category'] = 'Target families'
                results.append(pf_json)
        
        # find proteins
        ps = Protein.objects.filter(Q(name__icontains=q) | Q(entry_name__icontains=q) | Q(family__name__icontains=q),
            species__in=(species_list),
            source__in=(protein_source_list))[:10]
        for p in ps:
            p_json = {}
            p_json['id'] = p.id
            p_json['label'] = p.name + " [" + p.species.common_name + "]"
            p_json['slug'] = p.entry_name
            p_json['type'] = 'protein'
            p_json['category'] = 'Targets'
            results.append(p_json)

        # find protein aliases
        pas = ProteinAlias.objects.prefetch_related('protein').filter(name__icontains=q,
            protein__species__in=(species_list),
            protein__source__in=(protein_source_list))[:10]
        for pa in pas:
            pa_json = {}
            pa_json['id'] = pa.protein.id
            pa_json['label'] = pa.protein.name  + " [" + pa.protein.species.common_name + "]"
            pa_json['slug'] = pa.protein.entry_name
            pa_json['type'] = 'protein'
            pa_json['category'] = 'Targets'
            if pa_json not in results:
                results.append(pa_json)
        
        data = json.dumps(results)
    else:
        data = 'fail'
    mimetype = 'application/json'
    return HttpResponse(data, mimetype)