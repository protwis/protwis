from django.shortcuts import render
from django.views.generic import TemplateView, View
from django.http import HttpResponse
from django.views.decorators.cache import cache_page
from django.views.decorators.csrf import csrf_exempt

from common.diagrams_gpcr import DrawSnakePlot

from construct.models import *
from construct.functions import *
from protein.models import Protein, ProteinConformation
from structure.models import Structure
from mutation.models import Mutation

from datetime import datetime
import json
import copy
from collections import OrderedDict


# Create your views here.
@cache_page(60 * 60 * 24)
def detail(request, slug):

    # get constructs
    c = Construct.objects.filter(name=slug).prefetch_related(
        'modifications', 'mutations', 'deletions','insertions','insertions__insert_type',
        'expression','solubilization','purification','crystallization','crystal').all()[0]

    # get residues
    residues = Residue.objects.filter(protein_conformation__protein=c.protein).order_by('sequence_number').prefetch_related(
        'protein_segment', 'generic_number', 'display_generic_number')

    residues_lookup = {}
    for r in residues:
        residues_lookup[r.sequence_number] = r

    schematics = c.schematic()

    snake = cache.get(c.name+'_snake')
    if snake==None:
        print(c.name+'_snake no cache')
        snake = cache.set(c.name+'_snake', DrawSnakePlot(residues,c.protein.get_protein_class(),str(c.protein),nobuttons = True), 60*60*24*2) #two days
        snake = cache.get(c.name+'_snake')
    else:
        print(c.name+'_snake used cache')

    chunk_size = 10
    context = {'c':c, 'chunk_size': chunk_size, 'snake': snake, 'annotations': json.dumps(schematics['annotations']), 'schematics': schematics, 'residues_lookup': residues_lookup}
    return render(request,'construct/construct_detail.html',context)

def fetch_all_pdb(request):

    structures = Structure.objects.all()

    for s in structures:
        pdbname = str(s)
        print(pdbname)
        try:
            protein = Protein.objects.filter(entry_name=pdbname.lower()).get()
            d = fetch_pdb_info(pdbname,protein)

            #delete before adding new
            Construct.objects.filter(name=d['construct_crystal']['pdb_name']).delete()
            add_construct(d)
        except:
            print(pdbname,'failed')


    # d = fetch_pdb_info(slug,protein)

    # #delete before adding new
    # Construct.objects.filter(name=d['construct_crystal']['pdb_name']).delete()
    # add_construct(d)

    context = {'d':d}

    return render(request,'pdb_fetch.html',context)

def fetch_pdb(request, slug):

    protein = Protein.objects.filter(entry_name=slug.lower()).get()

    d = fetch_pdb_info(slug,protein)

    #delete before adding new
    print(d['construct_crystal']['pdb_name'])
    Construct.objects.filter(name__iexact=d['construct_crystal']['pdb_name']).delete()
    add_construct(d)
    context = {'d':d}

    return render(request,'pdb_fetch.html',context)

def fetch_pdb_for_webform(request, slug, **response_kwargs):

    slug = slug.lower()
    protein = Protein.objects.filter(entry_name=slug).get()

    d = fetch_pdb_info(slug,protein)
    d = convert_ordered_to_disordered_annotation(d)

    jsondata = json.dumps(d)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)

class ConstructBrowser(TemplateView):
    """
    Fetching construct data for browser
    """

    template_name = "construct_browser.html"

    def get_context_data (self, **kwargs):

        context = super(ConstructBrowser, self).get_context_data(**kwargs)
        try:
            cons = Construct.objects.all().prefetch_related(
                "crystal","mutations","purification","protein__family__parent__parent__parent", "insertions", "modifications", "deletions", "crystallization__chemical_lists",
                "protein__species","structure__pdb_code","structure__publication__web_link", "contributor")

            context['constructs'] = cache.get('construct_browser')
            if context['constructs']==None:
                context['constructs'] = []
                for c in cons:
                    c.schematics = c.schematic()
                    context['constructs'].append(c)

                cache.set('construct_browser', context['constructs'], 60*60*24*2) #two days
            else:
                print('construct_browser used cache')

        except Construct.DoesNotExist as e:
            pass

        return context

@csrf_exempt #jquery send post, so no csrf
def align(request):

    ids = json.loads(request.POST.get('ids'))

    c_ids = []
    s_ids = []
    for i in ids:
        if i.startswith('align'):
            s_ids.append(i.split('_')[1])
        else:
            c_ids.append(i)

    cons = Construct.objects.filter(pk__in=c_ids).prefetch_related(
                "crystal","mutations","purification","protein__family__parent__parent__parent", "insertions", "modifications", "deletions", "crystallization__chemical_lists",
                "protein__species","structure__pdb_code","structure__publication__web_link", "contributor")
    proteins = []
    constructs = OrderedDict()
    annotations = {}
    for c in cons:
        # print(c)
        proteins.append(c.protein)
        constructs[c.name] = c.protein.entry_name
        annotations[c.name] = c.schematic()['annotations']

    print(annotations)

    if len(s_ids):
        rs = Residue.objects.filter(protein_conformation__protein__in=proteins, protein_segment__slug__in=s_ids).prefetch_related(
        'protein_conformation__protein', 'protein_conformation__state', 'protein_segment',
        'generic_number__scheme', 'display_generic_number__scheme')
    else:
        s_ids = ['N-term','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','ICL4','H8','C-term']
        rs = Residue.objects.filter(protein_conformation__protein__in=proteins).prefetch_related(
        'protein_conformation__protein', 'protein_conformation__state', 'protein_segment',
        'generic_number__scheme', 'display_generic_number__scheme')

    print("residues",len(rs))

    distinct_gn = []
    ordered_gn = OrderedDict()
    distinct_segments = []
    overview = OrderedDict()
    segment_length = OrderedDict()
    for s in s_ids:
        overview[s] = OrderedDict()
        segment_length[s] = {'aligned':0, 'before':0,'after':0,'total':0}

    protein_lookup = {}
    print('build stuff')

    segment = ''
    protein = ''

    track_unaligned = {}

    #Find all possible generic numbers, to ensure gaps
    for r in rs.order_by('protein_conformation__id','sequence_number'):
        if segment!=r.protein_segment.slug or protein!=r.protein_conformation.protein.entry_name:
            no_encountered_gn = True
            length = 0
            length_before = 0
            length_after = 0

        segment = r.protein_segment.slug
        protein = r.protein_conformation.protein.entry_name

        if protein not in protein_lookup:
            protein_lookup[protein] = {}
            track_unaligned[protein] = {}

        if segment not in track_unaligned[protein]:
            track_unaligned[protein][segment] = {'before':[],'after':[]}

        if segment not in distinct_segments:
            distinct_segments.append(segment)
            overview[segment] = OrderedDict()

        if r.generic_number:
            no_encountered_gn = False
            gn = r.generic_number.label
            protein_lookup[protein][gn] = {'aa':r.amino_acid,'pos':r.sequence_number}
            gn_sort = gn.split('x')[1]
            gn_sort = float("0."+gn_sort)
            if gn not in distinct_gn:
                distinct_gn.append(gn)
                overview[segment][gn_sort] = [gn,{'aa':'-','pos':''}]
            length += 1
        else:
            if no_encountered_gn:
                track_unaligned[protein][segment]['before'].append({'aa':r.amino_acid,'pos':r.sequence_number})
                length_before += 1
            else:
                track_unaligned[protein][segment]['after'].append({'aa':r.amino_acid,'pos':r.sequence_number})
                length_after += 1

        if len(overview[segment])>segment_length[segment]['aligned']:
            segment_length[segment]['aligned'] = len(overview[segment])
        if length_before>segment_length[segment]['before']:
            segment_length[segment]['before'] = length_before
        if length_after>segment_length[segment]['after']:
            segment_length[segment]['after'] = length_after
        if segment_length[segment]['aligned']+segment_length[segment]['before']+segment_length[segment]['after']>segment_length[segment]['total']:
            segment_length[segment]['total'] = segment_length[segment]['aligned']+segment_length[segment]['before']+segment_length[segment]['after']

    # SORT generic residues to ensure correct order
    gn_list = ""
    ordered_summary = OrderedDict()
    for seg,gns in overview.items():
        ordered_summary[seg] = OrderedDict()
        #GN LIST
        gn_list += """<td class="ali-td ali-residue res-color-X">&nbsp;</td>"""
        if seg!='C-term':
            for _ in range(segment_length[seg]['before']):
                gn_list += """<td class="ali-td">&nbsp;</td>"""
        for gn in sorted(gns):
            ordered_summary[seg][gns[gn][0]] = {'aa':'-','pos':''}
            gn_list += """<td class="ali-td-generic-num">{}</td>""".format("x"+gns[gn][0].split("x")[1])

        if seg=='C-term':
            for _ in range(segment_length[seg]['before']):
                gn_list += """<td class="ali-td">&nbsp;</td>"""

        for _ in range(segment_length[seg]['after']):
            gn_list += """<td class="ali-td">&nbsp;</td>"""

    alignment = OrderedDict()
    alignment_print_sequence = ""
    for c,p in constructs.items():
        alignment[c] = copy.deepcopy(ordered_summary)
        alignment_print_sequence += '<tr>'
        for seg,gns in alignment[c].items():

            if p not in track_unaligned:
                track_unaligned[p] = {seg: {'before':[],'after':[]}}

            if p not in protein_lookup:
                protein_lookup[p] = {}

            if seg not in track_unaligned[p]:
                track_unaligned[p][seg] = {'before':[],'after':[]}

            alignment_print_sequence += """<td class="ali-td ali-residue res-color-_">&nbsp;</td>"""

            if seg!='C-term':
                for _ in range(segment_length[seg]['before']-len(track_unaligned[p][seg]['before'])):
                    alignment_print_sequence += """<td class="ali-td ali-residue res-color-">
                                                    -</td>"""

            for aa in track_unaligned[p][seg]['before']:
                if aa['pos'] in annotations[c]:
                    annotation = annotations[c][aa['pos']][0]
                    annotation_text = "<br>"+annotations[c][aa['pos']][1]
                else:
                    annotation = ''
                    annotation_text = ''

                alignment_print_sequence += """<td class="ali-td ali-residue res-color-{}">
                                                <div data-toggle="tooltip" data-placement="top" data-html="true"
                                                title="{}{}{}">{}</div></td>""".format(annotation,aa['aa'],aa['pos'],annotation_text,aa['aa'])
            for gn, aa in gns.items():
                if gn in protein_lookup[p]:
                    aa = protein_lookup[p][gn]
                    alignment[c][seg][gn] = aa

                if aa['pos'] in annotations[c]:
                    annotation = annotations[c][aa['pos']][0]
                    annotation_text = "<br>"+annotations[c][aa['pos']][1]
                else:
                    annotation = ''
                    annotation_text = ''
                alignment_print_sequence += """<td class="ali-td ali-residue res-color-{}">
                                                <div data-toggle="tooltip" data-placement="top" data-html="true"
                                                title="{}{}<br>SCHEME: {}{}">{}</div></td>""".format(annotation,aa['aa'],aa['pos'],gn,annotation_text,aa['aa'])

            for aa in track_unaligned[p][seg]['after']:
                if aa['pos'] in annotations[c]:
                    annotation = annotations[c][aa['pos']][0]
                    annotation_text = "<br>"+annotations[c][aa['pos']][1]
                else:
                    annotation = ''
                    annotation_text = ''
                alignment_print_sequence += """<td class="ali-td ali-residue res-color-{}">
                                                <div data-toggle="tooltip" data-placement="top" data-html="true"
                                                title="{}{}{}">{}</div></td>""".format(annotation,aa['aa'],aa['pos'],annotation_text,aa['aa'])

            for _ in range(segment_length[seg]['after']-len(track_unaligned[p][seg]['after'])):
                alignment_print_sequence += """<td class="ali-td ali-residue res-color-">
                                                -</td>"""
            if seg=='C-term':
                for _ in range(segment_length[seg]['before']-len(track_unaligned[p][seg]['before'])):
                    alignment_print_sequence += """<td class="ali-td ali-residue res-color-">
                                                    -</td>"""

        alignment_print_sequence += '</tr>'

    print('done',len(alignment_print_sequence))
    context = {'constructs': constructs,'alignment_print_sequence': alignment_print_sequence, 'segment_length' : segment_length, 'gn_list' : gn_list, 'segments': s_ids, 'c_ids': json.dumps(c_ids)} #, 'alignment_print_sequence': alignment_print_sequence

    return render(request,'align.html',context)
