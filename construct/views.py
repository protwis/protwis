from django.shortcuts import render
from django.views.generic import TemplateView, View
from django.http import HttpResponse
from django.views.decorators.cache import cache_page

from common.diagrams_gpcr import DrawSnakePlot

from construct.models import *
from construct.functions import *
from protein.models import Protein, ProteinConformation
from structure.models import Structure
from mutation.models import Mutation

from datetime import datetime
import json


# Create your views here.
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
            context['constructs'] = []
            for c in cons:
                c.schematics = c.schematic()
                context['constructs'].append(c)
        except Construct.DoesNotExist as e:
            pass

        return context
