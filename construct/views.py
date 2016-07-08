from django.shortcuts import render
from django.views.generic import TemplateView, View

from common.diagrams_gpcr import DrawSnakePlot

from construct.models import *
from protein.models import Protein, ProteinConformation
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
        snake = cache.get_or_set(c.name+'_snake', DrawSnakePlot(residues,c.protein.get_protein_class(),str(c.protein),nobuttons = True), 60*60*24*2) #two days
    else:
        print(c.name+'_snake used cache')

    chunk_size = 10
    context = {'c':c, 'chunk_size': chunk_size, 'snake': snake, 'annotations': json.dumps(schematics['annotations']), 'schematics': schematics, 'residues_lookup': residues_lookup}
    return render(request,'construct/construct_detail.html',context)

class ConstructBrowser(TemplateView):
    """
    Fetching construct data for browser
    """

    template_name = "construct_browser.html"

    def get_context_data (self, **kwargs):

        context = super(ConstructBrowser, self).get_context_data(**kwargs)
        try:
            cons = Construct.objects.all().prefetch_related(
                "crystal","mutations","purification")
            context['constructs'] = []
            for c in cons:
                c.schematics = c.schematic()
                context['constructs'].append(c)
        except Construct.DoesNotExist as e:
            pass

        return context
