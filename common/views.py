from django.http import HttpResponse
from django.shortcuts import render
from django.views.generic import TemplateView
from django.conf import settings

from common.selection import SimpleSelection, Selection, SelectionItem
from protein.models import Protein, ProteinFamily, ProteinSegment, Species, ProteinSource, ProteinSet
from residue.models import ResidueGenericNumber, ResidueNumberingScheme

import inspect
from collections import OrderedDict


class AbsTargetSelection(TemplateView):
    """An abstract class for the target selection page used in many apps. To use it in another app, create a class 
    based view for that app that extends this class"""
    template_name = 'common/targetselection.html'

    type_of_selection = 'targets'
    step = 1
    number_of_steps = 2
    title = 'SELECT TARGETS'
    description = 'Select targets by searching or browsing in the middle column. You can select entire target families or individual targets.\n\nSelected targets will appear in the right column, where you can edit the list.\n\nOnce you have selected all your targets, click the green button.'
    docs = False
    filters = True
    numbering_schemes = False
    search = True
    family_tree = True
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '#',
            'color': 'success',
        },
    }
    # OrderedDict to preserve the order of the boxes
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])

    # proteins and families
    #try - except block prevents manage.py from crashing - circular dependencies between protein - common 
    try:
        ppf = ProteinFamily.objects.get(slug='000')
        pfs = ProteinFamily.objects.filter(parent=ppf.id)
        ps = Protein.objects.filter(family=ppf)
        psets = ProteinSet.objects.all().prefetch_related('proteins')
        tree_indent_level = []
        action = 'expand'
        # remove the parent family (for all other families than the root of the tree, the parent should be shown)
        del ppf
    except Exception as e:
        pass

    # species
    sps = Species.objects.all()

    # numbering schemes
    gns = ResidueNumberingScheme.objects.all()

    def get_context_data(self, **kwargs):
        """get context from parent class (really only relevant for children of this class, as TemplateView does
        not have any context variables)"""
        context = super().get_context_data(**kwargs)

        # get selection from session and add to context
        # get simple selection from session
        simple_selection = self.request.session.get('selection', False)
        
        # create full selection and import simple selection (if it exists)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)

        # on the first page of a workflow, clear the selection
        if self.step == 1:
            selection.clear('reference')
            selection.clear('targets')
            selection.clear('segments')
            simple_selection = selection.exporter()
            self.request.session['selection'] = simple_selection

        context['selection'] = {}
        for selection_box, include in self.selection_boxes.items():
            if include:
                context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]

        if self.filters:
            context['selection']['species'] = selection.species
            context['selection']['annotation'] = selection.annotation

        # get attributes of this class and add them to the context
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]
        return context

class AbsReferenceSelection(AbsTargetSelection):
    type_of_selection = 'reference'
    step = 1
    number_of_steps = 3
    title = 'SELECT A REFERENCE TARGET'
    description = 'Select a reference target by searching or browsing in the right column.\n\nThe reference will be compared to the targets you select later in the workflow.\n\nOnce you have selected your reference target, you will be redirected to the next step.'
    selection_boxes = {}
    psets = [] # protein sets not applicable for this selection

class AbsBrowseSelection(AbsTargetSelection):
    type_of_selection = 'browse'
    step = 1
    number_of_steps = 1
    title = 'SELECT A TARGET OR FAMILY'
    description = 'Select a target or family by searching or browsing in the right column.'
    selection_boxes = {}
    psets = [] # protein sets not applicable for this selection

class AbsSegmentSelection(TemplateView):
    """An abstract class for the segment selection page used in many apps. To use it in another app, create a class 
    based view for that app that extends this class"""
    template_name = 'common/segmentselection.html'

    step = 2
    number_of_steps = 2
    title = 'SELECT SEQUENCE SEGMENTS'
    description = 'Select sequence segments in the middle column. You can expand helices and select individual residues by clicking on the down arrows next to each helix.\n\nSelected segments will appear in the right column, where you can edit the list.\n\nOnce you have selected all your segments, click the green button.'
    docs = '/docs/protein'
    segment_list = True
    buttons = {
        'continue': {
            'label': 'Show alignment',
            'url': '/alignment/render',
            'color': 'success',
        },
    }
    # OrderedDict to preserve the order of the boxes
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', True),
    ])

    ss = ProteinSegment.objects.filter(partial=False).prefetch_related('generic_numbers')
    ss_cats = ProteinSegment.objects.values_list('category').order_by('category').distinct('category')
    action = 'expand'

    def get_context_data(self, **kwargs):
        """get context from parent class (really only relevant for child classes of this class, as TemplateView does
        not have any context variables)"""
        context = super().get_context_data(**kwargs)

        # get selection from session and add to context
        # get simple selection from session
        simple_selection = self.request.session.get('selection', False)

        # create full selection and import simple selection (if it exists)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)

        context['selection'] = {}
        for selection_box, include in self.selection_boxes.items():
            if include:
                context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]

        # get attributes of this class and add them to the context
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]
        return context


def AddToSelection(request):
    """Receives a selection request, adds the selected item to session, and returns the updated selection"""
    selection_type = request.GET['selection_type']
    selection_subtype = request.GET['selection_subtype']
    selection_id = request.GET['selection_id']
    
    o = []
    if selection_type == 'reference' or selection_type == 'targets':
        if selection_subtype == 'protein':
            o.append(Protein.objects.get(pk=selection_id))
        elif selection_subtype == 'protein_set':
            selection_subtype = 'protein'
            pset = ProteinSet.objects.get(pk=selection_id)
            for protein in pset.proteins.all():
                o.append(protein)
        elif selection_subtype == 'family':
            o.append(ProteinFamily.objects.get(pk=selection_id))
        elif selection_subtype == 'set':
            o.append(ProteinSet.objects.get(pk=selection_id))
        elif selection_subtype == 'structure':
            o.append(Protein.objects.get(entry_name=selection_id.lower()))
    elif selection_type == 'segments':
        if selection_subtype == 'residue':
            o.append(ResidueGenericNumber.objects.get(pk=selection_id))
        else:
            o.append(ProteinSegment.objects.get(pk=selection_id))

    # get simple selection from session
    simple_selection = request.session.get('selection', False)
    
    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    for obj in o:
        # add the selected item to the selection
        selection_object = SelectionItem(selection_subtype, obj)
        selection.add(selection_type, selection_subtype, selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection
    
    return render(request, 'common/selection_lists.html', selection.dict(selection_type))

def RemoveFromSelection(request):
    """Removes one selected item from the session"""
    selection_type = request.GET['selection_type']
    selection_subtype = request.GET['selection_subtype']
    selection_id = request.GET['selection_id']
    
    # get simple selection from session
    simple_selection = request.session.get('selection', False)
    
    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # remove the selected item to the selection
    selection.remove(selection_type, selection_subtype, selection_id)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection
    
    return render(request, 'common/selection_lists.html', selection.dict(selection_type))

def ClearSelection(request):
    """Clears all selected items of the selected type from the session"""
    selection_type = request.GET['selection_type']
    
    # get simple selection from session
    simple_selection = request.session.get('selection', False)
    
    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # remove the selected item to the selection
    selection.clear(selection_type)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection
    
    return render(request, 'common/selection_lists.html', selection.dict(selection_type))

def SelectFullSequence(request):
    """Adds all segments to the selection"""
    selection_type = request.GET['selection_type']

    # get simple selection from session
    simple_selection = request.session.get('selection', False)
    
    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # get all segments
    segments = ProteinSegment.objects.filter(partial=False)
    for segment in segments:
        selection_object = SelectionItem(segment.category, segment)
        # add the selected item to the selection
        selection.add(selection_type, segment.category, selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection
    
    return render(request, 'common/selection_lists.html', selection.dict(selection_type))

def SelectAlignableSegments(request):
    """Adds all alignable segments to the selection"""
    selection_type = request.GET['selection_type']

    # get simple selection from session
    simple_selection = request.session.get('selection', False)
    
    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # get all segments
    segments = ProteinSegment.objects.filter(partial=False, slug__startswith='TM')
    for segment in segments:
        selection_object = SelectionItem(segment.category, segment)
        # add the selected item to the selection
        selection.add(selection_type, segment.category, selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection
    
    return render(request, 'common/selection_lists.html', selection.dict(selection_type))

def ToggleFamilyTreeNode(request):
    """Opens/closes a node in the family selection tree"""
    action = request.GET['action']
    type_of_selection = request.GET['type_of_selection']
    node_id = request.GET['node_id']
    parent_tree_indent_level = int(request.GET['tree_indent_level'])
    tree_indent_level = []
    for i in range(parent_tree_indent_level+1):
        tree_indent_level.append(0)
    parent_tree_indent_level = tree_indent_level[:]
    del parent_tree_indent_level[-1]

    # session
    simple_selection = request.session.get('selection')
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    ppf = ProteinFamily.objects.get(pk=node_id)
    if action == 'expand':
        pfs = ProteinFamily.objects.filter(parent=node_id)

        # species filter
        species_list = []
        for species in selection.species:
            species_list.append(species.item)

        # annotation filter
        protein_source_list = []
        for protein_source in selection.annotation:
            protein_source_list.append(protein_source.item)

        ps = Protein.objects.order_by('id').filter(family=ppf,
            species__in=(species_list),
            source__in=(protein_source_list))
        action = 'collapse'
    else:
        pfs = ps = {}
        action = 'expand'
    
    return render(request, 'common/selection_tree.html', {
        'action': action,
        'type_of_selection': type_of_selection,
        'ppf': ppf,
        'pfs': pfs,
        'ps': ps,
        'parent_tree_indent_level': parent_tree_indent_level,
        'tree_indent_level': tree_indent_level,
    })

def SelectionAnnotation(request):
    """Updates the selected level of annotation"""
    protein_source = request.GET['protein_source']
    
    if protein_source == 'All':
        pss = ProteinSource.objects.all()
    else:
        pss = ProteinSource.objects.filter(name=protein_source)

    # get simple selection from session
    simple_selection = request.session.get('selection', False)
    
    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # reset the annotation selection
    selection.clear('annotation')

    # add the selected items to the selection
    for ps in pss:
        selection_object = SelectionItem('annotation', ps)
        selection.add('annotation', 'annotation', selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection
    
    return render(request, 'common/selection_filters_annotation.html', selection.dict('annotation'))

def SelectionSpeciesPredefined(request):
    """Updates the selected species to predefined sets (Human and all)"""
    species = request.GET['species']

    # get simple selection from session
    simple_selection = request.session.get('selection', False)
    
    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)
    
    all_sps = Species.objects.all()
    sps = False
    if species == 'All':
        sps = all_sps
    elif species:
        sps = Species.objects.filter(common_name=species)
    elif not selection.species:
        sps = Species.objects.filter(common_name='Human') # if nothing is selected, select human

    if sps:
        # reset the species selection
        selection.clear('species')

        # add the selected items to the selection
        for sp in sps:
            selection_object = SelectionItem('species', sp)
            selection.add('species', 'species', selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    # add all species objects to context (for comparison to selected species)
    context = selection.dict('species')
    context['sps'] = all_sps
    
    return render(request, 'common/selection_filters_species.html', context)

def SelectionSpeciesToggle(request):
    """Updates the selected species arbitrary selections"""
    species_id = request.GET['species_id']

    all_sps = Species.objects.all()
    sps = Species.objects.filter(pk=species_id)

    # get simple selection from session
    simple_selection = request.session.get('selection', False)
    
    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # add the selected items to the selection
    for sp in sps:
        exists = selection.remove('species', 'species', species_id)
        if not exists:
            selection_object = SelectionItem('species', sp)
            selection.add('species', 'species', selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    # add all species objects to context (for comparison to selected species)
    context = selection.dict('species')
    context['sps'] = Species.objects.all()
    
    return render(request, 'common/selection_filters_species_selector.html', context)

def ExpandSegment(request):
    """Expands a segment to show it's generic numbers"""
    segment_id = request.GET['segment_id']
    
    # fetch the generic numbers
    context = {}
    context['generic_numbers'] = ResidueGenericNumber.objects.filter(protein_segment__id=segment_id,
        scheme__slug=settings.DEFAULT_NUMBERING_SCHEME)
    
    return render(request, 'common/segment_generic_numbers.html', context)

def SelectionSchemesPredefined(request):
    """Updates the selected numbering_schemes to predefined sets (GPCRdb and All)"""
    numbering_schemes = request.GET['numbering_schemes']

    # get simple selection from session
    simple_selection = request.session.get('selection', False)
    
    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)
    
    all_gns = ResidueNumberingScheme.objects.all()
    gns = False
    if numbering_schemes == 'All':
        gns = all_gns
    elif numbering_schemes:
        gns = ResidueNumberingScheme.objects.filter(slug=numbering_schemes)
    elif not selection.numbering_schemes:
        gns = ResidueNumberingScheme.objects.filter(slug='gpcrdb') # if nothing is selected, select gpcrdb

    if gns:
        # reset the species selection
        selection.clear('numbering_schemes')

        # add the selected items to the selection
        for gn in gns:
            selection_object = SelectionItem('numbering_schemes', gn)
            selection.add('numbering_schemes', 'numbering_schemes', selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    # add all species objects to context (for comparison to selected species)
    context = selection.dict('numbering_schemes')
    context['gns'] = all_gns
    
    return render(request, 'common/selection_filters_numbering_schemes.html', context)

def SelectionSchemesToggle(request):
    """Updates the selected numbering schemes arbitrary selections"""
    numbering_scheme_id = request.GET['numbering_scheme_id']
    print(numbering_scheme_id)
    all_gns = ResidueNumberingScheme.objects.all()
    gns = ResidueNumberingScheme.objects.filter(pk=numbering_scheme_id)

    # get simple selection from session
    simple_selection = request.session.get('selection', False)
    
    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # add the selected items to the selection
    for gn in gns:
        exists = selection.remove('numbering_schemes', 'numbering_schemes', numbering_scheme_id)
        if not exists:
            selection_object = SelectionItem('numbering_schemes', gn)
            selection.add('numbering_schemes', 'numbering_schemes', selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    # add all species objects to context (for comparison to selected species)
    context = selection.dict('numbering_schemes')
    context['gns'] = ResidueNumberingScheme.objects.all()
    
    return render(request, 'common/selection_filters_numbering_schemes.html', context)