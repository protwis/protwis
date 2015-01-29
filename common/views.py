from django.http import HttpResponse
from django.shortcuts import render
from django.views.generic import TemplateView

from common.classes import SimpleSelection
from common.classes import Selection
from common.classes import SelectionItem
from protein.models import Protein
from protein.models import ProteinFamily

import inspect
from collections import OrderedDict


class AbsTargetSelection(TemplateView):
    """An abstract class for the target selection page used in many apps. To use it in another app, create a class 
    based view for that app that extends this class"""
    template_name = 'common/targetselection.html'

    step = 1
    number_of_steps = 2
    title = 'Select targets'
    description = 'Select receptors by searching or browsing in the middle column. You can select entire receptor families or individual receptors.\n\nSelected receptors will appear in the right column, where you can edit the list.\n\nOnce you have selected all your receptors, click the green button.'
    docs = '/docs/protein'
    filters = True
    search = True
    family_tree = True
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/protein/segmentselection',
            'color': 'success',
        },
    }
    # OrderedDict to preserve the order of the boxes
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])

    ppf = ProteinFamily.objects.get(slug='000')
    pfs = ProteinFamily.objects.order_by('id').filter(parent=ppf.id) # FIXME move order_by to model
    ps = Protein.objects.filter(family=ppf)
    tree_indent_level = []
    action = 'expand'
    # remove the parent family (for all other families than the root of the tree, the parent should be shown)
    del ppf


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
    
    if selection_subtype == 'protein':
        p = Protein.objects.get(pk=selection_id)
        selection_object = SelectionItem('protein', p)
    elif selection_subtype == 'family':
        pf = ProteinFamily.objects.get(pk=selection_id)
        selection_object = SelectionItem('family', pf)
    elif selection_subtype == 'set':
        ps = ProteinSet.objects.get(pk=selection_id)
        selection_object = SelectionItem('set', ps)

    # get simple selection from session
    simple_selection = request.session.get('selection', False)
    
    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    # add the selected item to the selection
    sel_type = getattr(selection, selection_type)
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
    sel_type = getattr(selection, selection_type)
    selection.remove(selection_type, selection_subtype, selection_id)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection
    
    return render(request, 'common/selection_lists.html', selection.dict(selection_type))

def ClearSelection(request):
    """Clears all selected items of the selected type from the session"""
    selection_type = request.GET['selection_type']
    
    # create empty selections
    selection = Selection()
    simple_selection = SimpleSelection()

    # add simple selection to session
    request.session['selection'] = simple_selection
    
    return render(request, 'common/selection_lists.html', selection.dict(selection_type))

def ToggleFamilyTreeNode(request):
    """WRITEME"""
    action = request.GET['action']
    node_id = request.GET['node_id']
    parent_tree_indent_level = int(request.GET['tree_indent_level'])
    tree_indent_level = []
    for i in range(parent_tree_indent_level+1):
        tree_indent_level.append(0)
    # if action == 'collapse':
    #     del tree_indent_level[-1]
    parent_tree_indent_level = tree_indent_level[:]
    del parent_tree_indent_level[-1]

    ppf = ProteinFamily.objects.get(pk=node_id)
    if action == 'expand':
        pfs = ProteinFamily.objects.order_by('id').filter(parent=node_id) # FIXME move order_by to model
        ps = Protein.objects.order_by('id').filter(family=ppf)
        action = 'collapse'
    else:
        pfs = ps = {}
        action = 'expand'
    
    return render(request, 'common/selection_tree.html', {
        'action': action,
        'ppf': ppf,
        'pfs': pfs,
        'ps': ps,
        'parent_tree_indent_level': parent_tree_indent_level,
        'tree_indent_level': tree_indent_level,
    })