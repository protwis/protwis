from django.shortcuts import render

from common.classes import SimpleSelection
from common.classes import Selection
from common.classes import SelectionItem
from protein.models import Protein
from protein.models import ProteinFamily


def AddToSelection(request):
    """Receives a selection request, adds the selected item to session, and returns the current selection"""
    selection_type = request.GET['selection_type']
    selection_subtype = request.GET['selection_subtype']
    selection_id = request.GET['selection_id']
    
    if selection_subtype == 'protein':
        p = Protein.objects.get(entry_name=selection_id)
        selection_object = SelectionItem('protein', p)
    elif selection_subtype == 'family':
        pf = ProteinFamily.objects.get(slug=selection_id)
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

    # add the selected item to the selection (if it's not already selected)
    sel_type = getattr(selection, selection_type)
    if not selection_object in sel_type:
        selection.add(selection_type, selection_subtype, selection_object)

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection
    
    return render(request, 'common/selected_data.html', selection.render())

def ClearSelection(request):
    # create a blank selection
    selection = Selection()

    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection
    
    return render(request, 'common/selected_data.html', selection.render())