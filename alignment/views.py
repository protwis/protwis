from django.shortcuts import render
from django.conf import settings

from common.views import AbsTargetSelection
from common.views import AbsSegmentSelection
# from common.alignment_SITE_NAME import Alignment
Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')

from collections import OrderedDict


class TargetSelection(AbsTargetSelection):
    step = 1
    number_of_steps = 2
    docs = '/docs/alignment'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/alignment/segmentselection',
            'color': 'success',
        },
    }


class SegmentSelection(AbsSegmentSelection):
    step = 2
    number_of_steps = 2
    docs = '/docs/alignment'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', True),
    ])
    buttons = {
        'continue': {
            'label': 'Show alignment',
            'url': '/alignment/render',
            'color': 'success',
        },
    }


def render_alignment(request):
    # get the user selection from session
    simple_selection = request.session.get('selection', False)
    
    # create an alignment object
    a = Alignment()

    # load data from selection into the alignment
    a.load_proteins_from_selection(simple_selection)
    a.load_positions_from_selection(simple_selection)

    # build the alignment data matrix
    a.build_alignment_matrix()

    return render(request, 'alignment/alignment.html', {'a': a})