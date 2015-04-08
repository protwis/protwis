from django.shortcuts import render
from django.conf import settings

from common.views import AbsReferenceSelection
from common.views import AbsSegmentSelection
from common.views import AbsTargetSelection
# from common.alignment_SITE_NAME import Alignment
Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')

from collections import OrderedDict


class ReferenceSelection(AbsReferenceSelection):
    step = 1
    number_of_steps = 3
    docs = '/docs/similaritysearch'
    selection_boxes = OrderedDict([
        ('reference', True),
        ('segments', False),
        ('targets', False),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/similaritysearch/segmentselection',
            'color': 'success',
        },
    }


class SegmentSelection(AbsSegmentSelection):
    step = 2
    number_of_steps = 3
    docs = '/docs/similaritysearch'
    selection_boxes = OrderedDict([
        ('reference', True),
        ('segments', True),
        ('targets', False),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/similaritysearch/targetselection',
            'color': 'success',
        },
    }


class TargetSelection(AbsTargetSelection):
    step = 3
    number_of_steps = 3
    docs = '/docs/similaritysearch'
    selection_boxes = OrderedDict([
        ('reference', True),
        ('segments', True),
        ('targets', True),
    ])
    buttons = {
        'continue': {
            'label': 'Show similarity',
            'url': '/similaritysearch/render',
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

    # calculate identity and similarity of each row compared to the reference
    a.calculate_similarity()

    return render(request, 'similaritysearch/alignment.html', {'a': a})