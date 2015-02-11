from common.views import AbsReferenceSelection
from common.views import AbsSegmentSelection
from common.views import AbsTargetSelection

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
            'url': '/similaritysearch/results',
            'color': 'success',
        },
    }