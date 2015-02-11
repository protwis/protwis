from common.views import AbsTargetSelection
from common.views import AbsSegmentSelection

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
            'url': '/alignment/results',
            'color': 'success',
        },
    }