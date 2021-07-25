from django.shortcuts import render
from django.conf import settings

from common.views import AbsSegmentSelection
#from common.views import AbsTargetSelection
from common.views import AbsTargetSelectionTable
# from common.alignment_SITE_NAME import Alignment
from protwis.context_processors import site_title

Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')

from collections import OrderedDict


class TargetSelection(AbsTargetSelectionTable):
    step = 1
    number_of_steps = 2
    title = "SELECT RECEPTORS"
    description = "Select receptors in the table (below) or browse the classification tree (right). You can select entire" \
        + " families or individual receptors.\n\nOnce you have selected all your receptors, click the green button."
    docs = "sequences.html#similarity-matrix"
    selection_boxes = OrderedDict([
        ("reference", False),
        ("targets", True),
        ("segments", False),
    ])
    buttons = {
        "continue": {
            "label": "Next",
            "onclick": "submitSelection('/similaritymatrix/segmentselection');",
            "color": "success",
        },
    }

# class TargetSelection(AbsTargetSelection):
#     step = 1
#     number_of_steps = 2
#     docs = 'sequences.html#similarity-matrix'
#     selection_boxes = OrderedDict([
#         ('reference', False),
#         ('targets', True),
#         ('segments', False),
#     ])
#     buttons = {
#         'continue': {
#             'label': 'Continue to next step',
#             'url': '/similaritymatrix/segmentselection',
#             'color': 'success',
#         },
#     }


class SegmentSelection(AbsSegmentSelection):
    step = 2
    number_of_steps = 2
    docs = 'sequences.html#similarity-matrix'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', False),
        ('segments', True),
    ])
    buttons = {
        'continue': {
            'label': 'Show matrix',
            'url': '/similaritymatrix/render',
            'color': 'success',
        },
    }


def render_matrix(request):
    # get the user selection from session
    simple_selection = request.session.get('selection', False)

    # create an alignment object
    a = Alignment()

    # load data from selection into the alignment
    a.load_proteins_from_selection(simple_selection)
    a.load_segments_from_selection(simple_selection)

    # build the alignment data matrix
    a.build_alignment()

    # NOTE: NOT necessary for similarity matrix
    # calculate consensus sequence + amino acid and feature frequency
    # a.calculate_statistics()

    # calculate identity and similarity of each row compared to the reference
    a.calculate_similarity_matrix()

    return render(request, 'similaritymatrix/matrix.html', {'p': a.proteins, 'm': a.similarity_matrix})

def render_csv_matrix(request):
    # get the user selection from session
    simple_selection = request.session.get('selection', False)

    # create an alignment object
    a = Alignment()
    a.show_padding = False

    # load data from selection into the alignment
    a.load_proteins_from_selection(simple_selection)
    a.load_segments_from_selection(simple_selection)

    # build the alignment data matrix
    a.build_alignment()

    # calculate consensus sequence + amino acid and feature frequency
    # NOTE: NOT necessary for similarity matrix
    # a.calculate_statistics()

    # calculate identity and similarity of each row compared to the reference
    a.calculate_similarity_matrix()

    response = render(request, 'similaritymatrix/matrix_csv.html', {'p': a.proteins, 'm': a.similarity_matrix})
    response['Content-Disposition'] = "attachment; filename=" + site_title(request)["site_title"] + "_similaritymatrix.csv"
    return response
