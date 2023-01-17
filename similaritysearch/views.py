from django.shortcuts import render, redirect
from django.conf import settings

from common.views import AbsReferenceSelection
from common.views import AbsSegmentSelection
#from common.views import AbsTargetSelection
from common.views import AbsTargetSelectionTable
from protwis.context_processors import site_title

# from common.alignment_SITE_NAME import Alignment
Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')

from collections import OrderedDict


class ReferenceSelection(AbsReferenceSelection):
    step = 1
    number_of_steps = 3
    target_input = False
    docs = 'sequences.html#similarity-search-gpcrdb'
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
    docs = 'sequences.html#similarity-search-gpcrdb'
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


class TargetSelection(AbsTargetSelectionTable):
    step = 3
    number_of_steps = 3
    docs = "sequences.html#similarity-search-gpcrdb"
    title = "SELECT RECEPTORS"
    description = "Select receptors in the table (below) or browse the classification tree (right). You can select entire" \
        + " families or individual receptors.\n\nOnce you have selected all your receptors, click the green button."
    selection_boxes = OrderedDict([
        ("reference", True),
        ("segments", True),
        ("targets", True),
    ])
    buttons = {
        "continue": {
            "label": "Next",
            "onclick": "submitSelection('/similaritysearch/render');",
            "color": "success",
        },
    }


# class TargetSelection(AbsTargetSelection):
#     step = 3
#     number_of_steps = 3
#     docs = 'sequences.html#similarity-search-gpcrdb'
#     selection_boxes = OrderedDict([
#         ('reference', True),
#         ('segments', True),
#         ('targets', True),
#     ])
#     buttons = {
#         'continue': {
#             'label': 'Show similarity',
#             'url': '/similaritysearch/render',
#             'color': 'success',
#         },
#     }

def render_alignment(request):
    # get the user selection from session
    simple_selection = request.session.get('selection', False)

    if simple_selection == False or not simple_selection.targets or not simple_selection.reference:
        return redirect("/similaritysearch/referenceselection")

    # create an alignment object
    a = Alignment()

    # load data from selection into the alignment
    a.load_reference_protein_from_selection(simple_selection)
    # only show wildtype protein entries if selection type is family
    if len([t for t in simple_selection.targets if t.type=='family'])>0:
        a.load_proteins_from_selection(simple_selection, only_wildtype=True)
    else:
        a.load_proteins_from_selection(simple_selection)
    a.load_segments_from_selection(simple_selection)

    # build the alignment data matrix
    check = a.build_alignment()
    if check == 'Too large':
        return render(request, 'similaritysearch/error.html', {'proteins': len(a.proteins), 'residues':a.number_of_residues_total})

    # calculate consensus sequence + amino acid and feature frequency
    a.calculate_statistics()

    # calculate identity and similarity of each row compared to the reference
    a.calculate_similarity()

    num_of_sequences = len(a.proteins)
    num_residue_columns = len(a.positions) + len(a.segments)

    return render(request, 'similaritysearch/alignment.html', {'a': a, 'num_of_sequences': num_of_sequences,
        'num_residue_columns': num_residue_columns})

def render_fasta_alignment(request):
    # get the user selection from session
    simple_selection = request.session.get('selection', False)

    # create an alignment object
    a = Alignment()
    a.show_padding = False

    # load data from selection into the alignment
    a.load_reference_protein_from_selection(simple_selection)
    a.load_proteins_from_selection(simple_selection)
    a.load_segments_from_selection(simple_selection)

    # build the alignment data matrix
    a.build_alignment()

    # calculate identity and similarity of each row compared to the reference
    a.calculate_similarity()

    num_of_sequences = len(a.proteins)
    num_residue_columns = len(a.positions) + len(a.segments)

    response = render(request, 'alignment/alignment_fasta.html', {'a': a, 'num_of_sequences': num_of_sequences,
        'num_residue_columns': num_residue_columns}, content_type='text/fasta')
    response['Content-Disposition'] = "attachment; filename=" + site_title(request)["site_title"] + "_alignment.fasta"
    return response

def render_csv_alignment(request):
    # get the user selection from session
    simple_selection = request.session.get('selection', False)

    # create an alignment object
    a = Alignment()
    a.show_padding = False

    # load data from selection into the alignment
    a.load_reference_protein_from_selection(simple_selection)
    a.load_proteins_from_selection(simple_selection)
    a.load_segments_from_selection(simple_selection)

    # build the alignment data matrix
    a.build_alignment()

    # calculate consensus sequence + amino acid and feature frequency
    a.calculate_statistics()

    # calculate identity and similarity of each row compared to the reference
    a.calculate_similarity()

    num_of_sequences = len(a.proteins)
    num_residue_columns = len(a.positions) + len(a.segments)

    response = render(request, 'alignment/alignment_csv.html', {'a': a, 'num_of_sequences': num_of_sequences,
        'num_residue_columns': num_residue_columns}, content_type='text/fasta')
    response['Content-Disposition'] = "attachment; filename=" + site_title(request)["site_title"] + "_alignment.csv"
    return response
