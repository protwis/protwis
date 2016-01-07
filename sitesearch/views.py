from django.shortcuts import render
from django.conf import settings

from common.views import AbsReferenceSelection
from common.views import AbsSegmentSelection
from common.views import AbsTargetSelection
# from common.alignment_SITE_NAME import Alignment
Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')
from protein.models import ProteinSegment

from collections import OrderedDict


class TargetSelection(AbsTargetSelection):
    step = 1
    number_of_steps = 2
    docs = 'sites.html#site-search-manual'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/sitesearch/segmentselection',
            'color': 'success',
        },
    }


class SegmentSelection(AbsSegmentSelection):
    step = 2
    number_of_steps = 2
    docs = 'sites.html#site-search-manual'
    position_type = 'site_residue'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', True),
    ])
    buttons = {
        'continue': {
            'label': 'Show site',
            'url': '/sitesearch/render',
            'color': 'success',
        },
    }

    title = 'SELECT SITE RESIDUES'
    description = 'Select site residues. Click the down arrow next to each helix to expand the available residues' \
        + ' within that helix.\n\nSelected residues will appear in the right column, where you can edit the' \
        + ' list.\n\nIn the right column, select a chemical feature for each residue. When a feature has been' \
        + ' selected, a list of amino acids that match the feature will appear to the right of the residue.\n\n' \
        + ' The selected residues can be organised into separate interactions. An interaction can contain one or' \
        + ' more residues. To add an interaction, click the "Add interaction" button. Selected residues will be' \
        + ' added to the currently active interaction (shown in bold text). To change the active interaction, click' \
        + ' on the name of the interaction. Within an interaction, the number of residues required to match can be' \
        + ' specified in the "Min. match" selection box.\n\n Once you have selected your site residues, click the' \
        + ' green button.'

    ss = ProteinSegment.objects.filter(slug__in=settings.REFERENCE_POSITIONS, partial=False).prefetch_related(
        'generic_numbers')
    ss_cats = ss.values_list('category').order_by('category').distinct('category')


class TargetSelectionPdb(TargetSelection):
    step = 1
    number_of_steps = 3
    docs = 'sites.html#site-search-from-pdb-complex'
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/sitesearch/structureupload',
            'color': 'success',
        },
    }


class StructureUpload(AbsSegmentSelection):
    step = 2
    number_of_steps = 3
    docs = 'sites.html#site-search-from-pdb-complex'
    segment_list = False
    structure_upload = True

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'onclick': 'document.getElementById("form").submit()',
            'color': 'success',
        },
    }

    title = 'SELECT OR UPLOAD A STRUCTURE'
    description = 'Enter a PDB code or upload your own PDB file in the middle column, and click the green button.' \
        + '\n\nProtein-ligand interactions from the complex will be analyzed and shown on the next page.'


class SegmentSelectionPdb(SegmentSelection):
    step = 3
    number_of_steps = 3
    docs = 'sites.html#site-search-from-pdb-complex'

    ss = ProteinSegment.objects.filter(slug__in=settings.REFERENCE_POSITIONS, partial=False).prefetch_related(
        'generic_numbers')
    ss_cats = ss.values_list('category').order_by('category').distinct('category')


def render_alignment(request):
    # get the user selection from session
    simple_selection = request.session.get('selection', False)
    
    # create an alignment object
    a = Alignment()

    # group residue selection (for site search)
    a.use_residue_groups = True

    # load data from selection into the alignment
    a.load_reference_protein_from_selection(simple_selection)
    a.load_proteins_from_selection(simple_selection)
    a.load_segments_from_selection(simple_selection)

    # build the alignment data matrix
    a.build_alignment()

    # evaluate sites
    a.evaluate_sites(request)

    num_of_sequences = len(a.proteins)
    num_of_non_matching_sequences = len(a.non_matching_proteins)
    num_residue_columns = len(a.positions) + len(a.segments)

    context = {'a': a, 'num_of_sequences': num_of_sequences,
        'num_of_non_matching_sequences': num_of_non_matching_sequences, 'num_residue_columns': num_residue_columns}

    return render(request, 'sitesearch/alignment.html', context)

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
    num_of_non_matching_sequences = len(a.non_matching_proteins)
    num_residue_columns = len(a.positions) + len(a.segments)
    
    context = {'a': a, 'num_of_sequences': num_of_sequences,
        'num_of_non_matching_sequences': num_of_non_matching_sequences, 'num_residue_columns': num_residue_columns}

    response = render(request, 'alignment/alignment_fasta.html', context, content_type='text/fasta')
    response['Content-Disposition'] = "attachment; filename=" + settings.SITE_TITLE + "_alignment.fasta"
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
    response['Content-Disposition'] = "attachment; filename=" + settings.SITE_TITLE + "_alignment.csv"
    return response
