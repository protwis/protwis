from django.shortcuts import render
from django.conf import settings
from django.http import HttpResponse
from django.views.decorators.csrf import csrf_exempt

from common import definitions
from common.selection import SimpleSelection, Selection, SelectionItem

from common.views import AbsReferenceSelection
from common.views import AbsSegmentSelection
from common.views import AbsTargetSelection
# from common.alignment_SITE_NAME import Alignment
Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')
from protein.models import ProteinSegment
from residue.models import ResidueGenericNumberEquivalent

import os
from collections import OrderedDict
from io import BytesIO
import xlsxwriter, xlrd


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

    rsets = False
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
    rsets = False

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


def site_download(request):

    simple_selection = request.session.get('selection', False)
    outstream = BytesIO()
    wb = xlsxwriter.Workbook(outstream, {'in_memory': True})
    worksheet = wb.add_worksheet()
    row_count = 0

    for position in simple_selection.segments:
        if position.type == 'site_residue' and position.properties['site_residue_group']:
            group_id = position.properties['site_residue_group']
            #Values saved into the file come in the following order:
            # id    min_match   generic_number  numbering_scheme    feature     amino_acids
            #List of amino acids is not necessary, it is just for easy visual identification
            worksheet.write_row(row_count, 0, [group_id, simple_selection.site_residue_groups[group_id -1][0], position.item.label, position.item.scheme.slug, position.properties['feature'], ','.join(definitions.AMINO_ACID_GROUPS_OLD[position.properties['feature']])])
            row_count += 1
    wb.close()
    outstream.seek(0)
    response = HttpResponse(outstream.read(), content_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    response['Content-Disposition'] = "attachment; filename=site_definition.xlsx"

    return response


def site_upload(request):

    # get simple selection from session
    simple_selection = request.session.get('selection', False)

    # create full selection and import simple selection (if it exists)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)

    selection_type = 'segments'
    selection_subtype = 'site_residue'

    if request.FILES == {}:
        return render(request, 'common/selection_lists_sitesearch.html', '')

    #Overwriting the existing selection
    selection.clear(selection_type)

    workbook = xlrd.open_workbook(file_contents=request.FILES['xml_file'].read())
    worksheets = workbook.sheet_names()
    for worksheet_name in worksheets:
        worksheet = workbook.sheet_by_name(worksheet_name)
        for row in worksheet.get_rows():

    #for row in ws.rows:
            if len(row) < 5:
                continue
            group_id = int(row[0].value)
            min_match = int(row[1].value)
            try:
                position = ResidueGenericNumberEquivalent.objects.get(label=row[2].value, scheme__slug=row[3].value)
            except e:
                print(e)
                continue
            feature = row[4].value

            # update the selected group
            if not selection.active_site_residue_group:
                selection.active_site_residue_group = group_id
            if not selection.site_residue_groups:
                selection.site_residue_groups = [[]]
            if len(selection.site_residue_groups) < group_id:
                for x in group_id - len(selection.site_residue_groups):
                    selection.site_residue_groups[x].append([])
#            selection.site_residue_groups[group_id - 1].append(1)
            properties = {}
            properties['feature'] = feature
            properties['site_residue_group'] = selection.active_site_residue_group
            properties['amino_acids'] = ','.join(definitions.AMINO_ACID_GROUPS_OLD[feature]) if feature != 'custom' else row[5].value
            selection_object = SelectionItem(selection_subtype, position, properties)
            selection.add(selection_type, selection_subtype, selection_object)
            selection.site_residue_groups[group_id - 1][0] = min_match
    # export simple selection that can be serialized
    simple_selection = selection.exporter()
    # add simple selection to session
    request.session['selection'] = simple_selection

    context = selection.dict(selection_type)

    # amino acid groups
    amino_acid_groups = {
        'amino_acid_groups_old': definitions.AMINO_ACID_GROUPS_OLD,
        'amino_acid_group_names_old': definitions.AMINO_ACID_GROUP_NAMES_OLD }

    return render(request, 'common/selection_lists_sitesearch.html', context)


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
