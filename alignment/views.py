from django.shortcuts import render
from django.conf import settings
from django.http import HttpResponse
from django.views.generic import TemplateView
from django.db.models import Case, When
from django.core.cache import cache
from django.core.cache import caches
try:
    cache_alignment = caches['alignments']
except:
    cache_alignment = cache

from alignment import functions
from common import definitions
from common.selection import Selection
from common.views import AbsTargetSelection
from common.views import AbsSegmentSelection
from common.views import AbsMiscSelection
from structure.functions import BlastSearch

# from common.alignment_SITE_NAME import Alignment
Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')
from protein.models import Protein, ProteinSegment, ProteinFamily, ProteinSet
from residue.models import ResidueNumberingScheme, ResiduePositionSet

from collections import OrderedDict
from copy import deepcopy
import hashlib
import inspect
from io import BytesIO
import itertools
import json
import numpy as np
import os
import xlsxwriter
import xlrd


class TargetSelection(AbsTargetSelection):
    step = 1
    number_of_steps = 2
    docs = 'sequences.html#structure-based-alignments'
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

class PosTargetSelection(AbsTargetSelection):
    step = 1
    number_of_steps = 4
    docs = 'sequences.html#structure-based-alignments'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/alignment/negativegroupselection',
            'color': 'success',
        },
    }

class NegTargetSelection(AbsTargetSelection):

    step = 2
    number_of_steps = 4
    docs = 'sequences.html#structure-based-alignments'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/alignment/segmentselectionsignature',
            'color': 'success',
        },
    }

    def get_context_data(self, **kwargs):
        #A bit ugly solution to having two target sets without modifying half of common.selection
        context = super(NegTargetSelection, self).get_context_data(**kwargs)

        self.request.session['targets_pos'] = deepcopy(self.request.session.get('selection', False))
        del self.request.session['selection']

        return context

class TargetSelectionGprotein(AbsTargetSelection):
    step = 1
    number_of_steps = 2
    psets = False
    filters = True
    filter_gprotein = True

    docs = 'sequences.html#structure-based-alignments'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/alignment/segmentselectiongprot',
            'color': 'success',
        },
    }
    try:
        ppf = ProteinFamily.objects.get(slug="100_000")
        pfs = ProteinFamily.objects.filter(parent=ppf.id)
        ps = Protein.objects.filter(family=ppf)

        tree_indent_level = []
        action = 'expand'
        # remove the parent family (for all other families than the root of the tree, the parent should be shown)
        del ppf
    except Exception as e:
        pass

class TargetSelectionArrestin(AbsTargetSelection):
    step = 1
    number_of_steps = 2
    psets = False
    filters = True
    filter_gprotein = True

    docs = 'sequences.html#structure-based-alignments'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/alignment/segmentselectionsignature',
            'color': 'success',
        },
    }
    try:
        ppf = ProteinFamily.objects.get(slug="200_000")
        pfs = ProteinFamily.objects.filter(parent=ppf.id)
        ps = Protein.objects.filter(family=ppf)

        tree_indent_level = []
        action = 'expand'
        # remove the parent family (for all other families than the root of the tree, the parent should be shown)
        del ppf
    except Exception as e:
        pass

class SegmentSelection(AbsSegmentSelection):
    step = 2
    number_of_steps = 2
    docs = 'sequences.html#structure-based-alignments'
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

class SegmentSelectionGprotein(AbsSegmentSelection):
    step = 2
    number_of_steps = 2
    docs = 'sequences.html#structure-based-alignments'
    description = 'Select sequence segments in the middle column for G proteins. You can expand every structural element and select individual' \
        + ' residues by clicking on the down arrows next to each helix, sheet or loop.\n\n You can select the full sequence or show all structured regions at the same time.\n\nSelected segments will appear in the' \
        + ' right column, where you can edit the list.\n\nOnce you have selected all your segments, click the green' \
        + ' button.'

    template_name = 'common/segmentselection.html'

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

    position_type = 'gprotein'
    rsets = ResiduePositionSet.objects.filter(name__in=['Gprotein Barcode', 'YM binding site']).prefetch_related('residue_position')

    ss = ProteinSegment.objects.filter(partial=False, proteinfamily='Gprotein').prefetch_related('generic_numbers')
    ss_cats = ss.values_list('category').order_by('category').distinct('category')

class SegmentSelectionArrestin(AbsSegmentSelection):
    step = 2
    number_of_steps = 2
    docs = 'sequences.html#structure-based-alignments'
    description = 'Select sequence segments in the middle column for beta and visual arrestins. You can expand every structural element and select individual' \
        + ' residues by clicking on the down arrows next to each helix, sheet or loop.\n\n You can select the full sequence or show all structured regions at the same time.\n\nSelected segments will appear in the' \
        + ' right column, where you can edit the list.\n\nOnce you have selected all your segments, click the green' \
        + ' button.'

    template_name = 'common/segmentselection.html'

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

    position_type = 'arrestin'

    ## Add some Arrestin specific positions
    rsets = ResiduePositionSet.objects.filter(name__in=['Arrestin interface']).prefetch_related('residue_position')

    ## ProteinSegment for different proteins
    ss = ProteinSegment.objects.filter(partial=False, proteinfamily='Arrestin').prefetch_related('generic_numbers')
    ss_cats = ss.values_list('category').order_by('category').distinct('category')

class SegmentSelectionSignature(AbsSegmentSelection):
    step = 3
    number_of_steps = 4

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', False),
        ('segments', True),
    ])
    buttons = {
        'continue': {
            'label': 'Calculate sequence signature',
            'url': '/alignment/render_signature',
            'color': 'success',
        },
    }

class BlastSearchInput(AbsMiscSelection):
    step = 1
    number_of_steps = 1
    docs = 'sequences.html#similarity-search-blast'
    title = 'BLAST search'
    description = 'Enter a sequence into the text box and press the green button.'
    buttons = {
        'continue': {
            'label': 'BLAST',
            'onclick': 'document.getElementById("form").submit()',
            'color': 'success',
        },
    }
    selection_boxes = {}
    blast_input = True

class BlastSearchResults(TemplateView):
    """
    An interface for blast similarity search of the input sequence.
    """
    template_name="blast/blast_search_results.html"

    def post(self, request, *args, **kwargs):

        if 'human' in request.POST.keys():
            blast = BlastSearch(blastdb=os.sep.join([settings.STATICFILES_DIRS[0], 'blast', 'protwis_human_blastdb']), top_results=50)
            blast_out = blast.run(request.POST['input_seq'])
        else:
            blast = BlastSearch(top_results=50)
            blast_out = blast.run(request.POST['input_seq'])

        context = {}
        context['results'] = [(Protein.objects.get(pk=x[0]), x[1]) for x in blast_out]
        context["input"] = request.POST['input_seq']

        return render(request, self.template_name, context)

def render_alignment(request):
    # get the user selection from session
    simple_selection = request.session.get('selection', False)

    # create an alignment object
    a = Alignment()

    # load data from selection into the alignment
    a.load_proteins_from_selection(simple_selection)
    a.load_segments_from_selection(simple_selection)

    #create unique proteins_id
    protein_ids = []
    for p in a.proteins:
        protein_ids.append(p.pk)
    protein_list = ','.join(str(x) for x in sorted(protein_ids))

    #create unique proteins_id
    segments_ids = []
    for s in a.segments:
        segments_ids.append(s)
    segments_list = ','.join(str(x) for x in sorted(segments_ids))

    s = str(protein_list+"_"+segments_list)
    key = "ALIGNMENT_"+hashlib.md5(s.encode('utf-8')).hexdigest()
    return_html = cache_alignment.get(key)

    if return_html==None or 'Custom' in segments_ids:
        # build the alignment data matrix
        check = a.build_alignment()
        if check == 'Too large':
            return render(request, 'alignment/error.html', {'proteins': len(a.proteins), 'residues':a.number_of_residues_total})
        # calculate consensus sequence + amino acid and feature frequency
        a.calculate_statistics()

        num_of_sequences = len(a.proteins)
        num_residue_columns = len(a.positions) + len(a.segments)

        return_html = render(request, 'alignment/alignment.html', {'a': a, 'num_of_sequences': num_of_sequences,
            'num_residue_columns': num_residue_columns})
    if 'Custom' not in segments_ids:
        #update it if used
        cache_alignment.set(key,return_html, 60*60*24*7) #set alignment cache one week

    return return_html

def render_family_alignment(request, slug):
    # create an alignment object
    a = Alignment()

    # fetch proteins and segments
    proteins = Protein.objects.filter(family__slug__startswith=slug, sequence_type__slug='wt')

    if len(proteins)>50 and len(slug.split("_"))<4:
        # If alignment is going to be too big, only pick human.
        proteins = Protein.objects.filter(family__slug__startswith=slug, sequence_type__slug='wt', species__latin_name='Homo sapiens')

    if slug.startswith('100'):

        gsegments = definitions.G_PROTEIN_SEGMENTS

        preserved = Case(*[When(slug=pk, then=pos) for pos, pk in enumerate(gsegments['Full'])])
        segments = ProteinSegment.objects.filter(slug__in = gsegments['Full'], partial=False).order_by(preserved)
    else:
        segments = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR')
        if len(proteins)>50:
            # if a lot of proteins, exclude some segments
            segments = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR').exclude(slug__in=['N-term','C-term'])
        if len(proteins)>200:
            # if many more proteins exluclude more segments
            segments = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR').exclude(slug__in=['N-term','C-term']).exclude(category='loop')

    protein_ids = []
    for p in proteins:
        protein_ids.append(p.pk)
    protein_list = ','.join(str(x) for x in sorted(protein_ids))

    #create unique proteins_id
    segments_ids = []
    for s in segments:
        segments_ids.append(s.slug)
    segments_list = ','.join(str(x) for x in sorted(segments_ids))

    s = str(protein_list+"_"+segments_list)
    key = "ALIGNMENT_"+hashlib.md5(s.encode('utf-8')).hexdigest()
    return_html = cache_alignment.get(key)

    if return_html==None:
        # load data into the alignment
        a.load_proteins(proteins)
        a.load_segments(segments)

        # build the alignment data matrix
        a.build_alignment()

        # calculate consensus sequence + amino acid and feature frequency
        a.calculate_statistics()

        num_of_sequences = len(a.proteins)
        num_residue_columns = len(a.positions) + len(a.segments)

        return_html = render(request, 'alignment/alignment.html', {'a': a, 'num_of_sequences': num_of_sequences,
        'num_residue_columns': num_residue_columns})

    #update it if used
    cache_alignment.set(key,return_html, 60*60*24*7) #set alignment cache one week

    return return_html

def render_fasta_alignment(request):
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

    response = render(request, 'alignment/alignment_fasta.html', context={'a': a}, content_type='text/fasta')
    response['Content-Disposition'] = "attachment; filename=" + settings.SITE_TITLE + "_alignment.fasta"
    return response

def render_fasta_family_alignment(request, slug):
    # create an alignment object
    a = Alignment()
    a.show_padding = False

    # fetch proteins and segments
    proteins = Protein.objects.filter(family__slug__startswith=slug, sequence_type__slug='wt')
    segments = ProteinSegment.objects.filter(partial=False)

    # load data into the alignment
    a.load_proteins(proteins)
    a.load_segments(segments)

    # build the alignment data matrix
    a.build_alignment()

    response = render(request, 'alignment/alignment_fasta.html', context={'a': a}, content_type='text/fasta')
    response['Content-Disposition'] = "attachment; filename=" + settings.SITE_TITLE + "_alignment.fasta"
    return response

def render_csv_alignment(request):
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
    a.calculate_statistics()

    response = render(request, 'alignment/alignment_csv.html', context={'a': a}, content_type='text/csv')
    response['Content-Disposition'] = "attachment; filename=" + settings.SITE_TITLE + "_alignment.csv"
    return response

def render_reordered (request, group):

    #grab the selections from session data
    #targets set #1
    ss_pos = request.session.get('targets_pos', False)
    #targets set #2
    ss_neg = request.session.get('selection', False)

    aln = Alignment()

    if group == 'positive':
        aln.load_proteins_from_selection(ss_pos)
    elif group == 'negative':
        aln.load_proteins_from_selection(ss_neg)

    aln.load_segments_from_selection(ss_neg)
    aln.build_alignment()
    aln.calculate_statistics()

    return render(request, 'alignment/alignment_reordered.html', context={'aln': aln})

def render_signature (request):

    # Minimum difference between positive and negative set
    # TODO: put these params to settings?
    min_diff = 30
    # Minimum conservation of the feature
    min_cons = 0

    colors = {
        'neg': 'crimson',
        'pos': 'royalblue',
        'sma': 'gray',
        'pol': 'blueviolet',
        'hb': 'blueviolet',
        'hbd': 'blue',
        'hba': 'violet',
        'alhp': 'gold',
        'hp': 'gold',
        'ar': 'limegreen',
        'arhp': 'limegreen',
        'charge': 'blue',
        'lar': 'black'
        }
    # grab the selections from session data

    # targets set #1
    ss_pos = request.session.get('targets_pos', False)
    # targets set #2
    ss_neg = request.session.get('selection', False)

    # calculate the alignments
    aln_pos = Alignment()
    aln_pos.load_proteins_from_selection(ss_pos)
    aln_pos.load_segments_from_selection(ss_neg)
    aln_pos.build_alignment()
    aln_pos.calculate_statistics()

    aln_neg = Alignment()
    aln_neg.load_proteins_from_selection(ss_neg)
    aln_neg.load_segments_from_selection(ss_neg)
    aln_neg.build_alignment()
    aln_neg.calculate_statistics()

    # get common set of residues
    common_gn = deepcopy(aln_pos.generic_numbers)

    for scheme in aln_neg.numbering_schemes:
        for segment in aln_neg.segments:
            for pos in aln_neg.generic_numbers[scheme[0]][segment].items():
                if pos[0] not in common_gn[scheme[0]][segment].keys():
                    common_gn[scheme[0]][segment][pos[0]] = pos[1]
            common_gn[scheme[0]][segment] = OrderedDict(sorted(
                common_gn[scheme[0]][segment].items(),
                key=lambda x: x[0].split('x')
                ))

    # get fature and amino acid frequency normalized to the common set of residues
    feats_pos_norm = OrderedDict()
    feats_neg_norm = OrderedDict()
    aa_pos_norm = OrderedDict()
    aa_neg_norm = OrderedDict()

    for sid, segment in enumerate(aln_neg.segments):
        feats_pos_norm[segment] = np.array(
            [[x[0] for x in feat[sid]] for feat in aln_pos.feature_stats],
            dtype='int'
            )
        feats_neg_norm[segment] = np.array(
            [[x[0] for x in feat[sid]] for feat in aln_neg.feature_stats],
            dtype='int'
            )
        aa_pos_norm[segment] = np.array(
            [[x[0] for x in aa[sid]] for aa in aln_pos.amino_acid_stats],
            dtype='int'
            )
        aa_neg_norm[segment] = np.array(
            [[x[0] for x in aa[sid]] for aa in aln_neg.amino_acid_stats],
            dtype='int'
            )

    # normalize and calculate difference
    feats_delta = OrderedDict()
    aa_delta = OrderedDict()

    for segment in aln_neg.segments:
        #TODO: get the correct default numering scheme from settings
        for idx, res in enumerate(common_gn['gpcrdba'][segment].keys()):
            if res not in aln_pos.generic_numbers['gpcrdba'][segment].keys():
                feats_pos_norm[segment] = np.insert(feats_pos_norm[segment], idx, 0, axis=1)
                aa_pos_norm[segment] = np.insert(aa_pos_norm[segment], idx, 0, axis=1)
            elif res not in aln_neg.generic_numbers['gpcrdba'][segment].keys():
                feats_neg_norm[segment] = np.insert(feats_neg_norm[segment], idx, 0, axis=1)
                aa_neg_norm[segment] = np.insert(aa_neg_norm[segment], idx, 0, axis=1)

        # now the difference
        feats_delta[segment] = np.subtract(feats_pos_norm[segment], feats_neg_norm[segment])
        aa_delta[segment] = np.subtract(aa_pos_norm[segment], aa_neg_norm[segment])

    # adapt the feature stats to the format from Alignment objects
    feats_delta_disp = []
    for row, feat in enumerate(definitions.AMINO_ACID_GROUPS.keys()):
        tmp_row = []
        for segment in aln_neg.segments:
            #first item is the real value,
            # second is the assignmnent of color (via css)
            # 0 - red, 5 - yellow, 10 - green
            #third item is a tooltip
            tmp_row.append([[
                x,
                int(x/20)+5,
                "{} - {}".format(
                    feats_pos_norm[segment][row][y],
                    feats_neg_norm[segment][row][y]
                    )
                ] for y, x in enumerate(feats_delta[segment][row])])
        feats_delta_disp.append(tmp_row)

    # adapt the amino acids stats to the format from Alignment objects
    aa_delta_disp = []
    for row, aa in enumerate(definitions.AMINO_ACIDS.keys()):
        tmp_row = []
        for segment in aln_neg.segments:
            #first item is the real value,
            # second is the assignmnent of color (via css)
            # 0 - red, 5 - yellow, 10 - green
            #third item is a tooltip
            tmp_row.append([[
                x,
                int(x/20)+5,
                "{} - {}".format(
                    aa_pos_norm[segment][row][y],
                    aa_neg_norm[segment][row][y]
                    )
                ] for y, x in enumerate(aa_delta[segment][row])])
        aa_delta_disp.append(tmp_row)



    # save for later
    # signature_map = feats_delta.argmax(axis=0)


    # signature_data = []
    # for i, f_id in enumerate(signature_map):
    #     signature_data.append([
    #         i,
    #         colors[list(definitions.AMINO_ACID_GROUPS.keys())[f_id]],
    #         int(feats_delta[f_id][i])
    #         ])

    # options = {
    #     'xticks': residue_positions_common,
    #     'anchor': 'signature_plot'
    # }

    num_of_sequences_pos = len(aln_pos.proteins)
    num_residue_columns_pos = len(aln_pos.positions)
    num_of_sequences_neg = len(aln_neg.proteins)
    num_residue_columns_neg = len(aln_neg.positions)
    return_html = render(request, 'sequence_signature/sequence_signature.html', {
        # 'signature_data': json.dumps(signature_data),
        # 'signature_options': json.dumps(options),
        'num_residue_columns': len([x for x in common_gn['gpcrdba'][segment] for segment in aln_neg.segments]),
        'num_of_sequences_pos': num_of_sequences_pos,
        'num_residue_columns_pos': num_residue_columns_pos,
        'num_of_sequences_neg': num_of_sequences_neg,
        'num_residue_columns_neg': num_residue_columns_neg,

        'common_segments': aln_neg.segments,
        'common_generic_numbers': common_gn,
        'feats_signature': feats_delta_disp,
        'aa_signature': aa_delta_disp,
        'a_pos': aln_pos,
        'a_neg': aln_neg,
        })

    return return_html

def render_signature_excel (request):

    # version #1 - 8 sheets with separate pieces of signature outline

    definitions.AMINO_ACID_GROUPS['custom'] = ['A', 'Q']

    # step 1 - repeat the data preparation for a sequence signature

    # targets set #1
    ss_pos = request.session.get('targets_pos', False)
    # targets set #2
    ss_neg = request.session.get('selection', False)

    # calculate the alignments
    aln_pos = Alignment()
    aln_pos.load_proteins_from_selection(ss_pos)
    aln_pos.load_segments_from_selection(ss_neg)
    aln_pos.build_alignment()
    aln_pos.calculate_statistics()

    aln_neg = Alignment()
    aln_neg.load_proteins_from_selection(ss_neg)
    aln_neg.load_segments_from_selection(ss_neg)
    aln_neg.build_alignment()
    aln_neg.calculate_statistics()

    print(functions.get_numbering_schemes(aln_pos.proteins + aln_neg.proteins))
    # get common set of residues
    common_gn = deepcopy(aln_pos.generic_numbers)

    for scheme in aln_neg.numbering_schemes:
        for segment in aln_neg.segments:
            for pos in aln_neg.generic_numbers[scheme[0]][segment].items():
                if pos[0] not in common_gn[scheme[0]][segment].keys():
                    common_gn[scheme[0]][segment][pos[0]] = pos[1]
            common_gn[scheme[0]][segment] = OrderedDict(sorted(
                common_gn[scheme[0]][segment].items(),
                key=lambda x: x[0].split('x')
                ))

    # get fature and amino acid frequency normalized to the common set of residues
    feats_pos_norm = OrderedDict()
    feats_neg_norm = OrderedDict()
    aa_pos_norm = OrderedDict()
    aa_neg_norm = OrderedDict()

    for sid, segment in enumerate(aln_neg.segments):
        feats_pos_norm[segment] = np.array(
            [[x[0] for x in feat[sid]] for feat in aln_pos.feature_stats],
            dtype='int'
            )
        feats_neg_norm[segment] = np.array(
            [[x[0] for x in feat[sid]] for feat in aln_neg.feature_stats],
            dtype='int'
            )
        aa_pos_norm[segment] = np.array(
            [[x[0] for x in aa[sid]] for aa in aln_pos.amino_acid_stats],
            dtype='int'
            )
        aa_neg_norm[segment] = np.array(
            [[x[0] for x in aa[sid]] for aa in aln_neg.amino_acid_stats],
            dtype='int'
            )

    # normalize and calculate difference
    feats_delta = OrderedDict()
    aa_delta = OrderedDict()

    for segment in aln_neg.segments:
        #TODO: get the correct default numering scheme from settings
        for idx, res in enumerate(common_gn['gpcrdba'][segment].keys()):
            if res not in aln_pos.generic_numbers['gpcrdba'][segment].keys():
                feats_pos_norm[segment] = np.insert(feats_pos_norm[segment], idx, 0, axis=1)
                aa_pos_norm[segment] = np.insert(aa_pos_norm[segment], idx, 0, axis=1)
            elif res not in aln_neg.generic_numbers['gpcrdba'][segment].keys():
                feats_neg_norm[segment] = np.insert(feats_neg_norm[segment], idx, 0, axis=1)
                aa_neg_norm[segment] = np.insert(aa_neg_norm[segment], idx, 0, axis=1)

        # now the difference
        feats_delta[segment] = np.subtract(feats_pos_norm[segment], feats_neg_norm[segment])
        aa_delta[segment] = np.subtract(aa_pos_norm[segment], aa_neg_norm[segment])

    # adapt the feature stats to the format from Alignment objects
    feats_delta_disp = []
    for row, feat in enumerate(definitions.AMINO_ACID_GROUPS.keys()):
        tmp_row = []
        for segment in aln_neg.segments:
            #first item is the real value,
            # second is the assignmnent of color (via css)
            # 0 - red, 5 - yellow, 10 - green
            #third item is a tooltip
            tmp_row.append([[
                x,
                int(x/20)+5,
                "{} - {}".format(
                    feats_pos_norm[segment][row][y],
                    feats_neg_norm[segment][row][y]
                    )
                ] for y, x in enumerate(feats_delta[segment][row])])
        feats_delta_disp.append(tmp_row)

    # adapt the amino acids stats to the format from Alignment objects
    aa_delta_disp = []
    for row, aa in enumerate(definitions.AMINO_ACIDS.keys()):
        tmp_row = []
        for segment in aln_neg.segments:
            #first item is the real value,
            # second is the assignmnent of color (via css)
            # 0 - red, 5 - yellow, 10 - green
            #third item is a tooltip
            tmp_row.append([[
                x,
                int(x/20)+5,
                "{} - {}".format(
                    aa_pos_norm[segment][row][y],
                    aa_neg_norm[segment][row][y]
                    )
                ] for y, x in enumerate(aa_delta[segment][row])])
        aa_delta_disp.append(tmp_row)

    outstream = BytesIO()
    # wb = xlsxwriter.Workbook('excel_test.xlsx', {'in_memory': False})
    wb = xlsxwriter.Workbook(outstream, {'in_memory': True})
    # Feature stats for signature
    functions.prepare_excel_worksheet(
        wb,
        'signature_properties',
        aln_neg,
        common_gn,
        feats_delta_disp
    )
    # Amino acids stats for signature
    functions.prepare_excel_worksheet(
        wb,
        'signature_amino_acids',
        aln_neg,
        common_gn,
        aa_delta_disp,
        data_type='amino_acids'
    )
    # Feature stats for positive group alignment
    functions.prepare_excel_worksheet(
        wb,
        'positive_group_properties',
        aln_pos,
        data_block=aln_pos.feature_stats
    )
    # Amino acids stats for positive group alignment
    functions.prepare_excel_worksheet(
        wb,
        'positive_group_amino_acids',
        aln_pos,
        data_block=aln_pos.amino_acid_stats,
        data_type='amino_acids'
    )
    # Positive group alignment
    functions.prepare_excel_worksheet(
        wb,
        'positive_group_aln',
        aln_pos
    )
    # Feature stats for negative group alignment
    functions.prepare_excel_worksheet(
        wb,
        'negative_group_properties',
        aln_neg,
        data_block=aln_neg.feature_stats
    )
    # Amino acids stats for negative group alignment
    functions.prepare_excel_worksheet(
        wb,
        'negative_group_amino_acids',
        aln_neg,
        data_block=aln_neg.amino_acid_stats,
        data_type='amino_acids'
    )
    # Negative group alignment
    functions.prepare_excel_worksheet(
        wb,
        'negative_group_aln',
        aln_neg
    )

    wb.close()
    outstream.seek(0)
    response = HttpResponse(
        outstream.read(),
        content_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )
    response['Content-Disposition'] = "attachment; filename=sequence_signature.xlsx"

    return response
