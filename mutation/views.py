from django.shortcuts import render
from django.template import loader, Context
from django.db.models import Count, Min, Sum, Avg, Q
from django.http import HttpResponse
from django.conf import settings
from django.core.cache import cache
from django.views.decorators.cache import cache_page
from mutation.functions import *
from mutation.models import *

from common.views import AbsTargetSelection
from common.views import AbsSegmentSelection
from common.diagrams_gpcr import DrawHelixBox, DrawSnakePlot
from common import definitions

from residue.models import Residue,ResidueNumberingScheme, ResidueGenericNumberEquivalent
from residue.views import ResidueTablesDisplay
from protein.models import Protein,ProteinSegment,ProteinFamily,ProteinConformation
from interaction.models import ResidueFragmentInteraction, StructureLigandInteraction
from interaction.views import calculate
from interaction.forms import PDBform

from datetime import datetime
from collections import OrderedDict
import json
import yaml
import os
import copy
#env/bin/python3 -m pip install xlrd
import csv
from io import BytesIO
import re
import math
import urllib
import xlsxwriter #sudo pip3 install XlsxWriter
import operator

Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')

class TargetSelection(AbsTargetSelection):
    step = 1
    number_of_steps = 2
    docs = 'mutations.html#mutation-browser'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/mutations/segmentselection',
            'color': 'success',
        },
    }
    default_species = False


class SegmentSelection(AbsSegmentSelection):
    step = 2
    number_of_steps = 2
    docs = 'mutations.html#mutation-browser'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', True),
    ])
    buttons = {
        'continue': {
            'label': 'Show mutants',
            'url': '/mutations/render',
            'color': 'success',
        },
    }

def render_mutations(request, protein = None, family = None, download = None, receptor_class = None, gn = None, aa = None, **response_kwargs):

    # get the user selection from session
    simple_selection = request.session.get('selection', False)

     # local protein list
    proteins = []
    alignment_proteins = []
    original_positions = []
    used_schemes = {}
    # print("receptor_class",receptor_class,family)
    if receptor_class==None:

        if protein: # if protein static page

            proteins.append(Protein.objects.get(entry_name = protein))
            segments_ids = ProteinSegment.objects.filter(proteinfamily='GPCR').values('id')
            original_segments = ProteinSegment.objects.all()

        elif family:

            family_proteins = Protein.objects.filter(family__slug__startswith=family, sequence_type__slug='wt')
            for fp in family_proteins:
                proteins.append(fp)
                if fp.residue_numbering_scheme.slug not in used_schemes:
                    used_schemes[fp.residue_numbering_scheme.slug] = 0
                used_schemes[fp.residue_numbering_scheme.slug] += 1
            segments_ids = ProteinSegment.objects.filter(proteinfamily='GPCR').values('id')
            original_segments = ProteinSegment.objects.filter(proteinfamily='GPCR')
            used_scheme = max(used_schemes, key=used_schemes.get)

        else:

            # flatten the selection into individual proteins
            for target in simple_selection.targets:
                if target.type == 'protein':
                    proteins.append(target.item)
                elif target.type == 'family':
                    # species filter
                    species_list = []
                    for species in simple_selection.species:
                        species_list.append(species.item)

                    # annotation filter
                    protein_source_list = []
                    for protein_source in simple_selection.annotation:
                        protein_source_list.append(protein_source.item)

                    if species_list:
                        family_proteins = Protein.objects.filter(family__slug__startswith=target.item.slug,
                            species__in=(species_list),
                            source__in=(protein_source_list)).select_related('residue_numbering_scheme', 'species')
                    else:
                        family_proteins = Protein.objects.filter(family__slug__startswith=target.item.slug,
                            source__in=(protein_source_list)).select_related('residue_numbering_scheme', 'species')

                    for fp in family_proteins:
                        proteins.append(fp)

            original_segments = []
            original_positions = []
            segments_ids = []
            for segment in simple_selection.segments:
                if segment.type=='residue':
                    original_positions.append(segment.item.default_generic_number.label)
                else:
                    original_segments.append(segment.item)
                segments_ids.append(segment.item.id)

        #scheme
        used_schemes = {}
        species_list = {}
        longest_name = 0
        for entry in proteins:
            if entry.species.common_name=='Human':
                alignment_proteins.append(entry)
            if entry.residue_numbering_scheme.slug not in used_schemes:
                used_schemes[entry.residue_numbering_scheme.slug] = 0
            used_schemes[entry.residue_numbering_scheme.slug] += 1
            if entry.species.common_name not in species_list:
                if len(entry.species.common_name)>10 and len(entry.species.common_name.split())>1:
                    name = entry.species.common_name.split()[0][0]+". "+" ".join(entry.species.common_name.split()[1:])
                    if len(" ".join(entry.species.common_name.split()[1:]))>11:
                        name = entry.species.common_name.split()[0][0]+". "+" ".join(entry.species.common_name.split()[1:])[:8]+".."
                else:
                    name = entry.species.common_name
                species_list[entry.species.common_name] = name

                if len(re.sub('<[^>]*>', '', entry.name)+" "+name)>longest_name:
                    longest_name = len(re.sub('<[^>]*>', '', entry.name)+" "+name)


        if len(alignment_proteins)==0:
            alignment_proteins = proteins

        if len(proteins)==1:
            protein = proteins[0]

        used_scheme = max(used_schemes, key=used_schemes.get)
        mutations = MutationExperiment.objects.filter(
                                Q(protein__in=proteins),
                                Q(residue__protein_segment__in=original_segments) | Q(residue__generic_number__label__in=original_positions)
                                ).prefetch_related('residue__display_generic_number',
                                'residue__protein_segment','residue__generic_number','exp_func','exp_qual',
                                'exp_measure', 'exp_type', 'ligand_role', 'ligand','refs','raw',
                                'ligand__properities', 'refs__web_link', 'refs__web_link__web_resource', 'review__web_link__web_resource','protein','mutation__protein')
    else:
        # print(gn,receptor_class,aa)
        protein_ids = ''
        HelixBox = ''
        SnakePlot = ''
        numbering_schemes = ''
        mutations_pos_list = []
        family_proteins = Protein.objects.filter(family__slug__startswith=receptor_class, sequence_type__slug='wt').all()[0]
        used_schemes[family_proteins.residue_numbering_scheme.slug] = 1
        used_scheme = max(used_schemes, key=used_schemes.get)
        mutations = MutationExperiment.objects.filter(protein__family__slug__startswith=receptor_class,
                                residue__generic_number__label=gn, residue__amino_acid=aa).prefetch_related('residue__generic_number',
                                'residue__protein_segment','residue__generic_number','exp_func','exp_qual',
                                'exp_measure', 'exp_type', 'ligand_role', 'ligand','refs','raw',
                                'ligand__properities', 'refs__web_link', 'refs__web_link__web_resource', 'review__web_link__web_resource','protein','mutation__protein')

    mutations_list = {}
    mutations_list_seq = {}
    mutations_generic_number = {}
    mutations_display_generic_number = {}
    mutations_class_generic_number = {}
    context = {}
    context['data'] = ''
    context['header'] = ''
    context['segments'] = ''
    context['number_of_schemes'] = ''
    context['longest_name'] = {'div' : 0, 'height': 0}

    gn_lookup = {}

    residue_table_list = []

    gn_labels = mutations.values_list('mutation__residue__generic_number__label', flat=True)
    class_gns = ResidueGenericNumberEquivalent.objects.filter(default_generic_number__label__in = gn_labels, scheme__slug = used_scheme).prefetch_related('default_generic_number').all()

    for class_gn in class_gns:
        label = class_gn.default_generic_number
        gn_lookup[label] = class_gn

    import urllib.parse
    mutation_tables = ''
    for mutation in mutations:
        residue_table_list.append(mutation.residue.generic_number)

        # if not mutation.residue.generic_number: continue #cant map those without display numbers
        if mutation.residue.generic_number and mutation.residue.generic_number.label not in mutations_list: mutations_list[mutation.residue.generic_number.label] = []
        if mutation.residue.sequence_number not in mutations_list_seq: mutations_list_seq[mutation.residue.sequence_number] = [[]]

        exp_type = "N/A"
        if mutation.exp_func and mutation.exp_type:
            exp_type = mutation.exp_type.type+" ("+mutation.exp_func.func+")"
        elif mutation.exp_func:
            exp_type = "("+mutation.exp_func.func+")"


        if exp_type=='N/A' and (mutation.opt_basal_activity or mutation.opt_receptor_expression):
            # If N/A and optional data something is irrevalent.
            continue

        if mutation.ligand:
            ligand = mutation.ligand.name
        else:
            ligand = ''
        if mutation.exp_qual:
            qual = mutation.exp_qual.qual
        else:
            qual = ''

        mutations_list_seq[mutation.residue.sequence_number][0].append([mutation.foldchange,ligand.replace('\xe2', "").replace('\'', ""),qual])
        if mutation.residue.generic_number:
            mutations_list[mutation.residue.generic_number.label].append([mutation.foldchange,ligand.replace('\xe2', "").replace('\'', ""),qual])
            if mutation.residue.generic_number not in gn_lookup:
                class_gn = ResidueGenericNumberEquivalent.objects.filter(default_generic_number = mutation.residue.generic_number, scheme__slug = used_scheme).get()
                gn_lookup[mutation.residue.generic_number] = class_gn
            else:
                class_gn = gn_lookup[mutation.residue.generic_number]
            mutations_display_generic_number[mutation.raw.id] = mutation.residue.display_generic_number.label
            mutations_generic_number[mutation.raw.id] = mutation.residue.generic_number.label
            mutations_class_generic_number[mutation.raw.id] = class_gn.label
            gn_display = mutation.residue.display_generic_number.label
        else:
            gn_display = ''
        mutation.refs_link = ""
        mutation.refs_title = ""
        mutation.refs_main = ""
        mutation.review_link = ""
        mutation.review_title = ""
        mutation.review_main = ""
        if mutation.refs:
            mutation.refs_link = mutation.refs.web_link
            mutation.refs_title = mutation.refs.title
            mutation.refs_main = mutation.citation()
        if mutation.review:
            mutation.review_link = mutation.review.web_link
            mutation.review_title = mutation.review.title
            mutation.review_main = mutation.review_citation()

        smiles = ""
        if mutation.ligand and mutation.ligand.properities and mutation.ligand.properities.smiles:
            smiles = urllib.parse.quote_plus(mutation.ligand.properities.smiles)

        lig_name = ""
        lig_role_name = ""
        if mutation.ligand:
            lig_name = mutation.ligand.name
        if mutation.ligand_role:
            lig_role_name = mutation.ligand_role.name

        if not mutation.review_title:
            mutation.review_title = ''

        if not mutation.refs_title:
            mutation.refs_title = ''

        row = '''
                <tr>
                <td><a href="/protein/%s">%s</a></td>
                <td>%s</td>
                <td>%s</td>
                <td>%s</td>
                <td>%s => %s</td>

                <td><a class="foldchange" href="#" data-toogle="tooltip"  data-html="true" data-original-title="%s" data-placement="top"> %s</a></td>

                <td>%s</td>
                <td>
                <a class="smiles-tooltip" data-toggle="tooltip"  data-html="true" data-original-title="<img src='http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/%s/PNG'><br>%s<br>%s" data-placement="right">%s</a></td>
                <td>
                <a class="citation-tooltip" target="_blank" href="%s" data-toggle="tooltip"  data-container="body" data-html="true" data-original-title="%s" data-placement="left" >%s</a></td>
                <td>
                <a class="citation-tooltip" target="_blank" href="%s" data-toggle="tooltip"  data-container="body" data-html="true" data-original-title="%s" data-placement="left" >%s</a></td>
                 </tr>
        ''' % (mutation.protein.entry_name,mutation.protein.entry_name,gn_display,mutation.residue.sequence_number,mutation.residue.protein_segment.slug,
               mutation.residue.amino_acid, mutation.mutation.amino_acid,mutation.getCalculation(),mutation.getFoldorQual(),exp_type,smiles,lig_name,
               lig_role_name,lig_name,mutation.refs_link,mutation.refs_title,mutation.refs_main,mutation.review_link,mutation.review_title,mutation.review_main)
        mutation_tables += row
    # mutation_tables = ''.join(mutation_tables)
        #print(mutation.refs.web_link)

    if receptor_class==None and not download: #if not a small lookup
        # create an alignment object
        #print(proteins)
        #alignment_proteins = Protein.objects.filter(protein__in=proteins)
        protein_hash = hash(tuple(sorted(ProteinConformation.objects.filter(protein__in=alignment_proteins).values_list('id',flat=True))))

        excluded_segment = ['C-term','N-term']
        excluded_segment = []

        protein_ids = list(set([x.id for x in proteins]))
        if protein:
            #if protein do something else
            SnakePlot = proteins[0].get_snake_plot_no_buttons()
            HelixBox = proteins[0].get_helical_box_no_buttons()

            # Fix for plots: convert generic numbering of mutations_list_seq to protein positions
            if len(proteins)>1:
                mutations_pos_list = {}
                if len(mutations_list) > 0:
                    residuelist = Residue.objects.filter(protein_conformation__protein__entry_name=str(proteins[0]), generic_number__label__in=mutations_list.keys()).prefetch_related('display_generic_number','generic_number')
                    for residue in residuelist:
                        if residue.generic_number and residue.generic_number.label in mutations_list:
                            mutations_pos_list[residue.sequence_number] = mutations_list[residue.generic_number.label]
            else:
                mutations_pos_list = mutations_list_seq

            segments = ProteinSegment.objects.filter(proteinfamily='GPCR', pk__in=segments_ids,category='helix')
        else:
            segments = ProteinSegment.objects.filter(proteinfamily='GPCR').exclude(slug__in = excluded_segment).prefetch_related()
            segment_hash = hash(tuple(sorted(segments.values_list('id',flat=True))))

            consensus = cache.get(str(protein_hash)+"&"+str(segment_hash)+"&consensus")
            generic_number_objs = cache.get(str(protein_hash)+"&"+str(segment_hash)+"&generic_number_objs")
            if generic_number_objs == None or consensus == None or consensus == []:
                a = Alignment()

                a.load_proteins(alignment_proteins)

                a.load_segments(segments) #get all segments to make correct diagrams

                # build the alignment data matrix
                a.build_alignment()

                # calculate consensus sequence + amino acid and feature frequency
                a.calculate_statistics()

                consensus = a.full_consensus
                generic_number_objs = a.generic_number_objs
                cache.set(str(protein_hash)+"&"+str(segment_hash)+"&consensus",consensus)
                cache.set(str(protein_hash)+"&"+str(segment_hash)+"&generic_number_objs",generic_number_objs)

            residue_list = []
            generic_numbers = []
            reference_generic_numbers = {}
            count = 1 #build sequence_number

            mutations_pos_list = {}
            for aa in consensus:
                # for aa,v in a.full_consensus[seg].items():
                r = Residue()
                r.sequence_number =  aa.sequence_number #FIXME is this certain to be correct that the position in consensus is seq position?
                #print(aa,aa.family_generic_number,aa.generic_number)
                if aa.family_generic_number and aa.family_generic_number in generic_number_objs:
                    r.generic_number = generic_number_objs[aa.family_generic_number] #FIXME
                    if aa.family_generic_number in mutations_list:
                        if r.sequence_number not in mutations_pos_list:
                            mutations_pos_list[r.sequence_number] = []
                        mutations_pos_list[r.sequence_number].append(mutations_list[aa.family_generic_number])
                    r.segment_slug = aa.segment_slug
                    r.family_generic_number = aa.family_generic_number
                else:
                    r.segment_slug = aa.segment_slug
                    r.family_generic_number = aa.family_generic_number
                r.amino_acid = aa.amino_acid
                r.frequency = aa.frequency #Grab consensus information
                residue_list.append(r)

                count += 1

            HelixBox = DrawHelixBox(consensus,'Class A',str(protein_ids), nobuttons = 1)
            SnakePlot = DrawSnakePlot(consensus,'Class A',str(protein_ids), nobuttons = 1)

        numbering_schemes_selection = ['gpcrdba','gpcrdbb','gpcrdbc','gpcrdbf'] #there is a residue_numbering_scheme attribute on the protein model, so it's easy to find out
        numbering_schemes_selection = list(used_schemes.keys())
        numbering_schemes_selection = ['gpcrdb'] + numbering_schemes_selection #always use A for reference
        numbering_schemes = ResidueNumberingScheme.objects.filter(slug__in=numbering_schemes_selection).all()

        if ResidueNumberingScheme.objects.get(slug=settings.DEFAULT_NUMBERING_SCHEME) in numbering_schemes:
            default_scheme = ResidueNumberingScheme.objects.get(slug=settings.DEFAULT_NUMBERING_SCHEME)
        else:
            default_scheme = numbering_schemes[0]

    # prepare the dictionary
    # each helix has a dictionary of positions
    # default_generic_number or first scheme on the list is the key
    # value is a dictionary of other gn positions and residues from selected proteins

        if len(protein_ids)<20 and receptor_class==None: #too many to make meaningful residuetable / not run when download
            data = OrderedDict()
            for segment in segments:
                data[segment.slug] = OrderedDict()
                residues = Residue.objects.filter(protein_segment=segment, protein_conformation__protein__in=proteins,
                                                generic_number__in=residue_table_list).prefetch_related('protein_conformation__protein',
                                                'protein_conformation__state', 'protein_segment',
                                                'generic_number','display_generic_number','generic_number__scheme', 'alternative_generic_numbers__scheme')
                pos_list = residues.values_list('generic_number__label',flat=True)
                for scheme in numbering_schemes:
                    if scheme == default_scheme and scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                        for pos in pos_list:
                            data[segment.slug][pos] = {scheme.slug : pos, 'seq' : ['-']*len(proteins)}
                    elif scheme == default_scheme:
                        for pos in pos_list:
                            data[segment.slug][pos] = {scheme.slug : pos, 'seq' : ['-']*len(proteins)}

                for residue in residues:
                    alternatives = residue.alternative_generic_numbers.all()
                    pos = residue.generic_number
                    for alternative in alternatives:
                        if alternative.scheme not in numbering_schemes:
                            continue
                        scheme = alternative.scheme
                        if default_scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                            pos = residue.generic_number
                            if scheme == pos.scheme:
                                data[segment.slug][pos.label]['seq'][proteins.index(residue.protein_conformation.protein)] = str(residue)
                            else:
                                if scheme.slug not in data[segment.slug][pos.label].keys():
                                    data[segment.slug][pos.label][scheme.slug] = alternative.label
                                if alternative.label not in data[segment.slug][pos.label][scheme.slug]:
                                    data[segment.slug][pos.label][scheme.slug] += " "+alternative.label
                                data[segment.slug][pos.label]['seq'][proteins.index(residue.protein_conformation.protein)] = str(residue)
                        else:
                            if scheme.slug not in data[segment.slug][pos.label].keys():
                                data[segment.slug][pos.label][scheme.slug] = alternative.label
                            if alternative.label not in data[segment.slug][pos.label][scheme.slug]:
                                data[segment.slug][pos.label][scheme.slug] += " "+alternative.label
                            data[segment.slug][pos.label]['seq'][proteins.index(residue.protein_conformation.protein)] = str(residue)

            # Preparing the dictionary of list of lists. Dealing with tripple nested dictionary in django templates is a nightmare
            flattened_data = OrderedDict.fromkeys([x.slug for x in segments], [])

            for s in iter(flattened_data):
                flattened_data[s] = [[data[s][x][y.slug] for y in numbering_schemes]+data[s][x]['seq'] for x in sorted(data[s])]

            if len(species_list)>1:
                context['header'] = zip([x.short_name for x in numbering_schemes] + [x.name+" "+species_list[x.species.common_name] for x in proteins], [x.name for x in numbering_schemes] + [x.name for x in proteins],[x.name for x in numbering_schemes] + [x.entry_name for x in proteins])
            else:
                context['header'] = zip([x.short_name for x in numbering_schemes] + [x.name for x in proteins], [x.name for x in numbering_schemes] + [x.name for x in proteins],[x.name for x in numbering_schemes] + [x.entry_name for x in proteins])
            context['segments'] = [x.slug for x in segments if len(data[x.slug])]
            context['data'] = flattened_data
            context['number_of_schemes'] = len(numbering_schemes)
            context['longest_name'] = {'div' : longest_name*2, 'height': longest_name*2+75}
        else:
            context['data'] = ''
            context['header'] = ''
            context['segments'] = ''
            context['number_of_schemes'] = ''
            context['longest_name'] = {'div' : 0, 'height': 0}

    if download:
        raws = mutations.values('raw')
        rawmutations = MutationRaw.objects.filter(pk__in = raws).all()

        data = []
        for r in rawmutations:
            headers = []
            values = {}
            values = r.__dict__ #print(r.__dict__)
            # for field, val in r:
            #     headers.append(field)
            #     values[field] = val
            # print(values)
            if values['id'] in mutations_class_generic_number:
                values['generic'] = mutations_class_generic_number[values['id']]
            else:
                values['generic'] = ''
            data.append(values)
        headers = ['submitting_group','reference','review','data_container','data_container_number', 'protein', 'mutation_pos', 'generic', 'mutation_from', 'mutation_to',
        'ligand_name', 'ligand_idtype', 'ligand_id', 'ligand_class',
        'exp_type', 'exp_func',  'exp_wt_value',  'exp_wt_unit','exp_mu_effect_sign', 'exp_mu_effect_type', 'exp_mu_effect_value',
        'exp_fold_change',
        'exp_mu_effect_qual', 'exp_mu_effect_ligand_prop',  'exp_mu_ligand_ref', 'opt_receptor_expression', 'opt_basal_activity',
        'opt_gain_of_activity', 'opt_ligand_emax', 'opt_agonist', 'added_date'
         ] #'added_by',

        #EXCEL SOLUTION
        output = BytesIO()
        workbook = xlsxwriter.Workbook(output)
        worksheet = workbook.add_worksheet()

        col = 0
        for h in headers:
            worksheet.write(0, col, h)
            col += 1
        row = 1
        for d in data:
            col = 0
            for h in headers:
                worksheet.write(row, col, d[h])
                col += 1
            row += 1
        workbook.close()
        output.seek(0)
        xlsx_data = output.read()

        response = HttpResponse(xlsx_data,content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
        response['Content-Disposition'] = 'attachment; filename=GPCRdb_mutant_data.xlsx' #% 'mutations'
        return response

    else:
        return render(request, 'mutation/list.html', {'mutation_tables':mutation_tables,'mutations': mutations, 'HelixBox':HelixBox, 'SnakePlot':SnakePlot, 'data':context['data'],
            'header':context['header'], 'longest_name':context['longest_name'], 'segments':context['segments'], 'number_of_schemes':len(numbering_schemes), 'mutations_pos_list' : json.dumps(mutations_pos_list), 'protein_ids':str(protein_ids)})

# Create your views here.
def index(request):
    return HttpResponse("Hello, world. You're at the polls index.")

class designPDB(AbsTargetSelection):

    # Left panel
    step = 1
    number_of_steps = 1
    # docs = 'generic_numbering.html'  # FIXME

    # description = 'Select receptors to index by searching or browsing in the middle column. You can select entire' \
    #     + ' receptor families and/or individual receptors.\n\nSelected receptors will appear in the right column,' \
    #     + ' where you can edit the list.\n\nSelect which numbering schemes to use in the middle column.\n\nOnce you' \
    #     + ' have selected all your receptors, click the green button.'

    description = 'Upload a structure of a receptor and ligand bound. The tool will then deduce the interactions and suggest mutations to verify these.'

    # Middle section
    numbering_schemes = False
    filters = False
    search = True
    title = "Select PDB code or upload PDB file"

    template_name = 'mutation/designselectionpdb.html'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])

    # Buttons
    buttons = {
        'continue': {
            'label': 'Show results',
            'onclick': 'submitupload()',
            'color': 'success',
            #'url': 'calculate/'
        }
    }

    redirect_on_select = False

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        context['structures'] = ResidueFragmentInteraction.objects.values('structure_ligand_pair__structure__pdb_code__index', 'structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name').annotate(
            num_ligands=Count('structure_ligand_pair', distinct=True), num_interactions=Count('pk', distinct=True)).order_by('structure_ligand_pair__structure__pdb_code__index')
        context['form'] = PDBform()
        return context

class design(AbsTargetSelection):

    # Left panel
    step = 1
    number_of_steps = 1
    # docs = 'generic_numbering.html'  # FIXME

    # description = 'Select receptors to index by searching or browsing in the middle column. You can select entire' \
    #     + ' receptor families and/or individual receptors.\n\nSelected receptors will appear in the right column,' \
    #     + ' where you can edit the list.\n\nSelect which numbering schemes to use in the middle column.\n\nOnce you' \
    #     + ' have selected all your receptors, click the green button.'

    description = 'Get mutation suggestions based on interaction data from structures and mutagenesis experimental data.'

    # Middle section
    numbering_schemes = False
    filters = False
    search = True
    title = "Select a receptor"

    template_name = 'mutation/designselection.html'
    selection_only_receptors = True

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])

    # Buttons
    buttons = {
        'continue': {
            'label': 'Show results',
            'onclick': 'submitupload()',
            'color': 'success',
            #'url': 'calculate/'
        }
    }

    redirect_on_select = False

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        context['structures'] = ResidueFragmentInteraction.objects.values('structure_ligand_pair__structure__pdb_code__index', 'structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name').annotate(
            num_ligands=Count('structure_ligand_pair', distinct=True), num_interactions=Count('pk', distinct=True)).order_by('structure_ligand_pair__structure__pdb_code__index')
        context['form'] = PDBform()
        return context


def showcalculationPDB(request):
    if request.method == 'POST':
        form = PDBform(request.POST, request.FILES)

        if 'file' in request.FILES: #uploaded file
            print('user upload')
        else:
            print('pdb code entered')

        context = calculate(request)

        #print(context['residues'])
        matrix = definitions.DESIGN_SUBSTITUTION_MATRIX
        newresidues = []
        for r in context['residues']:
            if r['slug']=='acc':
                continue
            if r['slug'][:5]=='polar':
                scoretype = 'polar'
            elif r['slug'][:3]=='aro':
                scoretype = 'aromatic'
            elif r['slug'][:3]=='hyd':
                scoretype = 'hydrophobic'
            else:
                scoretype = 'unknown'
            possible_subs = ''
            if r['aa'] in matrix[scoretype]:
                possible_subs = matrix[scoretype][r['aa']][0]
                possible_rea = matrix[scoretype][r['aa']][1]

            i = 0
            s = ''
            for p in possible_subs:
                s += ', '.join(p) + " : "+ possible_rea[i] + "<br>"
                i += 1

            r['suggested'] = s
            newresidues.append(r)

        context['residues'] = newresidues

        return render(request, 'mutation/designpdb.html', context)


@cache_page(60 * 60 * 24 *7)
def coverage(request):

    context = {}

    #gpcr_class = '004' #class a

    families = ProteinFamily.objects.all()
    lookup = {}
    for f in families:
        lookup[f.slug] = f.name.replace("receptors","")

    class_proteins = Protein.objects.filter(family__slug__startswith="00", source__name='SWISSPROT').prefetch_related('family').order_by('family__slug')
    print("time 1")

    coverage = OrderedDict()

    temp = OrderedDict([
                        ('name',''),
                        ('interactions', 0),
                        ('receptor_i', 0) ,
                        ('mutations' , 0),
                        ('receptor_m', 0),
                        ('mutations_an' , 0),
                        ('receptor_m_an', 0),
                        ('receptor_t',0),
                        ('children', OrderedDict()) ,
                        ('fraction_i',0),
                        ('fraction_m',0),
                        ('fraction_m_an',0)
                        ])

    for p in class_proteins:
        #print(p,p.family.slug)
        fid = p.family.slug.split("_")
        if fid[0] not in coverage:
            coverage[fid[0]] = copy.deepcopy(temp)
            coverage[fid[0]]['name'] = lookup[fid[0]]
        if fid[1] not in coverage[fid[0]]['children']:
            coverage[fid[0]]['children'][fid[1]] = copy.deepcopy(temp)
            coverage[fid[0]]['children'][fid[1]]['name'] = lookup[fid[0]+"_"+fid[1]]
        if fid[2] not in coverage[fid[0]]['children'][fid[1]]['children']:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]] = copy.deepcopy(temp)
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['name'] = lookup[fid[0]+"_"+fid[1]+"_"+fid[2]][:28]
        if fid[3] not in coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children']:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]] = copy.deepcopy(temp)
            coverage[fid[0]]['receptor_t'] += 1
            coverage[fid[0]]['children'][fid[1]]['receptor_t'] += 1
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['receptor_t'] += 1
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['name'] = p.entry_name.split("_")[0] #[:10]
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['receptor_t'] = 1

    coverage3 = copy.deepcopy(coverage)
    print("time 2")


    if 1==1:
        class_interactions = ResidueFragmentInteraction.objects.filter(structure_ligand_pair__annotated=True).prefetch_related(
            'rotamer__residue__display_generic_number','interaction_type',
            'structure_ligand_pair__structure__protein_conformation__protein__parent__family',
            'structure_ligand_pair__ligand__properities',
            )

        class_mutations = MutationExperiment.objects.all().prefetch_related('protein__family','protein__parent__family','exp_func','residue__display_generic_number','mutation','refs__web_link', 'exp_qual','ligand').order_by('foldchange','exp_qual')

        generic = {}

        score_copy = {'score': {'a':0,'i':0,'i_weight':0,'m':0,'m_weight':0,'s':0,'s_weight':0} , 'interaction' : {},'mutation': {}}

        dump_interactions = []
        for i in class_interactions:
            fid = i.structure_ligand_pair.structure.protein_conformation.protein.parent.family.slug.split("_")
            interaction_type = i.interaction_type.slug
            interaction_type_class = i.interaction_type.type
            if i.rotamer.residue.display_generic_number:
                dgn = i.rotamer.residue.display_generic_number.label
            else:
                dgn = 'N/A'
            #dump_interactions.append([dgn,i.rotamer.residue.sequence_number, i.rotamer.residue.amino_acid,i.structure_ligand_pair.structure.pdb_code.index,interaction_type,interaction_type_class,i.interaction_type.name,i.structure_ligand_pair.structure.protein_conformation.protein.parent.entry_name])
            if interaction_type=='polar_backbone':
                continue
            if interaction_type=='acc':
                continue
            coverage[fid[0]]['interactions'] += 1
            coverage[fid[0]]['children'][fid[1]]['interactions'] += 1
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['interactions'] += 1
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['interactions'] += 1

            if coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['interactions']==1: #if first time receptor gets a point
                coverage[fid[0]]['receptor_i'] += 1
                coverage[fid[0]]['children'][fid[1]]['receptor_i'] += 1
                coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['receptor_i'] += 1
                coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['receptor_i'] = 1

                coverage[fid[0]]['fraction_i'] = coverage[fid[0]]['receptor_i']/coverage[fid[0]]['receptor_t']
                coverage[fid[0]]['children'][fid[1]]['fraction_i'] = coverage[fid[0]]['children'][fid[1]]['receptor_i']/coverage[fid[0]]['children'][fid[1]]['receptor_t']
                coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['fraction_i'] = coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['receptor_i']/coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['receptor_t']

        # FOR ALEX
        # print(dump_interactions)
        # with open('interactions_dump.csv', 'w', newline='') as myfile:
        #     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        #     wr.writerows(dump_interactions)
        total_r = 0
        total_r_un = 0 #unannotated
        total_m = 0 #annotated
        total_m_un = 0 #unannotated
        for m in class_mutations:
            # break
            fid = m.protein.family.slug.split("_")
            total_m_un += 1
            coverage[fid[0]]['mutations'] += 1
            coverage[fid[0]]['children'][fid[1]]['mutations'] += 1
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['mutations'] += 1
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['mutations'] += 1

            if coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['mutations']==1: #if first time receptor gets a point
                total_r_un += 1
                coverage[fid[0]]['receptor_m'] += 1
                coverage[fid[0]]['children'][fid[1]]['receptor_m'] += 1
                coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['receptor_m'] += 1

                coverage[fid[0]]['fraction_m'] = coverage[fid[0]]['receptor_m']/coverage[fid[0]]['receptor_t']
                coverage[fid[0]]['children'][fid[1]]['fraction_m'] = coverage[fid[0]]['children'][fid[1]]['receptor_m']/coverage[fid[0]]['children'][fid[1]]['receptor_t']
                coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['fraction_m'] = coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['receptor_m']/coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['receptor_t']

            if m.exp_func is not None or m.foldchange!=0 or m.exp_qual is not None or m.ligand is not None: #if exp with data
                total_m += 1
                coverage[fid[0]]['mutations_an'] += 1
                coverage[fid[0]]['children'][fid[1]]['mutations_an'] += 1
                coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['mutations_an'] += 1
                coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['mutations_an'] += 1

                if coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['mutations_an']==1: #if first time receptor gets a point
                    total_r += 1
                    coverage[fid[0]]['receptor_m_an'] += 1
                    coverage[fid[0]]['children'][fid[1]]['receptor_m_an'] += 1
                    coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['receptor_m_an'] += 1
                    coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['receptor_m_an'] += 1

                    coverage[fid[0]]['fraction_m_an'] = coverage[fid[0]]['receptor_m_an']/coverage[fid[0]]['receptor_t']
                    coverage[fid[0]]['children'][fid[1]]['fraction_m_an'] = coverage[fid[0]]['children'][fid[1]]['receptor_m_an']/coverage[fid[0]]['children'][fid[1]]['receptor_t']
                    coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['fraction_m_an'] = coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['receptor_m_an']/coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['receptor_t']

        print("Total R",total_r,"Total M",total_m," <-- annotated || unannotated -->","Total R",total_r_un,"Total M",total_m_un)
        context['totals'] = {'total_r':total_r,'total_r_un':total_r_un, 'total_m':total_m, 'total_m_un':total_m_un}
        generic = OrderedDict(sorted(generic.items(), key=lambda x: x[1]['score']['s_weight'], reverse=True))

    CSS_COLOR_NAMES = ["AliceBlue","AntiqueWhite","Aqua","Aquamarine","Azure","Beige","Bisque","Black","BlanchedAlmond","Blue","BlueViolet","Brown","BurlyWood","CadetBlue","Chartreuse","Chocolate","Coral","CornflowerBlue","Cornsilk","Crimson","Cyan","DarkBlue","DarkCyan","DarkGoldenRod","DarkGray","DarkGrey","DarkGreen","DarkKhaki","DarkMagenta","DarkOliveGreen","Darkorange","DarkOrchid","DarkRed","DarkSalmon","DarkSeaGreen","DarkSlateBlue","DarkSlateGray","DarkSlateGrey","DarkTurquoise","DarkViolet","DeepPink","DeepSkyBlue","DimGray","DimGrey","DodgerBlue","FireBrick","FloralWhite","ForestGreen","Fuchsia","Gainsboro","GhostWhite","Gold","GoldenRod","Gray","Grey","Green","GreenYellow","HoneyDew","HotPink","IndianRed","Indigo","Ivory","Khaki","Lavender","LavenderBlush","LawnGreen","LemonChiffon","LightBlue","LightCoral","LightCyan","LightGoldenRodYellow","LightGray","LightGrey","LightGreen","LightPink","LightSalmon","LightSeaGreen","LightSkyBlue","LightSlateGray","LightSlateGrey","LightSteelBlue","LightYellow","Lime","LimeGreen","Linen","Magenta","Maroon","MediumAquaMarine","MediumBlue","MediumOrchid","MediumPurple","MediumSeaGreen","MediumSlateBlue","MediumSpringGreen","MediumTurquoise","MediumVioletRed","MidnightBlue","MintCream","MistyRose","Moccasin","NavajoWhite","Navy","OldLace","Olive","OliveDrab","Orange","OrangeRed","Orchid","PaleGoldenRod","PaleGreen","PaleTurquoise","PaleVioletRed","PapayaWhip","PeachPuff","Peru","Pink","Plum","PowderBlue","Purple","Red","RosyBrown","RoyalBlue","SaddleBrown","Salmon","SandyBrown","SeaGreen","SeaShell","Sienna","Silver","SkyBlue","SlateBlue","SlateGray","SlateGrey","Snow","SpringGreen","SteelBlue","Tan","Teal","Thistle","Tomato","Turquoise","Violet","Wheat","White","WhiteSmoke","Yellow","YellowGreen"];

#LightBlue
    CSS_COLOR_NAMES = ["SteelBlue","SlateBlue","LightCoral","Orange","LightGreen","LightGray","PeachPuff","PaleGoldenRod"]

    print("time 3")
    coverage2 = copy.deepcopy(coverage)

    print("time 4")
    tree = OrderedDict({'name':'GPCRs','children':[]})
    i = 0
    n = 0
    for c,c_v in coverage.items():
        c_v['name'] = c_v['name'].split("(")[0]
        if c_v['name'].strip() in ['Other GPCRs','Class T (Taste 2)','Class B2']:
            # i += 1
            continue
            # pass
        children = []
        for lt,lt_v in c_v['children'].items():
            if lt_v['name'].strip() == 'Orphan' and c_v['name'].strip()=="Class A":
                # $pass
                continue
            children_rf = []
            for rf,rf_v in lt_v['children'].items():
                rf_v['name'] = rf_v['name'].split("<")[0]
                if rf_v['name'].strip() == 'Class T (Taste 2)':
                    continue
                children_r = []
                for r,r_v in rf_v['children'].items():
                    r_v['color'] = CSS_COLOR_NAMES[i]
                    r_v['sort'] = n
                    children_r.append(r_v)
                    n += 1
                rf_v['children'] = children_r
                rf_v['sort'] = n
                rf_v['color'] = CSS_COLOR_NAMES[i]
                children_rf.append(rf_v)
            lt_v['children'] = children_rf
            lt_v['sort'] = n
            lt_v['color'] = CSS_COLOR_NAMES[i]
            children.append(lt_v)
        c_v['children'] = children
        c_v['sort'] = n
        c_v['color'] = CSS_COLOR_NAMES[i]
        tree['children'].append(c_v)
        #tree = c_v
        #break
        i += 1

    print("time 5")
    tree2 = OrderedDict({'name':'GPCRs','children':[]})
    i = 0
    n = 0
    for c,c_v in coverage2.items():
        children = []
        for lt,lt_v in c_v['children'].items():
            # if lt_v['name'].strip() == 'Orphan':
            #     continue
            children_rf = []
            for rf,rf_v in lt_v['children'].items():
                children_r = []
                for r,r_v in rf_v['children'].items():
                    r_v['color'] = CSS_COLOR_NAMES[i]
                    r_v['sort'] = n
                    children_r.append(r_v)
                    n += 1
                rf_v['children'] = [{'name':'', 'color':CSS_COLOR_NAMES[i], 'children': children_r }]
                rf_v['sort'] = n
                rf_v['color'] = CSS_COLOR_NAMES[i]
                children_rf.append(rf_v)
            lt_v['children'] = [{'name':'', 'color':CSS_COLOR_NAMES[i], 'children': children_rf }]
            lt_v['sort'] = n
            lt_v['color'] = CSS_COLOR_NAMES[i]
            children.append(lt_v)
        c_v['children'] = [{'name':'', 'color':CSS_COLOR_NAMES[i], 'children': children }]
        c_v['sort'] = n
        c_v['color'] = CSS_COLOR_NAMES[i]
        tree2['children'].append(c_v)
        i += 1
       # break
    #print(json.dumps(tree))
    print("time 6")
    context['coverage'] = coverage3 #coverage
    context['tree'] = json.dumps(tree)
    context['tree2'] = json.dumps(tree2)
    print("time 7")
    # return render(request, 'mutation/coverage.html', context)
    return render(request, 'mutation/statistics.html', context)


def pocket(request):

    context = {}

    gpcr_class = '006' #class a 1 , c 4, f 6

    class_interactions = ResidueFragmentInteraction.objects.filter(
        structure_ligand_pair__structure__protein_conformation__protein__family__slug__startswith=gpcr_class, structure_ligand_pair__annotated=True).prefetch_related(
        'rotamer__residue__display_generic_number','interaction_type',
        'structure_ligand_pair__structure__protein_conformation__protein__parent',
        'structure_ligand_pair__ligand__properities')

    class_mutations = MutationExperiment.objects.filter(
        protein__family__slug__startswith=gpcr_class).prefetch_related('protein','residue__display_generic_number','mutation','refs__web_link', 'exp_qual','ligand').order_by('foldchange','exp_qual')

    generic = {}

    score_copy = {'score': {'a':0,'i':0,'i_weight':0,'m':0,'m_weight':0,'s':0,'s_weight':0, 'pdbs' : [], 'pos' : [] } , 'interaction' : {},'mutation': {}}

    for i in class_interactions:
        ligand = i.structure_ligand_pair.ligand.name
        receptor = i.structure_ligand_pair.structure.protein_conformation.protein.parent.entry_name
        receptor = receptor.split("_")[0]
        interaction_type = i.interaction_type.slug
        interaction_type_class = i.interaction_type.type
        pdb = i.structure_ligand_pair.structure.pdb_code.index
        pos = i.rotamer.residue.sequence_number
        aa = i.rotamer.residue.amino_acid
        if i.rotamer.residue.display_generic_number:
            gn = i.rotamer.residue.display_generic_number.label
            gn2 = i.rotamer.residue.generic_number.label
            gn = gn +" - "+gn2
        else:
            continue
        if gn not in generic:
            generic[gn] = copy.deepcopy(score_copy)
        if receptor not in generic[gn]['interaction']:
            generic[gn]['interaction'][receptor] = {}

        if ligand not in generic[gn]['interaction'][receptor]:
            generic[gn]['interaction'][receptor][ligand] = {}
        if interaction_type=='acc':
            if 'a' not in generic[gn]['interaction'][receptor][ligand]:
                generic[gn]['score']['a'] += 1
                generic[gn]['score']['s'] += 1
                generic[gn]['score']['pdbs'].append([pdb,pos,aa])
                generic[gn]['interaction'][receptor][ligand]['a'] = 1
        elif interaction_type!='acc':
            if 'i' not in generic[gn]['interaction'][receptor][ligand]:
                generic[gn]['interaction'][receptor][ligand]['i'] = 1
                generic[gn]['score']['i'] += 1
                generic[gn]['score']['s'] += 1
                generic[gn]['score']['s_weight'] += 1
    for m in class_mutations:
        if not m.ligand: #ignore non ligand
            continue
        receptor = m.protein.entry_name
        receptor = receptor.split("_")[0]
        ligand = m.ligand.name
        if m.residue.display_generic_number:
            gn = m.residue.display_generic_number.label
        else:
            continue
        if gn not in generic:
            generic[gn] = copy.deepcopy(score_copy)
        if receptor not in generic[gn]['mutation']:
            generic[gn]['mutation'][receptor] = {}
        if ligand not in generic[gn]['mutation'][receptor]:
            generic[gn]['score']['m'] += 1
            generic[gn]['score']['s'] += 1

            if m.foldchange>5:
                generic[gn]['score']['m_weight'] += 1
                generic[gn]['score']['s_weight'] += 1
                generic[gn]['mutation'][receptor][ligand] = {}
            elif m.exp_qual:
                if m.exp_qual.qual=='Abolish' or m.exp_qual.qual.find('abolish')!=-1 or m.exp_qual.qual.find('Abolish')!=-1:
                    generic[gn]['score']['m_weight'] += 1
                    generic[gn]['score']['s_weight'] += 1
                    generic[gn]['mutation'][receptor][ligand] = {}

    generic = OrderedDict(sorted(generic.items(), key=lambda x: x[1]['score']['a'], reverse=True))
    #print(generic)
    context['gn'] = generic
    return render(request, 'mutation/pocket.html', context)

def showcalculation(request):
    if request.method == 'POST':
        form = PDBform(request.POST, request.FILES)
        context = calculate(request)

    else:
        simple_selection = request.session.get('selection', False)
        proteins = []
        for target in simple_selection.targets:
            if target.type == 'protein':
                proteins.append(target.item)

        context = {}
        context['proteins'] = proteins

    print(context['proteins'])

    protein_ids = []
    family_ids = []
    parent_ids = []
    class_ids = []

    for p in context['proteins']:
        protein_ids.append(p.family) #first level is receptor across speciest, parent is then "family"
        family_ids.append(p.family.parent) #first level is receptor across speciest, parent is then "family"
        parent_ids.append(p.family.parent.parent) #see above, go to parent.parent to get true parent.
        class_ids.append(p.residue_numbering_scheme)
        family = p.family
        while family.parent.parent is not None:
            family = family.parent


    family_level_ids = protein_ids[0].slug.split("_")

    if len(context['proteins'])>1:
        return HttpResponse("Only pick one protein")


    residues = Residue.objects.filter(protein_conformation__protein=context['proteins'][0]).prefetch_related('protein_segment','display_generic_number','generic_number')

    HelixBox = DrawHelixBox(
                residues, context['proteins'][0], str(p), nobuttons=1)
    SnakePlot = DrawSnakePlot(
                residues, context['proteins'][0], str(p), nobuttons=1)

    lookup = {}
    lookup_with_pos = {}
    lookup_pos = {}
    for r in residues:
        if r.generic_number:
            lookup[r.generic_number.label] = r.amino_acid
            lookup_with_pos[r.generic_number.label] = r.amino_acid+str(r.sequence_number)
            lookup_pos[r.generic_number.label] = str(r.sequence_number)

    gpcr_class = family
    class_interactions = ResidueFragmentInteraction.objects.filter(
        structure_ligand_pair__structure__protein_conformation__protein__family__slug__startswith=gpcr_class.slug, structure_ligand_pair__annotated=True).prefetch_related(
        'rotamer__residue__generic_number','interaction_type',
        'structure_ligand_pair__structure__protein_conformation__protein__parent__family',
        'structure_ligand_pair__ligand__properities')

    class_mutations = MutationExperiment.objects.filter(
        protein__family__slug__startswith=gpcr_class.slug).prefetch_related('protein__family','residue__generic_number','mutation','refs__web_link', 'exp_qual','ligand').order_by('foldchange','exp_qual')

    generic = {}

    score_copy = {'score': {'a':0,'i':0,'i_aa':0,'i_weight':0,'m':0,'m_aa':0,'m_weight':0,'s':0,'s_weight':0,'s_aa':0}, 'homology' : {0:0,1:0,2:0,3:0} , 'interaction_aa' : {},'mutation_aa': {}, 'interaction' : {},'mutation': {}}

    for i in class_interactions:
        #continue
        ligand = i.structure_ligand_pair.ligand.name
        receptor = i.structure_ligand_pair.structure.protein_conformation.protein.parent.entry_name
        receptor = receptor.split("_")[0]
        family_id = i.structure_ligand_pair.structure.protein_conformation.protein.parent.family.slug.split("_")
        interaction_type = i.interaction_type.slug
        interaction_type_class = i.interaction_type.type

        if family_level_ids[0:3]==family_id[0:3]:
            pass
            #continue

        if interaction_type=='polar_backbone':
            continue
        if i.rotamer.residue.generic_number:
            gn = i.rotamer.residue.generic_number.label
        else:
            continue
        if gn not in generic:
            generic[gn] = copy.deepcopy(score_copy)
        if receptor not in generic[gn]['interaction']:
            generic[gn]['interaction'][receptor] = {}
        if ligand not in generic[gn]['interaction'][receptor]:
            generic[gn]['interaction'][receptor][ligand] = {}
            generic[gn]['score']['i'] += 1
            generic[gn]['score']['s'] += 1
            generic[gn]['score']['s_weight'] += 1
            #print(gn,receptor,ligand)
        if gn in lookup and interaction_type=='acc':
            generic[gn]['score']['a'] += 1
        if gn in lookup and interaction_type!='acc':
            if lookup[gn] == i.rotamer.residue.amino_acid:
                if receptor not in generic[gn]['interaction_aa']:
                    generic[gn]['interaction_aa'][receptor] = {}
                if ligand not in generic[gn]['interaction_aa'][receptor]:
                    generic[gn]['interaction_aa'][receptor][ligand] = {}
                    generic[gn]['score']['i_aa'] += 1
                    generic[gn]['score']['s_aa'] += 1
                    if family_id==family_level_ids:
                        generic[gn]['homology'][0] += 1
                    elif family_id[0:3]==family_level_ids[0:3]:
                        generic[gn]['homology'][1] += 1
                    elif family_id[0:2]==family_level_ids[0:2]:
                        generic[gn]['homology'][2] += 1
                    elif family_id[0:1]==family_level_ids[0:1]:
                        generic[gn]['homology'][3] += 1
                    else:
                        print("error",family_id,family_level_ids)

    for m in class_mutations:
        #continue
        if not m.ligand: #ignore non ligand
            continue
        receptor = m.protein.entry_name
        receptor = receptor.split("_")[0]
        family_id = m.protein.family.slug.split("_")
        ligand = m.ligand.name

        if family_level_ids[0:3]==family_id[0:3]:
            pass
            #continue

        if m.residue.generic_number:
            gn = m.residue.generic_number.label
        else:
            continue
        if gn not in generic:
            generic[gn] = copy.deepcopy(score_copy)
        if receptor not in generic[gn]['mutation']:
            generic[gn]['mutation'][receptor] = {}
        if ligand not in generic[gn]['mutation'][receptor]:
            generic[gn]['score']['m'] += 1
            generic[gn]['score']['s'] += 1

            if m.foldchange>5:
                generic[gn]['score']['m_weight'] += 1
                generic[gn]['score']['s_weight'] += 1

                generic[gn]['mutation'][receptor][ligand] = {}
                #print(gn,receptor,ligand)
            elif m.exp_qual:
                if m.exp_qual.qual=='Abolish' or m.exp_qual.qual.find('abolish')!=-1 or m.exp_qual.qual.find('Abolish')!=-1:
                   generic[gn]['score']['m_weight'] += 1
                   generic[gn]['score']['s_weight'] += 1
                   m.foldchange = 6

                   generic[gn]['mutation'][receptor][ligand] = {}
                   #print(gn,receptor,ligand)
        if gn in lookup:
            if lookup[gn] == m.residue.amino_acid and m.foldchange>5:
                if receptor not in generic[gn]['mutation_aa']:
                    generic[gn]['mutation_aa'][receptor] = {}
                if ligand not in generic[gn]['mutation_aa'][receptor]:
                    generic[gn]['mutation_aa'][receptor][ligand] = {}
                    generic[gn]['score']['m_aa'] += 1
                    generic[gn]['score']['s_aa'] += 1
                    if family_id==family_level_ids:
                        generic[gn]['homology'][0] += 1
                    elif family_id[0:3]==family_level_ids[0:3]:
                        generic[gn]['homology'][1] += 1
                    elif family_id[0:2]==family_level_ids[0:2]:
                        generic[gn]['homology'][2] += 1
                    elif family_id[0:1]==family_level_ids[0:1]:
                        generic[gn]['homology'][3] += 1
                    else:
                        print("error",family_id,family_level_ids)

    position_scores = OrderedDict(sorted(generic.items(), key=lambda x: (x[1]['homology'][0],x[1]['homology'][1],x[1]['homology'][2],x[1]['homology'][3]), reverse=True))
    pocket_scores = {}
    score = 21
    prev_homology = ''
    for g,s in position_scores.items():
        if s['homology']=={0: 0, 1: 0, 2: 0, 3: 0}:
            score = 0
        if s['homology']!=prev_homology: #give same score if same homology numbers
            score = max(0,score-1) #give less and less score value
            prev_homology = s['homology']
        pocket_scores[g] = score


    #NEW CLASS METHOD, then select closest
    class_interactions = ResidueFragmentInteraction.objects.filter(
        structure_ligand_pair__structure__protein_conformation__protein__family__slug__startswith=family.slug, structure_ligand_pair__annotated=True).prefetch_related('rotamer__residue__generic_number','interaction_type','structure_ligand_pair__structure__protein_conformation__protein__parent__family','structure_ligand_pair__structure__pdb_code','structure_ligand_pair__ligand__properities')

    class_mutations = MutationExperiment.objects.filter(
        protein__family__slug__startswith=family.slug).prefetch_related('protein__family','residue__generic_number','mutation','refs__web_link', 'exp_qual','ligand__properities').order_by('-foldchange','exp_qual')

    # class_proteins = Protein.objects.filter(family__slug__startswith=family.slug, source__name='SWISSPROT',species__common_name='Human').all()
    # family_proteins = Protein.objects.filter(family__parent__in=family_ids, source__name='SWISSPROT').all() #,species__common_name='Human'
    # ligand_proteins = Protein.objects.filter(family__parent__parent__in =parent_ids, source__name='SWISSPROT',species__common_name='Human').all()

    class_p = []
    for i in class_interactions:
        p = i.structure_ligand_pair.structure.protein_conformation.protein.parent
        if p not in class_p:
            class_p.append(p)
    for m in class_mutations:
        entry_name = m.protein.entry_name
        if (int(m.foldchange)!=0 or m.exp_qual):
            p = m.protein
            if p not in class_p:
                class_p.append(p)
    # for p in class_proteins: #remove these as they are never used.
    #     if p not in class_p:
    #         class_p.append(p)
    print("Proteins in alignment",len(class_p),'proteins with contributing data')

    # print(list(settings.REFERENCE_POSITIONS.keys()))
    segments = ProteinSegment.objects.filter(slug__in=list(settings.REFERENCE_POSITIONS.keys()))
    segments = ProteinSegment.objects.filter(category='helix') #only use helix for faster rendering.

    # a = Alignment()
    # a.load_reference_protein(context['proteins'][0])
    # a.load_proteins(family_proteins)
    # #segments = ProteinSegment.objects.filter(category='helix')
    # a.load_segments(segments)
    # a.build_alignment()
    # a.calculate_similarity()
    # a.calculate_statistics()
    # family_generic_aa_count = a.calculate_aa_count_per_generic_number()

    # print('alignment 1')


    #Consider caching result! Would be by protein since it compares protein to whole class.
    json_generic = str(context['proteins'][0])+'_generic.json'
    json_alternative = str(context['proteins'][0])+'_alternative.json'
    json_similarity_list = str(context['proteins'][0])+'_similarity_list.json'


    generic_aa_count = cache.get(json_generic)
    alternative_aa = cache.get(json_alternative)
    similarity_list = cache.get(json_similarity_list)

    # if os.path.isfile(json_generic) and os.path.isfile(json_alternative) and os.path.isfile(json_similarity_list) and 1==1: #DISABLE THIS AS IT MISFIRED WHEN NEW DATA
    #     generic_aa_count = json.load(open(json_generic, 'r'))
    #     alternative_aa = json.load(open(json_alternative, 'r'))
    #     similarity_list = json.load(open(json_similarity_list, 'r'))
    #     print('alignment (class) using cache')
    if generic_aa_count==None or alternative_aa==None or similarity_list==None:
        a = Alignment()
        a.load_reference_protein(context['proteins'][0])
        a.load_proteins(class_p)
        segments = ProteinSegment.objects.filter(category='helix') #only use helix for faster rendering.
        #segments = ProteinSegment.objects.all()
        a.load_segments(segments) #get all segments to make correct diagrams

        # build the alignment data matrix
        a.build_alignment()

        # calculate similarity
        a.calculate_similarity()

        a.calculate_statistics()
        generic_aa_count = a.calculate_aa_count_per_generic_number()
        alternative_aa = a.aa_count_with_protein

        print('alignment (class)')

        similarity_list = {}
        for p in a.proteins:
            similarity_list[p.protein.entry_name] = [int(p.identity),int(p.similarity),p.similarity_score]
            if (p.protein.entry_name==context['proteins'][0].entry_name):
                similarity_list[p.protein.entry_name] = [int(100),int(100),1000]

        # json.dump(generic_aa_count, open(json_generic, 'w'))
        # json.dump(alternative_aa, open(json_alternative, 'w'))
        # json.dump(similarity_list, open(json_similarity_list, 'w'))

        cache.set(json_generic,generic_aa_count)
        cache.set(json_alternative,alternative_aa)
        cache.set(json_similarity_list,similarity_list)
    else:
        print('alignment (class) using cache')

    results = {}
    mutant_lookup = {}
    distinct_species = {}
    distinct_ligands = {}
    distinct_big_decrease = {}
    level = 0

    empty_result = {'homology':{0:0, 1:0, 2:0, 3:0},'interaction': {0:[], 1:[], 2:[]}, 'mutant': {0:[], 1:[], 2:[] }, 'closestinteraction' : { 'similarity' : 0},
                'interactions': { },
                'bestmutation':{ 'species' : '', 'similarity' : 0, 'foldchange' : 0, 'qual' : '', 'allmut' : [], 'counts':OrderedDict(), 'counts_close':{}, 'bigdecrease_distinct':0,'bigdecrease':0,'decrease':0,'bigincrease':0,'nonsignificant':0,'nodata':0}}


    for i in class_interactions:
        #continue
        interaction_type = i.interaction_type.slug
        interaction_type_class = i.interaction_type.type

        if interaction_type_class=='hidden':
            interaction_type_class = 'accessible'

        if interaction_type=='polar_backbone':
            continue

        if interaction_type_class=='accessible':
            continue

        if i.rotamer.residue.generic_number:
            generic = i.rotamer.residue.generic_number.label

            entry_name = i.structure_ligand_pair.structure.protein_conformation.protein.parent.entry_name
            family_id = i.structure_ligand_pair.structure.protein_conformation.protein.parent.family.slug.split("_")
            pdbcode = i.structure_ligand_pair.structure.pdb_code.index
            ligand = i.structure_ligand_pair.ligand.name
            if i.structure_ligand_pair.ligand.properities:
                smiles = i.structure_ligand_pair.ligand.properities.smiles
            else:
                smiles = ""

            if family_level_ids[0:3]==family_id[0:3]:
                pass
                #continue

            if generic in lookup:
                if generic not in distinct_species:
                    distinct_species[generic] = []
                    distinct_ligands[generic] = []
                    distinct_big_decrease[generic] = []

                if lookup[generic] == i.rotamer.residue.amino_acid:


                    if entry_name.split("_")[0] not in distinct_species[generic]:
                        distinct_species[generic].append(entry_name.split("_")[0])

                    if ligand not in distinct_ligands[generic]:
                        distinct_ligands[generic].append(ligand)

                    if generic  in results:

                        if interaction_type_class not in results[generic]['interactions']:
                            results[generic]['interactions'][interaction_type_class] = []

                        results[generic]['interactions'][interaction_type_class].append({ 'species' : entry_name, 'similarity' : similarity_list[entry_name][1], 'pdbcode' : pdbcode, 'ligand' : ligand,'smiles':smiles})



                        if family_id==family_level_ids:
                            results[generic]['homology'][0] += 1
                        elif family_id[0:3]==family_level_ids[0:3]:
                            results[generic]['homology'][1] += 1
                        elif family_id[0:2]==family_level_ids[0:2]:
                            results[generic]['homology'][2] += 1
                        elif family_id[0:1]==family_level_ids[0:1]:
                            results[generic]['homology'][3] += 1
                        else:
                            print("error",family_id,family_level_ids)

                        if similarity_list[entry_name][1]>results[generic]['closestinteraction']['similarity']:
                            results[generic]['closestinteraction']['species'] = entry_name
                            results[generic]['closestinteraction']['similarity'] = similarity_list[entry_name][1]
                            results[generic]['closestinteraction']['type'] = interaction_type
                            results[generic]['closestinteraction']['type_class'] = interaction_type_class
                            results[generic]['closestinteraction']['pdbcode'] = pdbcode
                        elif similarity_list[entry_name][1]==results[generic]['closestinteraction']['similarity'] and results[generic]['closestinteraction']['type_class']=='hydrophobic':
                            results[generic]['closestinteraction']['species'] = entry_name
                            results[generic]['closestinteraction']['similarity'] = similarity_list[entry_name][1]
                            results[generic]['closestinteraction']['type'] = interaction_type
                            results[generic]['closestinteraction']['type_class'] = interaction_type_class
                            results[generic]['closestinteraction']['pdbcode'] = pdbcode
                        #if similarity_list[entry_name][1]>results[generic]['interactions'][interaction_type_class]['similarity']:
                            #results[generic]['interactions'][interaction_type_class]['similarity'] = similarity_list[entry_name][1]

                    else:
                        results[generic] = copy.deepcopy(empty_result)
                        mutant_lookup[generic] = []
                        if generic in lookup:
                            if (lookup[generic] == i.rotamer.residue.amino_acid): #only for same aa (FIXME substitution)
                                if interaction_type_class=='accessible':
                                        continue

                                if family_id==family_level_ids:
                                    results[generic]['homology'][0] += 1
                                elif family_id[0:3]==family_level_ids[0:3]:
                                    results[generic]['homology'][1] += 1
                                elif family_id[0:2]==family_level_ids[0:2]:
                                    results[generic]['homology'][2] += 1
                                elif family_id[0:1]==family_level_ids[0:1]:
                                    results[generic]['homology'][3] += 1
                                else:
                                    print("error",family_id,family_level_ids)
                                results[generic]['interactions'][interaction_type_class] = [{ 'species' : entry_name, 'similarity' : similarity_list[entry_name][1], 'pdbcode' : pdbcode, 'ligand' : ligand,'smiles': smiles}]

                                results[generic]['closestinteraction']['species'] = entry_name
                                results[generic]['closestinteraction']['similarity'] = similarity_list[entry_name][1]
                                results[generic]['closestinteraction']['type'] = interaction_type
                                results[generic]['closestinteraction']['type_class'] = interaction_type_class
                                results[generic]['closestinteraction']['pdbcode'] = pdbcode
        else:
            pass
            #print('no generic number',interaction_type)

    print('parsed interaction data',len(results))

    for m in class_mutations:
        #continue
        if m.residue.generic_number:
            generic = m.residue.generic_number.label
            entry_name = m.protein.entry_name
            family_id = m.protein.family.slug.split("_")
            if m.ligand:
                ligand = m.ligand.name
                if m.ligand.properities:
                    smiles = m.ligand.properities.smiles
                else:
                    smiles = ""
            else:
                ligand = "N/A"
                smiles = ""

            if family_level_ids[0:3]==family_id[0:3]:
                pass
                #continue

            if m.exp_qual:
                qual = m.exp_qual.qual +" "+m.exp_qual.prop
            else:
                qual = ''
            #only select positions where interaction data is present and mutant has real data
            if generic in lookup: # or similarity_list[entry_name][1]>60 (or is closely related.)
                if generic not in distinct_species:
                    distinct_species[generic] = []
                    distinct_ligands[generic] = []
                    distinct_big_decrease[generic] = []
                #skip data that is far away / disable this for now
                #if similarity_list[entry_name][1]<50:
                #   continue


                #Only look at same residues (Expand with substitution possibilities) FIXME
                if lookup[generic] == m.residue.amino_acid:

                    #if row is allowed due to mutant data, create entry if it isnt there.
                    if generic not in results and (int(m.foldchange)!=0 or qual!=''):
                        results[generic] = copy.deepcopy(empty_result)
                    elif not (int(m.foldchange)!=0 or qual!=''): #skip no data on non-interesting positions / potentially miss a bit of data if datamutant comes later.. risk! FIXME
                    #should be fixed with order by
                    # generic not in results and
                        continue

                    if m.foldchange>20:
                        results[generic]['bestmutation']['bigdecrease'] += 1
                        if entry_name.split("_")[0] not in distinct_big_decrease[generic]:
                            results[generic]['bestmutation']['bigdecrease_distinct'] += 1
                            distinct_big_decrease[generic].append(entry_name.split("_")[0])
                    elif m.foldchange>5:
                        results[generic]['bestmutation']['decrease'] += 1
                    elif m.foldchange<-5:
                        results[generic]['bestmutation']['bigincrease'] += 1
                    elif (m.foldchange<5 or m.foldchange>-5) and m.foldchange!=0:
                        results[generic]['bestmutation']['nonsignificant'] += 1
                    else:
                        if m.exp_qual:
                            #print( m.exp_qual.qual.find('abolish'))
                            #print(m.exp_qual.qual)
                            if m.exp_qual.qual=='Abolish' or m.exp_qual.qual.find('abolish')!=-1 or m.exp_qual.qual.find('Abolish')!=-1:
                                results[generic]['bestmutation']['bigdecrease'] += 1
                                if entry_name.split("_")[0] not in distinct_big_decrease[generic]:
                                    results[generic]['bestmutation']['bigdecrease_distinct'] += 1
                                    distinct_big_decrease[generic].append(entry_name.split("_")[0])

                                m.foldchange = 20 #insert a 'fake' foldchange to make it count
                            elif m.exp_qual.qual=='Gain of':
                                results[generic]['bestmutation']['bigincrease'] += 1
                            elif m.exp_qual.qual=='Increase':
                                results[generic]['bestmutation']['nonsignificant'] += 1
                            elif m.exp_qual.qual=='Decrease':
                                results[generic]['bestmutation']['nonsignificant'] += 1
                            else:
                                results[generic]['bestmutation']['nonsignificant'] += 1 #non-abolish qual
                        else:
                            results[generic]['bestmutation']['nodata'] += 1

                    if m.foldchange>5:
                        if entry_name.split("_")[0] not in distinct_species[generic]:
                            distinct_species[generic].append(entry_name.split("_")[0])
                        if ligand not in distinct_ligands[generic]:
                            distinct_ligands[generic].append(ligand)


                    #If next is closer in similarity replace "closest"
                    if int(m.foldchange)!=0 or qual!='': #FIXME qual values need a corresponding foldchange value to outrank other values
                        if ((similarity_list[entry_name][1]>=results[generic]['bestmutation']['similarity'] and
                                m.foldchange>results[generic]['bestmutation']['foldchange']) and lookup[generic] == m.residue.amino_acid):
                            results[generic]['bestmutation']['species'] = entry_name
                            results[generic]['bestmutation']['similarity'] = similarity_list[entry_name][1]
                            results[generic]['bestmutation']['foldchange'] = m.foldchange
                            results[generic]['bestmutation']['qual'] = qual
                            results[generic]['bestmutation']['aa'] = m.mutation.amino_acid

                        results[generic]['bestmutation']['allmut'].append([entry_name,m.foldchange,qual,m.mutation.amino_acid,similarity_list[entry_name][1],ligand,smiles])

                    if int(m.foldchange>5):
                        if family_id==family_level_ids:
                            results[generic]['homology'][0] += 1
                        elif family_id[0:3]==family_level_ids[0:3]:
                            results[generic]['homology'][1] += 1
                        elif family_id[0:2]==family_level_ids[0:2]:
                            results[generic]['homology'][2] += 1
                        elif family_id[0:1]==family_level_ids[0:1]:
                            results[generic]['homology'][3] += 1
                        else:
                            print("error",family_id,family_level_ids)

                    if m.mutation.amino_acid in results[generic]['bestmutation']['counts']:
                        results[generic]['bestmutation']['counts'][m.mutation.amino_acid] += 1
                    else:
                        results[generic]['bestmutation']['counts'][m.mutation.amino_acid] = 1

                    if m.mutation.amino_acid in results[generic]['bestmutation']['counts_close'] and similarity_list[entry_name][1]>60:
                        results[generic]['bestmutation']['counts_close'][m.mutation.amino_acid] += 1
                    elif similarity_list[entry_name][1]>60:
                        results[generic]['bestmutation']['counts_close'][m.mutation.amino_acid] = 1


                mutant_lookup[generic] = []
        else:
            pass
            #print('no generic number')
    print('parsed mutant data',len(results))

    #Fetch defined subsitution matrix for mutant design tool
    matrix = definitions.DESIGN_SUBSTITUTION_MATRIX

    summary = {}
    summary_score = []
    for res,values in results.items():

        if position_scores[res]['score']['a']==0: #skip those who have no evidence of being in pocket
            #pass
            continue

        if res in lookup:
            summary[res] = {}
            summary[res]['aa'] = lookup[res]
            summary[res]['aa_pos'] = lookup_with_pos[res]
            summary[res]['pos'] = lookup_pos[res]
        else: #skip those that are not present in reference
            continue

        summary[res]['receptor_class'] =  family.slug
        summary[res]['closestinteraction'] = values['closestinteraction']
        summary[res]['interactions'] = values['interactions']
        summary[res]['bestmutation'] = values['bestmutation']
        summary[res]['homology'] = values['homology']
        summary[res]['homology_color'] = ''
        summary[res]['alternatives'] = []
        summary[res]['interest_score'] = 0
        summary[res]['score_text'] = ''
        summary[res]['pocket'] = 0
        summary[res]['freq'] = 0
        summary[res]['freq_r'] = []
        if res in distinct_species:
            summary[res]['freq_r'] = distinct_species[res]
        summary[res]['freq_l'] = []
        if res in distinct_ligands:
            summary[res]['freq_l'] = distinct_ligands[res]
        summary[res]['pocket_rank'] = 0
        summary[res]['pocket_sum'] = position_scores[res]['score']['s_aa']
        summary[res]['pocket_acc'] = position_scores[res]['score']['a']
        summary[res]['pocket_homology'] = position_scores[res]['homology']
        summary[res]['distinct_interaction'] = {}
        summary[res]['ranking_score'] = 0
        #reverse counts

        #homology and position freq rank
        if summary[res]['homology'][0]>0:
            summary[res]['ranking_score']=4000
            summary[res]['homology_color'] = '#1e6a00'
        elif summary[res]['homology'][1]>0:
            summary[res]['ranking_score']=3000
            summary[res]['homology_color'] = '#519d32'
        elif summary[res]['homology'][2]>0:
            summary[res]['ranking_score']=2000
            summary[res]['homology_color'] = '#92c27f'
        elif summary[res]['homology'][3]>0:
            summary[res]['ranking_score']=1000
            summary[res]['homology_color'] = '#d3e6cc'


        if 'counts' in summary[res]['bestmutation']:
            summary[res]['bestmutation']['counts'] = OrderedDict(sorted(summary[res]['bestmutation']['counts'].items(), key=lambda x: x[1], reverse=True))

        summary[res]['interest_score'] += len(summary[res]['bestmutation']['allmut'])/10 #add a small value for each mutation for position
        summary[res]['score_text'] += '# mutants: '+str(len(summary[res]['bestmutation']['allmut'])/10)+'<br>'

        if res in pocket_scores:
            summary[res]['pocket_rank'] = 21-pocket_scores[res]
            summary[res]['interest_score'] += pocket_scores[res]
            summary[res]['ranking_score'] += pocket_scores[res]
            summary[res]['score_text'] += '# freq of interactions/mutants within class: '+ str(pocket_scores[res])+'<br>'
            summary[res]['pocket'] = pocket_scores[res]

        #Find alternatives for position (useful for specificity investigation)
        temp = 0
        if res in alternative_aa:
            for aalist in alternative_aa[res].items():
                if aalist[0] !=summary[res]['aa']:
                    for p in aalist[1]:
                        if similarity_list[p][1]>60: #close enough to suggest
                             summary[res]['alternatives'].append([p,similarity_list[p][1],aalist[0]])
                             summary[res]['interest_score'] += 1 #add a small value for these hits
                             temp += 1
            if temp: #if alternatives found
                summary[res]['score_text'] += '# alternatives: '+str(temp)+'<br>'

        # if res in generic_aa_count:
        #     #removed ligand_generic_aa_count[res][lookup[res]] -- not used
        #     summary[res]['conservation'] = [family_generic_aa_count[res][lookup[res]],0,generic_aa_count[res][lookup[res]], round(family_generic_aa_count[res][lookup[res]] / generic_aa_count[res][lookup[res]],1), round(family_generic_aa_count[res][lookup[res]] - generic_aa_count[res][lookup[res]],1)]

        #     summary[res]['interest_score'] += summary[res]['conservation'][4]/10
        #     summary[res]['score_text'] += '# cons span: '+str(summary[res]['conservation'][4]/10)+'<br>'

        scores = {'hyd':0,'aromatic':0,'polar':0,'unknown':0} #dict to keep track of scores to select subs.

        secondary_interaction = ''
        if 'type_class' in summary[res]['closestinteraction']:
            temp = 0
            summary[res]['freq'] = 1
            if summary[res]['closestinteraction']['type_class']=='hydrophobic':
                summary[res]['interest_score'] += 10
                temp = 10
                if 'polar' in summary[res]['interactions']: #if secondary there is polar
                    summary[res]['score_text'] += '# secondary interaction: '+str(25)+'<br>'
                    summary[res]['interest_score'] += 25
                    secondary_interaction = 'polar'
                elif 'aromatic' in summary[res]['interactions']: #if secondary there is polar
                    summary[res]['score_text'] += '# secondary interaction: '+str(10)+'<br>'
                    summary[res]['interest_score'] += 10
                    secondary_interaction = 'aromatic'
            elif summary[res]['closestinteraction']['type_class']=='polar':
                summary[res]['interest_score'] += 50
                temp = 50
                if 'hydrophobic' in summary[res]['interactions']:
                    secondary_interaction = 'hydrophobic'
            elif summary[res]['closestinteraction']['type_class']=='aromatic':
                summary[res]['interest_score'] += 20
                temp = 20
                if 'hydrophobic' in summary[res]['interactions']:
                    secondary_interaction = 'hydrophobic'

            summary[res]['score_text'] += '# closest interaction: '+str(temp)+'<br>'

        if 'type_class' in summary[res]['closestinteraction']: #if there is a closest interaction, calculate supporting data
            summary[res]['distinct_interaction'][summary[res]['closestinteraction']['type_class']] = 0
            distinct_species_i = [summary[res]['closestinteraction']['species'].split("_")[0]]
            similarity_sum = 0
            for interaction in summary[res]['interactions'][summary[res]['closestinteraction']['type_class']]:
                if interaction['species'].split("_")[0] not in distinct_species and interaction['species']!=summary[res]['closestinteraction']['species']:
                    similarity_sum += interaction['similarity']
                    distinct_species_i.append(interaction['species'].split("_")[0])
                    summary[res]['distinct_interaction'][summary[res]['closestinteraction']['type_class']] += 1
            if len(distinct_species)>1: #if there was supporting
                summary[res]['interest_score'] += min((len(distinct_species_i)-1)*5,temp) #can't add more than the actual closest interaction
                summary[res]['score_text'] += '# supporting interactions: '+str(min((len(distinct_species)-1)*5,temp))+'<br>'
                summary[res]['freq'] += len(distinct_species_i)

            #summary[res]['freq_r'] += list(set(distinct_species_i) - set(summary[res]['freq_r']))

        summary[res]['freq_r'] = sorted(summary[res]['freq_r'])

        #print(len(summary[res]['freq_r']),summary[res]['freq_r'])

        summary[res]['freq'] = len(summary[res]['freq_r'])

        temp = 0
        if summary[res]['bestmutation']['foldchange']>50:
            summary[res]['interest_score'] += 70
            temp = 70
        elif summary[res]['bestmutation']['foldchange']>10:
            summary[res]['interest_score'] += 50
            temp = 50
        elif summary[res]['bestmutation']['foldchange']>5:
            summary[res]['interest_score'] += 20
            temp = 20
        elif summary[res]['bestmutation']['foldchange']>1:
            summary[res]['interest_score'] += 0
            temp = 0

        wildcard_interaction = 0
        if temp>=20 and 'type_class' not in summary[res]['closestinteraction']: #if big foldchange but no interaction
            wildcard_interaction = 1
        elif summary[res]['bestmutation']['bigdecrease_distinct']>=1 and 'type_class' not in summary[res]['closestinteraction']: #if big decrease in class
            wildcard_interaction = 1

        if summary[res]['bestmutation']['foldchange']>5:
            summary[res]['score_text'] += '# closest foldchange: '+str(temp)+'<br>'
            #summary[res]['freq'] += summary[res]['bestmutation']['bigdecrease_distinct']

        if not summary[res]['bestmutation']['foldchange']>5 and summary[res]['bestmutation']['bigdecrease']>1:
            summary[res]['score_text'] += '# class decrease: '+str(10)+'<br>'
            summary[res]['interest_score'] += 10
            #summary[res]['freq'] += summary[res]['bestmutation']['bigdecrease_distinct']


        #summary[res]['scores'] = scores

        summary[res]['sub'] = ''
        summary[res]['existing_mutants_protein'] = ''
        summary[res]['existing_mutants_family'] = '' #Change to number of
        #summary[res]['existing_mutants_protein'] = 0
        #summary[res]['existing_mutants_family'] = 0


        summary[res]['suggestion'] = OrderedDict()
        if 'type_class' in summary[res]['closestinteraction']:
            interaction_type = summary[res]['closestinteraction']['type_class']
            if interaction_type in matrix:
                if summary[res]['aa'] in matrix[interaction_type]:
                    possible_subs = copy.deepcopy(matrix[interaction_type][summary[res]['aa']][0])
                    for i,sub in enumerate(possible_subs):
                        if len(sub)>1:
                            possible_subs[i] = '/'.join(sub)
                    possible_subs_text = matrix[interaction_type][summary[res]['aa']][1]
                    summary[res]['sub'] = possible_subs
                    summary[res]['subtext'] = possible_subs_text
                    summary[res]['subtextfull']  = list(zip(possible_subs, possible_subs_text))
                    summary[res]['suggestion'][interaction_type] = list(zip(possible_subs, possible_subs_text))

                    #Make list of references of existing mutants -- can be depreciated or moved to "popover"
                    if res in mutant_lookup:
                        for subs in possible_subs:
                            if len(subs)==1:
                                for m in mutant_lookup[res]:
                                    if summary[res]['aa']==m[0] : #and subs==m[1]
                                        if m[3]==0:
                                            summary[res]['existing_mutants_protein'] += "<br>"+m[0]+"=>"+m[1]+" "+m[2]
                                            #summary[res]['existing_mutants_protein'] += 1
                                        else:
                                            summary[res]['existing_mutants_family'] += "<br>"+m[0]+"=>"+m[1]+" "+m[2]
                                            #summary[res]['existing_mutants_family'] += 1
                            else:
                                for sub in subs:
                                    #print(res,summary[res]['aa'],sub,'mutant data?')
                                    for m in mutant_lookup[res]:
                                        if summary[res]['aa']==m[0] : #and sub==m[1]
                                            if m[3]==0:
                                                summary[res]['existing_mutants_protein'] += "<br>"+m[0]+"=>"+m[1]+" "+m[2]
                                                #summary[res]['existing_mutants_protein'] += 1
                                            else:
                                                summary[res]['existing_mutants_family'] += "<br>"+m[0]+"=>"+m[1]+" "+m[2]
                                                #summary[res]['existing_mutants_family'] += 1
                else: #if no aa in that type (ODD!) allow wildtype
                    wildcard_interaction = 1

            else:
                print('error',interaction_type,summary[res]['aa'])

        if secondary_interaction:
            if summary[res]['aa'] in matrix[secondary_interaction]:
                possible_subs = copy.deepcopy(matrix[secondary_interaction][summary[res]['aa']][0])
                for i,sub in enumerate(possible_subs):
                    if len(sub)>1:
                        possible_subs[i] = '/'.join(sub)
                possible_subs_text = matrix[secondary_interaction][summary[res]['aa']][1]
                summary[res]['suggestion'][secondary_interaction +' (secondary)'] = list(zip(possible_subs, possible_subs_text))

        if wildcard_interaction:
            ignore_hyd = 0
            if summary[res]['aa'] in list(matrix['aromatic'].keys())+list(matrix['polar'].keys()):
                #ignore hydrophobic sug if there is from arom and polar
                ignore_hyd=1
            for interaction_type in matrix:
                if interaction_type=='hydrophobic' and ignore_hyd:
                    continue
                if summary[res]['aa'] in matrix[interaction_type]:
                    possible_subs = copy.deepcopy(matrix[interaction_type][summary[res]['aa']][0])
                    for i,sub in enumerate(possible_subs):
                        if len(sub)>1:
                            possible_subs[i] = '/'.join(sub)
                    possible_subs_text = matrix[interaction_type][summary[res]['aa']][1]
                    summary[res]['suggestion'][interaction_type] = list(zip(possible_subs, possible_subs_text))

        summary[res]['final_suggestions'] = {}
        for interaction_type,sugs in summary[res]['suggestion'].items():
            for sug in sugs:
                if sug[0] in summary[res]['final_suggestions']:
                    summary[res]['final_suggestions'][sug[0]][interaction_type] = sug[1]
                else:
                    summary[res]['final_suggestions'][sug[0]] = OrderedDict()
                    summary[res]['final_suggestions'][sug[0]][interaction_type] = sug[1]

        #print(summary[res]['final_suggestions'])
        summary[res]['suggestion'] = summary[res]['final_suggestions']

        #summary_score.append([summary[res]['interest_score'],summary[res],res])
        summary[res]['ranking_score'] += summary[res]['interest_score']*0.00001 #add a tiny ordering by magnitude if all else is  the same
        summary_score.append([summary[res]['ranking_score'],summary[res],res])

    print('made summary')


    sorted_summary = sorted(summary_score,key=lambda x: x[0],reverse=True)
    new_summary = OrderedDict();
    diagram_summary = []
    for res in sorted_summary:
        if res[0]>10:
            new_summary[res[2]] = res[1]
            diagram_summary.append({'pos':res[1]['pos'],'homology':res[1]['homology']})

    #print(sorted_summary)

    #summary = sorted_summary

    context['results'] = new_summary
    context['results_json'] = list(diagram_summary)
    context['family_ids'] = family_ids
    context['parent_ids'] = parent_ids
    context['HelixBox'] = HelixBox
    context['SnakePlot'] = SnakePlot

    print('sending to render')
    return render(request, 'mutation/design.html', context)


# Create your views here.
def ajax(request, slug, **response_kwargs):
    if '[' in slug:
        x = ast.literal_eval(urllib.parse.unquote(slug))
        mutations = MutationExperiment.objects.filter(protein__pk__in=x).order_by('residue__sequence_number').prefetch_related('residue')

        mutations_list = {}
        mutations_generic_number = {}
        residue_table_list = []
        for mutation in mutations:
            if not mutation.residue.generic_number: continue #cant map those without display numbers
            if mutation.residue.generic_number.label not in mutations_list: mutations_list[mutation.residue.generic_number.label] = []
            if mutation.ligand:
                ligand = mutation.ligand.name
            else:
                ligand = ''
            if mutation.exp_qual:
                qual = mutation.exp_qual.qual
            else:
                qual = ''
            mutations_list[mutation.residue.generic_number.label].append([mutation.foldchange,ligand,qual])

        a = Alignment()
        proteins = Protein.objects.filter(pk__in=x).all()
        a.load_proteins(proteins)
        segments = ProteinSegment.objects.all().filter().prefetch_related()
        a.load_segments(segments) #get all segments to make correct diagrams

        # build the alignment data matrix
        a.build_alignment()

        # calculate consensus sequence + amino acid and feature frequency
        a.calculate_statistics()

        jsondata = {}
        for aa in a.full_consensus:
            if aa.family_generic_number and aa.family_generic_number in a.generic_number_objs:
                if aa.family_generic_number in mutations_list:
                    jsondata[aa.sequence_number] = mutations_list[aa.family_generic_number]

    else:
        slug_without_species = slug.split('_')[0]
        mutations = MutationExperiment.objects.filter(protein__entry_name=slug).order_by('residue__sequence_number').prefetch_related('residue')
        # mutations = MutationExperiment.objects.filter(protein__entry_name__startswith=slug_without_species).order_by('residue__sequence_number').prefetch_related('residue')
        jsondata = {}
        for mutation in mutations:
            if mutation.residue.sequence_number not in jsondata: jsondata[mutation.residue.sequence_number] = []
            if mutation.ligand:
                ligand = mutation.ligand.name
            else:
                ligand = ''
            if mutation.exp_qual:
                qual = mutation.exp_qual.qual
            else:
                qual = ''
            jsondata[mutation.residue.sequence_number].append([mutation.foldchange,ligand,qual])

    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)

def ajaxSegments(request, slug, segments, **response_kwargs):
    if '[' in slug:
        x = ast.literal_eval(urllib.parse.unquote(slug))
        segments = ast.literal_eval(urllib.parse.unquote(segments))
        mutations = MutationExperiment.objects.filter(protein__pk__in=x,residue__protein_segment__slug__in=segments).order_by('residue__sequence_number').prefetch_related('residue')
    else:
        mutations = MutationExperiment.objects.filter(protein__entry_name=slug).order_by('residue__sequence_number').prefetch_related('residue')
    jsondata = {}
    for mutation in mutations:
        if mutation.residue.sequence_number not in jsondata: jsondata[mutation.residue.sequence_number] = []
        jsondata[mutation.residue.sequence_number].append(mutation.foldchange)

    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)

def importmutation(request):

    rows = loaddatafromexcel('/vagrant/protwis/mutation/import.xlsx')

    rows = analyse_rows(rows)

    whattoreturn = []
    c = 0
    skipped = 0
    inserted = 0
    for r in rows:
        #print(r)
        raw_id = insert_raw(r)
        ref_id = check_reference(r['reference'])
        lig_id = get_ligand(r)

        if (ref_id==None):
            whattoreturn.append(['Skipped due to no working reference',r])
            skipped += 1
            continue

        protein_id = 0
        residue_id = 0

        check=Protein.objects.filter(entry_name=r['protein'])
        if check.exists():
            check=Protein.objects.get(entry_name=r['protein'])
            protein_id = check
        else:
            whattoreturn.append(['Skipped due to no protein',r['protein']])
            skipped += 1
            continue

        check=Residue.objects.filter(protein_conformation__protein=protein_id,sequence_number=r['mutation_pos'])
        if check.exists():
            check=Residue.objects.get(protein_conformation__protein=protein_id,sequence_number=r['mutation_pos'])
            residue_id = check
        else:
            whattoreturn.append(['Skipped due to no residue',r['protein'],r['mutation_pos']])
            skipped += 1
            continue
            # residue_id = Residue()
            # residue_id.protein = protein_id
            # residue_id.sequence_number = r['mutation_pos']
            # residue_id.amino_acid = r['mutation_from']
            # residue_id.save()

        obj, created = MutationLigandClass.objects.get_or_create(classname=r['ligand_class'])
        ligclass_id = obj

        obj, created = MutationLigandRef.objects.get_or_create(reference=r['exp_mu_ligand_ref'])
        ligref_id = obj

        obj, created = MutationExperimentalType.objects.get_or_create(type=r['exp_type'])
        exp_type_id = obj

        obj, created = MutationFunc.objects.get_or_create(func=r['exp_func'])
        exp_func_id = obj

        obj, created = MutationMeasure.objects.get_or_create(measure=r['exp_mu_effect_type'])
        exp_measure_id = obj

        obj, created = MutationQual.objects.get_or_create(qual=r['exp_mu_effect_qual'], prop=r['exp_mu_effect_ligand_prop'])
        exp_qual_id = obj

        obj, created = MutationLigandRef.objects.get_or_create(reference=r['exp_mu_ligand_ref'])
        effect_ligand_reference_id = obj

        obj, created =  MutationOptional.objects.get_or_create(type=r['opt_type'], wt=r['opt_wt'], mu=r['opt_mu'], sign=r['opt_sign'], percentage=r['opt_percentage'], qual=r['opt_qual'], agonist=r['opt_agonist'])
        exp_opt_id = obj

        obj, created =  Mutation.objects.get_or_create(amino_acid=r['mutation_to'],protein=protein_id, residue=residue_id)
        mutation_id = obj



        logtypes = ['pEC50','pIC50','pK']


        foldchange = 0
        typefold = ''
        if r['exp_mu_effect_type']=='Activity/affinity' and r['exp_wt_value']!=0:

            if re.match("(" + ")|(".join(logtypes) + ")", r['exp_type']):  #-log values!
                foldchange = round(math.pow(10,-r['exp_mu_value_raw'])/pow(10,-r['exp_wt_value']),3);
                typefold = r['exp_type']+"_log"
            else:
                foldchange = round(r['exp_mu_value_raw']/r['exp_wt_value'],3);
                typefold = r['exp_type']+"_not_log"


            if foldchange<1 and foldchange!=0:
                foldchange = -round((1/foldchange),3)
            elif r['exp_mu_effect_type'] =='Fold effect (mut/wt)':
                foldchange = round(r['exp_mu_value_raw'],3);
                if foldchange<1: foldchange = -round((1/foldchange),3);


        obj, created = MutationExperiment.objects.get_or_create(
        refs=ref_id,
        protein=protein_id,
        residue=residue_id, #MISSING
        ligand=lig_id,
        ligand_class=ligclass_id,
        ligand_ref = ligref_id,
        raw = raw_id,
        optional = exp_opt_id,
        exp_type=exp_type_id,
        exp_func=exp_func_id,
        exp_measure = exp_measure_id,
        exp_qual = exp_qual_id,

        mutation=mutation_id,
        wt_value=r['exp_wt_value'], #
        wt_unit=r['exp_wt_unit'],

        mu_value = r['exp_mu_value_raw'],
        mu_sign = r['exp_mu_effect_sign'],
        foldchange = foldchange
        #foldchange = 1

        #added_by='munk',
        #added_date=datetime.now()
        )
        #print(foldchange)
        mut_id = obj.id


        #whattoreturn.append([protein_id,residue_id,raw_id,ref_id,lig_id,ligclass_id,exp_type_id,exp_func_id,exp_measure_id,exp_qual_id,typefold,foldchange,mut_id])
        whattoreturn.append(['Inserted',protein_id.entry_name,residue_id.sequence_number,foldchange])
        inserted += 1
        c += 1
        #if c>10: break

    print('DONE!')
    context = {'rows': whattoreturn, 'skipped' : skipped, 'inserted' : inserted}

    #return HttpResponse(ref_id)
    return render(request,'mutation/index.html',context)
