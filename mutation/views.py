from django.shortcuts import render
from django.template import loader, Context
from django.db.models import Count, Min, Sum, Avg, Q
from django.http import HttpResponse
from django.conf import settings
from mutation.functions import *
from mutation.models import *

from common.views import AbsTargetSelection
from common.views import AbsSegmentSelection
from common.diagrams_gpcr import DrawHelixBox, DrawSnakePlot
from common import definitions

from residue.models import Residue,ResidueNumberingScheme
from residue.views import ResidueTablesDisplay
from protein.models import Protein,ProteinSegment
from interaction.models import ResidueFragmentInteraction, StructureLigandInteraction
from interaction.views import calculate
from interaction.forms import PDBform

from datetime import datetime
from collections import OrderedDict
import json
#env/bin/python3 -m pip install xlrd

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

def render_mutations(request, protein = None, family = None, download = None, **response_kwargs):  

    # get the user selection from session
    simple_selection = request.session.get('selection', False)
    
     # local protein list
    proteins = []

    if protein: # if protein static page

        proteins.append(Protein.objects.get(entry_name = protein))
        segments_ids = ProteinSegment.objects.all().values('id')
        original_segments = ProteinSegment.objects.all()

    elif family:

        family_proteins = Protein.objects.filter(family__slug__startswith=family, sequence_type__slug='wt')
        for fp in family_proteins:
            proteins.append(fp)
        segments_ids = ProteinSegment.objects.all().values('id')
        original_segments = ProteinSegment.objects.all()

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
        segments_ids = []
        for segment in simple_selection.segments:
            original_segments.append(segment.item)
            segments_ids.append(segment.item.id)

       #scheme
    used_schemes = {}
    for protein in proteins:
        if protein.residue_numbering_scheme.slug not in used_schemes:
            used_schemes[protein.residue_numbering_scheme.slug] = 0
        used_schemes[protein.residue_numbering_scheme.slug] += 1

    used_scheme = max(used_schemes, key=used_schemes.get)

    mutations = MutationExperiment.objects.filter(protein__in=proteins, 
                            residue__protein_segment__in=original_segments).prefetch_related('residue__display_generic_number',
                            'residue__protein_segment','residue__generic_number','exp_func','exp_qual',
                            'exp_measure', 'exp_type', 'ligand_role', 'ligand','refs','raw',
                            'ligand__properities', 'refs__web_link', 'refs__web_link__web_resource')
     
    mutations_list = {}
    mutations_display_generic_number = {}

    residue_table_list = []
    for mutation in mutations:
        residue_table_list.append(mutation.residue.generic_number)
        if not mutation.residue.display_generic_number: continue #cant map those without display numbers
        if mutation.residue.display_generic_number.label not in mutations_list: mutations_list[mutation.residue.display_generic_number.label] = []
        if mutation.ligand:
            ligand = mutation.ligand.name
        else:
            ligand = ''
        if mutation.exp_qual:
            qual = mutation.exp_qual.qual
        else:
            qual = ''
        mutations_list[mutation.residue.display_generic_number.label].append([mutation.foldchange,ligand,qual])

        mutations_display_generic_number[mutation.raw.id] = mutation.residue.display_generic_number.label

    # create an alignment object
    a = Alignment()

    a.load_proteins(proteins)
    segments = ProteinSegment.objects.all().filter().prefetch_related()
    a.load_segments(segments) #get all segments to make correct diagrams

    # build the alignment data matrix
    a.build_alignment()

    # calculate consensus sequence + amino acid and feature frequency
    a.calculate_statistics()

    residue_list = []
    generic_numbers = []
    reference_generic_numbers = {}
    count = 1 #build sequence_number

    mutations_pos_list = {}
    for aa in a.full_consensus:
        # for aa,v in a.full_consensus[seg].items():
        r = Residue()
        r.sequence_number =  aa.sequence_number #FIXME is this certain to be correct that the position in consensus is seq position? 
        if aa.family_generic_number and aa.display_generic_number:
            r.display_generic_number = aa.display_generic_number #FIXME
            if r.display_generic_number.label in mutations_list:
                if r.sequence_number not in mutations_pos_list: 
                    mutations_pos_list[r.sequence_number] = []
                mutations_pos_list[r.sequence_number].append(mutations_list[r.display_generic_number.label])
            r.segment_slug = aa.segment_slug
            r.family_generic_number = aa.family_generic_number
        else:
            r.segment_slug = aa.segment_slug
            r.family_generic_number = aa.family_generic_number
        r.amino_acid = aa.amino_acid
        r.frequency = aa.frequency #Grab consensus information
        residue_list.append(r)

        count += 1         

    protein_ids = list(set([x.id for x in proteins]))
    HelixBox = DrawHelixBox(a.full_consensus,'Class A',str(protein_ids), nobuttons = 1)
    SnakePlot = DrawSnakePlot(a.full_consensus,'Class A',str(protein_ids), nobuttons = 1)

    context = {}
    numbering_schemes_selection = ['gpcrdba','gpcrdbb','gpcrdbc','gpcrdbf'] #there is a residue_numbering_scheme attribute on the protein model, so it's easy to find out
    numbering_schemes_selection = list(used_schemes.keys())
    numbering_schemes_selection = ['gpcrdb'] + numbering_schemes_selection #always use A for reference
    numbering_schemes = ResidueNumberingScheme.objects.filter(slug__in=numbering_schemes_selection).all()

    segments = ProteinSegment.objects.filter(pk__in=segments_ids,category='helix')

    if ResidueNumberingScheme.objects.get(slug=settings.DEFAULT_NUMBERING_SCHEME) in numbering_schemes:
        default_scheme = ResidueNumberingScheme.objects.get(slug=settings.DEFAULT_NUMBERING_SCHEME)
    else:
        default_scheme = numbering_schemes[0]

    # prepare the dictionary
    # each helix has a dictionary of positions
    # default_generic_number or first scheme on the list is the key
    # value is a dictionary of other gn positions and residues from selected proteins 
    if len(protein_ids)<100 and not download: #too many to make meaningful residuetable / not run when download
        data = OrderedDict()
        for segment in segments:
            data[segment.slug] = OrderedDict()
            residues = Residue.objects.filter(protein_segment=segment, protein_conformation__protein__in=proteins, 
                                            generic_number__in=residue_table_list).prefetch_related('protein_conformation__protein', 
                                            'protein_conformation__state', 'protein_segment',
                                            'generic_number','display_generic_number','generic_number__scheme', 'alternative_generic_numbers__scheme')
            for scheme in numbering_schemes:
                if scheme == default_scheme and scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                    for pos in list(set([x.generic_number.label for x in residues if x.protein_segment == segment])):
                        data[segment.slug][pos] = {scheme.slug : pos, 'seq' : ['-']*len(proteins)}
                elif scheme == default_scheme:
                    for pos in list(set([x.generic_number.label for x in residues if x.protein_segment == segment])):
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
        
        context['header'] = zip([x.short_name for x in numbering_schemes] + [x.name for x in proteins], [x.name for x in numbering_schemes] + [x.name for x in proteins],[x.name for x in numbering_schemes] + [x.entry_name for x in proteins])
        context['segments'] = [x.slug for x in segments if len(data[x.slug])]
        context['data'] = flattened_data
        context['number_of_schemes'] = len(numbering_schemes)
    else:
        context['data'] = ''
        context['header'] = ''
        context['segments'] = ''
        context['number_of_schemes'] = ''

    if download:
        raws = mutations.values('raw')
        rawmutations = MutationRaw.objects.filter(pk__in = raws).all()
        
        data = []
        for r in rawmutations:
            headers = []
            values = {}
            for field, val in r:
                headers.append(field)
                values[field] = val
            if values['id'] in mutations_display_generic_number:
                values['generic'] = mutations_display_generic_number[values['id']]
            else:
                values['generic'] = ''
            data.append(values)
        headers = ['reference', 'protein', 'mutation_pos', 'generic', 'mutation_from', 'mutation_to', 
        'ligand_name', 'ligand_idtype', 'ligand_id', 'ligand_class',
        'exp_type', 'exp_func',  'exp_wt_value',  'exp_wt_unit','exp_mu_effect_sign', 'exp_mu_effect_type', 'exp_mu_effect_value', 
        'exp_mu_effect_qual', 'exp_mu_effect_ligand_prop',  'exp_mu_ligand_ref', 'opt_type', 'opt_wt',
        'opt_mu', 'opt_sign', 'opt_percentage', 'opt_qual','opt_agonist', 'added_date'                
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
        return render(request, 'mutation/list.html', {'mutations': mutations, 'HelixBox':HelixBox, 'SnakePlot':SnakePlot, 'data':context['data'], 
            'header':context['header'], 'segments':context['segments'], 'number_of_schemes':len(numbering_schemes), 'mutations_pos_list' : json.dumps(mutations_pos_list), 'protein_ids':str(protein_ids)})

# Create your views here.
def index(request):
    return HttpResponse("Hello, world. You're at the polls index.")

class designPDB(AbsTargetSelection):

    # Left panel
    step = 1
    number_of_steps = 1
    docs = 'generic_numbering.html'  # FIXME

    # description = 'Select receptors to index by searching or browsing in the middle column. You can select entire' \
    #     + ' receptor families and/or individual receptors.\n\nSelected receptors will appear in the right column,' \
    #     + ' where you can edit the list.\n\nSelect which numbering schemes to use in the middle column.\n\nOnce you' \
    #     + ' have selected all your receptors, click the green button.'

    description = 'Mutant Design Tool'

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
    docs = 'generic_numbering.html'  # FIXME

    # description = 'Select receptors to index by searching or browsing in the middle column. You can select entire' \
    #     + ' receptor families and/or individual receptors.\n\nSelected receptors will appear in the right column,' \
    #     + ' where you can edit the list.\n\nSelect which numbering schemes to use in the middle column.\n\nOnce you' \
    #     + ' have selected all your receptors, click the green button.'

    description = 'Mutant Design Tool'

    # Middle section
    numbering_schemes = False
    filters = False
    search = True
    title = "Select annotated receptor interactions, PDB code or upload PDB file"

    template_name = 'mutation/designselection.html'

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

            if r['slug'][:5]=='polar':
                scoretype = 'polar'
            elif r['slug'][:3]=='aro':
                scoretype = 'aromatic'
            elif r['slug'][:3]=='hyd':
                scoretype = 'hyd'
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

def showcalculation(request):
    print(request.method)
    if request.method == 'POST':
        form = PDBform(request.POST, request.FILES)

        if 'file' in request.FILES: #uploaded file
            print('user upload')
        else:
            print('pdb code entered')

        context = calculate(request)

    else:
        print('protein selected')
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





    if len(context['proteins'])>1:
        return HttpResponse("Only pick one protein")

    #print(family_ids,parent_ids)

    residues = Residue.objects.filter(protein_conformation__protein=context['proteins'][0]).prefetch_related('display_generic_number')

    lookup = {}
    for r in residues:
        if r.display_generic_number: 
            lookup[r.display_generic_number.label] = r.amino_acid

    protein_interaction_pairs = StructureLigandInteraction.objects.filter(structure__protein_conformation__protein__parent__family__in=protein_ids,annotated=True)
  
    family_interaction_pairs = StructureLigandInteraction.objects.filter(structure__protein_conformation__protein__parent__family__parent__in=family_ids,annotated=True)
    
    parent_interaction_pairs = StructureLigandInteraction.objects.filter(structure__protein_conformation__protein__parent__family__parent__parent__in=parent_ids,annotated=True)

    protein_interactions = ResidueFragmentInteraction.objects.filter(
        structure_ligand_pair__structure__protein_conformation__protein__parent__family__in=protein_ids, structure_ligand_pair__annotated=True).order_by('rotamer__residue__sequence_number').prefetch_related('rotamer__residue__display_generic_number','interaction_type')

    family_interactions = ResidueFragmentInteraction.objects.filter(
        structure_ligand_pair__structure__protein_conformation__protein__parent__family__parent__in=family_ids, structure_ligand_pair__annotated=True).order_by('rotamer__residue__sequence_number').prefetch_related('rotamer__residue__display_generic_number','interaction_type')
    
    family_interactions = family_interactions.exclude(structure_ligand_pair__structure__protein_conformation__protein__parent__family__in=protein_ids)

    parent_interactions = ResidueFragmentInteraction.objects.filter(
        structure_ligand_pair__structure__protein_conformation__protein__parent__family__parent__parent__in=parent_ids, structure_ligand_pair__annotated=True).prefetch_related('rotamer__residue__display_generic_number','interaction_type')

    parent_interactions = parent_interactions.exclude(structure_ligand_pair__structure__protein_conformation__protein__parent__family__parent__in=family_ids)

    protein_mutations = MutationExperiment.objects.filter(protein__family__in=protein_ids).order_by('refs__year').prefetch_related('residue__display_generic_number','mutation','refs__web_link')

    family_mutations = MutationExperiment.objects.filter(protein__family__parent__in=family_ids).order_by('refs__year').prefetch_related('residue__display_generic_number','mutation','refs__web_link')

    family_mutations = family_mutations.exclude(protein__family__in=protein_ids)

    parent_mutations = MutationExperiment.objects.filter(protein__family__parent__parent__in=parent_ids).order_by('refs__year').prefetch_related('residue__display_generic_number','mutation','refs__web_link')

    parent_mutations = parent_mutations.exclude(protein__family__parent__in=family_ids)

    #NEW CLASS METHOD, then select closest
    class_interactions = ResidueFragmentInteraction.objects.filter(
        structure_ligand_pair__structure__protein_conformation__protein__family__slug__startswith=family.slug, structure_ligand_pair__annotated=True).prefetch_related('rotamer__residue__display_generic_number','interaction_type','structure_ligand_pair__structure__protein_conformation__protein__parent')

    class_mutations = MutationExperiment.objects.filter(
        protein__family__slug__startswith=family.slug).prefetch_related('protein','residue__display_generic_number','mutation','refs__web_link', 'exp_qual')

    class_proteins = Protein.objects.filter(family__slug__startswith=family.slug, source__name='SWISSPROT',species__common_name='Human').all()

    family_proteins = Protein.objects.filter(family__parent__in=family_ids, source__name='SWISSPROT',species__common_name='Human').all()
    ligand_proteins = Protein.objects.filter(family__parent__parent__in =parent_ids, source__name='SWISSPROT',species__common_name='Human').all()

    class_p = []

    for i in class_interactions:
        p = i.structure_ligand_pair.structure.protein_conformation.protein.parent
        if p not in class_p:
            class_p.append(p)
    for m in class_mutations:
        p = m.protein
        if p not in class_p:
            class_p.append(p)
    for p in class_proteins:
        if p not in class_p:
            class_p.append(p)
    print("Proteins in alignment",len(class_p),'proteins in family',len(family_proteins),'proteins in ligand class',len(ligand_proteins))

    #print(list(settings.REFERENCE_POSITIONS.keys()))
    segments = ProteinSegment.objects.filter(slug__in=list(settings.REFERENCE_POSITIONS.keys()))

    a = Alignment()
    a.load_reference_protein(context['proteins'][0])
    a.load_proteins(family_proteins)
    #segments = ProteinSegment.objects.filter(category='helix')
    a.load_segments(segments)
    a.build_alignment()
    a.calculate_similarity()
    a.calculate_statistics()
    family_generic_aa_count = a.calculate_aa_count_per_generic_number()

    a = Alignment()
    a.load_reference_protein(context['proteins'][0])
    a.load_proteins(ligand_proteins)
    #segments = ProteinSegment.objects.filter(category='helix')
    a.load_segments(segments)
    a.build_alignment()
    a.calculate_similarity()
    a.calculate_statistics()
    ligand_generic_aa_count = a.calculate_aa_count_per_generic_number()


    a = Alignment()
    a.load_reference_protein(context['proteins'][0])
    a.load_proteins(class_p)
    #segments = ProteinSegment.objects.filter(category='helix')
    #segments = ProteinSegment.objects.all()
    a.load_segments(segments) #get all segments to make correct diagrams

    # build the alignment data matrix
    a.build_alignment()

    # calculate similarity
    a.calculate_similarity()

    a.calculate_statistics()
    generic_aa_count = a.calculate_aa_count_per_generic_number()
    alternative_aa = a.aa_count_with_protein

    similarity_list = {}
    for p in a.proteins:
        similarity_list[p.protein.entry_name] = [int(p.identity),int(p.similarity),p.similarity_score]
        if (p.protein.entry_name==context['proteins'][0].entry_name):
            similarity_list[p.protein.entry_name] = [int(100),int(100),1000]
        #print(p.protein.entry_name,similarity_list[p.protein.entry_name])
    #print(similarity_list)

    results = {}
    mutant_lookup = {}
    level = 0
    for data in [protein_interactions,family_interactions,parent_interactions]:

        for i in data:
            interaction_type = i.interaction_type.slug
            if i.rotamer.residue.display_generic_number:
                generic = i.rotamer.residue.display_generic_number.label
                if generic in results:
                    results[generic]['interaction'][level].append([interaction_type,i.rotamer.residue.amino_acid])
                else:
                    results[generic] = {'interaction': {0:[], 1:[], 2:[]}, 'mutant': {0:[], 1:[], 2:[] }, 
                    'bestinteraction': { 'closest' : { 'similarity' : 0} }, 
                    'bestmutation':{ 'species' : '', 'similarity' : 0, 'foldchange' : '', 'qual' : '', 'allmut' : [], 'counts':{}, 'counts_close':{}, 'bigdecrease':0,'bigincrease':0,'nonsignificant':0,'nodata':0} }
                    results[generic]['interaction'][level].append([interaction_type,i.rotamer.residue.amino_acid])
                    mutant_lookup[generic] = []
            else:
                pass
                #print('no generic number',interaction_type)
        level += 1

    for i in class_interactions:
        interaction_type = i.interaction_type.slug
        interaction_type_class = i.interaction_type.type
        if i.rotamer.residue.display_generic_number:
            generic = i.rotamer.residue.display_generic_number.label
            entry_name = i.structure_ligand_pair.structure.protein_conformation.protein.parent.entry_name
            if generic in lookup:
                if generic in results and lookup[generic] == i.rotamer.residue.amino_acid:
                    if interaction_type_class not in results[generic]['bestinteraction']:
                        results[generic]['bestinteraction'][interaction_type_class] = { 'species' : '', 'similarity' : 0}
                    if similarity_list[entry_name][1]>results[generic]['bestinteraction']['closest']['similarity']:
                        results[generic]['bestinteraction']['closest']['species'] = entry_name
                        results[generic]['bestinteraction']['closest']['similarity'] = similarity_list[entry_name][1]
                        results[generic]['bestinteraction']['closest']['type'] = interaction_type
                    if similarity_list[entry_name][1]>results[generic]['bestinteraction'][interaction_type_class]['similarity']:
                        results[generic]['bestinteraction'][interaction_type_class]['species'] = entry_name
                        results[generic]['bestinteraction'][interaction_type_class]['similarity'] = similarity_list[entry_name][1]

            else:
                results[generic] = {'interaction': {0:[], 1:[], 2:[]}, 'mutant': {0:[], 1:[], 2:[] }, 
                'bestinteraction': { 'closest' : { 'similarity' : 0} }, 
                'bestmutation':{ 'species' : '', 'similarity' : 0, 'foldchange' : '', 'qual' : '', 'allmut' : [], 'counts':{}, 'counts_close':{}, 'bigdecrease':0,'bigincrease':0,'nonsignificant':0,'nodata':0}}
                if generic in lookup:
                    if (lookup[generic] == i.rotamer.residue.amino_acid): #only for same aa (FIXME substitution)
                        results[generic]['bestinteraction']['closest']['species'] = entry_name
                        results[generic]['bestinteraction']['closest']['similarity'] = similarity_list[entry_name][1]
                        results[generic]['bestinteraction']['closest']['type'] = interaction_type
                        results[generic]['bestinteraction'][interaction_type_class] = { 'species' : entry_name, 'similarity' : similarity_list[entry_name][1]}
                
                mutant_lookup[generic] = []
        else:
            pass
            #print('no generic number',interaction_type)

    level = 0
    for data in [protein_mutations,family_mutations,parent_mutations]: 
        for m in data:
            if m.residue.display_generic_number:
                generic = m.residue.display_generic_number.label
                # if m.foldchange==0:
                #     continue #skip null foldchange (should add qualitive later)
                if generic in results:
                    results[generic]['mutant'][level].append([m.foldchange,m.residue.amino_acid]) #add more info
                    if level==0 or level==1: #save mutants that will be interesting for user
                        mutant_lookup[generic].append([m.residue.amino_acid,m.mutation.amino_acid,str(m.refs.web_link)+" "+str(m.refs.year),level])
                else:
                    #results[generic] = {'interaction': {0:[], 1:[], 2:[] }, 'mutant': {0:[], 1:[], 2:[] } }
                    #results[generic]['mutant'][level].append([m.foldchange,m.residue.amino_acid,m.mutation.amino_acid,m.refs.web_link])#add more info
                    #mutant_lookup[generic] = []
                    #mutant_lookup[generic].append([m.residue.amino_acid,m.mutation.amino_acid,str(m.refs.web_link)])
                    pass ### DO NOT calc on postions without interaction data

            else:
                pass
                #print(m.residue.sequence_number,m.foldchange)  
        level += 1

    for m in class_mutations:
        if m.residue.display_generic_number:
            generic = m.residue.display_generic_number.label
            entry_name = m.protein.entry_name
            if m.exp_qual:
                qual = m.exp_qual.qual +" "+m.exp_qual.prop
            else:
                qual = ''
            if generic in lookup:
                if generic in results:
                    if (similarity_list[entry_name][1]>=results[generic]['bestmutation']['similarity'] and lookup[generic] == m.residue.amino_acid):
                        results[generic]['bestmutation']['species'] = entry_name
                        #if results[generic]['bestmutation']['similarity']!=similarity_list[entry_name][1]:#reset allmut
                        #    results[generic]['bestmutation']['allmut'] = []
                        results[generic]['bestmutation']['similarity'] = similarity_list[entry_name][1]
                        results[generic]['bestmutation']['foldchange'] = m.foldchange
                        results[generic]['bestmutation']['qual'] = qual
                        results[generic]['bestmutation']['aa'] = m.mutation.amino_acid

                    if lookup[generic] == m.residue.amino_acid:

                        results[generic]['bestmutation']['allmut'].append([entry_name,m.foldchange,qual,m.mutation.amino_acid,results[generic]['bestmutation']['similarity']])

                        if m.foldchange>5:
                            results[generic]['bestmutation']['bigdecrease'] += 1
                        elif m.foldchange<-5:
                            results[generic]['bestmutation']['bigincrease'] += 1
                        elif (m.foldchange<5 or m.foldchange>-5) and m.foldchange!=0:
                            results[generic]['bestmutation']['nonsignificant'] += 1
                        else:
                            if m.exp_qual:
                                if m.exp_qual.qual=='Abolish':
                                    results[generic]['bestmutation']['bigdecrease'] += 1
                                elif m.exp_qual.qual=='Gain of':
                                    results[generic]['bestmutation']['bigincrease'] += 1
                                elif m.exp_qual.qual=='Increase':
                                    results[generic]['bestmutation']['nonsignificant'] += 1
                                elif m.exp_qual.qual=='Decrease':
                                    results[generic]['bestmutation']['nonsignificant'] += 1
                            else:
                                results[generic]['bestmutation']['nodata'] += 1

                        if m.mutation.amino_acid in results[generic]['bestmutation']['counts']:
                            results[generic]['bestmutation']['counts'][m.mutation.amino_acid] += 1
                        else:
                            results[generic]['bestmutation']['counts'][m.mutation.amino_acid] = 1

                        if m.mutation.amino_acid in results[generic]['bestmutation']['counts_close'] and similarity_list[entry_name][1]>60:
                            results[generic]['bestmutation']['counts_close'][m.mutation.amino_acid] += 1
                        elif similarity_list[entry_name][1]>60:
                            results[generic]['bestmutation']['counts_close'][m.mutation.amino_acid] = 1



            else:
                results[generic] = {'interaction': {0:[], 1:[], 2:[]}, 'mutant': {0:[], 1:[], 2:[] }, 
                'bestinteraction': { 'closest' : { 'similarity' : 0} }, 
                'bestmutation':{ 'species' : '', 'similarity' : 0, 'foldchange' : '', 'qual' : '', 'allmut' : [], 'counts':{}, 'counts_close':{}, 'bigdecrease':0,'bigincrease':0,'nonsignificant':0,'nodata':0 }}
                if generic in lookup:
                    if (lookup[generic] == m.residue.amino_acid): #only for same aa (FIXME substitution)
                        results[generic]['bestmutation']['species'] = entry_name
                        results[generic]['bestmutation']['similarity'] = similarity_list[entry_name][1]
                        results[generic]['bestmutation']['foldchange'] = m.foldchange
                        results[generic]['bestmutation']['qual'] = qual
                        results[generic]['bestmutation']['aa'] = m.mutation.amino_acid
                        results[generic]['bestmutation']['allmut'].append([entry_name,m.foldchange,qual,m.mutation.amino_acid,results[generic]['bestmutation']['similarity']])
                        results[generic]['bestmutation']['counts'][m.mutation.amino_acid] = 1

                        if similarity_list[entry_name][1]>60:
                            results[generic]['bestmutation']['counts_close'][m.mutation.amino_acid] = 1

                        if m.foldchange>5:
                            results[generic]['bestmutation']['bigdecrease'] += 1
                        elif m.foldchange<-5:
                            results[generic]['bestmutation']['bigincrease'] += 1
                        elif (m.foldchange<5 or m.foldchange>-5) and m.foldchange!=0:
                            results[generic]['bestmutation']['nonsignificant'] += 1
                        else:
                            if m.exp_qual:
                                if m.exp_qual.qual=='Abolish':
                                    results[generic]['bestmutation']['bigdecrease'] += 1
                                elif m.exp_qual.qual=='Gain of':
                                    results[generic]['bestmutation']['bigincrease'] += 1
                                elif m.exp_qual.qual=='Increase':
                                    results[generic]['bestmutation']['nonsignificant'] += 1
                                elif m.exp_qual.qual=='Decrease':
                                    results[generic]['bestmutation']['nonsignificant'] += 1
                            else:
                                results[generic]['bestmutation']['nodata'] += 1

                
                mutant_lookup[generic] = []
        else:
            pass
            #print('no generic number')


    matrix = definitions.DESIGN_SUBSTITUTION_MATRIX

    summary = {}
    for res,values in results.items():
        #print(res)


        if res in lookup:
            summary[res] = {}
            summary[res]['aa'] = lookup[res]
        else: #skip those that are not present in reference
            continue
        summary[res]['bestinteraction'] = values['bestinteraction']
        summary[res]['bestmutation'] = values['bestmutation']

        for aalist in alternative_aa[res].items():
            if aalist[0] !=summary[res]['aa']:
                for p in aalist[1]:
                    if similarity_list[p][1]>60: #close enough to suggest
                        print (res,p,similarity_list[p][1],aalist[0])

        if res in generic_aa_count:
            summary[res]['conservation'] = [family_generic_aa_count[res][lookup[res]],ligand_generic_aa_count[res][lookup[res]],generic_aa_count[res][lookup[res]]]
        scores = {'hyd':0,'aromatic':0,'polar':0,'unknown':0} #dict to keep track of scores to select subs.
        for type, level in values.items():
            #print(type)
            if type!='interaction': continue
            summary[res][type] = {}
            for level, values in level.items():

                if level==0: #weight for scoring
                    weight = 5
                elif level==1:
                    weight = 4
                elif level==2:
                    weight = 1

                #print(level)
                #print(values) 
                temp = {}
                temp_same_aa = {}
                for value in values:
                    if value[0][:2]=='HB':
                        scoretype = 'polar'
                    elif value[0][:3]=='aro':
                        scoretype = 'aromatic'
                    elif value[0][:3]=='HYD':
                        scoretype = 'hyd'
                    elif value[0][:5]=='polar':
                        scoretype = 'polar'
                    elif value[0][:3]=='aro':
                        scoretype = 'aromatic'
                    elif value[0][:3]=='hyd':
                        scoretype = 'hyd'
                    else:
                        scoretype = 'unknown'

                    if value[0] in temp:
                        temp[value[0]] += 1
                    else:
                        temp[value[0]] = 1

                    if value[1]==summary[res]['aa']: #keep seperat those with same AA (should be weighted higher)
                        if value[0] in temp_same_aa:
                            temp_same_aa[value[0]] += 1
                        else:
                            temp_same_aa[value[0]] = 1

                        scores[scoretype] += 1*weight



                temp = sorted(temp.items(), key=operator.itemgetter(1),reverse=True)
                temp_same_aa = sorted(temp_same_aa.items(), key=operator.itemgetter(1),reverse=True)
                s = ''
                for t in temp:
                    s += t[0] + '('+str(t[1])+')'

                s += '<br><font color=red>'

                for t in temp_same_aa:
                    s += t[0] + '('+str(t[1])+')'

                s += '</font>'
                summary[res][type][level]= s
        scores = sorted(scores.items(), key=operator.itemgetter(1),reverse=True)
        summary[res]['scores'] = scores

        summary[res]['sub'] = ''
        summary[res]['existing_mutants_protein'] = ''
        summary[res]['existing_mutants_family'] = '' #Change to number of
        #summary[res]['existing_mutants_protein'] = 0
        #summary[res]['existing_mutants_family'] = 0
        if scores[0][1]>=1: #minimum to look at the top score
            interaction_type = scores[0][0]
            if summary[res]['aa'] in matrix[interaction_type]:
                possible_subs = matrix[interaction_type][summary[res]['aa']][0]
                summary[res]['sub'] = possible_subs
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

                                    #print(m)
                        else:
                            for sub in subs:
                                print(res,summary[res]['aa'],sub,'mutant data?')
                                for m in mutant_lookup[res]:
                                    if summary[res]['aa']==m[0] : #and sub==m[1]
                                        if m[3]==0:
                                            summary[res]['existing_mutants_protein'] += "<br>"+m[0]+"=>"+m[1]+" "+m[2]
                                            #summary[res]['existing_mutants_protein'] += 1
                                        else:
                                            summary[res]['existing_mutants_family'] += "<br>"+m[0]+"=>"+m[1]+" "+m[2]
                                            #summary[res]['existing_mutants_family'] += 1

                                        #print(m)
            else:
                print('error',interaction_type,summary[res]['aa'])
    #print(mutant_lookup)

    context['results'] = summary
    context['family_ids'] = family_ids
    context['parent_ids'] = parent_ids

    context['protein_interaction_pairs'] = len(protein_interaction_pairs)
    context['family_interaction_pairs'] = len(family_interaction_pairs)
    context['parent_interaction_pairs'] = len(parent_interaction_pairs)

    context['protein_interactions'] = protein_interactions
    context['family_interactions'] = family_interactions
    context['parent_interactions'] = parent_interactions

    context['protein_mutations'] = protein_mutations
    context['family_mutations'] = family_mutations
    context['parent_mutations'] = parent_mutations


    return render(request, 'mutation/design.html', context)


# Create your views here.
def ajax(request, slug, **response_kwargs):
    if '[' in slug:
        x = ast.literal_eval(urllib.parse.unquote(slug))
        mutations = MutationExperiment.objects.filter(protein__pk__in=x).order_by('residue__sequence_number').prefetch_related('residue')
    else:
        mutations = MutationExperiment.objects.filter(protein__entry_name=slug).order_by('residue__sequence_number').prefetch_related('residue')
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