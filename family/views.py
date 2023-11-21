from django.shortcuts import get_object_or_404, render
from django.conf import settings
from django.http import HttpResponse
from django.views import generic
from django.core.cache import cache
from django.views.decorators.cache import cache_page

from common.diagrams_gpcr import DrawHelixBox, DrawSnakePlot

from protein.models import Protein, ProteinFamily, ProteinSegment, ProteinConformation
from residue.models import Residue,ResidueGenericNumber
from mutation.models import MutationExperiment
from structure.models import Structure
from interaction.models import ResidueFragmentInteraction
from mutational_landscape.models import NaturalMutations, PTMs


Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')

import json

def hex_to_RGB(hex):
  ''' "#FFFFFF" -> [255,255,255] '''
  # Pass 16 to the integer function for change of base
  return [int(hex[i:i+2], 16) for i in range(1,6,2)]


def RGB_to_hex(RGB):
  ''' [255,255,255] -> "#FFFFFF" '''
  # Components need to be integers for hex to make sense
  RGB = [int(x) for x in RGB]
  return "#"+"".join(["0{0:x}".format(v) if v < 16 else
            "{0:x}".format(v) for v in RGB])

def color_dict(gradient):
  ''' Takes in a list of RGB sub-lists and returns dictionary of
    colors in RGB and hex form for use in a graphing function
    defined later on '''
  return [RGB_to_hex(RGB) for RGB in gradient]

def linear_gradient(start_hex="#4682B4", finish_hex="#FFB347", n=10):
    # http://bsou.io/posts/color-gradients-with-python
  ''' returns a gradient list of (n) colors between
    two hex colors. start_hex and finish_hex
    should be the full six-digit color string,
    inlcuding the number sign ("#FFFFFF") '''
  # Starting and ending colors in RGB form
  s = hex_to_RGB(start_hex)
  f = hex_to_RGB(finish_hex)
  # Initilize a list of the output colors with the starting color
  RGB_list = [s]
  # Calcuate a color at each evenly spaced value of t from 1 to n
  for t in range(1, n):
    # Interpolate RGB vector for color at the current value of t
    curr_vector = [
      int(s[j] + (float(t)/(n-1))*(f[j]-s[j]))
      for j in range(3)
    ]
    # Add it to our list of output colors
    RGB_list.append(curr_vector)

  return color_dict(RGB_list)

@cache_page(60 * 60 * 24 * 7)
def detail(request, slug):
    # FULL class A is too big for consensus
    if slug == "001":
        return HttpResponse("Displaying a consensus of all class A receptors is currently not supported.")

    # get family
    pf = ProteinFamily.objects.get(slug=slug)

    # get family list
    ppf = pf
    families = [ppf.name]
    while ppf.parent.parent:
        families.append(ppf.parent.name)
        ppf = ppf.parent
    families.reverse()

    # number of proteins
    proteins = Protein.objects.filter(family__slug__startswith=pf.slug, sequence_type__slug='wt')
    no_of_proteins = proteins.count()
    no_of_human_proteins = Protein.objects.filter(family__slug__startswith=pf.slug, species__id=1,
        sequence_type__slug='wt').count()
    list_proteins = list(proteins.values_list('pk',flat=True))


    # get structures of this family
    structures = Structure.objects.filter(protein_conformation__protein__parent__family__slug__startswith=slug
        ).exclude(structure_type__slug__startswith='af-').order_by('-representative', 'resolution').prefetch_related('pdb_code__web_resource')

    mutations = MutationExperiment.objects.filter(protein__in=proteins).prefetch_related('residue__generic_number',
                                'exp_qual', 'ligand')


    interactions = ResidueFragmentInteraction.objects.filter(
        structure_ligand_pair__structure__protein_conformation__protein__parent__in=proteins,
        structure_ligand_pair__annotated=True).exclude(interaction_type__type='hidden').prefetch_related(
            'rotamer__residue__generic_number', 'interaction_type'
            ).order_by('rotamer__residue__sequence_number')
    interaction_list = {}
    for interaction in interactions:
        if interaction.rotamer.residue.generic_number:
            gn = interaction.rotamer.residue.generic_number.label
            aa = interaction.rotamer.residue.amino_acid
            interactiontype = interaction.interaction_type.name
            if gn not in interaction_list:
                interaction_list[gn] = []
            interaction_list[gn].append([aa, interactiontype])

    ## Variants
    NMs = NaturalMutations.objects.filter(
        protein__in=proteins).prefetch_related('residue__generic_number')

    natural_mutation_list = {}
    max_snp_pos = 1
    for NM in NMs:
        if NM.residue.generic_number:
            if NM.residue.generic_number.label in natural_mutation_list:
                natural_mutation_list[NM.residue.generic_number.label]['val'] += 1
                if not str(NM.amino_acid) in natural_mutation_list[NM.residue.generic_number.label]['AA']:
                    natural_mutation_list[NM.residue.generic_number.label]['AA'] = natural_mutation_list[NM.residue.generic_number.label]['AA'] + str(NM.amino_acid) + ' '

                if natural_mutation_list[NM.residue.generic_number.label]['val'] > max_snp_pos:
                    max_snp_pos = natural_mutation_list[NM.residue.generic_number.label]['val']
            else:
                natural_mutation_list[NM.residue.generic_number.label] = {'val':1, 'AA': NM.amino_acid + ' '}

    ## PTMs
    ptms = PTMs.objects.filter(
        protein__in=proteins).prefetch_related('residue__generic_number')

    ptm_list = {}
    for ptm in ptms:
        if ptm.residue.generic_number:
            if ptm.residue.generic_number.label in ptm_list:
                ptm_list[ptm.residue.generic_number.label]['val'] += 1
                if not str(ptm.modification) in ptm_list[ptm.residue.generic_number.label]['mod']:
                    ptm_list[ptm.residue.generic_number.label]['mod'] = ptm_list[ptm.residue.generic_number.label]['mod'] + ', ' + str(ptm.modification)
            else:
                ptm_list[ptm.residue.generic_number.label] = {'val':1, 'mod': ptm.modification}

    # CMs = CancerMutations.objects.filter(
    #     protein__in=proteins).prefetch_related('residue__generic_number')
    #
    # cancer_mutation_list = {}
    # max_cancer_pos = 1
    # for CM in CMs:
    #     if CM.residue.generic_number:
    #         if CM.residue.generic_number.label in cancer_mutation_list:
    #             cancer_mutation_list[CM.residue.generic_number.label]['val'] += 1
    #             if not str(CM.amino_acid) in cancer_mutation_list[CM.residue.generic_number.label]['AA']:
    #                 cancer_mutation_list[CM.residue.generic_number.label]['AA'] = cancer_mutation_list[CM.residue.generic_number.label]['AA'] + str(CM.amino_acid)
    #
    #             if cancer_mutation_list[CM.residue.generic_number.label]['val'] > max_cancer_pos:
    #                 max_cancer_pos = cancer_mutation_list[CM.residue.generic_number.label]['val']
    #         else:
    #             cancer_mutation_list[CM.residue.generic_number.label] = {'val':1, 'AA': CM.amino_acid}
    #
    # DMs = DiseaseMutations.objects.filter(
    #     protein__in=proteins).prefetch_related('residue__generic_number')

    # disease_mutation_list = {}
    # max_disease_pos = 1
    # for DM in DMs:
    #     if DM.residue.generic_number:
    #         if DM.residue.generic_number.label in disease_mutation_list:
    #             disease_mutation_list[DM.residue.generic_number.label]['val'] += 1
    #             if not str(DM.amino_acid) in disease_mutation_list[DM.residue.generic_number.label]['AA']:
    #                 disease_mutation_list[DM.residue.generic_number.label]['AA'] = disease_mutation_list[DM.residue.generic_number.label]['AA'] + str(DM.amino_acid)
    #
    #             if disease_mutation_list[DM.residue.generic_number.label]['val'] > max_cancer_pos:
    #                 max_cancer_pos = disease_mutation_list[DM.residue.generic_number.label]['val']
    #         else:
    #             disease_mutation_list[DM.residue.generic_number.label] = {'val':1, 'AA': DM.amino_acid}

    mutations_list = {}
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
        mutations_list[mutation.residue.generic_number.label].append([mutation.foldchange,ligand.replace("'", "\\'"),qual])

    try:
        pc = ProteinConformation.objects.get(protein__family__slug=slug, protein__sequence_type__slug='consensus')
    except ProteinConformation.DoesNotExist:
        try:
            # In case of single members, not all families have a set consensus - grab the consensus of that single member
            pc = ProteinConformation.objects.get(protein__family__slug__startswith=slug, protein__sequence_type__slug='consensus')
        except ProteinConformation.DoesNotExist:
            try:
                pc = ProteinConformation.objects.get(protein__family__slug=slug, protein__species_id=1,
                    protein__sequence_type__slug='wt')
            except ProteinConformation.DoesNotExist:
                # In case of single members, not all families have a set consensus - grab the human representative of that single member
                pc = ProteinConformation.objects.get(protein__family__slug__startswith=slug, protein__species_id=1,
                    protein__sequence_type__slug='wt')

    residues = Residue.objects.filter(protein_conformation=pc).order_by('sequence_number').prefetch_related(
        'protein_segment', 'generic_number', 'display_generic_number')

    jsondata = {}
    jsondata_interaction = {}
    jsondata_natural_mutations = {}
    jsondata_ptms = {}
    # jsondata_cancer_mutations = {}
    # jsondata_disease_mutations = {}
    for r in residues:
        if r.generic_number:
            if r.generic_number.label in mutations_list:
                jsondata[r.sequence_number] = [mutations_list[r.generic_number.label]]
            if r.generic_number.label in interaction_list:
                jsondata_interaction[r.sequence_number] = interaction_list[r.generic_number.label]
            if r.generic_number.label in natural_mutation_list:
                jsondata_natural_mutations[r.sequence_number] = natural_mutation_list[r.generic_number.label]
            if r.generic_number.label in ptm_list:
                jsondata_ptms[r.sequence_number] = ptm_list[r.generic_number.label]
            # if r.generic_number.label in cancer_mutation_list:
                # jsondata_cancer_mutations[r.sequence_number] = cancer_mutation_list[r.generic_number.label]
            # if r.generic_number.label in disease_mutation_list:
                # jsondata_disease_mutations[r.sequence_number] = disease_mutation_list[r.generic_number.label]

    jsondata_natural_mutations['color'] = linear_gradient(start_hex="#c79494", finish_hex="#c40100", n=max_snp_pos)
    # jsondata_cancer_mutations['color'] = linear_gradient(start_hex="#d8baff", finish_hex="#422d65", n=max_cancer_pos)
    # jsondata_disease_mutations['color'] = linear_gradient(start_hex="#ffa1b1", finish_hex="#6e000b", n=max_disease_pos)

    HelixBox = DrawHelixBox(residues, 'Class A', 'family_diagram_preloaded_data')
    SnakePlot = DrawSnakePlot(residues, 'Class A', 'family_diagram_preloaded_data')

    # process residues and return them in chunks of 10
    # this is done for easier scaling on smaller screens
    chunk_size = 10
    r_chunks = []
    r_buffer = []
    last_segment = False
    border = False
    title_cell_skip = 0
    for i, r in enumerate(residues):
        # title of segment to be written out for the first residue in each segment
        segment_title = False

        # keep track of last residues segment (for marking borders)
        if r.protein_segment.slug != last_segment:
            last_segment = r.protein_segment.slug
            border = True

        # if on a border, is there room to write out the title? If not, write title in next chunk
        if i == 0 or (border and len(last_segment) <= (chunk_size - i % chunk_size)):
            segment_title = True
            border = False
            title_cell_skip += len(last_segment) # skip cells following title (which has colspan > 1)

        if i and i % chunk_size == 0:
            r_chunks.append(r_buffer)
            r_buffer = []

        r_buffer.append((r, segment_title, title_cell_skip))

        # update cell skip counter
        if title_cell_skip > 0:
            title_cell_skip -= 1
    if r_buffer:
        r_chunks.append(r_buffer)

    context = {'pf': pf, 'families': families, 'structures': structures, 'no_of_proteins': no_of_proteins,
        'no_of_human_proteins': no_of_human_proteins, 'HelixBox':HelixBox, 'SnakePlot':SnakePlot,
        'mutations':mutations, 'r_chunks': r_chunks, 'chunk_size': chunk_size, 'mutations_pos_list' : json.dumps(jsondata),'interaction_pos_list' : json.dumps(jsondata_interaction), 'natural_mutations_pos_list': json.dumps(jsondata_natural_mutations), 'ptms_pos_list': json.dumps(jsondata_ptms)} # ,'cancer_mutations_pos_list': json.dumps(jsondata_cancer_mutations), 'disease_mutations_pos_list': json.dumps(jsondata_disease_mutations)

    return render(request, 'family/family_detail.html', context)
