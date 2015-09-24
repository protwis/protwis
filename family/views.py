from django.shortcuts import get_object_or_404, render
from django.conf import settings
from django.views import generic

from common.diagrams_gpcr import DrawHelixBox, DrawSnakePlot

from protein.models import Protein, ProteinFamily, ProteinSegment
from residue.models import Residue,ResidueGenericNumber
from mutation.models import MutationExperiment
Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')

def detail(request, slug):
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
    no_of_proteins = Protein.objects.filter(family__slug__startswith=pf.slug).count()
    no_of_human_proteins = Protein.objects.filter(family__slug__startswith=pf.slug, species__id=1).count()

        # create an alignment object
    a = Alignment()

    # fetch proteins and segments
    proteins = Protein.objects.filter(family__slug__startswith=slug, sequence_type__slug='wt')
    segments = ProteinSegment.objects.filter(partial=False)

    mutations = MutationExperiment.objects.filter(protein__in=proteins)
    
    # load data into the alignment
    a.load_proteins(proteins)
    a.load_segments(segments)

    # build the alignment data matrix
    a.build_alignment()

    # calculate consensus sequence + amino acid and feature frequency
    a.calculate_statistics()

    residue_list = []
    generic_numbers = []
    reference_generic_numbers = {}
    count = 0 #build sequence_number

    ################################################################################
    #FIXME -- getting matching display_generic_numbers kinda randomly.

    for seg in a.consensus: #Grab list of generic_numbers to lookup for their display numbers
        for aa in a.consensus[seg]:
            if "x" in aa:
                generic_numbers.append(aa)

    generic_ids = Residue.objects.filter(generic_number__label__in=generic_numbers).values('id').distinct('generic_number__label').order_by('generic_number__label')
    rs = Residue.objects.filter(id__in=generic_ids).prefetch_related('display_generic_number','generic_number')

    for r in rs: #make lookup dic.
        reference_generic_numbers[r.generic_number.label] = r
    ################################################################################
    
    for seg in a.consensus:
        for aa,v in a.consensus[seg].items():
            r = Residue()
            r.sequence_number =  count #FIXME is this certain to be correct that the position in consensus is seq position? 
            if "x" in aa:
                r.display_generic_number = reference_generic_numbers[aa].display_generic_number #FIXME
                r.segment_slug = seg
                r.family_generic_number = aa
            else:
                r.segment_slug = seg
                r.family_generic_number = aa
            r.amino_acid = v[0]
            r.extra = v[2] #Grab consensus information
            residue_list.append(r)

            count += 1         
    HelixBox = DrawHelixBox(residue_list,'Class A',str('test'))
    SnakePlot = DrawSnakePlot(residue_list,'Class A',str('test'))

    return render(request, 'family/family_detail.html', {'pf': pf, 'families': families,
        'no_of_proteins': no_of_proteins, 'no_of_human_proteins': no_of_human_proteins, 'a':a, 'HelixBox':HelixBox, 'SnakePlot':SnakePlot, 'mutations':mutations})