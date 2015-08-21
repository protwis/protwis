from django.shortcuts import get_object_or_404, render
from django.conf import settings
from django.views import generic

from common.diagrams_gpcr import DrawHelixBox, DrawSnakePlot

from protein.models import Protein, ProteinFamily, ProteinSegment
from residue.models import Residue,ResidueGenericNumber
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
    
    # load data into the alignment
    a.load_proteins(proteins)
    a.load_segments(segments)

    # build the alignment data matrix
    a.build_alignment()

    # calculate consensus sequence + amino acid and feature frequency
    a.calculate_statistics()

    #print(a.consensus)
    residue_list = []
    generic_numbers = []
    reference_segments = {}
    reference_generic_numbers = {}
    count = 0 #build sequence_number

    for seg in a.consensus:
        for aa in seg.items():
            if "x" in aa[0]:
                generic_numbers.append(aa[0])

    generic_ids = Residue.objects.filter(generic_number__label__in=generic_numbers).values('id').distinct('generic_number__label').order_by('generic_number__label')
    rs = Residue.objects.filter(id__in=generic_ids).prefetch_related('display_generic_number','generic_number')

    for r in rs:
        reference_generic_numbers[r.generic_number.label] = r

    for seg in a.consensus:
        for aa in seg.items():
            r = Residue()
            r.sequence_number =  count
            if "x" in aa[0]:
                #r.generic_number = reference_generic_numbers[aa[0]].generic_number
                r.display_generic_number = reference_generic_numbers[aa[0]].display_generic_number
                #r.protein_segment = reference_generic_numbers[aa[0]].protein_segment
                #print(aa[0][0])
                if int(aa[0][0])<8:
                    segment_slug = "TM"+str(aa[0][0])
                elif aa[0][0]=='8':
                    segment_slug = "H"+str(aa[0][0])
                r.segment_slug = segment_slug
                r.family_generic_number = aa[0]
            else:
                segment_slug = aa[0][:-5] #remove -0001 from slug
                r.segment_slug = segment_slug
                r.family_generic_number = aa[0]
                # if segment_slug not in reference_segments:
                #     reference_residue = Residue.objects.filter(protein_segment__slug=segment_slug).prefetch_related('protein_segment')[:1].get()
                #     reference_segments[segment_slug] = reference_residue
                # r.protein_segment = reference_segments[segment_slug].protein_segment
            r.amino_acid = aa[1][0]
            r.extra = aa[1][2]
            residue_list.append(r)

            count += 1

    HelixBox = DrawHelixBox(residue_list,'Class A',str('test'))
    SnakePlot = DrawSnakePlot(residue_list,'Class A',str('test'))

    return render(request, 'family/family_detail.html', {'pf': pf, 'families': families,
        'no_of_proteins': no_of_proteins, 'no_of_human_proteins': no_of_human_proteins, 'a':a, 'HelixBox':HelixBox, 'SnakePlot':SnakePlot})