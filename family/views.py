from django.shortcuts import get_object_or_404, render
from django.views import generic

from protein.models import Protein, ProteinFamily


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

    return render(request, 'family/family_detail.html', {'pf': pf, 'families': families,
        'no_of_proteins': no_of_proteins, 'no_of_human_proteins': no_of_human_proteins})