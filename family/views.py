from django.shortcuts import get_object_or_404, render
from django.views import generic

from protein.models import  ProteinFamily


def detail(request, slug):
    # get family
    pf = ProteinFamily.objects.get(slug=slug)

    return render(request, 'family/family_detail.html', {'pf': pf})