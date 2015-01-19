from django.shortcuts import get_object_or_404, render
from django.views import generic

from protein.models import Protein


class IndexView(generic.ListView):
    model = Protein

    def get_queryset(self):
        return Protein.objects.all()


class DetailView(generic.DetailView):
    model = Protein
    slug_field = 'entry_name'