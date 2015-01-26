from django.shortcuts import get_object_or_404, render
from django.views import generic

from protein.models import Protein
from protein.models import ProteinFamily

import inspect


class IndexView(generic.ListView):
    model = Protein

    def get_queryset(self):
        return Protein.objects.all()


class DetailView(generic.DetailView):
    model = Protein
    slug_field = 'entry_name'


class TargetSelection(generic.TemplateView):
    template_name = 'protein/targetselection.html'

    step = 1
    number_of_steps = 2
    title = 'Select targets'
    description = 'Select receptors by searching or browsing in the middle column. You can select entire receptor families or individual receptors.\n\nSelected receptors will appear in the right column, where you can edit the list.\n\nOnce you have selected all your receptors, click the green button.'
    docs = '/docs/protein'
    pfs = ProteinFamily.objects.all()
    ps = Protein.objects.all()

    def get_context_data(self, **kwargs):
        # get context from parent class (really only relevant for child classes of this class, as TemplateView does
        # not have any context variables)
        context = super().get_context_data(**kwargs)

        # get attributes of this class and add them to the context
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]
        return context
