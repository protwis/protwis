from django.shortcuts import render
from django.views.generic import TemplateView
from django.http import HttpResponse
from structure.models import Structure


class StructureBrowser(TemplateView):

    template_name = "structure_browser.html"

    def get_context_data(self, **kwargs):

        context = super(StructureBrowser, self).get_context_data(**kwargs)
        try:
            context['crystals'] = Structure.objects.all()
        except Structure.DoesNotExist as e:
            pass

        return context



class Statistics(TemplateView):

    template_name = 'statistics.html'
    pass


#==============================================================================