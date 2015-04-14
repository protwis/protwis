from django.shortcuts import render
from django.views.generic import TemplateView
from django.http import HttpResponse
from build_structures.models import Structure


class StructureBrowser(TemplateView):
    pass