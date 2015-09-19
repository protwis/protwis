from django.shortcuts import render
from django.conf import settings
from django.views import generic
from django.core import serializers

from documentation.models import Documentation

import json, os

def index(request):
    data = serializers.serialize( "python", Documentation.objects.all(), fields=('title','description','image'))
    return render(request, 'documentation/overview.html', {'documentation': data})

def crystals(request):
    title = Documentation.objects.get(title__icontains='crystal')
    return render(request, 'documentation/index.html', {'documentation': title})

def diagrams(request):
    title = Documentation.objects.get(title__icontains='diagrams')
    return render(request, 'documentation/index.html', {'documentation': title})

def numbering(request):
    title = Documentation.objects.get(title__icontains='numbering')
    return render(request, 'documentation/index.html', {'documentation': title})

def pharmacophore(request):
    title = Documentation.objects.get(title__icontains='pharmacophore')
    return render(request, 'documentation/index.html', {'documentation': title})

def sequences(request):
    title = Documentation.objects.get(title__icontains='sequences')
    return render(request, 'documentation/index.html', {'documentation': title})

def similarities(request):
    title = Documentation.objects.get(title__icontains='similarities')
    return render(request, 'documentation/index.html', {'documentation': title})
