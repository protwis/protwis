from django.shortcuts import render
from django.conf import settings
from django.views import generic
from django.views.static import serve
from django.core.files import File 

from pages.models import Pages
from common.models import ReleaseNotes
import os

def index(request):
	return render(request,'pages/index.html', {'site_title': settings.SITE_TITLE})

def releasenotes(request):
    context = {}
    context['release_notes'] = ReleaseNotes.objects.all()
    return render(request, 'pages/releasenotes.html', context)

def contact(request):
    title = Pages.objects.get(title__icontains='contact')
    return render(request, 'pages/index.html', {'pages': title})

def citing(request):
    title = Pages.objects.get(title__icontains='citing')
    return render(request, 'pages/index.html', {'pages': title})

def contribute(request):
    title = Pages.objects.get(title__icontains='contribute')
    return render(request, 'pages/index.html', {'pages': title})

def meetings(request):
    title = Pages.objects.get(title__icontains='meetings')
    return render(request, 'pages/index.html', {'pages': title})

def update(request):
    title = Pages.objects.get(title__icontains='update')
    return render(request, 'pages/index.html', {'pages': title})

def servers(request):
    title = Pages.objects.get(title__icontains='servers')
    return render(request, 'pages/index.html', {'pages': title})

def about(request):
    title = Pages.objects.get(title__icontains='about')
    return render(request, 'pages/index.html', {'pages': title})

def acknowledgements(request):
    title = Pages.objects.get(title__icontains='acknowledgements')
    return render(request, 'pages/index.html', {'pages': title})

def legalnotice(request):
	title = Pages.objects.get(title__icontains='legal notice')
	return render(request, 'pages/index.html', {'pages': title})

def linking(request):
    title = Pages.objects.get(title__icontains='linking')
    return render(request, 'pages/index.html', {'pages': title})