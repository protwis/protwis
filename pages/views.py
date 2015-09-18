from django.shortcuts import render
from django.conf import settings
from django.views import generic
from django.views.static import serve
from django.core.files import File 

from pages.models import Pages
import os

def index(request):
	return render(request,'pages/index.html', {'site_title': settings.SITE_TITLE})

def contact(request):
    title = Pages.objects.get(title__icontains='contact')
    return render(request, 'pages/index.html', {'pages': title})

def citing(request):
    title = Pages.objects.get(title__icontains='citing')
    return render(request, 'pages/index.html', {'pages': title})

def contributors(request):
    title = Pages.objects.get(title__icontains='contributors')
    return render(request, 'pages/index.html', {'pages': title})

def contribute(request):
    title = Pages.objects.get(title__icontains='contribute')
    return render(request, 'pages/index.html', {'pages': title})

def committees(request):
    title = Pages.objects.get(title__icontains='committees')
    return render(request, 'pages/index.html', {'pages': title})

def meetings(request):
    title = Pages.objects.get(title__icontains='meetings')
    return render(request, 'pages/index.html', {'pages': title})

def update(request):
    title = Pages.objects.get(title__icontains='update')
    return render(request, 'pages/index.html', {'pages': title})

def poster(request):
    filepath = open('/static/home/poster/GPCRDB_Poster.pdf',"r")
    django_file = File('/static/home/poster/GPCRDB_Poster.pdf') 
    t = loader.get_gemplate('pages/index.html') 
    c = Context({'file':django_file}) 
    return HttpResponse(t.render(c)) 

def servers(request):
	title = Pages.objects.get(title__icontains='servers')
	return render(request, 'pages/index.html', {'pages': title})