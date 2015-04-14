from django.shortcuts import render
from django.http import HttpResponse
from interaction.functions import *
# Create your views here.

def index(request):
    #return HttpResponse("Hello, world. You're at the polls index.")
	context = {}
	#return HttpResponse(ref_id)
	return render(request,'interaction/index.html',context)

def calculate(request):		
    return HttpResponse("You pressed!")
	#context = {}
	#return HttpResponse(ref_id)
	#return render(request,'interaction/index.html',context)
