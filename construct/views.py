from django.shortcuts import render

from structure.models import Construct

from datetime import datetime


# Create your views here.
def constructs(request):
    construct=Construct.objects.all()



    context = {'constructs': construct}

    return render(request,'structure/index_construct.html',context)

