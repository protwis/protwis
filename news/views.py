from django.shortcuts import render
from django.conf import settings
from django.views import generic
from django.core import serializers

from news.models import News

import json, os

def index(request):
    data = serializers.serialize( "python", News.objects.all().order_by('-date'), fields=('date','image','html'))
    return render(request, 'news/index.html', {'news': data})