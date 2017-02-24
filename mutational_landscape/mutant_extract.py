from django.shortcuts import render
from django.template import loader, Context
from django.db.models import Count, Min, Sum, Avg, Q
from django.http import HttpResponse
from django.conf import settings
from django.core.cache import cache
from django.views.decorators.cache import cache_page
from mutation.functions import *
from mutation.models import *

mutations = MutationExperiment.objects.filter(protein__entry_name='adrb2_human').order_by('residue__sequence_number').prefetch_related('residue')

print(mutations)