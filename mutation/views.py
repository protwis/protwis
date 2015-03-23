from django.shortcuts import render
from django.http import HttpResponse
from mutation.testing import *

from mutation.models import *
from datetime import datetime
#env/bin/python3 -m pip install xlrd


# Create your views here.
def index(request):
    return HttpResponse("Hello, world. You're at the polls index.")

def importmutation(request):

	rows = loaddatafromexcel('/vagrant/shared/protwis/mutation/import.xlsx')

	rows = analyse_rows(rows)

	whattoreturn = []

	for r in rows:
		raw_id = insert_raw(r)



	return HttpResponse(raw_id)