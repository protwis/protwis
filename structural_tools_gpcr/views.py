from django.http import HttpResponse
from django.shortcuts import render
from django.views.generic import TemplateView
from django import forms

#Class for starting page of generic numbers assignment
class GenericNumberingStart(TemplateView):
    template_name = 'upload_file_form'
    upload_form_data = {
        "pdb_file": forms.FileField(),
        }
    upload_file_form = forms.Form(upload_form_data)


class FileUploadForm(Form):
    
    pass


def UploadFile(request):
    
    pass