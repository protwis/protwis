from django.http import HttpResponse
from django.shortcuts import render
from django.views.generic import TemplateView
from django import forms

#Class for starting page of generic numbers assignment
class GenericNumberingStart(TemplateView):
    template_name = 'generic_numbering.html'
    
    #Left panel
    step = 1
    number_of_steps = 2
    title = "UPLOAD A PDB FILE"
    description = """
    Upload a pdb file you want to be annotated with generic numbers. Note that "CA" atoms will be assigned a number in GPCRdb notation, and "N" atoms will be annotated with Ballesteros-Weinstein scheme.
    
    Once you have selected all your targets, click the green button.
        """

    #Input file form data
    header = "Upload your pdb file"
    upload_form_data = {
        "pdb_file": forms.FileField(),
        }
    upload_file_form = forms.Form(upload_form_data)
    form_id = 'gn_pdb_file'
    url = ''
    
    #Buttons
