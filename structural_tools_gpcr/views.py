from django.views.generic import TemplateView
from django import forms
from django.core.context_processors import csrf


import inspect

#Class for starting page of generic numbers assignment
class GenericNumberingStart(TemplateView):
    template_name = 'common_structural_tools.html'
    
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
    form_code = forms.Form()
    form_code.fields = upload_form_data
    form_id = 'gn_pdb_file'
    url = '/structural_tools_gpcr/gn_results'
    mid_section = "upload_file_form.html"

    #Buttons
    buttons = {
        'continue' : {
            'label' : 'Assign generic numbers',
            #'url' : '',#'/structural_tools_gpcr/gn_results',
            'color' : 'success',
            },
        }

    def get_context_data(self, **kwargs):
        context = super(GenericNumberingStart, self).get_context_data(**kwargs)
        # get attributes of this class and add them to the context
        context['form_code'] = str(self.form_code)
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]
        return context

#Class rendering results from generic numbers assignment
class GenericNumberingResults(TemplateView):
    template_name='common_structural_tools.html'

    #Left panel - blank
    #Mid section
    mid_section = 'gn_results.html'
    #Buttons - none

    def get_context_data(self, **kwargs):
        context =  csrf(super(GenericNumberingResults, self).get_context_data(**kwargs))

        if self.request.method == 'POST':
            pdb_upload = kwargs['pdb_file']
        context['debug_data'] = kwargs['pdb_file']

        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]
        return context
