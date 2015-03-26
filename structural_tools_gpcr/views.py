from django.views.generic import TemplateView
from django.http import HttpResponse
from django import forms
from django.shortcuts import render
from django.core.context_processors import csrf
import structural_tools_gpcr.assign_generic_numbers as gn
import inspect, os
from io import StringIO
from Bio.PDB import PDBIO


#===================Assignment of generic numbers==============================
#Class for starting page of generic numbers assignment
class GenericNumberingIndex(TemplateView):

    template_name = 'common_structural_tools.html'
    
    #Left panel
    step = 1
    number_of_steps = 1
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
            'color' : 'success',
            },
        }


    def get_context_data(self, **kwargs):

        context = super(GenericNumberingIndex, self).get_context_data(**kwargs)
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


    def post(self, request, *args, **kwargs):

        generic_numbering = gn.GenericNumbering(StringIO(request.FILES['pdb_file'].file.read().decode('UTF-8')))
        out_struct = generic_numbering.assign_generic_numbers()
        out_stream = StringIO()
        io = PDBIO()
        io.set_structure(out_struct)
        io.save(out_stream)
        request.session['outfile'] = { request.FILES['pdb_file'].name : out_stream, }
        self.input_file = request.FILES['pdb_file'].name

        context =  super(GenericNumberingResults, self).get_context_data(**kwargs)
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        return render(request, self.template_name, context)


    def get_context_data(self, **kwargs):

        context =  super(GenericNumberingResults, self).get_context_data(**kwargs)
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        return context

#==============================================================================

#========================Superposition of structures===========================
#Class for starting page of superposition workflow
class SuperpositionWorkflowIndex(TemplateView):
    
    template_name = "common_structural_tools.html"

    pass


#==============================================================================

def ServeOutfile(request, outfile):
    
    root, ext = os.path.splitext(outfile)
    out_stream = request.session['outfile'][outfile]

    response = HttpResponse(mimetype="chemical/x-pdb")
    response['Content-Disposition'] = 'attachment; filename="{}_GPCRDB.pdb"'.format(root)
    response.write(out_stream.get_value())

    return response