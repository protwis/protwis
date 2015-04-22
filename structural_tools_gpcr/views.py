from django.views.generic import TemplateView
from django.http import HttpResponse
from django import forms
from django.shortcuts import render
from django.core.context_processors import csrf

from structural_tools_gpcr.assign_generic_numbers import GenericNumbering
from structural_tools_gpcr.structural_superposition import ProteinSuperpose
from common.views import AbsSegmentSelection
from common.selection import Selection

import inspect, os
from io import StringIO
from collections import OrderedDict

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
    header = "Upload your pdb file:"
    upload_form_data = {
        "pdb_file": forms.FileField(),
        }
    form_code = forms.Form()
    form_code.fields = upload_form_data
    form_id = 'gn_pdb_file'
    url = '/structural_tools_gpcr/generic_numbering_results'
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

        generic_numbering = GenericNumbering(StringIO(request.FILES['pdb_file'].file.read().decode('UTF-8')))
        out_struct = generic_numbering.assign_generic_numbers()
        out_stream = StringIO()
        io = PDBIO()
        io.set_structure(out_struct)
        io.save(out_stream)
        if len(out_stream.getvalue()) > 0:
            request.session['outfile'] = { request.FILES['pdb_file'].name : out_stream, }
            self.success = True
            self.outfile = request.FILES['pdb_file'].name
        else:
            self.input_file = request.FILES['pdb_file'].name
            self.success = False

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

    #Left panel
    step = 1
    number_of_steps = 3
    title = "UPLOAD YOUR FILES"
    description = """
    Upload a pdb file for reference structure, and one or more files that will be superposed. You can also select the structures from crystal structure browser.

    Once you have uploaded/selected all your targets, click the green button.
    """

    header = "Upload or select your structures:"
    upload_form_data = {
        'ref_file' : forms.FileField(label="Reference structure"),
        'alt_files' : forms.FileField(label="Structure(s) to superpose"),
        }
    form_code = forms.Form()
    form_code.fields = upload_form_data
    form_id = 'superpose_files'
    url = '/structural_tools_gpcr/superposition_workflow_selection'
    mid_section = 'upload_file_form.html'

    #Buttons
    buttons = {
        'continue' : {
            'label' : 'Select segments',
            'color' : 'success',
            }
        }


    def get_context_data(self, **kwargs):

        context = super(SuperpositionWorkflowIndex, self).get_context_data(**kwargs)
        # get attributes of this class and add them to the context
        context['form_code'] = str(self.form_code)
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        return context



#Class rendering selection box for sequence segments
class SuperpositionWorkflowSelection(AbsSegmentSelection):

    template_name='common_structural_tools.html'

    #Left panel
    step = 2
    number_of_steps = 3

    #Mid section
    mid_section = 'common/segmentselection.html'

    #Right panel
    segment_list = True
    buttons = {
        'continue': {
            'label': 'Superpose proteins',
            'url': '/structural_tools_gpcr/superposition_workflow_results',
            'color': 'success',
        },
    }
    # OrderedDict to preserve the order of the boxes
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', False),
        ('segments', True),
    ])

    def post(self, request, *args, **kwargs):

        request.session['ref_file'] = request.FILES['ref_file']
        request.session['alt_files'] = request.FILES['alt_files']
        simple_selection = request.session.get('selection', False)

        # create full selection and import simple selection (if it exists)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)
        
        context =  super(SuperpositionWorkflowSelection, self).get_context_data(**kwargs)
        context['selection'] = {}
        for selection_box, include in self.selection_boxes.items():
            if include:
                context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]


        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        return render(request, self.template_name, context)






#Class rendering results from superposition workflow
class SuperpositionWorkflowResults(TemplateView):

    template_name='common_structural_tools.html'

    #Left panel - blank
    #Mid section
    mid_section = 'superposition_workflow_results.html'
    #Buttons - none


    #def post(self, request, *args, **kwargs):
    def get_context_data(self, **kwargs):

        context =  super(SuperpositionWorkflowResults, self).get_context_data(**kwargs)
        
        simple_selection = self.request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)

        superposition = ProteinSuperpose(StringIO(self.request.session['ref_file'].file.read().decode('UTF-8')),[StringIO(self.request.session['alt_files'].file.read().decode('UTF-8'))], selection)

        out_structs = superposition.run()
        if len(out_structs) == 0:
            self.success = False
        elif len(out_structs) == 1:
            out_stream = StringIO()
            io = PDBIO()
            io.set_structure(out_structs[0])
            io.save(out_stream)
            if len(out_stream.getvalue()) > 0:
                self.session['outfiles'] = { request.FILES['alt_files'].name : out_stream, }
                #self.input_file = request.FILES['pdb_file'].name
                self.success = True
                self.outfile = request.FILES['alt_files'].name

        
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        #return render(self.request, self.template_name, context) 
        return context

       

#==============================================================================

def ServePdbOutfile(request, outfile):
    
    root, ext = os.path.splitext(outfile)
    out_stream = request.session['outfile'][outfile]
    print(request.session['outfile'][outfile])

    response = HttpResponse(content_type="chemical/x-pdb")
    response['Content-Disposition'] = 'attachment; filename="{}_GPCRDB.pdb"'.format(root)
    response.write(out_stream.getvalue())

    return response