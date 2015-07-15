from django.shortcuts import render
from django.views.generic import TemplateView
from django.http import HttpResponse, JsonResponse
from django import forms

from protein.models import Gene, ProteinSegment
from structure.models import Structure
from structure.functions import CASelector, SelectionParser, GenericNumbersSelector
from structure.assign_generic_numbers_gpcr import GenericNumbering
from structure.structural_superposition import ProteinSuperpose,FragmentSuperpose
from common.views import AbsSegmentSelection,AbsReferenceSelection
from common.selection import Selection
from common.extensions import MultiFileField
from common.alignment import Alignment

import inspect
import os
import zipfile
from copy import deepcopy
from io import StringIO, BytesIO
from collections import OrderedDict
from Bio.PDB import PDBIO, PDBParser

class StructureBrowser(TemplateView):
    """
    Fetching Structure data for browser
    """

    template_name = "structure_browser.html"

    def get_context_data (self, **kwargs):

        context = super(StructureBrowser, self).get_context_data(**kwargs)
        try:
            context['crystals'] = Structure.objects.prefetch_related()
        except Structure.DoesNotExist as e:
            pass

        return context



class StructureStatistics(TemplateView):
    """
    So not ready that EA wanted to publish it.
    """

    template_name = 'structure_statistics.html'

    def get_context_data (self, **kwargs):
        context = super(StructureStatistics, self).get_context_data(**kwargs)

        #Prepare chart with unique crystallized receptors by year
        all_structs = list(Structure.objects.all())
        years = list(set([x.publication_date.year for x in all_structs]))
        unique_structs = self.get_unique_structures(all_structs)
        families = list(set([x.protein_conformation.protein.get_protein_family() for x in unique_structs]))
        
        extra = {
            'x_axis_format': '',
            'y_axis_format': 'f',
            }
        context['charttype'] = "multiBarChart"
        context['chartdata'] = self.get_per_family_cumulative_data_series(years, families, unique_structs)
        context['extra'] = extra

        context['charttype2'] = "multiBarChart"
        context['chartdata2'] = self.get_per_family_data_series(years, families, unique_structs)
        context['extra2'] = extra

        return context


    def get_unique_structures(self, structures):
        """
        Prepare a list of unique crystallized receptors. Uniqueness is evaluated by Gene object.
        """
        uniques = []
        genes = Gene.objects.all()
        for gene in genes:
            coded_proteins = [x.entry_name for x in gene.proteins.all()]
            for struct in structures:
                if struct.protein_conformation.protein.parent.entry_name in coded_proteins and struct not in uniques:
                    uniques.append(struct)
        return uniques


    def get_per_family_data_series(self, years, families, structures):
        """
        Prepare data for multiBarGraph of unique crystallized receptors. Returns data series for django-nvd3 wrapper.
        """
        series = {'x' : years,}
        data = {}
        for year in years:
            for family in families:
                if family not in data.keys():
                    data[family] = []
                count = 0
                for structure in structures:
                    if structure.protein_conformation.protein.get_protein_family() == family and structure.publication_date.year == year:
                        count += 1
                data[family].append(count)
        for idx, family in enumerate(data.keys()):
            series['name{:n}'.format(idx+1)] = family
            series['y{:n}'.format(idx+1)] = data[family]
        return series


    def get_per_family_cumulative_data_series(self, years, families, structures):
        """
        Prepare data for multiBarGraph of unique crystallized receptors. Returns data series for django-nvd3 wrapper.
        """
        series = {'x' : years,}
        data = {}
        for year in years:
            for family in families:
                if family not in data.keys():
                    data[family] = []
                count = 0
                for structure in structures:
                    if structure.protein_conformation.protein.get_protein_family() == family and structure.publication_date.year == year:
                        count += 1
                if len(data[family]) > 0:
                    data[family].append(count + data[family][-1])
                else:
                    data[family].append(count)

        for idx, family in enumerate(data.keys()):
            series['name{:n}'.format(idx+1)] = family
            series['y{:n}'.format(idx+1)] = data[family]
        return series


class GenericNumberingIndex(TemplateView):
    """
    Starting page of generic numbering assignment workflow.
    """
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
    header = "Select a file to upload:"
    upload_form_data = {
        "pdb_file": forms.FileField(),
        }
    form_code = forms.Form()
    form_code.fields = upload_form_data
    form_id = 'gn_pdb_file'
    url = '/structure/generic_numbering_results'
    mid_section = "upload_file_form.html"
    form_height = 200
    #Buttons
    buttons = {
        'continue' : {
            'label' : 'Assign generic numbers',
            'color' : 'success',
            },
        }


    def get_context_data (self, **kwargs):

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

    template_name = 'common_structural_tools.html'

    #Left panel - blank
    #Mid section
    mid_section = 'gn_results.html'
    #Buttons - none


    def post (self, request, *args, **kwargs):

        generic_numbering = GenericNumbering(StringIO(request.FILES['pdb_file'].file.read().decode('UTF-8',"ignore")))
        out_struct = generic_numbering.assign_generic_numbers()
        out_stream = StringIO()
        io = PDBIO()
        io.set_structure(out_struct)
        io.save(out_stream)
        if len(out_stream.getvalue()) > 0:
            request.session['outfile'] = { request.FILES['pdb_file'].name : out_stream, }
            self.success = True
            self.outfile = request.FILES['pdb_file'].name
            self.replacement_tag = 'GPCRDB'
        else:
            self.input_file = request.FILES['pdb_file'].name
            self.success = False

        context = super(GenericNumberingResults, self).get_context_data(**kwargs)
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        return render(request, self.template_name, context)


    def get_context_data (self, **kwargs):

        context = super(GenericNumberingResults, self).get_context_data(**kwargs)
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
    #
    upload_form_data = OrderedDict([
        ('ref_file', forms.FileField(label="Reference structure")),
        ('alt_files', MultiFileField(label="Structure(s) to superpose", max_num=10, min_num=1)),
        ('exclusive', forms.BooleanField(label='Download only superposed subset of atoms', widget=forms.CheckboxInput())),
        ])
    form_code = forms.Form()
    form_code.fields = upload_form_data
    form_id = 'superpose_files'
    url = '/structure/superposition_workflow_selection'
    mid_section = 'superposition_workflow_upload_file_form.html'

    #Buttons
    buttons = {
        'continue' : {
            'label' : 'Select segments',
            'color' : 'success',
            }
        }

    # OrderedDict to preserve the order of the boxes
    selection_boxes = OrderedDict([('reference', True),
        ('targets', True),
        ('segments', False)])

    def get_context_data (self, **kwargs):

        context = super(SuperpositionWorkflowIndex, self).get_context_data(**kwargs)

        # get selection from session and add to context
        # get simple selection from session
        simple_selection = self.request.session.get('selection', False)

        # create full selection and import simple selection (if it exists)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)

        #Clearing selections for fresh run
        print(self.kwargs)
        if 'clear' in self.kwargs.keys():
            selection.clear('reference')
            selection.clear('targets')
            selection.clear('segments')

        context['selection'] = {}
        for selection_box, include in self.selection_boxes.items():
            if include:
                context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]

        # get attributes of this class and add them to the context
        context['form_code'] = str(self.form_code)
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        return context



#Class rendering selection box for sequence segments
class SuperpositionWorkflowSelection(AbsSegmentSelection):

    #Left panel
    step = 2
    number_of_steps = 3

    #Mid section
    #mid_section = 'segment_selection.html'

    #Right panel
    segment_list = True
    buttons = {
        'continue': {
            'label': 'Superpose proteins',
            'url': '/structure/superposition_workflow_results',
            'color': 'success',
        },
    }
    # OrderedDict to preserve the order of the boxes
    selection_boxes = OrderedDict([('reference', False),
        ('targets', False),
        ('segments', True),])


    def post (self, request, *args, **kwargs):

        # create full selection and import simple selection (if it exists)
        simple_selection = request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)
        if 'exclusive' in request.POST:
            request.session['exclusive'] = True
        else:
            request.session['exclusive'] = False
        if 'ref_file' in request.FILES:
            request.session['ref_file'] = request.FILES['ref_file']
        elif selection.reference != []:
            print(selection.reference[0].item)
        if 'alt_files' in request.FILES:
            request.session['alt_files'] = request.FILES.getlist('alt_files')
        elif selection.targets != []:
            print(selection.targets)

        
        context = super(SuperpositionWorkflowSelection, self).get_context_data(**kwargs)
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

    template_name = 'common_structural_tools.html'

    #Left panel - blank
    #Mid section
    mid_section = 'superposition_results.html'
    #Buttons - none


    def get_context_data (self, **kwargs):

        context = super(SuperpositionWorkflowResults, self).get_context_data(**kwargs)
        
        simple_selection = self.request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)
        if 'ref_file' in self.request.session.keys():
            ref_file = StringIO(self.request.session['ref_file'].file.read().decode('UTF-8'))
        elif selection.reference != []:
            ref_file = StringIO(Structure.objects.get(protein_conformation__protein = selection.reference[0].item).pdb_data.pdb)
        if 'alt_files' in self.request.session.keys():
            alt_files = [StringIO(alt_file.file.read().decode('UTF-8')) for alt_file in self.request.session['alt_files']]
        elif selection.targets != []:
            alt_files = [StringIO(Structure.objects.get(protein_conformation__protein = x.item).pdb_data.pdb) for x in selection.targets]
        superposition = ProteinSuperpose(deepcopy(ref_file),alt_files, selection)
        out_structs = superposition.run()
        if len(out_structs) == 0:
            self.success = False
        elif len(out_structs) >= 1:
            io = PDBIO()            
            out_stream = BytesIO()
            zipf = zipfile.ZipFile(out_stream, 'w')

            if 'ref_file' in self.request.session.keys():
                ref_struct = PDBParser().get_structure('ref', ref_file)[0]
                ref_name = self.request.session['ref_file'].name
            elif selection.reference != []:
                ref_struct = PDBParser().get_structure('ref', StringIO(Structure.objects.get(protein_conformation__protein = selection.reference[0].item).pdb_data.pdb))[0]
                ref_name = "{!s}.pdb".format(str(selection.reference[0].item).upper())

            if 'alt_files' in self.request.session.keys():
                alt_file_names = [x.name for x in self.request.session['alt_files']]
            elif selection.targets != []:
                alt_file_names = ["{!s}_superposed.pdb".format(str(y)) for y in [Structure.objects.get(protein_conformation__protein = x.item) for x in selection.targets]]

            if self.request.session['exclusive']:
                consensus_gn_set = CASelector(SelectionParser(selection), ref_struct, out_structs).get_consensus_gn_set()
                io.set_structure(ref_struct)
                tmp = StringIO()
                io.save(tmp, GenericNumbersSelector(consensus_gn_set))
                zipf.writestr(ref_name, tmp.getvalue())
                for alt_struct, alt_file_name in zip(out_structs, alt_file_names):
                    tmp = StringIO()
                    io.set_structure(alt_struct)
                    io.save(tmp, GenericNumbersSelector(consensus_gn_set))
                    zipf.writestr(alt_file_name, tmp.getvalue())
                zipf.close()
                if len(out_stream.getvalue()) > 0:
                    self.request.session['outfile'] = { "Superposed_substructures.zip" : out_stream, }
                    self.outfile = "Superposed_substructures.zip"
                    self.success = True
                    self.zip = 'zip'
            else:
                for alt_struct, alt_file_name in zip(out_structs, alt_file_names):
                    tmp = StringIO()
                    io.set_structure(alt_struct)
                    io.save(tmp)
                    zipf.writestr(alt_file_name, tmp.getvalue())
                zipf.close()
                if len(out_stream.getvalue()) > 0:
                    self.request.session['outfile'] = { "Superposed_structures.zip" : out_stream, }
                    self.outfile = "Superposed_structures.zip"
                    self.success = True
                    self.zip = 'zip'
        
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        return context



class FragmentSuperpositionIndex(TemplateView):

    template_name = 'common_structural_tools.html'

    #Left panel
    step = 1
    number_of_steps = 1
    title = "SUPERPOSE FRAGMENTS OF CRYSTAL STRUCTURES"
    description = """
    The tool implements a fragment-based pharmacophore method, as published in <a href='http://www.ncbi.nlm.nih.gov/pubmed/25286328'>Fidom K, et al (2015)</a>. Interacting ligand moiety - residue pairs extracted from selected crystal structures of GPCRs are superposed onto the input pdb file based on gpcrdb generic residue numbers. Resulting aligned ligand fragments can be used for placement of pharmacophore features.

    Upload a pdb file you want to superpose the interacting moiety - residue pairs.
    
    Once you have selected all your targets, click the green button.
        """

    #Input file form data
    header = "Select a file to upload:"
    upload_form_data = OrderedDict([
        ("pdb_file", forms.FileField()),
        ("similarity", forms.ChoiceField(choices=(('identical','Use fragments with identical residues'),
                     ('similar','Use fragments with residues of similar properties')),
            widget=forms.RadioSelect())),
        ("representative", forms.ChoiceField(choices=(('closest','Use fragments from the evolutionary closest crystal structure'),
                     ('any','Use all available fragments')), widget=forms.RadioSelect())),
        ])
    form_code = forms.Form()
    form_code.fields = upload_form_data
    form_code.initial={'similarity': 'similar', 'representative': 'closest'}
    form_id = 'fragments'
    url = '/structure/fragment_superposition_results'
    mid_section = "upload_file_form.html"
    form_height = 350

    #Buttons
    buttons = {
        'continue' : {
            'label' : 'Retrieve fragments',
            'color' : 'success',
            },
        }


    def get_context_data (self, **kwargs):

        context = super(FragmentSuperpositionIndex, self).get_context_data(**kwargs)
        # get attributes of this class and add them to the context
        context['form_code'] = str(self.form_code)
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        return context



class FragmentSuperpositionResults(TemplateView):

    template_name = "common_structural_tools.html"

    #Left panel - blank
    #Mid section
    mid_section = 'superposition_results.html'
    #Buttons - none

    def post (self, request, *args, **kwargs):
        
        frag_sp = FragmentSuperpose(StringIO(request.FILES['pdb_file'].file.read().decode('UTF-8', 'ignore')),request.FILES['pdb_file'].name)
        superposed_fragments = []
        superposed_fragments_repr = []
        print(request.POST)
        if request.POST['similarity'] == 'identical':
            if request.POST['representative'] == 'any':
                superposed_fragments = frag_sp.superpose_fragments()
            else:
                superposed_fragments_repr = frag_sp.superpose_fragments(representative=True)
                superposed_fragments = frag_sp.superpose_fragments()
        else:
            if request.POST['representative'] == 'any':
                superposed_fragments = frag_sp.superpose_fragments(use_similar=True)
            else:
                superposed_fragments_repr = frag_sp.superpose_fragments(representative=True, use_similar=True)
                superposed_fragments = frag_sp.superpose_fragments(use_similar=True)
        if superposed_fragments == []:
            self.message = "No fragments were aligned."
        else:
            io = PDBIO()
            out_stream = BytesIO()
            zipf = zipfile.ZipFile(out_stream, 'a')
            for fragment, pdb_data in superposed_fragments:
                io.set_structure(pdb_data)
                tmp = StringIO()
                io.save(tmp)
                if request.POST['representative'] == 'any':
                    zipf.writestr(fragment.generate_filename(), tmp.getvalue())
                else:
                    zipf.writestr("all_fragments//{!s}".format(fragment.generate_filename()), tmp.getvalue())
            if superposed_fragments_repr != []:
                for fragment, pdb_data in superposed_fragments_repr:
                    io.set_structure(pdb_data)
                    tmp = StringIO()
                    io.save(tmp)
                    zipf.writestr("representative_fragments//{!s}".format(fragment.generate_filename()), tmp.getvalue())
            zipf.close()
            if len(out_stream.getvalue()) > 0:
                request.session['outfile'] = { 'interacting_moiety_residue_fragments.zip' : out_stream, }
                self.outfile = 'interacting_moiety_residue_fragments.zip'
                self.success = True
                self.zip = 'zip'
                self.message = '{:n} fragments were superposed.'.format(len(superposed_fragments))

        context = super(FragmentSuperpositionResults, self).get_context_data(**kwargs)
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        return render(request, self.template_name, context)
       

#==============================================================================
class TemplateTargetSelection(AbsReferenceSelection):
    """
    Starting point for template selection workflow. Target selection. 
    """
    
    type_of_selection = 'reference'
    #Left panel
    description = 'Select targets by searching or browsing in the middle column. You can select entire target families or individual targets.\n\nSelected targets will appear in the right column, where you can edit the list.\n\nOnce you have selected all your targets, either proceed with all TMs alignment ("Find template" button) or specify the sequence segments manualy ("Advanced segment selection" button).'
    step = 1
    number_of_steps = 2

    #Mid section

    #Right panel
    buttons = {
        'continue': {
            'label': 'Find template',
            'url': '/structure/template_browser',
            'color': 'success',
        },
        'segments' : {
            'label' : 'Advanced segment selection',
            'url' : '/structure/template_segment_selection',
            'color' : 'info',
        },
    }



#==============================================================================
class TemplateSegmentSelection(AbsSegmentSelection):
    """
    Advanced selection of sequence segments for template search.
    """
   #Left panel
    step = 2
    number_of_steps = 2

    #Mid section
    #mid_section = 'segment_selection.html'

    #Right panel
    segment_list = True
    buttons = {
        'continue': {
            'label': 'Find template',
            'url': '/structure/template_browser',
            'color': 'success',
        },
    }
    # OrderedDict to preserve the order of the boxes
    selection_boxes = OrderedDict([('reference', True),
        ('targets', False),
        ('segments', True),])    



#==============================================================================
class TemplateBrowser(TemplateView):
    """
    Fetching Structure data and ordering by similarity
    """

    template_name = "structure_browser.html"

    def get_context_data (self, **kwargs):

        context = super(TemplateBrowser, self).get_context_data(**kwargs)

        # get simple selection from session
        simple_selection = self.request.session.get('selection', False)
        a = Alignment()
        a.load_reference_protein_from_selection(simple_selection)
        if simple_selection.segments != []:
            a.load_segments_from_selection(simple_selection)
        else:
            a.load_segments(ProteinSegment.objects.filter(slug__in=['TM1', 'TM2', 'TM3', 'TM4','TM5','TM6', 'TM7']))
        a.load_proteins([x.protein_conformation.protein for x in list(Structure.objects.all())])
        a.build_alignment()
        a.calculate_similarity()
        print(len(a.proteins[1:]))
        context['crystals'] = []
        for prot in a.proteins[1:]:
            context['crystals'].append(Structure.objects.get(protein_conformation__protein__entry_name=prot.protein.entry_name))
        print(context['crystals'])
        #try:
        #    context['crystals'] = Structure.objects.all()
        #except Structure.DoesNotExist as e:
        #    pass

        return context


#==============================================================================
def ServePdbOutfile (request, outfile, replacement_tag):
    
    root, ext = os.path.splitext(outfile)
    out_stream = request.session['outfile'][outfile]
    response = HttpResponse(content_type="chemical/x-pdb")
    response['Content-Disposition'] = 'attachment; filename="{}_{}.pdb"'.format(root, replacement_tag)
    response.write(out_stream.getvalue())

    return response


def ServeZipOutfile (request, outfile):
    
    out_stream = request.session['outfile'][outfile]
    response = HttpResponse(content_type="application/zip")
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(outfile)
    response.write(out_stream.getvalue())

    return response