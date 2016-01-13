from django.shortcuts import render
from django.conf import settings
from django.views.generic import TemplateView, View
from django.http import HttpResponse, JsonResponse, HttpResponseRedirect
from django.db.models import Count, Q, Prefetch
from django import forms

from protein.models import Gene, ProteinSegment
from structure.models import Structure
from structure.functions import CASelector, SelectionParser, GenericNumbersSelector, SubstructureSelector, check_gn
from structure.assign_generic_numbers_gpcr import GenericNumbering
from structure.structural_superposition import ProteinSuperpose,FragmentSuperpose
from interaction.models import ResidueFragmentInteraction,StructureLigandInteraction
from protein.models import Protein
from common.views import AbsSegmentSelection,AbsReferenceSelection
from common.selection import Selection, SelectionItem
from common.extensions import MultiFileField
Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')

import inspect
import os
import zipfile
import math
import json
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
            context['structures'] = Structure.objects.all().select_related(
                "pdb_code__web_resource",
                "protein_conformation__protein__species",
                "protein_conformation__protein__source",
                "protein_conformation__protein__family__parent__parent__parent",
                "publication__web_link__web_resource").prefetch_related(
                "stabilizing_agents",
                "protein_conformation__protein__parent__endogenous_ligands__properities__ligand_type",
                Prefetch("ligands", queryset=StructureLigandInteraction.objects.filter(
                annotated=True).prefetch_related('ligand__properities__ligand_type', 'ligand_role')))
        except Structure.DoesNotExist as e:
            pass

        return context

def StructureDetails(request, pdbname):
    """
    Show structure details
    """
    pdbname = pdbname
    structures = ResidueFragmentInteraction.objects.values('structure_ligand_pair__ligand__name','structure_ligand_pair__pdb_reference','structure_ligand_pair__annotated').filter(structure_ligand_pair__structure__pdb_code__index=pdbname).annotate(numRes = Count('pk', distinct = True)).order_by('-numRes')
    resn_list = ''

    for structure in structures:
        if structure['structure_ligand_pair__annotated']:
            resn_list += ",\""+structure['structure_ligand_pair__pdb_reference']+"\""
    print(resn_list)

    crystal = Structure.objects.get(pdb_code__index=pdbname)
    p = Protein.objects.get(protein=crystal.protein_conformation.protein)
    residues = ResidueFragmentInteraction.objects.filter(structure_ligand_pair__structure__pdb_code__index=pdbname, structure_ligand_pair__annotated=True).order_by('rotamer__residue__sequence_number')
    return render(request,'structure_details.html',{'pdbname': pdbname, 'structures': structures, 'crystal': crystal, 'protein':p, 'residues':residues, 'annotated_resn': resn_list})

def ServePdbDiagram(request, pdbname):       
    structure=Structure.objects.filter(pdb_code__index=pdbname) 
    if structure.exists():
        structure=structure.get()
    else:
         quit() #quit!

    if structure.pdb_data is None:
        quit()

    response = HttpResponse(structure.pdb_data.pdb, content_type='text/plain')
    return response

    
def ServePdbLigandDiagram(request,pdbname,ligand):      
    pair = StructureLigandInteraction.objects.filter(structure__pdb_code__index=pdbname).filter(Q(ligand__properities__inchikey=ligand) | Q(ligand__name=ligand)).exclude(pdb_file__isnull=True).get()
    response = HttpResponse(pair.pdb_file.pdb, content_type='text/plain')
    return response

class StructureStatistics(TemplateView):
    """
    So not ready that EA wanted to publish it.
    """
    template_name = 'structure_statistics.html'

    def get_context_data (self, **kwargs):
        context = super().get_context_data(**kwargs)

        #Prepare chart with unique crystallized receptors by year
        all_structs = list(Structure.objects.all().prefetch_related('protein_conformation__protein'))
        years = self.get_years_range(list(set([x.publication_date.year for x in all_structs])))
        unique_structs = list(Structure.objects.order_by('protein_conformation__protein__parent', 'state',
            'publication_date', 'resolution').distinct('protein_conformation__protein__parent').prefetch_related('protein_conformation__protein'))
        families = list(set([x.protein_conformation.protein.get_protein_family() for x in unique_structs]))
        classes = [x.protein_conformation.protein.get_protein_class() for x in unique_structs]
        
        tmp = OrderedDict()
        for x in sorted(list(set(classes))):
            tmp[x] = classes.count(x)
        #Basic stats
        context['all_structures'] = len(all_structs)
        context['unique_structures'] = len(unique_structs)
        context['unique_by_class'] = tmp
        context['unique_complexes'] = len(StructureLigandInteraction.objects.filter(annotated=True).distinct('structure__protein_conformation__protein__family__name', 'ligand__name'))

        context['chartdata'] = self.get_per_family_cumulative_data_series(years, families, unique_structs)
        context['chartdata_y'] = self.get_per_family_data_series(years, families, unique_structs)
        context['chartdata_all'] = self.get_per_family_cumulative_data_series(years, families, all_structs)
        context['chartdata_reso'] = self.get_resolution_coverage_data_series(all_structs)

        return context


    def get_years_range(self, years_list):

        min_y = min(years_list)
        max_y = max(years_list)
        return range(min_y, max_y+1)


    def get_per_family_data_series(self, years, families, structures):
        """
        Prepare data for multiBarGraph of unique crystallized receptors. Returns data series for django-nvd3 wrapper.
        """
        series = []
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
        for family in families:
            series.append({"values": 
                [{
                    'x': years[i],
                    'y': j
                    } for i, j in enumerate(data[family])],
                "key": family,
                "yAxis": "1"})
        return json.dumps(series)


    def get_per_family_cumulative_data_series(self, years, families, structures):
        """
        Prepare data for multiBarGraph of unique crystallized receptors. Returns data series for django-nvd3 wrapper.
        """
        series = []
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
        for family in families:
            series.append({"values": 
                [{
                    'x': years[i],
                    'y': j
                    } for i, j in enumerate(data[family])],
                "key": family,
                "yAxis": "1"})
        return json.dumps(series)


    def get_resolution_coverage_data_series(self, structures):
        """
        Prepare data for multiBarGraph of resolution coverage of available crystal structures.
        """
        #Resolutions boundaries
        reso_min = float(min([round(x.resolution, 1) for x in structures]))
        reso_max = float(max([round(x.resolution, 1) for x in structures]))
        step = (reso_max - reso_min)/10

        brackets = [reso_min + step*x for x in range(10)] + [reso_max]

        reso_count = []
        bracket_labels = []
        for idx, bracket in enumerate(brackets):
            if idx == 0:
                reso_count.append(len([x for x in structures if x.resolution <= bracket]))
                bracket_labels.append('< {:.1f}'.format(bracket))
            else:
                reso_count.append(len([x for x in structures if bracket-step < x.resolution <= bracket]))
                bracket_labels.append('{:.1f}-{:.1f}'.format(brackets[idx-1],bracket))
        
        return json.dumps([{"values": [{
                    'x': bracket_labels[i],
                    'y': j
                    } for i, j in enumerate(reso_count)],
                "key": 'Resolution coverage',
                "yAxis": "1"}])
         
            



class GenericNumberingIndex(TemplateView):
    """
    Starting page of generic numbering assignment workflow.
    """
    template_name = 'common_structural_tools.html'
    
    #Left panel
    step = 1
    number_of_steps = 2
    documentation_url = settings.DOCUMENTATION_URL
    docs = 'structures.html#pdb-file-residue-numbering'
    title = "UPLOAD A PDB FILE"
    description = """
    Upload a pdb file to be annotated with generic numbers from GPCRdb.

    The numbers can be visualized in molecular viewers such as PyMOL, with scripts available with the output files.
    
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

    #Left panel
    step = 1
    number_of_steps = 2
    title = "SELECT SUBSTRUCTURE"
    description = 'Download the desired substructures.'
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
            request.session['gn_outfile'] = out_stream
            request.session['gn_outfname'] = request.FILES['pdb_file'].name
            self.success = True
        else:
            self.input_file = request.FILES['pdb_file'].name
            self.success = False

        context = super(GenericNumberingResults, self).get_context_data(**kwargs)
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        return render(request, self.template_name, context)



class GenericNumberingSelection(AbsSegmentSelection):
    """
    Segment selection for download of annotated substructure.
    """

    step = 2
    number_of_steps = 2

    docs = 'structures.html#pdb-file-residue-numbering'

    #Mid section
    #mid_section = 'segment_selection.html'

    #Right panel
    segment_list = True
    buttons = {
        'continue': {
            'label': 'Download substructure',
            'url': '/structure/generic_numbering_results/substr',
            'color': 'success',
        },
    }
    # OrderedDict to preserve the order of the boxes
    selection_boxes = OrderedDict([('reference', False),
        ('targets', False),
        ('segments', True),])


    def get_context_data(self, **kwargs):
        
        context = super(GenericNumberingSelection, self).get_context_data(**kwargs)

        simple_selection = self.request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)

        context['selection'] = {}
        context['selection']['site_residue_groups'] = selection.site_residue_groups
        context['selection']['active_site_residue_group'] = selection.active_site_residue_group
        for selection_box, include in self.selection_boxes.items():
            if include:
                context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]

        # get attributes of this class and add them to the context
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]
        return context


class GenericNumberingDownload(View):
    """
    Serve the (sub)structure depending on user's choice.
    """
    def get(self, request, *args, **kwargs):

        if self.kwargs['substructure'] == 'custom':
            return HttpResponseRedirect('/structure/generic_numbering_selection')

        simple_selection = self.request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)
        out_stream = StringIO()
        io = PDBIO()
        request.session['gn_outfile'].seek(0)
        gn_struct = PDBParser(PERMISSIVE=True, QUIET=True).get_structure(request.session['gn_outfname'], request.session['gn_outfile'])[0]

        if self.kwargs['substructure'] == 'full':
            io.set_structure(gn_struct)
            io.save(out_stream)

        if self.kwargs['substructure'] == 'substr':
            io.set_structure(gn_struct)
            io.save(out_stream, GenericNumbersSelector(parsed_selection=SelectionParser(selection)))

        root, ext = os.path.splitext(request.session['gn_outfname'])
        response = HttpResponse(content_type="chemical/x-pdb")
        response['Content-Disposition'] = 'attachment; filename="{}_GPCRDB.pdb"'.format(root)
        response.write(out_stream.getvalue())

        return response


#==============================================================================

#========================Superposition of structures===========================
#Class for starting page of superposition workflow
class SuperpositionWorkflowIndex(TemplateView):
    
    template_name = "common_structural_tools.html"

    #Left panel
    step = 1
    number_of_steps = 3
    documentation_url = settings.DOCUMENTATION_URL
    docs = 'structures.html#structure-superposition'
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
        #('exclusive', forms.BooleanField(label='Download only superposed subset of atoms', widget=forms.CheckboxInput())),
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
        if 'clear' in self.kwargs.keys():
            selection.clear('reference')
            selection.clear('targets')
            selection.clear('segments')
            if 'alt_files' in self.request.session.keys():
                del self.request.session['alt_files']
            if 'ref_file' in self.request.session.keys():
                del self.request.session['ref_file']
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

    docs = 'structures.html#structure-superposition'

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

        if 'ref_file' in request.FILES:
            request.session['ref_file'] = request.FILES['ref_file']
        if 'alt_files' in request.FILES:
            request.session['alt_files'] = request.FILES.getlist('alt_files')
        
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

    def get_context_data(self, **kwargs):

        self.buttons = {
            'continue': {
                'label': 'Download substructures',
                'url': '/structure/superposition_workflow_results/custom',
                'color': 'success',
            },
        }
        context = super(SuperpositionWorkflowSelection, self).get_context_data(**kwargs)

        simple_selection = self.request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)

        context['selection'] = {}
        context['selection']['site_residue_groups'] = selection.site_residue_groups
        context['selection']['active_site_residue_group'] = selection.active_site_residue_group
        for selection_box, include in self.selection_boxes.items():
            if include:
                context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]

        # get attributes of this class and add them to the context
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]
        return context


#Class rendering results from superposition workflow
class SuperpositionWorkflowResults(TemplateView):
    """
    Select download mode for the superposed structures. Full structures, superposed fragments only, select substructure.
    """

    template_name = 'common_structural_tools.html'

    #Left panel
    step = 3
    number_of_steps = 3
    title = "SELECT SUBSTRUCTURE"
    description = 'Download the desired substructures.'

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
            ref_file = StringIO(selection.reference[0].item.get_cleaned_pdb())
        if 'alt_files' in self.request.session.keys():
            alt_files = [StringIO(alt_file.file.read().decode('UTF-8')) for alt_file in self.request.session['alt_files']]
        elif selection.targets != []:
            alt_files = [StringIO(x.item.get_cleaned_pdb()) for x in selection.targets if x.type == 'structure']
        superposition = ProteinSuperpose(deepcopy(ref_file),alt_files, selection)
        out_structs = superposition.run()
        if 'alt_files' in self.request.session.keys():
            alt_file_names = [x.name for x in self.request.session['alt_files']]
        else:
            alt_file_names = ['{}_{}.pdb'.format(x.item.protein_conformation.protein.parent.entry_name, x.item.pdb_code.index) for x in selection.targets if x.type == 'structure']
        if len(out_structs) == 0:
            self.success = False
        elif len(out_structs) >= 1:
            io = PDBIO()            
            self.request.session['alt_structs'] = {}
            for alt_struct, alt_file_name in zip(out_structs, alt_file_names):
                tmp = StringIO()
                io.set_structure(alt_struct)
                io.save(tmp)
                self.request.session['alt_structs'][alt_file_name] = tmp

            self.success = True
        
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        return context

class SuperpositionWorkflowDownload(View):
    """
    Serve the (sub)structures depending on user's choice.
    """

    def get(self, request, *args, **kwargs):

        if self.kwargs['substructure'] == 'select':
            return HttpResponseRedirect('/structure/superposition_workflow_selection')

        io = PDBIO()            
        out_stream = BytesIO()
        zipf = zipfile.ZipFile(out_stream, 'w')
        simple_selection = self.request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)
        self.alt_substructure_mapping = {}
        #reference
        if 'ref_file' in request.session.keys():
            self.request.session['ref_file'].file.seek(0)
            ref_struct = PDBParser(PERMISSIVE=True, QUIET=True).get_structure('ref', StringIO(self.request.session['ref_file'].file.read().decode('UTF-8')))[0]
            gn_assigner = GenericNumbering(structure=ref_struct)
            gn_assigner.assign_generic_numbers()
            self.ref_substructure_mapping = gn_assigner.get_substructure_mapping_dict()
            ref_name = self.request.session['ref_file'].name
        elif selection.reference != []:
            ref_struct = PDBParser(PERMISSIVE=True, QUIET=True).get_structure('ref', StringIO(selection.reference[0].item.get_cleaned_pdb()))[0]
            gn_assigner = GenericNumbering(structure=ref_struct)
            gn_assigner.assign_generic_numbers()
            self.ref_substructure_mapping = gn_assigner.get_substructure_mapping_dict()
            ref_name = '{}_{}_ref.pdb'.format(selection.reference[0].item.protein_conformation.protein.parent.entry_name, selection.reference[0].item.pdb_code.index)

        alt_structs = {}
        for alt_id, st in self.request.session['alt_structs'].items():
            st.seek(0)
            alt_structs[alt_id] = PDBParser(PERMISSIVE=True, QUIET=True).get_structure(alt_id, st)[0]
            gn_assigner = GenericNumbering(structure=alt_structs[alt_id])
            gn_assigner.assign_generic_numbers()
            self.alt_substructure_mapping[alt_id] = gn_assigner.get_substructure_mapping_dict()

        if self.kwargs['substructure'] == 'full':

            io.set_structure(ref_struct)
            tmp = StringIO()
            io.save(tmp)
            zipf.writestr(ref_name, tmp.getvalue())

            for alt_name in self.request.session['alt_structs']:
                tmp = StringIO()
                io.set_structure(alt_structs[alt_name])
                io.save(tmp)
                zipf.writestr(alt_name, tmp.getvalue())

        elif self.kwargs['substructure'] == 'substr':

            consensus_gn_set = CASelector(SelectionParser(selection), ref_struct, alt_structs.values()).get_consensus_gn_set()
            io.set_structure(ref_struct)
            tmp = StringIO()
            io.save(tmp, GenericNumbersSelector(consensus_gn_set))
            zipf.writestr(ref_name, tmp.getvalue())
            for alt_name in self.request.session['alt_structs']:
                tmp = StringIO()
                io.set_structure(alt_structs[alt_name])
                io.save(tmp, GenericNumbersSelector(consensus_gn_set))
                zipf.writestr(alt_name, tmp.getvalue())

        elif self.kwargs['substructure'] == 'custom':

            io.set_structure(ref_struct)
            tmp = StringIO()
            io.save(tmp, SubstructureSelector(self.ref_substructure_mapping, parsed_selection=SelectionParser(selection)))
            zipf.writestr(ref_name, tmp.getvalue())
            for alt_name in self.request.session['alt_structs']:
                tmp = StringIO()
                io.set_structure(alt_structs[alt_name])
                io.save(tmp, SubstructureSelector(self.alt_substructure_mapping[alt_name], parsed_selection=SelectionParser(selection)))
                zipf.writestr(alt_name, tmp.getvalue())

        zipf.close()
        if len(out_stream.getvalue()) > 0:
            response = HttpResponse(content_type="application/zip")
            response['Content-Disposition'] = 'attachment; filename="Superposed_structures.zip"'
            response.write(out_stream.getvalue())

        if 'ref_file' in request.FILES:
            request.session['ref_file'] = request.FILES['ref_file']
        if 'alt_files' in request.FILES:
            request.session['alt_files'] = request.FILES.getlist('alt_files')


        return response


class FragmentSuperpositionIndex(TemplateView):

    template_name = 'common_structural_tools.html'

    #Left panel
    step = 1
    number_of_steps = 1

    documentation_url = settings.DOCUMENTATION_URL
    docs = 'sites.html#pharmacophore-generation'
    
    title = "SUPERPOSE FRAGMENTS OF CRYSTAL STRUCTURES"
    description = """
    The tool implements a fragment-based pharmacophore method, as published in <a href='http://www.ncbi.nlm.nih.gov/pubmed/25286328'>Fidom K, et al (2015)</a>. Interacting ligand moiety - residue pairs extracted from selected crystal structures of GPCRs are superposed onto the input pdb file based on gpcrdb generic residue numbers. Resulting aligned ligand fragments can be used for placement of pharmacophore features.

    Upload a pdb file you want to superpose the interacting moiety - residue pairs.
    
    Once you have selected all your targets, click the green button.
        """

    #Input file form data
    header = "Select a file to upload:"
    #Can't control the class properly - staying with the dirty explicit html code
    form_id='fragments'
    form_code = """
    Pdb file:<input id="id_pdb_file" name="pdb_file" type="file" /></br>
    Similarity:</br>
    <input id="similarity" name="similarity" type="radio" value="identical" /> Use fragments with identical residues</br>
    <input checked="checked" id="similarity" name="similarity" type="radio" value="similar" /> Use fragments with residues of similar properties</br>

    Fragments:</br>
    <input checked="checked" id="representative" name="representative" type="radio" value="closest" /> Use fragments from the evolutionary closest crystal structure</br>
    <input id="representative" name="representative" type="radio" value="any" /> Use all available fragments</br></br>
    State:<select id="id_state" name="state">
    <option value="active">Antagonist-bound structures</option>
    <option value="inactive">Agonist-bound structures</option>
    </select>
    """
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
    mid_section = 'fragment_superposition_results.html'
    #Buttons - none

    def post (self, request, *args, **kwargs):
        
        frag_sp = FragmentSuperpose(StringIO(request.FILES['pdb_file'].file.read().decode('UTF-8', 'ignore')),request.FILES['pdb_file'].name)
        superposed_fragments = []
        superposed_fragments_repr = []
        if request.POST['similarity'] == 'identical':
            if request.POST['representative'] == 'any':
                superposed_fragments = frag_sp.superpose_fragments()
            else:
                superposed_fragments_repr = frag_sp.superpose_fragments(representative=True, state=request.POST['state'])
                superposed_fragments = frag_sp.superpose_fragments()
        else:
            if request.POST['representative'] == 'any':
                superposed_fragments = frag_sp.superpose_fragments(use_similar=True)
            else:
                superposed_fragments_repr = frag_sp.superpose_fragments(representative=True, use_similar=True, state=request.POST['state'])
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
    # Left panel
    description = 'Select a reference target by searching or browsing in the middle column.' \
        + '\n\nThe selected reference target will appear in the right column.' \
        + '\n\nOnce you have selected your reference target, either proceed with all TMs alignment ("Find template"' \
        + 'button) or specify the sequence segments manualy ("Advanced segment selection" button).'
    step = 1
    number_of_steps = 2
    redirect_on_select = False

    docs = 'structures.html#template-selection'

    # Mid section

    # Right panel
    buttons = OrderedDict()
    buttons['continue'] = {
        'label': 'Find template',
        'url': '/structure/template_browser',
        'color': 'success',
    }
    buttons['segments'] = {
        'label' : 'Advanced segment selection',
        'url' : '/structure/template_segment_selection',
        'color' : 'info',
    }

    selection_boxes = OrderedDict([('reference', True),
        ('targets', False),
        ('segments', False),])



#==============================================================================
class TemplateSegmentSelection(AbsSegmentSelection):
    """
    Advanced selection of sequence segments for template search.
    """
   #Left panel
    step = 2
    number_of_steps = 2

    docs = 'structures.html#template-selection'

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

    template_name = "template_browser.html"

    def get_context_data (self, **kwargs):

        context = super(TemplateBrowser, self).get_context_data(**kwargs)

        # get simple selection from session
        simple_selection = self.request.session.get('selection', False)
        
        # make an alignment
        a = Alignment()
        a.ignore_alternative_residue_numbering_schemes = True

        # load the selected reference into the alignment
        a.load_reference_protein_from_selection(simple_selection)
        
        # fetch 
        qs = Structure.objects.all().select_related(
            "pdb_code__web_resource",
            "protein_conformation__protein__species",
            "protein_conformation__protein__source",
            "protein_conformation__protein__family__parent__parent__parent",
            "publication__web_link__web_resource").prefetch_related(
            "stabilizing_agents",
            "protein_conformation__protein__parent__endogenous_ligands__properities__ligand_type",
            Prefetch("ligands", queryset=StructureLigandInteraction.objects.filter(
            annotated=True).prefetch_related('ligand__properities__ligand_type', 'ligand_role')))

        # Dirty but fast
        qsd = {}
        for st in qs:
            qsd[st.protein_conformation.protein.id] = st
        
        # add proteins to the alignment
        a.load_proteins([x.protein_conformation.protein for x in qs])
        
        if simple_selection.segments != []:
            a.load_segments_from_selection(simple_selection)
        else:
            a.load_segments(ProteinSegment.objects.filter(slug__in=['TM1', 'TM2', 'TM3', 'TM4','TM5','TM6', 'TM7']))
        
        a.build_alignment()
        a.calculate_similarity()

        context['structures'] = []
        for prot in a.proteins[1:]:
            try:
                context['structures'].append([prot.similarity, prot.identity, qsd[prot.protein.id]])
                del qsd[prot.protein.id]
            except KeyError:
                pass
        return context

#==============================================================================
class PDBClean(TemplateView):
    """
    Extraction, packing and serving out the pdb records selected via structure/template browser.
    """
    template_name = "pdb_download.html"

    def post(self, request, *args, **kwargs):

        context = super(PDBClean, self).get_context_data(**kwargs)

        self.posted = True
        pref = True
        water = False
        hets = False
        
        if 'pref_chain' not in request.POST.keys():
            pref = False
        if 'water' in request.POST.keys():
            water = True
        if 'hets' in request.POST.keys():
            hets = True

        # get simple selection from session
        simple_selection = request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)
        out_stream = BytesIO()
        io = PDBIO()
        zipf = zipfile.ZipFile(out_stream, 'w', zipfile.ZIP_DEFLATED)
        if selection.targets != []:
            for selected_struct in [x for x in selection.targets if x.type == 'structure']:
                struct_name = '{}_{}.pdb'.format(selected_struct.item.protein_conformation.protein.parent.entry_name, selected_struct.item.pdb_code.index)
                if hets:
                    lig_names = [x.pdb_reference for x in StructureLigandInteraction.objects.filter(structure=selected_struct.item, annotated=True)]
                else:
                    lig_names = None
                gn_assigner = GenericNumbering(structure=PDBParser(QUIET=True).get_structure(struct_name, StringIO(selected_struct.item.get_cleaned_pdb(pref, water, lig_names)))[0])
                tmp = StringIO()
                io.set_structure(gn_assigner.assign_generic_numbers())
                request.session['substructure_mapping'] = gn_assigner.get_substructure_mapping_dict()
                io.save(tmp)
                zipf.writestr(struct_name, tmp.getvalue())
                del gn_assigner, tmp
            for struct in selection.targets:
                selection.remove('targets', 'structure', struct.item.id)
            # export simple selection that can be serialized
            simple_selection = selection.exporter()

            request.session['selection'] = simple_selection
            request.session['cleaned_structures'] = out_stream

        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]
        return render(request, self.template_name, context)
        

    def get_context_data (self, **kwargs):

        context = super(PDBClean, self).get_context_data(**kwargs)
        self.success = False
        self.posted = False

        # get simple selection from session
        simple_selection = self.request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)
        if selection.targets != []:
            self.success = True

        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]

        return context


class PDBSegmentSelection(AbsSegmentSelection):

    #Left panel
    step = 2
    number_of_steps = 2

    #Mid section
    #mid_section = 'segment_selection.html'

    #Right panel
    segment_list = True

    # OrderedDict to preserve the order of the boxes
    selection_boxes = OrderedDict([('reference', False),
        ('targets', False),
        ('segments', True),])

    def get_context_data(self, **kwargs):

        self.buttons = {
            'continue': {
                'label': 'Download substructures',
                'url': '/structure/pdb_download/custom',
                'color': 'success',
            },
        }
        context = super(PDBSegmentSelection, self).get_context_data(**kwargs)

        simple_selection = self.request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)

        context['selection'] = {}
        context['selection']['site_residue_groups'] = selection.site_residue_groups
        context['selection']['active_site_residue_group'] = selection.active_site_residue_group
        for selection_box, include in self.selection_boxes.items():
            if include:
                context['selection'][selection_box] = selection.dict(selection_box)['selection'][selection_box]

        # get attributes of this class and add them to the context
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        for a in attributes:
            if not(a[0].startswith('__') and a[0].endswith('__')):
                context[a[0]] = a[1]
        return context

#==============================================================================
class PDBDownload(View):
    """
    Serve the PDB (sub)structures depending on user's choice.
    """
    
    def get(self, request, *args, **kwargs):

        if self.kwargs['substructure'] == 'select':
            return HttpResponseRedirect('/structure/pdb_segment_selection')

        if self.kwargs['substructure'] == 'full':
            out_stream = request.session['cleaned_structures']

        elif self.kwargs['substructure'] == 'custom':
            simple_selection = request.session.get('selection', False)
            selection = Selection()
            if simple_selection:
                selection.importer(simple_selection)
            io = PDBIO()
            zipf_in = zipfile.ZipFile(request.session['cleaned_structures'], 'r')
            out_stream = BytesIO()
            zipf_out = zipfile.ZipFile(out_stream, 'w', zipfile.ZIP_DEFLATED)
            for name in zipf_in.namelist():
                tmp = StringIO()
                io.set_structure(PDBParser(QUIET=True).get_structure(name, StringIO(zipf_in.read(name).decode('utf-8')))[0])
                io.save(tmp, SubstructureSelector(request.session['substructure_mapping'], parsed_selection=SelectionParser(selection)))
                zipf_out.writestr(name, tmp.getvalue())

            zipf_in.close()
            zipf_out.close()
            del request.session['substructure_mapping']
        if len(out_stream.getvalue()) > 0:
            response = HttpResponse(content_type="application/zip")
            response['Content-Disposition'] = 'attachment; filename="pdb_structures.zip"'
            response.write(out_stream.getvalue())

        return response

#==============================================================================
def ConvertStructuresToProteins(request):
    "For alignment from structure browser"

    simple_selection = request.session.get('selection', False)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)
    if selection.targets != []:
        for struct in selection.targets:
            prot = struct.item.protein_conformation.protein.parent
            selection.remove('targets', 'structure', struct.item.id)
            selection.add('targets', 'protein', SelectionItem('protein', prot))
    # export simple selection that can be serialized
    simple_selection = selection.exporter()

    # add simple selection to session
    request.session['selection'] = simple_selection

    return HttpResponseRedirect('/alignment/segmentselection')


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

def RenderTrees(request):
    number = request.GET['number']
    tree = open(settings.STATICFILES_DIRS[0] +'/home/images/00'+number+'_tree.xml').read()
    legend = open(settings.STATICFILES_DIRS[0] +'/home/images/00'+number+'_legend.svg').read()
    context = {'tree':tree, 'leg':legend, 'num':number}
    return render(request, 'phylogenetic_trees.html', context)
