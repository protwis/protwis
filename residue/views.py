from django.conf import settings
from django.shortcuts import render
from django.views.generic import TemplateView

from common.views import AbsTargetSelection
from common.selection import Selection
from protein.models import ProteinSegment
from residue.models import Residue,ResidueNumberingScheme

class TargetSelection(AbsTargetSelection):
    pass

class ResidueTablesSelection(AbsTargetSelection):

    # Left panel
    step = 1
    number_of_steps = 2
    
    description = 'Select receptors to index by searching or browsing in the middle column. You can select entire receptor families and/or individual receptors.\n\nSelected receptors will appear in the right column, where you can edit the list.\n\nSelect which numbering schemes to use in the middle column. By default, only the GPCRDB numbering scheme is selected.\n\nOnce you have selected all your receptors, click the green button.'


    # Middle section
    numbering_schemes = True


    # Buttons
    buttons = {
        'continue' : {
            'label' : 'Select segments',
            'url' : '/residue/residuetabledisplay',
            'color' : 'success',
            }
        }


class ResidueTablesDisplay(TemplateView):
    """
    A class rendering the residue numbering table.
    """
    template_name = 'residue_table.html'

    def get_context_data(self, **kwargs):
        """
        Get the selection data (proteins and numbering schemes) and prepare it for display.
        """
        context = super().get_context_data(**kwargs)

        # get the selection from session
        simple_selection = self.request.session.get('selection', False)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)
        # extract numbering schemes and proteins
        numbering_schemes = [x.item for x in selection.numbering_schemes]
        print(numbering_schemes)
        proteins = [x.item for x in selection.targets]
        print(proteins)
        
        # get the helices (TMs only at first)
        segments = ProteinSegment.objects.filter(category='helix')

        if ResidueNumberingScheme.objects.get(slug=settings.DEFAULT_NUMBERING_SCHEME) in numbering_schemes:
            default_scheme = ResidueNumberingScheme.objects.get(slug=settings.DEFAULT_NUMBERING_SCHEME)
        else:
            default_scheme = numbering_schemes[0]

        # prepare the dictionary
        # each helix has a dictionary of positions
        # default_generic_number or first scheme on the list is the key
        # value is a dictionary of other gn positions and residues from selected proteins 
        data = {}
        for segment in segments:
            data[segment.slug] = {}
            residues = Residue.objects.filter(protein_segment=segment, protein_conformation__protein__in=proteins).prefetch_related('generic_number','alternative_generic_numbers')
            for scheme in numbering_schemes:
                if scheme == default_scheme and scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                    for pos in list(set([x.generic_number.label for x in residues if x.protein_segment == segment])):
                        data[segment.slug][pos] = {scheme.slug : pos, 'seq' : ['-']*len(proteins)}
                        #data[segment.slug][pos][] = pos
                elif scheme == default_scheme:
                    for pos in list(set([x.alternative_generic_numbers.filter(scheme__slug=scheme.slug).label for x in residues if x.protein_segment == segment])):
                        data[segment.slug][pos] = {scheme.slug : pos, 'seq' : ['-']*len(proteins)}
                        #data[segment.slug][pos][scheme.slug] = pos
            
            for residue in residues:
                for scheme in numbering_schemes:
                    #if scheme == default_scheme:
                    #    continue
                    if default_scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                        pos = residue.generic_number
                        if scheme == pos.scheme:
                            data[segment.slug][pos.label]['seq'][proteins.index(residue.protein_conformation.protein)] = str(residue)
                        else:
                            if scheme.slug not in data[segment.slug][pos.label].keys():
                                data[segment.slug][pos.label][scheme.slug] = residue.alternative_generic_numbers.get(scheme__slug=scheme.slug).label
                            data[segment.slug][pos.label]['seq'][proteins.index(residue.protein_conformation.protein)] = str(residue)
                    else:
                        pos = residue.alternative_generic_numbers.get(scheme__slug=default_scheme.slug)
                        if scheme.slug not in data[segment.slug][pos.label].keys():
                            data[segment.slug][pos.label][scheme.slug] = residue.alternative_generic_numbers.get(scheme__slug=scheme.slug).label
                        data[segment.slug][pos.label]['seq'][proteins.index(residue.protein_conformation.protein)] = str(residue)

        for y in segments:
            print(y)
            for pos in sorted([x for x in data[y.slug].keys()]):
                try:
                    print(data[y.slug][pos])
                except:
                    pass

        return context


def SelectionItems(request):

    simple_selection = request.session.get('selection', False)
    selection = Selection()
    if simple_selection:
        selection.importer(simple_selection)
    context = {}
    context['items'] = selection.numbering_schemes
    print(context['items'])
    return render(request, 'selection_view.html', context)