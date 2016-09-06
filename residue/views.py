from django.conf import settings
from django.shortcuts import render
from django.views.generic import TemplateView

from common.views import AbsTargetSelection
from common.selection import Selection
from protein.models import ProteinSegment, Protein
from residue.models import Residue,ResidueNumberingScheme

from collections import OrderedDict

import re

class TargetSelection(AbsTargetSelection):
    pass

class ResidueTablesSelection(AbsTargetSelection):

    # Left panel
    step = 1
    number_of_steps = 2
    docs = 'generic_numbering.html'
    
    description = 'Select receptors to index by searching or browsing in the middle column. You can select entire' \
        + ' receptor families and/or individual receptors.\n\nSelected receptors will appear in the right column,' \
        + ' where you can edit the list.\n\nSelect which numbering schemes to use in the middle column.\n\nOnce you' \
        + ' have selected all your receptors, click the green button.'


    # Middle section
    numbering_schemes = True


    # Buttons
    buttons = {
        'continue' : {
            'label' : 'Show residue numbers',
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

        # get the user selection from session
        simple_selection = self.request.session.get('selection', False)
        
         # local protein list
        proteins = []

        # flatten the selection into individual proteins
        for target in simple_selection.targets:
            if target.type == 'protein':
                proteins.append(target.item)
            elif target.type == 'family':
                # species filter
                species_list = []
                for species in simple_selection.species:
                    species_list.append(species.item)

                # annotation filter
                protein_source_list = []
                for protein_source in simple_selection.annotation:
                    protein_source_list.append(protein_source.item)

                if species_list:
                    family_proteins = Protein.objects.filter(family__slug__startswith=target.item.slug,
                        species__in=(species_list),
                        source__in=(protein_source_list)).select_related('residue_numbering_scheme', 'species')
                else:
                    family_proteins = Protein.objects.filter(family__slug__startswith=target.item.slug,
                        source__in=(protein_source_list)).select_related('residue_numbering_scheme', 'species')
                    
                for fp in family_proteins:
                    proteins.append(fp)

            longest_name = 0
            species_list = {}
            for protein in proteins:
                if protein.species.common_name not in species_list:
                    if len(protein.species.common_name)>10 and len(protein.species.common_name.split())>1:
                        name = protein.species.common_name.split()[0][0]+". "+" ".join(protein.species.common_name.split()[1:])
                        if len(" ".join(protein.species.common_name.split()[1:]))>11:
                            name = protein.species.common_name.split()[0][0]+". "+" ".join(protein.species.common_name.split()[1:])[:8]+".."
                    else:
                        name = protein.species.common_name
                    species_list[protein.species.common_name] = name

                    if len(re.sub('<[^>]*>', '', protein.name)+" "+name)>longest_name:
                        longest_name = len(re.sub('<[^>]*>', '', protein.name)+" "+name)

        # get the selection from session
        selection = Selection()
        if simple_selection:
             selection.importer(simple_selection)
        # # extract numbering schemes and proteins
        numbering_schemes = [x.item for x in selection.numbering_schemes]
        
        # # get the helices (TMs only at first)
        segments = ProteinSegment.objects.filter(category='helix')

        if ResidueNumberingScheme.objects.get(slug=settings.DEFAULT_NUMBERING_SCHEME) in numbering_schemes:
            default_scheme = ResidueNumberingScheme.objects.get(slug=settings.DEFAULT_NUMBERING_SCHEME)
        else:
            default_scheme = numbering_schemes[0]

        # prepare the dictionary
        # each helix has a dictionary of positions
        # default_generic_number or first scheme on the list is the key
        # value is a dictionary of other gn positions and residues from selected proteins 
        data = OrderedDict()
        print(segments)
        for segment in segments:
            data[segment.slug] = OrderedDict()
            residues = Residue.objects.filter(protein_segment=segment, protein_conformation__protein__in=proteins).prefetch_related('protein_conformation__protein', 'protein_conformation__state', 'protein_segment',
                'generic_number__scheme', 'display_generic_number__scheme', 'alternative_generic_numbers__scheme')
            for scheme in numbering_schemes:
                if scheme == default_scheme and scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                    for pos in list(set([x.generic_number.label for x in residues if x.protein_segment == segment])):
                        data[segment.slug][pos] = {scheme.slug : pos, 'seq' : ['-']*len(proteins)}
                elif scheme == default_scheme:
                    for pos in list(set([x.generic_number.label for x in residues if x.protein_segment == segment])):
                            data[segment.slug][pos] = {scheme.slug : pos, 'seq' : ['-']*len(proteins)}

            for residue in residues:
                alternatives = residue.alternative_generic_numbers.all()
                pos = residue.generic_number
                for alternative in alternatives:
                    scheme = alternative.scheme
                    if default_scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                        pos = residue.generic_number
                        if scheme == pos.scheme:
                            data[segment.slug][pos.label]['seq'][proteins.index(residue.protein_conformation.protein)] = str(residue)
                        else:
                            if scheme.slug not in data[segment.slug][pos.label].keys():
                                data[segment.slug][pos.label][scheme.slug] = alternative.label
                            if alternative.label not in data[segment.slug][pos.label][scheme.slug]:
                                data[segment.slug][pos.label][scheme.slug] += " "+alternative.label
                            data[segment.slug][pos.label]['seq'][proteins.index(residue.protein_conformation.protein)] = str(residue)
                    else:
                        if scheme.slug not in data[segment.slug][pos.label].keys():
                            data[segment.slug][pos.label][scheme.slug] = alternative.label
                        if alternative.label not in data[segment.slug][pos.label][scheme.slug]:
                            data[segment.slug][pos.label][scheme.slug] += " "+alternative.label
                        data[segment.slug][pos.label]['seq'][proteins.index(residue.protein_conformation.protein)] = str(residue)
        print(data)
        # Preparing the dictionary of list of lists. Dealing with tripple nested dictionary in django templates is a nightmare
        flattened_data = OrderedDict.fromkeys([x.slug for x in segments], [])
        print(flattened_data)
        for s in iter(flattened_data):
            flattened_data[s] = [[data[s][x][y.slug] for y in numbering_schemes]+data[s][x]['seq'] for x in sorted(data[s])]
        
        context['header'] = zip([x.short_name for x in numbering_schemes] + [x.name+" "+species_list[x.species.common_name] for x in proteins], [x.name for x in numbering_schemes] + [x.name for x in proteins],[x.name for x in numbering_schemes] + [x.entry_name for x in proteins])
        context['segments'] = [x.slug for x in segments]
        context['data'] = flattened_data
        context['number_of_schemes'] = len(numbering_schemes)
        context['longest_name'] = {'div' : longest_name*2, 'height': longest_name*2+75}

        return context