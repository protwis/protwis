from django.conf import settings

from common.selection import Selection
from protein.models import Protein
from residue.models import Residue
from residue.models import ResidueGenericNumber
from residue.models import ResidueNumberingScheme

from collections import OrderedDict
from copy import deepcopy


class Alignment:
    """A class representing a protein sequence alignment, with or without a reference sequence"""
    def __init__(self):
        self.reference = False
        self.proteins = []
        self.segments = OrderedDict()
        self.positions = []
        self.matrix = []
        self.numbering_scheme = ResidueNumberingScheme.objects.get(slug=settings.DEFAULT_NUMBERING_SCHEME)

    def __str__(self):
        return str(self.__dict__)

    def calculate_similarity(self):
        """Calculate the sequence identity/similarity of every selected protein compared to a selected reference"""
        pass

    def load_proteins_from_selection(self, simple_selection):
        """Read user selection and fetch selected proteins from the DB"""
        # create full selection and import simple selection (if it exists)
        selection = Selection()
        if simple_selection:
            selection.importer(simple_selection)

        # flatten the selection into individual proteins
        for target in selection.targets:
            if target.type == 'protein':
                self.proteins.append(target.item)
            elif target.type == 'family':
                # species filter
                species_list = []
                for species in selection.species:
                    species_list.append(species.item)

                # annotation filter
                protein_source_list = []
                for protein_source in selection.annotation:
                    protein_source_list.append(protein_source.item)
                    
                family_proteins = Protein.objects.filter(family__slug__startswith=target.item.slug,
                    species__in=(species_list),
                    source__in=(protein_source_list))
                for fp in family_proteins:
                    self.proteins.append(fp)

    def load_positions_from_selection(self, simple_selection):
        """Read user selection and add selected protein segments/residue positions"""
        for segment in simple_selection.segments:
            segment_residues = ResidueGenericNumber.objects.filter(protein_segment=segment.item,
                scheme=self.numbering_scheme).order_by('label')
            self.segments[segment.item.slug] = OrderedDict()
            for segment_residue in segment_residues:
                self.segments[segment.item.slug][segment_residue.label] = []

    def build_alignment_matrix(self):
        """Fetch selected residues from DB and build an alignment matrix"""
        for p in self.proteins:
            row = []
            for segment, positions in self.segments.items():
                s = []

                # loop all positions in this segment
                for pos, display_numbers in positions.items():
                    try:
                        # fetch the residue with this generic number
                        r = Residue.objects.get(generic_number__scheme=self.numbering_scheme,
                            generic_number__label=pos, protein=p)
                        
                        # add position to the list of positions that are not empty
                        if pos not in self.positions:
                            self.positions.append(pos)

                        # add display number to list of display numbers for this position
                        if r.display_generic_number.label not in self.segments[segment][pos]:
                            self.segments[segment][pos].append(r.display_generic_number.label)

                        # append the residue to the matrix
                        s.append([pos, r.display_generic_number.label, r.amino_acid])
                    except:
                        s.append([pos, False, '-'])
                row.append(s)
            self.matrix.append(row)
            self.positions.sort()
        self.merge_generic_numbers()
        self.clear_empty_positions()

    def clear_empty_positions(self):
        """Remove empty columns from the segments and matrix"""
        segments = deepcopy(self.segments) # deepcopy is required because the dictionary changes during the loop
        for segment, positions in segments.items():
            for pos, dns in positions.items():
                if pos not in self.positions:
                    del self.segments[segment][pos]

        matrix = deepcopy(self.matrix) # deepcopy is required because the list changes during the loop
        for i, row in enumerate(matrix):
            for j, s in enumerate(row):
                for p in s:
                    if p[0] not in self.positions:
                        self.matrix[i][j].remove(p)

    def merge_generic_numbers(self):
        """Check whether there are many display numbers for each position, and merge them if there are"""
        segments = deepcopy(self.segments) # deepcopy is required because the dictionary changes during the loop
        for segment, positions in segments.items():
            for pos, dns in positions.items():
                if len(dns) == 1:
                    self.segments[segment][pos] = self.format_generic_number(dns[0])
                else:
                    for dn in dns:
                        split_dn = dn.split('x')

    def format_generic_number(self, generic_number):
        """A placeholder for an instance specific function"""
        return generic_number