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
        self.numbering_scheme = ResidueNumberingScheme.objects.get(slug='oliveira')

    def __str__(self):
        return str(self.__dict__)

    def calculate_similarity(self):
        for position in self.segments:
            pass

    def load_proteins_from_selection(self, simple_selection):
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
        for segment in simple_selection.segments:
            segment_residues = ResidueGenericNumber.objects.filter(protein_segment=segment.item,
                scheme=self.numbering_scheme).order_by('label')
            self.segments[segment.item.slug] = []
            for segment_residue in segment_residues:
                # FIXME move this to a instance specific function
                split_label = segment_residue.label.split("x")
                formatted_label = split_label[0] + "<br>x" + split_label[1]
                ###############
                self.segments[segment.item.slug].append([segment_residue.label, formatted_label])

    def build_alignment_matrix(self):
        for p in self.proteins:
            row = []
            for segment, positions in self.segments.items():
                s = []
                for pos in positions:
                    try:
                        r = Residue.objects.get(generic_number__label=pos[0], protein=p)
                        gns = r.generic_number.all()
                        for gn in gns:
                            if gn.scheme.slug == 'oliveira':
                                if gn.label not in self.positions:
                                    self.positions.append(gn.label)
                                break
                        s.append([pos[0], r.amino_acid])
                    except:
                        s.append([pos[0], '-'])
                row.append(s)
            self.matrix.append(row)
            self.positions.sort()

    def clear_empty_positions(self):
        segments = deepcopy(self.segments) # deepcopy is required because the dictionary changes during the loop
        for title, s in segments.items():
            for p in s:
                if p[0] not in self.positions:
                    self.segments[title].remove(p)

        matrix = deepcopy(self.matrix) # deepcopy is required because the list changes during the loop
        for i, row in enumerate(matrix):
            for j, s in enumerate(row):
                for p in s:
                    if p[0] not in self.positions:
                        self.matrix[i][j].remove(p)
