from common.selection import Selection
from protein.models import Protein

from collections import OrderedDict

class Alignment:
    """A class representing a protein sequence alignment, with or without a reference sequence"""
    def __init__(self):
        self.reference = False
        self.proteins = []
        self.segments = OrderedDict()
        self.matrix = []

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
                family_proteins = Protein.objects.filter(family__slug__startswith=target.item.slug)
                for fp in family_proteins:
                    self.proteins.append(fp)

    def load_positions_from_selection(self, simple_selection):
        self.segments['TM1'] = [140, 141, 142, 143, 144, 145, 146, 147, 148, 149]
        self.segments['Custom'] = [243, 251, 252]

    def build_alignment_matrix(self):
        for p in self.proteins:
            r = []
            for segment,positions in self.segments.items():
                s = []
                for pos in positions:
                    s.append(p.sequence[pos])
                r.append(s)
            self.matrix.append(r)