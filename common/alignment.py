from django.conf import settings

from common.selection import Selection
from protein.models import Protein
from residue.models import Residue
from residue.models import ResidueGenericNumber
from residue.models import ResidueNumberingScheme

from collections import OrderedDict
from copy import deepcopy
from operator import itemgetter


class Alignment:
    """A class representing a protein sequence alignment, with or without a reference sequence"""
    def __init__(self):
        self.reference = False
        self.proteins = []
        self.segments = OrderedDict()
        self.numbering_schemes = {}
        self.generic_numbers = OrderedDict()
        self.positions = []
        self.matrix = []
        self.default_numbering_scheme = ResidueNumberingScheme.objects.get(slug=settings.DEFAULT_NUMBERING_SCHEME)

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
                    source__in=(protein_source_list)).select_related('residue_numbering_scheme', 'species')
                for fp in family_proteins:
                    self.proteins.append(fp)
        
        # update numbering scheme list
        for p in self.proteins:
            if p.residue_numbering_scheme.slug not in self.numbering_schemes:
                self.numbering_schemes[p.residue_numbering_scheme.slug] = p.residue_numbering_scheme.name
        # order and convert numbering scheme dict to tuple
        self.numbering_schemes = sorted(self.numbering_schemes.items(), key=itemgetter(0))

    def load_positions_from_selection(self, simple_selection):
        """Read user selection and add selected protein segments/residue positions"""
        for segment in simple_selection.segments:
            segment_residues = ResidueGenericNumber.objects.filter(protein_segment=segment.item,
                scheme=self.default_numbering_scheme).order_by('label')
            
            # generic numbers in the schemes of all selected proteins
            for ns in self.numbering_schemes:
                if ns[0] not in self.generic_numbers:
                    self.generic_numbers[ns[0]] = OrderedDict()
                self.generic_numbers[ns[0]][segment.item.slug] = OrderedDict()
                for segment_residue in segment_residues:
                    self.generic_numbers[ns[0]][segment.item.slug][segment_residue.label] = []

            # segments
            self.segments[segment.item.slug] = []
            for segment_residue in segment_residues:
                self.segments[segment.item.slug].append(segment_residue.label)


    def build_alignment_matrix(self):
        """Fetch selected residues from DB and build an alignment matrix"""
        for p in self.proteins:
            row = []
            for segment, positions in self.segments.items():
                s = []
                first_residue_found = False
                # counters to keep track of gaps at the end of a segment
                gap_counter = 0
                position_counter = 1

                # fetch all residues for this segment. Prefetch the generic numbers (a lot faster)
                rs = Residue.objects.filter(generic_number__scheme=self.default_numbering_scheme,
                    generic_number__label__in=positions, protein=p).prefetch_related('generic_number', 'display_generic_number__scheme', 'alternative_generic_number__scheme')

                # create a dict of residues with their generic number as key for easy lookup
                res = {}
                for r in rs:
                    res[r.generic_number.label] = r

                # loop all positions in this segment
                for pos in positions:
                    try:
                        # find the residue record from the dict defined above
                        r = res[pos]
                        
                        # add position to the list of positions that are not empty
                        if pos not in self.positions:
                            self.positions.append(pos)

                        # add display number to list of display numbers for this position
                        if r.display_generic_number.label not in self.generic_numbers[p.residue_numbering_scheme.slug][segment][pos]:
                            self.generic_numbers[p.residue_numbering_scheme.slug][segment][pos].append(r.display_generic_number.label)
                            pass

                        # add display numbers for other numbering schemes of selected proteins
                        if len(self.numbering_schemes) > 1:
                            for arn in r.alternative_generic_number.all():
                                for ns in self.numbering_schemes:
                                    if arn.scheme.slug == ns[0] and arn.scheme.slug != p.residue_numbering_scheme.slug:
                                        self.generic_numbers[arn.scheme.slug][segment][pos].append(arn.label)
                                        break

                        # append the residue to the matrix
                        s.append([pos, r.display_generic_number.label, r.amino_acid, r.display_generic_number.scheme.short_name])
                        first_residue_found = True

                        # reset gap counter
                        gap_counter = 0
                    except:
                        if first_residue_found:
                            s.append([pos, False, '-'])

                            # update gap counter
                            gap_counter += 1

                            # if this is the last residue and there are gaps and the end of the segment, update them to
                            # end gaps
                            if (position_counter) == len(positions):
                                for i in range(gap_counter):
                                    s[len(positions)-(i+1)][2] = '_'
                        else:
                            s.append([pos, False, '_'])
                    
                    # update position counter
                    position_counter += 1
                row.append(s)
            self.matrix.append(row)
            self.positions.sort()
        self.merge_generic_numbers()
        self.clear_empty_positions()

    def clear_empty_positions(self):
        """Remove empty columns from the segments and matrix"""
        # segments and generic_numbers
        # deepcopy is required because the dictionary changes during the loop
        generic_numbers = deepcopy(self.generic_numbers)
        for ns, segments in generic_numbers.items():
            for segment, positions in segments.items():
                for pos, dn in positions.items():
                    if pos not in self.positions:
                        del self.generic_numbers[ns][segment][pos]
                        
                        if pos in self.segments[segment]:
                            self.segments[segment].remove(pos)

        # matrix
        matrix = deepcopy(self.matrix) # deepcopy is required because the list changes during the loop
        for i, row in enumerate(matrix):
            for j, s in enumerate(row):
                for p in s:
                    if p[0] not in self.positions:
                        self.matrix[i][j].remove(p)

    def merge_generic_numbers(self):
        """Check whether there are many display numbers for each position, and merge them if there are"""
        # deepcopy is required because the dictionary changes during the loop
        generic_numbers = deepcopy(self.generic_numbers)
        for ns, segments in generic_numbers.items():
            for segment, positions in segments.items():
                for pos, dns in positions.items():
                    if len(dns) == 1:
                        self.generic_numbers[ns][segment][pos] = self.format_generic_number(dns[0])
                    else:
                        self.generic_numbers[ns][segment][pos] = '-'.join(dns)

    def format_generic_number(self, generic_number):
        """A placeholder for an instance specific function"""
        return generic_number