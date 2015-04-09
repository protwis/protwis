from django.conf import settings

from common.selection import Selection
from protein.models import Protein
from residue.models import Residue
from residue.models import ResidueGenericNumber
from residue.models import ResidueNumberingScheme

from collections import OrderedDict
from copy import deepcopy
from operator import itemgetter
from Bio.SubsMat import MatrixInfo


class Alignment:
    """A class representing a protein sequence alignment, with or without a reference sequence"""
    def __init__(self):
        self.reference = False
        self.proteins = []
        self.segments = OrderedDict()
        self.numbering_schemes = {}
        self.generic_numbers = OrderedDict()
        self.positions = []
        self.default_numbering_scheme = ResidueNumberingScheme.objects.get(slug=settings.DEFAULT_NUMBERING_SCHEME)
        
        # refers to which Protein attribute to order by (identity, similarity or similarity score)
        self.order_by = 'similarity_score' 

    def __str__(self):
        return str(self.__dict__)

    def load_reference_protein(self, protein):
        """Loads a protein into the alignment as a reference"""
        self.proteins.insert(0, protein)
        self.update_numbering_schemes()

    def load_reference_protein_from_selection(self, simple_selection):
        """Read user selection and add selected reference protein"""
        if simple_selection and simple_selection.reference:
            self.load_reference_protein(simple_selection.reference[0].item)

    def load_proteins(self, proteins):
        """Load a list of protein objects into the alignment"""
        for protein in proteins:
            if protein not in self.proteins:
                self.proteins.append(protein)
        self.update_numbering_schemes()

    def load_proteins_from_selection(self, simple_selection):
        """Read user selection and add selected proteins"""
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
                    
                family_proteins = Protein.objects.filter(family__slug__startswith=target.item.slug,
                    species__in=(species_list),
                    source__in=(protein_source_list)).select_related('residue_numbering_scheme', 'species')
                for fp in family_proteins:
                    proteins.append(fp)

        # load protein list
        self.load_proteins(proteins)

    def load_segments(self, segments):
        for segment in segments:
            segment_residues = ResidueGenericNumber.objects.filter(protein_segment=segment,
                scheme=self.default_numbering_scheme).order_by('label')
            
            # generic numbers in the schemes of all selected proteins
            for ns in self.numbering_schemes:
                if ns[0] not in self.generic_numbers:
                    self.generic_numbers[ns[0]] = OrderedDict()
                self.generic_numbers[ns[0]][segment.slug] = OrderedDict()
                for segment_residue in segment_residues:
                    self.generic_numbers[ns[0]][segment.slug][segment_residue.label] = []

            # segments
            self.segments[segment.slug] = []
            for segment_residue in segment_residues:
                self.segments[segment.slug].append(segment_residue.label)

    def load_segments_from_selection(self, simple_selection):
        """Read user selection and add selected protein segments/residue positions"""
        # local segment list
        segments = []

        # read selection
        for segment in simple_selection.segments:
            segments.append(segment.item)

        # load segment positions
        self.load_segments(segments)

    def update_numbering_schemes(self):
        """Update numbering scheme list"""
        self.numbering_schemes = {}
        for p in self.proteins:
            if p.residue_numbering_scheme.slug not in self.numbering_schemes:
                self.numbering_schemes[p.residue_numbering_scheme.slug] = p.residue_numbering_scheme.name
        
        # order and convert numbering scheme dict to tuple
        self.numbering_schemes = sorted(self.numbering_schemes.items(), key=itemgetter(0))

    def build_alignment(self):
        """Fetch selected residues from DB and build an alignment"""
        # make a list of generic numbers from all segments
        generic_numbers = []
        for segment, gns in self.segments.items():
            generic_numbers += gns

        # fetch all residues for the selected proteins. Prefetch the generic numbers (a lot faster)
        # only fetch the alternative generic numbers if more than one numbering scheme is needed
        if len(self.numbering_schemes) > 1:
            rs = Residue.objects.filter(generic_number__scheme=self.default_numbering_scheme,
                generic_number__label__in=generic_numbers, protein__in=self.proteins).prefetch_related(
                'protein', 'protein_segment', 'generic_number', 'display_generic_number__scheme',
                'alternative_generic_number__scheme')
        else:
            rs = Residue.objects.filter(generic_number__scheme=self.default_numbering_scheme,
                generic_number__label__in=generic_numbers, protein__in=self.proteins).prefetch_related(
                'protein', 'protein_segment', 'generic_number', 'display_generic_number__scheme')

        # create a dict of proteins, segments and residues
        proteins = {}
        for r in rs:
            if r.protein.entry_name not in proteins:
                proteins[r.protein.entry_name] = {}
            if r.protein_segment.slug not in proteins[r.protein.entry_name]:
                proteins[r.protein.entry_name][r.protein_segment.slug] = {}
            proteins[r.protein.entry_name][r.protein_segment.slug][r.generic_number.label] = r

        for p in self.proteins:
            row = []
            for segment, positions in self.segments.items():
                s = []
                first_residue_found = False
                # counters to keep track of gaps at the end of a segment
                gap_counter = 0
                position_counter = 1

                # loop all positions in this segment
                for pos in positions:
                    try:
                        # find the residue record from the dict defined above
                        r = proteins[p.entry_name][segment][pos]
                        
                        # add position to the list of positions that are not empty
                        if pos not in self.positions:
                            self.positions.append(pos)

                        # add display number to list of display numbers for this position
                        if r.display_generic_number.label not in (
                            self.generic_numbers[p.residue_numbering_scheme.slug][segment][pos]):
                            self.generic_numbers[p.residue_numbering_scheme.slug][segment][pos].append(
                                r.display_generic_number.label)
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
            p.alignment = row
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
                        # remove position from generic numbers dict
                        del self.generic_numbers[ns][segment][pos]
                        
                        # remove position from segment dict
                        if pos in self.segments[segment]:
                            self.segments[segment].remove(pos)

        # proteins
        proteins = deepcopy(self.proteins) # deepcopy is required because the list changes during the loop
        for i, protein in enumerate(proteins):
            for j, s in enumerate(protein.alignment):
                for p in s:
                    if p[0] not in self.positions:
                        self.proteins[i].alignment[j].remove(p)

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

    def calculate_similarity(self):
        """Calculate the sequence identity/similarity of every selected protein compared to a selected reference"""
        gaps = ['-', '_']
        for i, protein in enumerate(self.proteins):
            # skip the first row, as it is the reference
            if i != 0:
                identities = []
                similarities = []
                similarity_scores = []
                for j, s in enumerate(protein.alignment):
                    for k, p in enumerate(s):
                        reference_residue = self.proteins[0].alignment[j][k][2]
                        protein_residue = self.proteins[i].alignment[j][k][2]
                        if not (reference_residue in gaps and protein_residue in gaps):
                            # identity
                            if protein_residue == reference_residue:
                                identities.append(1)
                            else:
                                identities.append(0)

                            # similarity
                            if reference_residue in gaps or protein_residue in gaps:
                                similarities.append(0)
                                similarity_scores.append(0)
                            else:
                                pair = (protein_residue, reference_residue)
                                similarity = self.score_match(pair, MatrixInfo.blosum62)
                                if similarity > 0:
                                    similarities.append(1)
                                else:
                                    similarities.append(0)
                                similarity_scores.append(similarity)
                
                # update the protein
                if identities:
                    self.proteins[i].identity = "{:10.0f}".format(sum(identities) / len(identities) * 100)
                if similarities:
                    self.proteins[i].similarity = "{:10.0f}".format(sum(similarities) / len(similarities) * 100)
                    self.proteins[i].similarity_score = sum(similarity_scores)

        # order protein list by similarity score
        ref = self.proteins.pop(0)
        # self.proteins.sort(key=lambda x: x.similarity_score, reverse=True)
        self.proteins.sort(key=lambda x: getattr(x, self.order_by), reverse=True)
        self.proteins.insert(0, ref)

    def score_match(self, pair, matrix):
        if pair not in matrix:
            return matrix[(tuple(reversed(pair)))]
        else:
            return matrix[pair]