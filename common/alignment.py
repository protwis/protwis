from django.conf import settings

from common.selection import Selection
from protein.models import Protein, ProteinConformation, ProteinState, ProteinSegment, ProteinFusionProtein
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
        self.states = [settings.DEFAULT_PROTEIN_STATE] # inactive, active etc
        
        # refers to which ProteinConformation attribute to order by (identity, similarity or similarity score)
        self.order_by = 'similarity_score' 

    def __str__(self):
        return str(self.__dict__)

    def load_reference_protein(self, protein):
        """Loads a protein into the alignment as a reference"""
        self.reference = True

        # fetch the selected conformations of the protein
        # only one can be selected for similarity search, therefore state[0] (other states ignored, if defined)
        try:
            pconf = ProteinConformation.objects.get(protein=protein,
                state=ProteinState.objects.get(slug=self.states[0]))
        except ProteinConformation.DoesNotExist:
            raise Exception ('Protein conformation {} not found for protein {}'.format(self.states[0],
                protein.entry_name))

        self.proteins.insert(0, pconf)
        self.update_numbering_schemes()

    def load_reference_protein_from_selection(self, simple_selection):
        """Read user selection and add selected reference protein"""
        if simple_selection and simple_selection.reference:
            self.load_reference_protein(simple_selection.reference[0].item)

    def load_proteins(self, proteins):
        """Load a list of protein objects into the alignment"""
        for protein in proteins:
            for state in self.states:
                try:
                    pconf = ProteinConformation.objects.get(protein=protein,
                        state=ProteinState.objects.get(slug=state))
                except ProteinConformation.DoesNotExist:
                    raise Exception ('Protein conformation {} not found for protein {}'.format(state,
                        protein.entry_name))
                if pconf not in self.proteins:
                    self.proteins.append(pconf)
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
        for s in segments:
            # fetch split segments (e.g. ICL3_1 and ICL3_2)
            alt_segments = ProteinSegment.objects.filter(slug__startswith=s.slug)
            for segment in alt_segments:
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
        for pc in self.proteins:
            if pc.protein.residue_numbering_scheme.slug not in self.numbering_schemes:
                rnsn = pc.protein.residue_numbering_scheme.name
                self.numbering_schemes[pc.protein.residue_numbering_scheme.slug] = rnsn
        
        # order and convert numbering scheme dict to tuple
        self.numbering_schemes = sorted(self.numbering_schemes.items(), key=itemgetter(0))

    def build_alignment(self):
        """Fetch selected residues from DB and build an alignment"""
        if len(self.numbering_schemes) > 1:
            rs = Residue.objects.filter(
                protein_segment__slug__in=self.segments, protein_conformation__in=self.proteins).prefetch_related(
                'protein_conformation__protein', 'protein_segment', 'generic_number', 'display_generic_number__scheme',
                'alternative_generic_number__scheme')
        else:
            rs = Residue.objects.filter(
                protein_segment__slug__in=self.segments, protein_conformation__in=self.proteins).prefetch_related(
                'protein_conformation__protein', 'protein_segment', 'generic_number', 'display_generic_number__scheme')

        # create a dict of proteins, segments and residues
        proteins = {}
        segment_counters = {}
        fusion_protein_inserted = {}
        for r in rs:
            # handle split segments, slug_1 and slug_2 are merged into slug, with a fusion protein in between
            pss = r.protein_segment.slug.split("_")
            if len(pss) > 1:
                ps = pss[0]
                segment_part = int(pss[1])
            else:
                ps = r.protein_segment.slug
                segment_part = 1

            # identifiers for protein/state/segment
            pcid = r.protein_conformation.protein.entry_name + "-" + r.protein_conformation.state.slug
            
            if pcid not in proteins:
                proteins[pcid] = {}
            if ps not in proteins[pcid]:
                proteins[pcid][ps] = {}
            if pcid not in segment_counters:
                segment_counters[pcid] = {}
            if ps not in segment_counters[pcid]:
                segment_counters[pcid][ps] = 1
            else:
                segment_counters[pcid][ps] += 1
            if pcid not in fusion_protein_inserted:
                fusion_protein_inserted[pcid] = {}
            if ps not in fusion_protein_inserted[pcid]:
                fusion_protein_inserted[pcid][ps] = False

            # user generic numbers as keys for aligned segments
            if r.generic_number:
                proteins[pcid][ps][r.generic_number.label] = r
            # use residue numbers in segment as keys for non-aligned segments
            else:
                # position label
                pos_label = ps + "-" + str("%04d" % (segment_counters[pcid][ps],))

                # insert fusion protein
                if not fusion_protein_inserted[pcid][ps] and segment_part == 2:
                    fp = ProteinFusionProtein.objects.get(protein=r.protein_conformation.protein,
                        segment_after=r.protein_segment)
                    fusion_pos_label = ps + "-" + str("%04d" % (segment_counters[pcid][ps]-1,)) + "-fusion"
                    proteins[pcid][ps][fusion_pos_label] = Residue(amino_acid=fp.protein_fusion.name)
                    if fusion_pos_label not in self.segments[ps]:
                        self.segments[ps].append(fusion_pos_label)
                    fusion_protein_inserted[pcid][ps] = True

                # residue
                proteins[pcid][ps][pos_label] = r
                if pos_label not in self.segments[ps]:
                    self.segments[ps].append(pos_label)

        # remove split segments from segment list and order segment positions
        for segment, positions in self.segments.items():
            s = segment.split("_")
            if len(s) > 1:
                del self.segments[segment]
            else:
                self.segments[segment].sort()

        # remove split segments from generic numbers list
        for ns, gs in self.generic_numbers.items():
            for segment, positions in gs.items():
                s = segment.split("_")
                if len(s) > 1:
                    del self.generic_numbers[ns][segment]

        for pc in self.proteins:
            row = []
            for segment, positions in self.segments.items():
                s = []
                first_residue_found = False
                
                # counters to keep track of gaps at the end of a segment
                gap_counter = 0
                position_counter = 1

                # numbering scheme
                ns = pc.protein.residue_numbering_scheme.slug

                # loop all positions in this segment
                for pos in positions:
                    try:
                        # find the residue record from the dict defined above
                        pcid = pc.protein.entry_name + "-" + pc.state.slug
                        r = proteins[pcid][segment][pos]
                        
                        # add position to the list of positions that are not empty
                        if pos not in self.positions:
                            self.positions.append(pos)

                        # add display number to list of display numbers for this position
                        if r.display_generic_number:
                            if r.display_generic_number.label not in self.generic_numbers[ns][segment][pos]:
                                self.generic_numbers[ns][segment][pos].append(r.display_generic_number.label)
                        else:
                            if pos not in self.generic_numbers[ns][segment]:
                                self.generic_numbers[ns][segment][pos] = []

                        # add display numbers for other numbering schemes of selected proteins
                        if r.generic_number and len(self.numbering_schemes) > 1:
                            for arn in r.alternative_generic_number.all():
                                for ns in self.numbering_schemes:
                                    if (arn.scheme.slug == ns[0] and
                                        arn.scheme.slug != ns):
                                        self.generic_numbers[arn.scheme.slug][segment][pos].append(arn.label)
                                        break

                        # append the residue to the matrix
                        if r.generic_number:
                            s.append([pos, r.display_generic_number.label, r.amino_acid,
                                r.display_generic_number.scheme.short_name])
                        else:
                            s.append([pos, "", r.amino_acid, ""])

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
            pc.alignment = row
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
                    if not dns: # don't format if there are no numbers
                        self.generic_numbers[ns][segment][pos] = ""
                    elif len(dns) == 1:
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
        self.proteins.sort(key=lambda x: getattr(x, self.order_by), reverse=True)
        self.proteins.insert(0, ref)

    def score_match(self, pair, matrix):
        if pair not in matrix:
            return matrix[(tuple(reversed(pair)))]
        else:
            return matrix[pair]