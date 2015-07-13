from django.conf import settings

from common.selection import Selection
from common.definitions import *
from protein.models import Protein, ProteinConformation, ProteinState, ProteinSegment, ProteinFusionProtein
from residue.models import Residue
from residue.models import ResidueGenericNumber
from residue.models import ResidueNumberingScheme
from structure.models import Structure

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
        self.consensus = []
        self.amino_acids = []
        self.amino_acid_stats = []
        self.features = []
        self.feature_stats = []
        self.default_numbering_scheme = ResidueNumberingScheme.objects.get(slug=settings.DEFAULT_NUMBERING_SCHEME)
        self.states = [settings.DEFAULT_PROTEIN_STATE] # inactive, active etc
        
        # refers to which ProteinConformation attribute to order by (identity, similarity or similarity score)
        self.order_by = 'similarity'

        # name of custom segment, where individually selected postions are collected
        self.custom_segment_label = 'Custom'

        # gap symbols
        self.gaps = ['-', '_']

        # when true, gaps at the beginning or end of a segment have a different symbol than other gaps
        self.show_padding = True

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
        # fetch all protein conformations
        protein_conformations = ProteinConformation.objects.order_by('protein__family__slug',
            'protein__entry_name').filter(protein__in=proteins,
            state__slug__in=self.states).select_related('protein__residue_numbering_scheme', 'protein__species',
            'state')
        pconfs = OrderedDict()
        for pconf in protein_conformations:
            pconf_label = pconf.__str__()
            if pconf_label not in pconfs:
                pconfs[pconf_label] = {}
            pconfs[pconf_label] = pconf

        for pconf_label, pconf in pconfs.items():
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
        selected_residue_positions = []
        for s in segments:
            # handle residue positions separately
            if s.__class__.__name__ == 'ResidueGenericNumber':
                selected_residue_positions.append(s)
                continue
            # fetch split segments (e.g. ICL3_1 and ICL3_2)
            alt_segments = ProteinSegment.objects.filter(slug__startswith=s.slug)
            for segment in alt_segments:
                segment_positions = ResidueGenericNumber.objects.filter(protein_segment=segment,
                    scheme=self.default_numbering_scheme).order_by('label')
                
                # generic numbers in the schemes of all selected proteins
                self.load_generic_numbers(segment.slug, segment_positions)

                # segments
                self.segments[segment.slug] = []
                for segment_residue in segment_positions:
                    self.segments[segment.slug].append(segment_residue.label)
        
        # combine individual residue positions into a custom segment
        if selected_residue_positions:
            if self.custom_segment_label not in self.segments:
                self.segments[self.custom_segment_label] = []
            for residue_position in selected_residue_positions:
                self.segments[self.custom_segment_label].append(residue_position.label)
            self.load_generic_numbers(self.custom_segment_label, selected_residue_positions)

    def load_segments_from_selection(self, simple_selection):
        """Read user selection and add selected protein segments/residue positions"""
        # local segment list
        segments = []

        # read selection
        for segment in simple_selection.segments:
            segments.append(segment.item)

        # load segment positions
        self.load_segments(segments)
    
    def load_generic_numbers(self, segment_slug, residues):
        """Loads generic numbers in the schemes of all selected proteins"""
        for ns in self.numbering_schemes:
            if ns[0] not in self.generic_numbers:
                self.generic_numbers[ns[0]] = OrderedDict()
            self.generic_numbers[ns[0]][segment_slug] = OrderedDict()
            for segment_residue in residues:
                self.generic_numbers[ns[0]][segment_slug][segment_residue.label] = []

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
        # fetch segment residues
        if len(self.numbering_schemes) > 1:
            rs = Residue.objects.filter(
                protein_segment__slug__in=self.segments, protein_conformation__in=self.proteins).prefetch_related(
                'protein_conformation__protein', 'protein_conformation__state', 'protein_segment',
                'generic_number__scheme', 'display_generic_number__scheme', 'alternative_generic_numbers__scheme')
        else:
            rs = Residue.objects.filter(
                protein_segment__slug__in=self.segments, protein_conformation__in=self.proteins).prefetch_related(
                'protein_conformation__protein', 'protein_conformation__state', 'protein_segment',
                'generic_number__scheme', 'display_generic_number__scheme')
        
        # fetch individually selected residues (Custom segment)
        if self.custom_segment_label in self.segments:
            crs = Residue.objects.filter(
                generic_number__label__in=self.segments[self.custom_segment_label],
                protein_conformation__in=self.proteins).prefetch_related('protein_conformation__protein',
                'protein_conformation__state', 'protein_segment', 'generic_number__scheme',
                'display_generic_number__scheme')

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

            # identifiers for protein/state
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
        
        # individually selected residues (Custom segment)
        if self.custom_segment_label in self.segments:
            for r in crs:
                ps = self.custom_segment_label
                pcid = r.protein_conformation.protein.entry_name + "-" + r.protein_conformation.state.slug
                if pcid not in proteins:
                    proteins[pcid] = {}
                if ps not in proteins[pcid]:
                    proteins[pcid][ps] = {}
                proteins[pcid][ps][r.generic_number.label] = r

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
                ns_slug = pc.protein.residue_numbering_scheme.slug

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
                            if r.display_generic_number.label not in self.generic_numbers[ns_slug][segment][pos]:
                                self.generic_numbers[ns_slug][segment][pos].append(r.display_generic_number.label)
                        else:
                            if pos not in self.generic_numbers[ns_slug][segment]:
                                self.generic_numbers[ns_slug][segment][pos] = []

                        # add display numbers for other numbering schemes of selected proteins
                        if r.generic_number and len(self.numbering_schemes) > 1:
                            for arn in r.alternative_generic_numbers.all():
                                for ns in self.numbering_schemes:
                                    if (arn.scheme.slug == ns[0] and arn.scheme.slug != ns):
                                        self.generic_numbers[arn.scheme.slug][segment][pos].append(arn.label)
                                        break

                        # append the residue to the matrix
                        if r.generic_number:
                            s.append([pos, r.display_generic_number.label, r.amino_acid,
                                r.display_generic_number.scheme.short_name, r.sequence_number])
                        else:
                            s.append([pos, "", r.amino_acid, "", r.sequence_number])

                        first_residue_found = True

                        # reset gap counter
                        gap_counter = 0
                    except:
                        if self.show_padding:
                            padding_symbol = '_'
                        else:
                            padding_symbol = '-'
                        if first_residue_found:
                            s.append([pos, False, '-', 0])

                            # update gap counter
                            gap_counter += 1

                            # if this is the last residue and there are gaps and the end of the segment, update them to
                            # end gaps
                            if self.show_padding:
                                if (position_counter) == len(positions):
                                    for i in range(gap_counter):
                                        s[len(positions)-(i+1)][2] = padding_symbol
                        else:
                            s.append([pos, False, padding_symbol, 0])
                    
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

    def calculate_statistics(self):
        """Calculate consesus sequence and amino acid and feature frequency"""
        aa_count = []
        feature_count = []
        most_freq_aa = []
        amino_acids = OrderedDict([(a, 0) for a in AMINO_ACIDS]) # from common.definitions
        self.amino_acids = AMINO_ACIDS.keys()
        features = OrderedDict([(a, 0) for a in AMINO_ACID_GROUPS])
        self.features = AMINO_ACID_GROUP_NAMES.values()
        for i, p in enumerate(self.proteins):
            for j, s in enumerate(p.alignment):
                if i == 0:
                    aa_count.append(OrderedDict())
                    feature_count.append(OrderedDict())
                    most_freq_aa.append(OrderedDict())
                for p in s:
                    generic_number = p[0]
                    amino_acid = p[2]
                    
                    # init counters
                    if generic_number not in aa_count[j]:
                        aa_count[j][generic_number] = amino_acids.copy()
                    if generic_number not in feature_count[j]:
                        feature_count[j][generic_number] = features.copy()
                    if generic_number not in most_freq_aa[j]:
                        most_freq_aa[j][generic_number] = [[], 0]

                    # stop here if this is gapped position (no need to collect stats on those)
                    if amino_acid in self.gaps:
                        continue

                    # update amino acid counter for this generic number
                    aa_count[j][generic_number][amino_acid] += 1

                    # update feature counter for this generic number
                    for feature, members in AMINO_ACID_GROUPS.items():
                        if amino_acid in members:
                            feature_count[j][generic_number][feature] += 1

                    # update most frequent amino_acids for this generic number
                    if aa_count[j][generic_number][amino_acid] > most_freq_aa[j][generic_number][1]:
                        most_freq_aa[j][generic_number] = [[amino_acid], aa_count[j][generic_number][amino_acid]]
                    elif aa_count[j][generic_number][amino_acid] == most_freq_aa[j][generic_number][1]:
                        if amino_acid not in most_freq_aa[j][generic_number][0]:
                            most_freq_aa[j][generic_number][0].append(amino_acid)


        # merge the amino acid counts into a consensus sequence
        num_proteins = len(self.proteins)
        for i, s in enumerate(most_freq_aa):
            self.consensus.append(OrderedDict())
            for p, r in s.items():
                conservation = str(round(r[1]/num_proteins*100))
                if len(conservation) == 1:
                    cons_interval = '0'
                else:
                    # the intervals are defined as 0-10, where 0 is 0-9, 1 is 10-19 etc. Used for colors.
                    cons_interval = conservation[:-1]
                num_freq_aa = len(r[0])
                if num_freq_aa == 1:
                    self.consensus[i][p] = [r[0][0], cons_interval,
                    r[0][0] + ' ' + str(round(r[1]/num_proteins*100)) + '%']
                elif num_freq_aa > 1:
                    self.consensus[i][p] = ['+', cons_interval,
                    '/'.join(r[0]) + ' ' + str(round(r[1]/num_proteins*100)) + '%']

        # process amino acid frequency
        for i, amino_acid in enumerate(AMINO_ACIDS):
            self.amino_acid_stats.append([])
            for j, segment_num in enumerate(aa_count):
                self.amino_acid_stats[i].append([])
                k = 0
                for gn, aas in segment_num.items():
                    self.amino_acid_stats[i][j].append([])
                    for aa, freq in aas.items():
                        if aa == amino_acid:
                            frequency = str(round(freq/num_proteins*100))
                            if len(frequency) == 1:
                                freq_interval = '0'
                            else:
                                # intervals defined in the same way as for the consensus sequence
                                freq_interval = frequency[:-1]
                            self.amino_acid_stats[i][j][k] = [frequency, freq_interval]
                    k += 1

        # process feature frequency
        for i, feature in enumerate(AMINO_ACID_GROUPS):
            self.feature_stats.append([])
            for j, segment_num in enumerate(feature_count):
                self.feature_stats[i].append([])
                k = 0
                for gn, fs in segment_num.items():
                    self.feature_stats[i][j].append([])
                    for f, freq in fs.items():
                        if f == feature:
                            frequency = str(round(freq/num_proteins*100))
                            if len(frequency) == 1:
                                freq_interval = '0'
                            else:
                                # intervals defined in the same way as for the consensus sequence
                                freq_interval = frequency[:-1]
                            self.feature_stats[i][j][k] = [frequency, freq_interval]
                    k += 1

    def calculate_similarity(self):
        """Calculate the sequence identity/similarity of every selected protein compared to a selected reference"""
        for i, protein in enumerate(self.proteins):
            # skip the first row, as it is the reference
            if i == 0:
                continue

            identities = []
            similarities = []
            similarity_scores = []
            for j, s in enumerate(protein.alignment):
                for k, p in enumerate(s):
                    reference_residue = self.proteins[0].alignment[j][k][2]
                    protein_residue = self.proteins[i].alignment[j][k][2]
                    if not (reference_residue in self.gaps and protein_residue in self.gaps):
                        # identity
                        if protein_residue == reference_residue:
                            identities.append(1)
                        else:
                            identities.append(0)

                        # similarity
                        if reference_residue in self.gaps or protein_residue in self.gaps:
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
        order_by_value = int(getattr(self.proteins[0], self.order_by))
        if order_by_value:
            self.proteins.sort(key=lambda x: getattr(x, self.order_by), reverse=True)
        self.proteins.insert(0, ref)

    def score_match(self, pair, matrix):
        if pair not in matrix:
            return matrix[(tuple(reversed(pair)))]
        else:
            return matrix[pair]

class AlignedReferenceTemplate(Alignment):
    ''' Creates a structure based alignment between reference protein and target proteins that are made up from the 
        best available unique structures. It marks the best match as the main template structure.

        @param reference_protein: Protein object of reference protein. \n
        @param segments: list of segment ids to be considered in the alignment, e.g. ['TM1','TM2']. \n
        @param query_states: list of activation sites considered. \n
        @param order_by: str of ordering the aligned proteins. Identity, similarity or simscore.
        @param provide_main_temlpate_structure: Structure object, use only when aligning loops and when the main 
        template is already known.
    '''
    def __init__(self, reference_protein, segments, query_states, order_by, provide_main_template_structure=None,
                 provide_similarity_table=None):
        super(AlignedReferenceTemplate, self).__init__()
        self.query_states = query_states
        self.order_by = order_by
        self.reference_protein = reference_protein
        self.load_reference_protein(reference_protein)
        self.load_proteins_by_structure()
        self.load_segments(ProteinSegment.objects.filter(slug__in=segments))
        self.build_alignment()
        self.calculate_similarity()
        self.reference_protein = self.proteins[0]
        self.main_template_protein = None
        self.ordered_proteins = []
        if provide_main_template_structure==None:
            self.main_template_structure = None
        else:
            self.main_template_structure = provide_main_template_structure
        segment_type = [str(x)[:2] for x in segments]
        if provide_similarity_table==None:
            self.provide_similarity_table = None
        else:
            self.provide_similarity_table = provide_similarity_table
        if 'TM' in segment_type and ('IC' not in segment_type or 'EC' not in segment_type):
            self.similarity_table = self.create_helix_similarity_table()
        elif 'IC' in segment_type or 'EC' in segment_type and 'TM' not in segment_type:
            self.loop_table = OrderedDict()            
            self.similarity_table = self.create_loop_similarity_table()
        if self.main_template_structure==None:
            self.main_template_structure = self.get_main_template()
        self.segment_labels = segments
        self.reference_dict = OrderedDict()
        self.template_dict = OrderedDict()
        self.alignment_dict = OrderedDict()
       
    def __repr__(self):
        return '<AlignedReferenceTemplate: Reference: {} ; Template: {}>'.format(self.reference_protein.protein.entry_name, 
                                                                                 self.main_template_structure)

    def load_proteins_by_structure(self):
        ''' Loads proteins into alignment based on available structures in the database.
        '''
        self.structures_data = Structure.objects.filter(state__name__in=self.query_states).order_by(
            'protein_conformation__protein__parent','resolution').distinct('protein_conformation__protein__parent')#.exclude(protein_conformation__protein__parent__family__parent=self.reference_protein.family.parent_id)
        self.load_proteins(
            [Protein.objects.get(id=target.protein_conformation.protein.parent.id) for target in self.structures_data])
  
    def get_main_template(self):
        ''' Returns main template structure after checking for matching helix start and end positions.
        '''
        try:
            for protein in self.ordered_proteins:
                if protein.protein==self.reference_protein.protein:
                    ref_positions = []
                    for segment in protein.alignment:
                        all_ref_positions = []
                        for res in segment:
                            if res[1]!=False and res[1]!='':
                               all_ref_positions.append(res[0])
                        ref_positions.append([all_ref_positions[0],all_ref_positions[-1]])
                else:
                    temp_positions = []
                    for segment in protein.alignment:
                        all_temp_positions = []
                        for res in segment:
                            if res[1]!=False and res[1]!='':
                               all_temp_positions.append(res[0])
                        temp_positions.append([all_temp_positions[0],all_temp_positions[-1]])
                    print(protein)
                    print(ref_positions)
                    print(temp_positions)
                    if ref_positions==temp_positions:
                        for struct in self.similarity_table:
                            if protein.protein==struct.protein_conformation.protein.parent:
                                self.main_template_protein = protein
                                return struct
        except:
            return None

    def create_helix_similarity_table(self):
        ''' Creates an ordered dictionary of structure objects, where templates are sorted by similarity and resolution.
        '''
        temp_list = []
        self.ordered_proteins = [self.proteins[0]]
        similarity_table = OrderedDict()
        for protein in self.proteins:
            if protein.protein!=self.reference_protein.protein:
                matches = self.structures_data.filter(protein_conformation__protein__parent__id=protein.protein.id)
                temp_list.append((list(matches)[0], int(protein.similarity), float(list(matches)[0].resolution), protein))
        sorted_list = sorted(temp_list, key=lambda x: (-x[1],x[2]))
        for i in sorted_list:
            similarity_table[i[0]] = i[1]
            self.ordered_proteins.append(i[3])
        return similarity_table

    def create_loop_similarity_table(self):
        ''' Creates an ordered dictionary of structure objects, where templates are sorted by similarity and resolution.
            Only templates that have the same loop length as the reference are considered. If the provided main 
            structure has the same length as the reference, no other structure is inserted to the table.
        '''
        temp_list = []
        similarity_table = OrderedDict()
        self.main_template_protein = [p for p in self.proteins if 
                                    p.protein==self.main_template_structure.protein_conformation.protein.parent][0]
        for protein in self.proteins:
            if protein.protein==self.reference_protein.protein:
                ref_length = []
                for res in protein.alignment[0]:
                    if res[1]!=False:
                        ref_length.append(res[0])
            elif protein.protein==self.main_template_protein.protein:
                main_temp_length = []
                main_struct_sim = int(protein.similarity)
                for res in protein.alignment[0]:
                    if res[1]!=False:
                        main_temp_length.append(res[0])
            else:
                temp_length = []
                for res in protein.alignment[0]:
                    if res[1]!=False:
                        temp_length.append(res[0])
                if self.provide_similarity_table==None:
                    match = self.structures_data.filter(protein_conformation__protein__parent__id=protein.protein.id)
                    temp_list.append((list(match)[0], temp_length, int(protein.similarity), 
                                      float(list(match)[0].resolution), protein))
                else:
                    for struct, score in self.provide_similarity_table.items():
                        if struct.protein_conformation.protein.parent.id==protein.protein.id:
                            match = struct
                            break
                    temp_list.append((match, temp_length, score, float(match.resolution), protein))
        if len(ref_length)!=len(main_temp_length):
            alt_temps = [entry for entry in temp_list if len(entry[1])==len(ref_length)]
            sorted_list = sorted(alt_temps, key=lambda x: (-x[2],x[3]))
            for i in sorted_list:
                similarity_table[i[0]] = i[2]
            try:
                self.main_template_protein = sorted_list[0][4]
                self.main_template_structure = sorted_list[0][0]
                self.loop_table = similarity_table
            except:
                self.main_template_protein = None
                self.main_template_structure = None
                self.loop_table = None
                return None            
        else:
            if self.provide_similarity_table==None:
                similarity_table[self.main_template_structure] = main_struct_sim
            else:
                similarity_table[self.main_template_structure] = self.provide_similarity_table[
                                                                                        self.main_template_structure]
            self.loop_table = similarity_table
        return similarity_table


    def enhance_best_alignment(self):
        ''' Creates an alignment between reference and main_template where matching residues are depicted with the 
            one-letter residue code, mismatches with '.', gaps with '-', gaps due to shorter sequences with 'x'.
        '''
        if not self.main_template_protein: raise AssertionError(
        "No main template with same helix endings. No homology model will be built for {}.".format(self.reference_protein))
        segment_count = 0
        for ref_segment, temp_segment, segment_label in zip(self.reference_protein.alignment,
                                                            self.main_template_protein.alignment, self.segment_labels):
            segment_count+=1
            ref_segment_dict,temp_segment_dict,align_segment_dict = OrderedDict(), OrderedDict(), OrderedDict()
            for ref_position, temp_position in zip(ref_segment,temp_segment):
                if ref_position[1]!=False and temp_position[1]!=False:
                    if ref_position[0]==temp_position[0]:
                        ref_segment_dict[ref_position[0]]=ref_position[2]
                        temp_segment_dict[temp_position[0]]=temp_position[2]
                        if ref_position[2]==temp_position[2]:
                            align_segment_dict[ref_position[0]]=ref_position[2]
                        else:
                            align_segment_dict[ref_position[0]]='.'
                    else:
                        print("Error: Generic numbers don't align")
                            
                elif ref_position[1]!=False and temp_position[1]==False:
                    ref_segment_dict[ref_position[0]]=ref_position[2]                    
                    if temp_position[2]=='-':
                        temp_segment_dict[temp_position[0]]='-'
                        align_segment_dict[temp_position[0]]='-'
                    elif temp_position[2]=='_':
                        temp_segment_dict[temp_position[0]]='x'
                        align_segment_dict[temp_position[0]]='x'
                        
                elif ref_position[2]=='-' and temp_position[1]!=False:
                    ref_segment_dict[ref_position[0]]='-'
                    temp_segment_dict[temp_position[0]]=temp_position[2]
                    align_segment_dict[ref_position[0]]='-'

            self.reference_dict[segment_label] = ref_segment_dict
            self.template_dict[segment_label] = temp_segment_dict
            self.alignment_dict[segment_label] = align_segment_dict
