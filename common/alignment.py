from django.conf import settings

from common.selection import Selection
from common.definitions import *
from protein.models import Protein, ProteinConformation, ProteinState, ProteinSegment, ProteinFusionProtein
from residue.models import Residue
from residue.models import ResidueGenericNumber, ResidueGenericNumberEquivalent
from residue.models import ResidueNumberingScheme
from structure.models import Structure, Rotamer, StructureSegment, StructureSegmentModeling, StructureCoordinates

from collections import OrderedDict
from copy import deepcopy
from operator import itemgetter
from Bio.SubsMat import MatrixInfo
import logging


class Alignment:
    """A class representing a protein sequence alignment, with or without a reference sequence"""
    def __init__(self):
        self.reference = False
        self.proteins = []
        self.non_matching_proteins = [] # proteins that do not match user specified site definitions
        self.segments = OrderedDict()
        self.numbering_schemes = {}
        self.generic_numbers = OrderedDict()
        self.generic_number_objs = {}
        self.positions = []
        self.consensus = OrderedDict()
        self.forced_consensus = OrderedDict() # consensus sequence where all conflicts are solved by rules
        self.full_consensus = [] # consensus sequence with full residue objects
        self.similarity_matrix = OrderedDict()
        self.amino_acids = []
        self.amino_acid_stats = []
        self.aa_count = OrderedDict()
        self.aa_count_with_protein = OrderedDict()
        self.features = []
        self.feature_stats = []
        self.default_numbering_scheme = ResidueNumberingScheme.objects.get(slug=settings.DEFAULT_NUMBERING_SCHEME)
        self.states = [settings.DEFAULT_PROTEIN_STATE] # inactive, active etc
        self.use_residue_groups = False
        self.ignore_alternative_residue_numbering_schemes = False # set to true if no numbering is to be displayed
        
        # refers to which ProteinConformation attribute to order by (identity, similarity or similarity score)
        self.order_by = 'similarity'

        # name of custom segment, where individually selected postions are collected
        self.custom_segment_label = 'Custom'
        self.interaction_segment_label = 'I'

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
        # FIXME take many conformational states into account
        try:
            pconf = ProteinConformation.objects.get(protein=protein)
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
        # fetch all conformations of selected proteins
        # FIXME only show inactive?
        protein_conformations = ProteinConformation.objects.order_by('protein__family__slug',
            'protein__entry_name').filter(protein__in=proteins).select_related('protein__residue_numbering_scheme',
            'protein__species', 'state')
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
                
                if species_list:
                    family_proteins = Protein.objects.filter(family__slug__startswith=target.item.slug,
                        species__in=(species_list), source__in=(protein_source_list)).select_related(
                        'residue_numbering_scheme', 'species')
                else:
                    family_proteins = Protein.objects.filter(family__slug__startswith=target.item.slug,
                        source__in=(protein_source_list)).select_related('residue_numbering_scheme', 'species')

                for fp in family_proteins:
                    proteins.append(fp)

        # load protein list
        self.load_proteins(proteins)

    def load_segments(self, selected_segments):
        selected_residue_positions = []
        for s in selected_segments:
            if hasattr(s, 'item'):
                selected_segment = s.item
            else:
                selected_segment = s
                
            # handle residue positions separately
            if selected_segment.__class__.__name__ == 'ResidueGenericNumberEquivalent':
                segment_residue = [selected_segment]
                
                # if in site search, add residue group to selected item
                if 'site_residue_group' in s.properties:
                    segment_residue.append(s.properties['site_residue_group'])
                selected_residue_positions.append(segment_residue)
                continue
            
            # fetch split segments (e.g. ECL2_before and ECL2_after)
            alt_segments = ProteinSegment.objects.filter(slug__startswith=selected_segment.slug)

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
            if self.use_residue_groups:
                interaction_groups = {}
                for residue_position in selected_residue_positions:
                    segment_slug = self.interaction_segment_label + str(residue_position[1])
                    if segment_slug not in self.segments:
                        self.segments[segment_slug] = []
                    if segment_slug not in interaction_groups:
                        interaction_groups[segment_slug] = []
                    self.segments[segment_slug].append(residue_position[0].default_generic_number.label)
                    interaction_groups[segment_slug].append(residue_position)
                for segment_slug, group_residues in interaction_groups.items():
                    self.load_generic_numbers(segment_slug, group_residues)
            else:
                segment_slug = self.custom_segment_label
                if segment_slug not in self.segments:
                    self.segments[segment_slug] = []
                for residue_position in selected_residue_positions:
                    self.segments[segment_slug].append(residue_position[0].default_generic_number.label)
                self.load_generic_numbers(segment_slug, selected_residue_positions)

    def load_segments_from_selection(self, simple_selection):
        """Read user selection and add selected protein segments/residue positions"""
        # local segment list
        segments = []

        # read selection
        for segment in simple_selection.segments:
            segments.append(segment)

        # load segment positions
        self.load_segments(segments)
    
    def load_generic_numbers(self, segment_slug, residues):
        """Loads generic numbers in the schemes of all selected proteins"""
        for ns in self.numbering_schemes:
            if ns[0] not in self.generic_numbers:
                self.generic_numbers[ns[0]] = OrderedDict()
            self.generic_numbers[ns[0]][segment_slug] = OrderedDict()
            for segment_residue in residues:
                if segment_slug == self.custom_segment_label or self.use_residue_groups:
                    residue_position = segment_residue[0].default_generic_number.label
                else:
                    residue_position = segment_residue.label
                self.generic_numbers[ns[0]][segment_slug][residue_position] = []

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
        if not self.ignore_alternative_residue_numbering_schemes and len(self.numbering_schemes) > 1:
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
        crs = {}
        for segment in self.segments:
            if segment == self.custom_segment_label or self.use_residue_groups:
                crs[segment] = Residue.objects.filter(
                    generic_number__label__in=self.segments[segment],
                    protein_conformation__in=self.proteins).prefetch_related('protein_conformation__protein',
                    'protein_conformation__state', 'protein_segment', 'generic_number__scheme',
                    'display_generic_number__scheme')

        # create a dict of proteins, segments and residues
        proteins = {}
        segment_counters = {}
        aligned_residue_encountered = {}
        fusion_protein_inserted = {}
        for r in rs:
            ps = r.protein_segment.slug

            # identifiers for protein/state
            pcid = r.protein_conformation.protein.entry_name + "-" + r.protein_conformation.state.slug
            
            # update protein dict
            if pcid not in proteins:
                proteins[pcid] = {}
            if ps not in proteins[pcid]:
                proteins[pcid][ps] = {}
            
            # update aligned residue tracker
            if pcid not in aligned_residue_encountered:
                aligned_residue_encountered[pcid] = {}
            if ps not in aligned_residue_encountered[pcid]:
                aligned_residue_encountered[pcid][ps] = False

            # what part of the segment is this? There are 4 possibilities:
            # 1. The aligned part (for both fully and partially aligned segments)
            # 2. The part before the aligned part in a partially aligned segment
            # 3. The part after the aligned part in a partially aligned segment
            # 4. An unaligned segment (then there is only one part)
            if r.generic_number:
                segment_part = 1
            elif ps in settings.REFERENCE_POSITIONS and not aligned_residue_encountered[pcid][ps]:
                segment_part = 2
            elif ps in settings.REFERENCE_POSITIONS and aligned_residue_encountered[pcid][ps]:
                segment_part = 3
            else:
                segment_part = 4

            # update segment counters
            if pcid not in segment_counters:
                segment_counters[pcid] = {}
            if segment_part == 3:
                part_ps = ps + '_after'
            else:
                part_ps = ps
            if part_ps not in segment_counters[pcid]:
                segment_counters[pcid][part_ps] = 1
            else:
                segment_counters[pcid][part_ps] += 1
            

            # update fusion protein tracker
            if pcid not in fusion_protein_inserted:
                fusion_protein_inserted[pcid] = {}
            if ps not in fusion_protein_inserted[pcid]:
                fusion_protein_inserted[pcid][ps] = False

            # user generic numbers as keys for aligned segments
            if r.generic_number:
                proteins[pcid][ps][r.generic_number.label] = r

                # register the presence of an aligned residue
                aligned_residue_encountered[pcid][ps] = True
            # use custom keys for non-aligned segments
            else:
                # label prefix + index
                # Unaligned segments should be split in the middle, with the first part "left aligned", and the second
                # "right aligned". If there is an aligned part of the segment, it goes in the middle.
                if segment_part == 2:
                    prefix = '00-'
                elif segment_part == 3:
                    prefix = 'zz-'
                else:
                    prefix = '01-'

                # Note that there is not enough information to assign correct indicies to "right aligned" residues, but
                # those are corrected below
                index = str("%04d" % (segment_counters[pcid][part_ps],))

                # position label
                pos_label =  prefix + ps + "-" + index

                # insert fusion protein FIXME add this
                # if not fusion_protein_inserted[pcid][ps] and aligned_residue_encountered[pcid][ps]:
                #     fp = ProteinFusionProtein.objects.get(protein=r.protein_conformation.protein,
                #         segment_after=r.protein_segment)
                #     fusion_pos_label = ps + "-" + str("%04d" % (segment_counters[pcid][ps]-1,)) + "-fusion"
                #     proteins[pcid][ps][fusion_pos_label] = Residue(amino_acid=fp.protein_fusion.name)
                #     if fusion_pos_label not in self.segments[ps]:
                #         self.segments[ps].append(fusion_pos_label)
                #     fusion_protein_inserted[pcid][ps] = True

                # residue
                proteins[pcid][ps][pos_label] = r

        # correct alignment of split segments
        for pcid, segments in proteins.items():
            for ps, positions in segments.items():
                pos_num = 1
                pos_num_after = 1
                for pos_label in sorted(positions):
                    res_obj = proteins[pcid][ps][pos_label]
                    right_align = False
                    # In a "normal", non split, unaligned segment, is this past the middle?
                    if (pos_label.startswith('01-')
                        and res_obj.protein_segment.category != 'terminus'
                        and pos_num > (segment_counters[pcid][ps] / 2 + 0.5)):
                        right_align = True
                    # In an partially aligned segment (prefixed with 00), where conserved residues are lacking, treat
                    # as an unaligned segment
                    elif (pos_label.startswith('00-')
                        and not aligned_residue_encountered[pcid][ps]
                        and pos_num > (segment_counters[pcid][ps] / 2 + 0.5)
                        or res_obj.protein_segment.slug == 'N-term'):
                        right_align = True
                    # In an N-terminus, always right align everything
                    elif pos_label.startswith('01-') and res_obj.protein_segment.slug == 'N-term':
                        right_align = True

                    if right_align:
                        # if so, "right align" from here using a zz prefixed label
                        updated_index = 'zz' + pos_label[2:]
                        proteins[pcid][ps][updated_index] = proteins[pcid][ps].pop(pos_label)
                        pos_label = updated_index

                    if pos_label.startswith('zz-'):
                        segment_label_after = ps + '_after' # parts after a partly aligned segment start with zz
                        if segment_label_after in segment_counters[pcid]:
                            segment_length = segment_counters[pcid][segment_label_after]
                            counter = pos_num_after
                        
                        # this might be the "second part" of an unaligned segment, e.g.
                        # AAAA----AAAAA
                        # AAAAAAAAAAAAA
                        else: 
                            segment_length = segment_counters[pcid][ps]
                            counter = pos_num

                        updated_index = pos_label[:-4] + str(9999 - (segment_length - counter))
                        proteins[pcid][ps][updated_index] = proteins[pcid][ps].pop(pos_label)
                        pos_label = updated_index
                        pos_num_after += 1
                    if pos_label not in self.segments[ps]:
                        self.segments[ps].append(pos_label)
                    pos_num += 1
        
        # individually selected residues (Custom segment)
        for segment in self.segments:
            if segment == self.custom_segment_label or self.use_residue_groups:
                for r in crs[segment]:
                    ps = segment
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

        for pc in self.proteins:
            row = OrderedDict()
            row_list = [] # FIXME redundant, remove when dependecies are removed
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
                        if (not self.ignore_alternative_residue_numbering_schemes and len(self.numbering_schemes) > 1):
                            if r.generic_number:
                                for arn in r.alternative_generic_numbers.all():
                                    for ns in self.numbering_schemes:
                                        if (arn.scheme.slug == ns[0] and arn.scheme.slug != ns_slug):
                                            self.generic_numbers[arn.scheme.slug][segment][pos].append(arn.label)
                                            break
                            else:
                                for ns in self.numbering_schemes:
                                    if pos not in self.generic_numbers[ns[0]][segment] and ns[0] != ns_slug:
                                        self.generic_numbers[ns[0]][segment][pos] = []

                        # append the residue to the matrix
                        if r.generic_number:
                            # s.append([pos, r.display_generic_number.label, r.amino_acid,
                            #   r.display_generic_number.scheme.short_name, r.sequence_number])

                            s.append([pos, r.display_generic_number.label, r.amino_acid,
                                r.display_generic_number.scheme.short_name, r.sequence_number, r.generic_number.label])


                            # update generic residue object dict
                            if pos not in self.generic_number_objs:
                                self.generic_number_objs[pos] = r.display_generic_number
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
                row[segment] = s
                row_list.append(s) # FIXME redundant, remove when dependecies are removed
            pc.alignment = row
            pc.alignment_list = row_list # FIXME redundant, remove when dependecies are removed
        self.sort_generic_numbers()
        self.merge_generic_numbers()
        self.clear_empty_positions()

    def clear_empty_positions(self):
        """Remove empty columns from the segments and matrix"""
        # segments and  
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
            for j, s in protein.alignment.items():
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

    def sort_generic_numbers(self):
        # remove split segments from generic numbers list
        for ns, gs in self.generic_numbers.items():
            for segment, positions in gs.items():
                s = segment.split("_")
                if len(s) > 1:
                    del self.generic_numbers[ns][segment]
                else:
                    ordered_generic_numbers = OrderedDict()
                    for gn in sorted(self.generic_numbers[ns][segment], key=lambda x: x.split('x')):
                        ordered_generic_numbers[gn] = self.generic_numbers[ns][segment][gn]
                    self.generic_numbers[ns][segment] = ordered_generic_numbers 

    def format_generic_number(self, generic_number):
        """A placeholder for an instance specific function"""
        return generic_number

    def calculate_statistics(self):
        """Calculate consesus sequence and amino acid and feature frequency"""
        feature_count = OrderedDict()
        most_freq_aa = OrderedDict()
        amino_acids = OrderedDict([(a, 0) for a in AMINO_ACIDS]) # from common.definitions
        self.amino_acids = AMINO_ACIDS.keys()
        features = OrderedDict([(a, 0) for a in AMINO_ACID_GROUPS])
        self.features = AMINO_ACID_GROUP_NAMES.values()
        for i, p in enumerate(self.proteins):
            entry_name = p.protein.entry_name
            for j, s in p.alignment.items():
                if i == 0:
                    self.aa_count[j] = OrderedDict()
                    feature_count[j] = OrderedDict()
                    most_freq_aa[j] = OrderedDict()
                for p in s:
                    generic_number = p[0]
                    amino_acid = p[2]

                    # stop here if this is gapped position (no need to collect stats on those)
                    if amino_acid in self.gaps:
                        continue
                    
                    # init counters
                    if generic_number not in self.aa_count[j]:
                        self.aa_count[j][generic_number] = amino_acids.copy()
                        if generic_number in self.generic_number_objs:
                            self.aa_count_with_protein[generic_number] = {}
                    if generic_number not in feature_count[j]:
                        feature_count[j][generic_number] = features.copy()
                    if generic_number not in most_freq_aa[j]:
                        most_freq_aa[j][generic_number] = [[], 0]

                    # update amino acid counter for this generic number
                    self.aa_count[j][generic_number][amino_acid] += 1
                    if generic_number in self.generic_number_objs:
                        if amino_acid not in self.aa_count_with_protein[generic_number]:
                            self.aa_count_with_protein[generic_number][amino_acid] = []
                        if entry_name not in self.aa_count_with_protein[generic_number][amino_acid]:
                            self.aa_count_with_protein[generic_number][amino_acid].append(entry_name)

                    # update feature counter for this generic number
                    for feature, members in AMINO_ACID_GROUPS.items():
                        if amino_acid in members:
                            feature_count[j][generic_number][feature] += 1

                    # update most frequent amino_acids for this generic number
                    if self.aa_count[j][generic_number][amino_acid] > most_freq_aa[j][generic_number][1]:
                        most_freq_aa[j][generic_number] = [[amino_acid], self.aa_count[j][generic_number][amino_acid]]
                    elif self.aa_count[j][generic_number][amino_acid] == most_freq_aa[j][generic_number][1]:
                        if amino_acid not in most_freq_aa[j][generic_number][0]:
                            most_freq_aa[j][generic_number][0].append(amino_acid)

        # merge the amino acid counts into a consensus sequence
        num_proteins = len(self.proteins)
        sequence_counter = 1
        for i, s in most_freq_aa.items():
            self.consensus[i] = OrderedDict()
            self.forced_consensus[i] = OrderedDict()
            for p, r in s.items():
                conservation = str(round(r[1]/num_proteins*100))
                if len(conservation) == 1:
                    cons_interval = '0'
                else:
                    # the intervals are defined as 0-10, where 0 is 0-9, 1 is 10-19 etc. Used for colors.
                    cons_interval = conservation[:-1]
                

                # forced consensus sequence uses the first residue to break ties
                self.forced_consensus[i][p] = r[0][0]

                # consensus sequence displays + in tie situations
                num_freq_aa = len(r[0])
                if num_freq_aa == 1:
                    self.consensus[i][p] = [r[0][0], cons_interval,
                    r[0][0] + ' ' + str(round(r[1]/num_proteins*100)) + '%']
                elif num_freq_aa > 1:
                    self.consensus[i][p] = ['+', cons_interval,
                    '/'.join(r[0]) + ' ' + str(round(r[1]/num_proteins*100)) + '%']

                # create a residue object full consensus
                res = Residue()
                res.sequence_number = sequence_counter
                if p in self.generic_number_objs:
                    res.display_generic_number = self.generic_number_objs[p]
                res.family_generic_number = p
                res.segment_slug = i
                res.amino_acid = r[0][0]
                res.frequency = self.consensus[i][p][2]
                self.full_consensus.append(res)

                # update sequence counter
                sequence_counter += 1

        # process amino acid frequency
        for i, amino_acid in enumerate(AMINO_ACIDS):
            self.amino_acid_stats.append([])
            j = 0
            for segment, segment_num in self.aa_count.items():
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
                j += 1

        # process feature frequency
        for i, feature in enumerate(AMINO_ACID_GROUPS):
            self.feature_stats.append([])
            j = 0
            for segment, segment_num in feature_count.items():
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
                j += 1

    def calculate_aa_count_per_generic_number(self):
        ''' Small function to return a dictionary of display_generic_number and the frequency of each AA '''
        generic_lookup_aa_freq = {}
        num_proteins = len(self.proteins)
        for j, a in self.aa_count.items():
            for g, p in a.items():
                for aa, c in p.items():
                    if g in self.generic_number_objs:
                        if self.generic_number_objs[g].label in generic_lookup_aa_freq:
                            generic_lookup_aa_freq[self.generic_number_objs[g].label][aa] = round(c/num_proteins*100)
                        else:
                            generic_lookup_aa_freq[self.generic_number_objs[g].label] = {aa: round(c/num_proteins*100) }
        return generic_lookup_aa_freq


    def calculate_similarity(self):
        """Calculate the sequence identity/similarity of every selected protein compared to a selected reference"""
        for i, protein in enumerate(self.proteins):
            # skip the first row, as it is the reference
            if i == 0:
                continue

            # calculate identity, similarity and similarity score to the reference
            calc_values = self.pairwise_similarity(self.proteins[0], self.proteins[i])
            
            # update the protein
            if calc_values:
                self.proteins[i].identity = calc_values[0]
                self.proteins[i].similarity = calc_values[1]
                self.proteins[i].similarity_score = calc_values[2]

        # order protein list by similarity score
        ref = self.proteins.pop(0)
        order_by_value = int(getattr(self.proteins[0], self.order_by))
        if order_by_value:
            self.proteins.sort(key=lambda x: getattr(x, self.order_by), reverse=True)
        self.proteins.insert(0, ref)

    def calculate_similarity_matrix(self):
        """Calculate a matrix of sequence identity/similarity for every selected protein"""
        self.similarity_matrix = OrderedDict()
        for i, protein in enumerate(self.proteins):
            protein_key = protein.protein.entry_name
            protein_name = "[" + protein.protein.species.common_name + "] " + protein.protein.name
            self.similarity_matrix[protein_key] = {'name': protein_name, 'values': []}
            for k, protein in enumerate(self.proteins):
                # calculate identity, similarity and similarity score to the reference
                calc_values = self.pairwise_similarity(self.proteins[i], self.proteins[k])
                if k == i:
                    value = '-'
                elif k < i:
                    value = calc_values[1].strip()
                elif k > i:
                    value = calc_values[0].strip()
                
                if value == '-':
                    color_class = "-"
                else:
                    if int(value) < 10:
                        color_class = 0
                    else:
                        color_class = str(value)[:-1]
                self.similarity_matrix[protein_key]['values'].append([value, color_class])

    def evaluate_sites(self, request):
        """Evaluate which user selected site definitions match each protein sequence"""
        # get simple selection from session
        simple_selection = request.session.get('selection', False)
        
        # format site definititions
        site_defs = {}
        for position in simple_selection.segments:
            if position.type == 'site_residue' and position.properties['site_residue_group']:
                group_id = position.properties['site_residue_group']
                if group_id not in site_defs:
                    # min match is the value of the first item in each groups list (hence the [0])
                    # site_defs example:
                    # {
                    #     1: {
                    #         'min_match': 2,
                    #         'positions': {
                    #             {'3x51': 'hbd', '6x50': 'pos'}
                    #         }
                    #     }
                    # }
                    site_defs[group_id] = {'min_match': simple_selection.site_residue_groups[group_id -1][0],
                        'positions': {}}

                site_defs[group_id]['positions'][position.item.label] = position.properties['feature']

        # go through all proteins and match against site definitions
        for protein in self.proteins:
            for k, segment in enumerate(protein.alignment.values(), start = 1):
                num_matched = 0
                min_match = site_defs[k]['min_match']
                for position in segment:
                    # position example: ['6x49', '6.49x49', 'L', 'GPCRdb(A)', 282, 282]
                    if position[2] in AMINO_ACID_GROUPS[site_defs[k]['positions'][position[0]]]:
                        num_matched += 1
                        if num_matched >= min_match:
                            break
                else:
                    # if the protein sequence does not match the definitions, store it in non_matching_proteins
                    self.non_matching_proteins.append(protein)
                    break

        # remove non-matching proteins from protein list
        self.proteins = [p for p in self.proteins if p not in self.non_matching_proteins]

    def pairwise_similarity(self, protein_1, protein_2):
        """Calculate the identity, similarity and similarity score between a pair of proteins"""
        identities = []
        similarities = []
        similarity_scores = []
        for j, s in protein_2.alignment.items():
            for k, p in enumerate(s):
                reference_residue = protein_1.alignment[j][k][2]
                protein_residue = protein_2.alignment[j][k][2]
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
        
        # format the calculated values
        if identities and similarities:
            identity = "{:10.0f}".format(sum(identities) / len(identities) * 100)
            similarity = "{:10.0f}".format(sum(similarities) / len(similarities) * 100)
            similarity_score = sum(similarity_scores)
            return identity, similarity, similarity_score
        else:
            return False

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
                 provide_similarity_table=None, main_pdb_array=None):
        super(AlignedReferenceTemplate, self).__init__()
        self.logger = logging.getLogger('homology_modeling')
        self.segment_labels = segments
        self.reference_protein = Protein.objects.get(entry_name=reference_protein)
        if provide_main_template_structure==None and provide_similarity_table==None:
            self.query_states = query_states
            self.order_by = order_by
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
            self.provide_main_template_structure = False
        else:
            self.main_template_structure = provide_main_template_structure
            self.provide_main_template_structure = True
        segment_type = [str(x)[:2] for x in segments]
        if provide_similarity_table==None:
            self.provide_similarity_table = None
        else:
            self.provide_similarity_table = provide_similarity_table
        if main_pdb_array!=None:
            self.main_pdb_array = main_pdb_array
        if 'TM' in segment_type: #and ('IC' not in segment_type or 'EC' not in segment_type):
            self.similarity_table = self.create_helix_similarity_table()
        elif 'IC' in segment_type or 'EC' in segment_type and 'TM' not in segment_type:
            self.loop_table = OrderedDict()            
            self.similarity_table = self.create_loop_similarity_table()
        if self.main_template_structure==None:
            self.main_template_structure = self.get_main_template()
        self.reference_dict = OrderedDict()
        self.template_dict = OrderedDict()
        self.alignment_dict = OrderedDict()
       
    def __repr__(self):
        return '<AlignedReferenceTemplate: Reference: {} ; Template: {}>'.format(self.reference_protein.protein.entry_name, 
                                                                                 self.main_template_structure)

    def load_proteins_by_structure(self):
        ''' Loads proteins into alignment based on available structures in the database.
        '''
        self.structures_data = Structure.objects.filter(
            state__name__in=self.query_states, protein_conformation__protein__parent__family__parent__parent__parent=
            self.reference_protein.family.parent.parent.parent, representative='t').order_by(
            'protein_conformation__protein__parent','resolution').distinct('protein_conformation__protein__parent')
        self.load_proteins(
            [Protein.objects.get(id=target.protein_conformation.protein.parent.id) for target in self.structures_data])
  
    def get_main_template(self):
        ''' Returns main template structure after checking for matching helix start and end positions.
        '''
        try:
            for st in self.similarity_table:
                if st.protein_conformation.protein.parent==self.ordered_proteins[1].protein:
                    self.main_template_protein = self.ordered_proteins[1]
                    return st
        except:
            pass

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
            Only templates that have the same loop length as the reference are considered.
        '''
        temp_list, temp_list1, temp_list2, temp_list_mid = [],[],[],[]
        similarity_table = OrderedDict()
        self.main_template_protein = self.main_template_structure.protein_conformation.protein.parent        
        ref_seq = Residue.objects.filter(protein_conformation__protein=self.reference_protein, 
                                         protein_segment__slug=self.segment_labels[0])                                        
        prot_conf = ProteinConformation.objects.get(protein=self.reference_protein)
        segment_order = []
        for i in list(Residue.objects.filter(protein_conformation=prot_conf)):
            if i.protein_segment.slug not in segment_order:
                segment_order.append(i.protein_segment.slug)
        prev_seg = segment_order[segment_order.index(self.segment_labels[0])-1]
        next_seg = segment_order[segment_order.index(self.segment_labels[0])+1]
        orig_before_gns = [i.replace('.','x') for i in list(self.main_pdb_array[prev_seg].keys())[-4:]]
        orig_after_gns = [j.replace('.','x') for j in list(self.main_pdb_array[next_seg].keys())[:4]]                                         
                                         
        last_before_gn = orig_before_gns[-1]
        first_after_gn = orig_after_gns[0]

        if self.segment_labels[0]=='ECL2':
            try:
                ref_ECL2 = self.ECL2_slicer(ref_seq)
            except:
                ref_ECL2 = None
        for struct, similarity in self.provide_similarity_table.items():
            if (self.segment_labels[0]=='ECL2' and ref_ECL2==None and StructureCoordinates.objects.get(structure=struct,
                                                        protein_segment__slug=self.segment_labels[0]).description.text!='Full'):
                continue
            elif (self.segment_labels[0]!='ECL2' and StructureCoordinates.objects.get(structure=struct,
                                                          protein_segment__slug=self.segment_labels[0]).description.text!='Full'):
                continue
            protein = struct.protein_conformation.protein.parent
            if protein==self.main_template_protein:
                main_template_mid_failed = False
                main_temp_seq = Residue.objects.filter(protein_conformation=struct.protein_conformation, 
                                         protein_segment__slug=self.segment_labels[0])
                try:
                    if self.segment_labels[0]=='ECL2' and ref_ECL2!=None:
                        main_temp_ECL2 = self.ECL2_slicer(main_temp_seq)
                        rota = [x for x in Rotamer.objects.filter(structure=struct, residue__in=main_temp_ECL2[1]) if x.pdbdata.pdb.startswith('COMPND')==False]
                        if len(rota)==3:
                            temp_list_mid.append((struct, 3, similarity, float(struct.resolution), protein))  
                        if len(ref_ECL2[0])==len(main_temp_ECL2[0]) and len(ref_ECL2[2])==len(main_temp_ECL2[2]):
                            temp_list1.append((struct, len(main_temp_ECL2[0]), similarity, float(struct.resolution),protein))
                            temp_list2.append((struct, len(main_temp_ECL2[2]), similarity, float(struct.resolution),protein))
                        elif len(ref_ECL2[0])==len(main_temp_ECL2[0]) and len(ref_ECL2[2])!=len(main_temp_ECL2[2]):
                            temp_list1.append((struct, len(main_temp_ECL2[0]), similarity, float(struct.resolution),protein))
                        elif len(ref_ECL2[0])!=len(main_temp_ECL2[0]) and len(ref_ECL2[2])==len(main_temp_ECL2[2]):
                            temp_list2.append((struct, len(main_temp_ECL2[2]), similarity, float(struct.resolution),protein))
                    else:
                        main_template_mid_failed = True
                        raise Exception()
                except:
                    if len(ref_seq)==len(main_temp_seq):
                        similarity_table[self.main_template_structure] = self.provide_similarity_table[
                                                                                            self.main_template_structure]
                        temp_list.append((struct, len(main_temp_seq), similarity, float(struct.resolution), protein))
            else:
                temp_length, temp_length1, temp_length2 = [],[],[]
                try:
                    alt_last_gn = Residue.objects.get(protein_conformation=struct.protein_conformation, 
                                                      generic_number__label=last_before_gn)
                    alt_first_gn= Residue.objects.get(protein_conformation=struct.protein_conformation, 
                                                      generic_number__label=first_after_gn)
                    temp_length = alt_first_gn.sequence_number-alt_last_gn.sequence_number-1
                    alt_seq = Residue.objects.filter(protein_conformation=struct.protein_conformation, 
                                           sequence_number__in=list(range(alt_last_gn.sequence_number+1,
                                                                          alt_first_gn.sequence_number)))
                    if self.segment_labels[0]=='ECL2' and ref_ECL2!=None and main_template_mid_failed==False:
                        alt_ECL2 = self.ECL2_slicer(alt_seq)
                        alt_rota = [x for x in Rotamer.objects.filter(structure=struct, residue__in=alt_ECL2[1]) if x.pdbdata.pdb.startswith('COMPND')==False]
                        if len(alt_rota)==3:
                            temp_list_mid.append((struct, 3, similarity, float(struct.resolution), protein))
                        if len(ref_ECL2[0])==len(alt_ECL2[0]) and len(ref_ECL2[2])==len(alt_ECL2[2]):
                            temp_length1 = len(alt_ECL2[0])
                            temp_length2 = len(alt_ECL2[2])
                        elif len(ref_ECL2[0])==len(alt_ECL2[0]) and len(ref_ECL2[2])!=len(alt_ECL2[2]):
                            temp_length1 = len(alt_ECL2[0])
                            temp_length2 = -1
                        elif len(ref_ECL2[0])!=len(alt_ECL2[0]) and len(ref_ECL2[2])==len(alt_ECL2[2]):
                            temp_length1 = -1
                            temp_length2 = len(alt_ECL2[2])
                        elif len(ref_ECL2[0])!=len(alt_ECL2[0]) and len(ref_ECL2[2])!=len(alt_ECL2[2]):
                            temp_length1 = -1
                            temp_length2 = -1
                    elif len(ref_seq)!=len(alt_seq):
                        continue
                    before_nums = list(range(alt_last_gn.sequence_number-3, alt_last_gn.sequence_number+1))
                    after_nums = list(range(alt_first_gn.sequence_number, alt_first_gn.sequence_number+4))
                    alt_before8 = Residue.objects.filter(protein_conformation__protein=protein,
                                                         sequence_number__in=before_nums)
                    alt_after8 = Residue.objects.filter(protein_conformation__protein=protein,
                                                         sequence_number__in=after_nums)
                    alt_before_gns = [r.generic_number.label for r in alt_before8]
                    alt_after_gns = [r.generic_number.label for r in alt_after8]
                    if orig_before_gns==alt_before_gns and orig_after_gns==alt_after_gns:
                        pass
                    else:
                        raise Exception()
                except:
                    temp_length, temp_length1, temp_length2 = -1,-1,-1
                temp_list.append((struct, temp_length, similarity, float(struct.resolution), protein))
                if self.segment_labels[0]=='ECL2' and ref_ECL2!=None:
                    temp_list1.append((struct, temp_length1, similarity, float(struct.resolution), protein))
                    temp_list2.append((struct, temp_length2, similarity, float(struct.resolution), protein))
        if self.segment_labels[0]=='ECL2' and ref_ECL2!=None:
            ECL2_1 = self.order_sim_table(temp_list1, ref_ECL2[0], OrderedDict())
            ECL2_mid = self.order_sim_table(temp_list_mid, ref_ECL2[1], OrderedDict())
            ECL2_2 = self.order_sim_table(temp_list2, ref_ECL2[2], OrderedDict())
            self.loop_table = OrderedDict([('ECL2_1',ECL2_1),('ECL2_mid',ECL2_mid),('ECL2_2',ECL2_2)])
            if len(ECL2_mid)==0:
                self.loop_table=None
            return self.loop_table
        else:
            return self.order_sim_table(temp_list, ref_seq, similarity_table)
                    
    def order_sim_table(self, temp_list, ref_seq, similarity_table):     
        alt_temps = [entry for entry in temp_list if entry[1]==len(ref_seq)]
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
        return similarity_table

    def ECL2_slicer(self, queryset):
        x50 = queryset.get(generic_number__label='45x50').sequence_number
        queryset_l = list(queryset)
        x50_i = x50-queryset_l[0].sequence_number
        ECL2_1 = queryset_l[:x50_i]
        ECL2_mid = queryset_l[x50_i:x50_i+3]
        ECL2_2 = queryset_l[x50_i+3:]
        if len(ECL2_mid)<3:
            raise AssertionError()
        return[ECL2_1,ECL2_mid,ECL2_2]
                
    def enhance_best_alignment(self):
        ''' Creates an alignment between reference and main_template where matching residues are depicted with the 
            one-letter residue code, mismatches with '.', gaps with '-', gaps due to shorter sequences with 'x'.
        '''
        if not self.main_template_protein: 
            self.logger.error(
            '''No main template with same helix endings. 
               No homology model will be built for {}.'''.format(self.reference_protein))
            return None
        
        for ref_seglab, temp_seglab in zip(self.reference_protein.alignment, self.main_template_protein.alignment):
            if 'TM' in ref_seglab or ref_seglab in ['ICL1','ECL1','ICL2','ECL2','H8']:
                ref_segment_dict,temp_segment_dict,align_segment_dict = OrderedDict(), OrderedDict(), OrderedDict()
                for ref_position, temp_position in zip(self.reference_protein.alignment[ref_seglab],
                                                       self.main_template_protein.alignment[temp_seglab]):
                    if ref_position[1]!=False and temp_position[1]!=False and ref_position[1]!='' and temp_position!='':
                        ref_segment_dict[ref_position[0]]=ref_position[2]
                        temp_segment_dict[temp_position[0]]=temp_position[2]
                        if ref_position[2]==temp_position[2]:
                            align_segment_dict[ref_position[0]]=ref_position[2]
                        else:
                            align_segment_dict[ref_position[0]]='.'
                    elif ref_position[1]=='' and temp_position[1]=='':                        
                        ref_segment_dict[str(ref_position[4])]=ref_position[2]
                        temp_segment_dict[str(temp_position[4])]=temp_position[2]
                        if ref_position[2]==temp_position[2]:
                            align_segment_dict[str(ref_position[4])]=ref_position[2]
                        else:
                            align_segment_dict[ref_position[4]]='.'                                
                    elif ref_position[1]!=False and temp_position[1]==False and ref_position[1]!='':
                        ref_segment_dict[ref_position[0]]=ref_position[2]                    
                        if temp_position[2]=='-':
                            temp_segment_dict[temp_position[0]]='-'
                            align_segment_dict[temp_position[0]]='-'
                        elif temp_position[2]=='_':
                            temp_segment_dict[temp_position[0]]='x'
                            align_segment_dict[temp_position[0]]='x'
                    elif ref_position[1]=='' and temp_position[1]==False:
                        ref_segment_dict[str(ref_position[4])]=ref_position[2]                    
                        if temp_position[2]=='-':
                            temp_segment_dict[temp_position[0]]='-'
                            align_segment_dict[temp_position[0]]='-'
                        elif temp_position[2]=='_':
                            temp_segment_dict[temp_position[0]]='x'
                            align_segment_dict[temp_position[0]]='x'
                    elif ref_position[2]=='-' and temp_position[1]!=False and temp_position[1]!='':
                        ref_segment_dict[ref_position[0]]='-'
                        temp_segment_dict[temp_position[0]]=temp_position[2]
                        align_segment_dict[ref_position[0]]='-'
                    elif (ref_position[2]=='-' or ref_position[2]=='_') and temp_position[1]=='':
                        ref_segment_dict[ref_position[0]]='-'
                        temp_segment_dict[str(temp_position[4])]=temp_position[2]
                        align_segment_dict[ref_position[0]]='-'
                    elif ref_position[2]=='_' and temp_position[1]!=False:
                        ref_segment_dict[ref_position[0]]='x'
                        temp_segment_dict[temp_position[0]]=temp_position[2]
                        align_segment_dict[ref_position[0]]='x'
    
                self.reference_dict[ref_seglab] = ref_segment_dict
                self.template_dict[ref_seglab] = temp_segment_dict
                self.alignment_dict[ref_seglab] = align_segment_dict

        for r_seglab, t_seglab, a_seglab in zip(self.reference_dict,self.template_dict,self.alignment_dict):
            if r_seglab in ['ICL1','ECL1','ICL2']:
                if len(list(self.reference_dict[r_seglab].keys()))==0:
                    well_aligned = False
                else:
                    well_aligned = True
                    for r, t, a in zip(self.reference_dict[r_seglab],self.template_dict[t_seglab],self.alignment_dict[a_seglab]):
                        if 'x' in r and 'x' in t and self.alignment_dict[a_seglab][a] in ['-','x']:
                            well_aligned = False
                        if self.reference_dict[r_seglab][r]=='-' and self.template_dict[t_seglab][t]!='-':
                            well_aligned = False
                    if well_aligned==False:
                        del self.reference_dict[r_seglab]
                        del self.template_dict[t_seglab]
                        del self.alignment_dict[a_seglab]
            elif r_seglab=='ECL2':
                del self.reference_dict[r_seglab]
                del self.template_dict[t_seglab]
                del self.alignment_dict[a_seglab]
        return self
