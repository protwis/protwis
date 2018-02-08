"""
A module for generating sequence signatures for the given two sets of proteins.
"""
from django.conf import settings

from alignment.functions import *
#from common.alignment import Alignment
Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')

from common.definitions import AMINO_ACID_GROUPS, AMINO_ACID_GROUP_NAMES


from collections import OrderedDict
from copy import deepcopy
import numpy as np
import re

class SequenceSignature:
    """
    A class handling the sequence signature.
    """

    def __init__(self):

        self.aln_pos = Alignment()
        self.aln_neg = Alignment()

        self.features_normalized_pos = OrderedDict()
        self.features_normalized_neg = OrderedDict()
        self.features_frequency_difference = OrderedDict()
        self.features_frequency_diff_display = []

        self.freq_cutoff = 30
        self.common_gn = OrderedDict()

    def setup_alignments_from_selection(self, positive_selection, negative_selection):
        """
        The function gathers necessary information from provided selections
        and runs the calculations of the sequence alignments independently for
        both protein sets. It also finds the common set of residue positions.

        Arguments:
            positive_selection {Selection} -- selection containing first group of proteins
            negative_selection {[type]} -- selction containing second group of proteins along with the user-selcted sequence segments for the alignment
        """

        self.aln_pos.load_proteins_from_selection(positive_selection)
        self.aln_neg.load_proteins_from_selection(negative_selection)

        # In case positive and negative sets come from different classes
        # unify the numbering schemes
        self.common_schemes = self.merge_numbering_schemes()
        self.aln_pos.numbering_schemes = self.common_schemes
        self.aln_neg.numbering_schemes = self.common_schemes
        # now load the segments and generic numbers
        self.aln_pos.load_segments_from_selection(negative_selection)
        self.aln_neg.load_segments_from_selection(negative_selection)

        self.aln_pos.build_alignment()
        self.aln_pos.calculate_statistics()
        self.aln_neg.build_alignment()
        self.aln_neg.calculate_statistics()

        self.common_gn = deepcopy(self.aln_pos.generic_numbers)
        for scheme in self.aln_neg.numbering_schemes:
            for segment in self.aln_neg.segments:
                for pos in self.aln_neg.generic_numbers[scheme[0]][segment].items():
                    if pos[0] not in self.common_gn[scheme[0]][segment].keys():
                        self.common_gn[scheme[0]][segment][pos[0]] = pos[1]
                self.common_gn[scheme[0]][segment] = OrderedDict(sorted(
                    self.common_gn[scheme[0]][segment].items(),
                    key=lambda x: x[0].split('x')
                    ))

    def calculate_signature(self):
        """
        Calculates the feature frequency difference between two protein sets.
        Generates the full differential matrix as well as maximum difference for a position (for scatter plot).
        """
        for sid, segment in enumerate(self.aln_neg.segments):
            self.features_normalized_pos[segment] = np.array(
                [[x[0] for x in feat[sid]] for feat in self.aln_pos.feature_stats],
                dtype='int'
                )
            self.features_normalized_neg[segment] = np.array(
                [[x[0] for x in feat[sid]] for feat in self.aln_neg.feature_stats],
                dtype='int'
                )

        for segment in self.aln_neg.segments:
            #TODO: get the correct default numering scheme from settings
            for idx, res in enumerate(self.common_gn['gpcrdba'][segment].keys()):
                if res not in self.aln_pos.generic_numbers['gpcrdba'][segment].keys():
                    self.features_normalized_pos[segment] = np.insert(self.features_normalized_pos[segment], idx, 0, axis=1)
                    # aa_pos_norm[segment] = np.insert(aa_pos_norm[segment], idx, 0, axis=1)
                elif res not in self.aln_neg.generic_numbers['gpcrdba'][segment].keys():
                    self.features_normalized_neg[segment] = np.insert(self.features_normalized_neg[segment], idx, 0, axis=1)
                    # aa_neg_norm[segment] = np.insert(aa_neg_norm[segment], idx, 0, axis=1)

            # now the difference
            self.features_frequency_difference[segment] = np.subtract(
                self.features_normalized_pos[segment],
                self.features_normalized_neg[segment]
                )
        # Version with display data
        for row, feat in enumerate(AMINO_ACID_GROUPS.keys()):
            tmp_row = []
            for segment in self.aln_neg.segments:
                #first item is the real value,
                # second is the assignmnent of color (via css)
                # 0 - red, 5 - yellow, 10 - green
                #third item is a tooltip
                tmp_row.append([[
                    x,
                    int(x/20)+5,
                    "{} - {}".format(
                        self.features_normalized_pos[segment][row][y],
                        self.features_normalized_neg[segment][row][y]
                        )
                    ] for y, x in enumerate(self.features_frequency_difference[segment][row])])
            self.features_frequency_diff_display.append(tmp_row)

    def prepare_display_data(self):

        options = {
            'num_residue_columns': len(sum([[x for x in self.common_gn['gpcrdba'][segment]] for segment in self.aln_neg.segments], [])),
            'num_of_sequences_pos': len(self.aln_pos.proteins),
            'num_residue_columns_pos': len(self.aln_pos.positions),
            'num_of_sequences_neg': len(self.aln_neg.proteins),
            'num_residue_columns_neg': len(self.aln_neg.positions),
            'common_segments': self.aln_neg.segments,
            'common_generic_numbers': self.common_gn,
            'feats_signature': self.features_frequency_diff_display,
            'a_pos': self.aln_pos,
            'a_neg': self.aln_neg,
        }
        return options

    def merge_numbering_schemes(self):
        """
        Extract all of the numbering schemes used for a set of proteins.

        Arguments:
            proteins {selection} -- A set of proteins to analyze
        """

        numbering_schemes = {}
        for prot in self.aln_pos.proteins + self.aln_neg.proteins:
            if prot.protein.residue_numbering_scheme.slug not in numbering_schemes:
                rnsn = prot.protein.residue_numbering_scheme.name
                numbering_schemes[prot.protein.residue_numbering_scheme.slug] = rnsn
        # order and convert numbering scheme dict to tuple
        return sorted(numbering_schemes.items(), key=lambda x: x[0])

    def prepare_excel_worksheet(self, workbook, worksheet_name, aln='positive', data='alignment'):
        """
        A function saving alignment data subset into the excel spreadsheet.
        It adds a worksheet to an existing workbook and saves only a selected subset of alignment data.
        For a complete save of the alignment it needs to be wrapped with additional code.

        The outline of the excel worksheet is similar to the one of html page.
        First column shows nunbering schemes, protein list, etc
        The frequency data start from column B

        Arguments:
            workbook {xlrsxwriter.Workbook} -- object to add workseet to
            worksheet_name {string} -- name for the new workseet

        Keyword Arguments:
            alignment {string} -- alignment to extract data from.
                                    Possible choices: positive, negative, signature
            data {string} -- data type to save to workshet: 'alignment' or 'features' frequencies
        """

        props = AMINO_ACID_GROUP_NAMES.values()
        worksheet = workbook.add_worksheet(worksheet_name)

        if aln == 'positive':
            numbering_schemes = self.aln_pos.numbering_schemes
            generic_numbers_set = self.aln_pos.generic_numbers
            alignment = self.aln_pos
            if data == 'features':
                data_block = self.aln_pos.feature_stats
        elif aln == 'negative':
            numbering_schemes = self.aln_neg.numbering_schemes
            generic_numbers_set = self.aln_neg.generic_numbers
            alignment = self.aln_neg
            if data == 'features':
                data_block = self.aln_neg.feature_stats
        else:
            numbering_schemes = self.common_schemes
            generic_numbers_set = self.common_gn
            if data == 'features':
                data_block = self.features_frequency_diff_display

        # First column, numbering schemes
        for row, scheme in enumerate(numbering_schemes):
            worksheet.write(1 + row, 0, scheme[1])

        # First column, stats
        if data == 'features':
            for offset, prop in enumerate(props):
                worksheet.write(3 + len(numbering_schemes) + offset, 0, prop)

        # First column, protein list (for alignment) and line for consensus sequence
        else:
            for offset, prot in enumerate(alignment.proteins):
                worksheet.write(
                    3 + len(numbering_schemes) + offset,
                    0,
                    prot.protein.entry_name
                )
            worksheet.write(
                3 + len(numbering_schemes) + len(alignment.proteins),
                0,
                'CONSENSUS'
                )

        # Second column and on
        # Segments
        offset = 0
        for segment in generic_numbers_set['gpcrdba'].keys():
            worksheet.merge_range(
                0,
                1 + offset,
                0,
                len(generic_numbers_set['gpcrdba'][segment]) + offset - 1,
                segment
            )
            offset += len(generic_numbers_set['gpcrdba'][segment])

        # Generic numbers
        offset = 1
        for row, item in enumerate(generic_numbers_set.items()):
            scheme = item[0]
            segment = item[1]
            for sn, gn_list in segment.items():
                for col, gn_pair in enumerate(gn_list.items()):
                    tm, bw, gpcrdb = re.split('\.|x', strip_html_tags(gn_pair[1]))
                    worksheet.write(
                        1 + 3 * row,
                        col + offset,
                        tm
                    )
                    worksheet.write(
                        2 + 3 * row,
                        col + offset,
                        bw
                    )
                    worksheet.write(
                        3 + 3*row,
                        col + offset,
                        gpcrdb
                    )
                offset += len(gn_list.items())

        # Stats
        if data == 'features':
            offset = 1 + 3 * len(numbering_schemes)

            for row, prop in enumerate(data_block):
                col_offset = 0
                for segment in prop:
                    for col, freq in enumerate(segment):
                        cell_format = workbook.add_format(get_format_props(freq[1]))
                        worksheet.write(
                            offset + row,
                            1 + col + col_offset,
                            freq[0] if isinstance(freq[0], int) else int(freq[0]),
                            cell_format
                        )
                    col_offset += len(segment)
        # Alignment
        else:
            offset = 1 + 3 * len(alignment.numbering_schemes)

            for row, data in enumerate(alignment.proteins):
                col_offset = 0
                for segment, sequence in data.alignment.items():
                    for col, res in enumerate(sequence):
                        cell_format = workbook.add_format(get_format_props(res=res[2]))
                        worksheet.write(
                            offset + row,
                            1 + col + col_offset,
                            res[2],
                            cell_format
                        )
                    col_offset += len(sequence)
            # Consensus sequence
            row = 1 + 3 * len(alignment.numbering_schemes) + len(alignment.proteins)
            col_offset = 0
            for segment, sequence in alignment.consensus.items():
                for col, data in enumerate(sequence.items()):
                    res = data[1]
                    cell_format = workbook.add_format(get_format_props(res=res[0]))
                    worksheet.write(
                        row,
                        1 + col + col_offset,
                        res[0],
                        cell_format
                    )
                    cell_format = workbook.add_format(get_format_props(res[1]))
                    worksheet.write(
                        row + 1,
                        1 + col + col_offset,
                        res[2],
                        cell_format
                    )

                col_offset += len(sequence.items())