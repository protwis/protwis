"""
A set of utility functions for alignment processing.
"""
import re
from collections import OrderedDict
from common import definitions


def get_numbering_schemes(proteins):
    """
    Extract all of the numbering schemes used for a set of proteins.

    Arguments:
        proteins {selection} -- A set of proteins to analyze
    """

    numbering_schemes = {}
    for prot in proteins:
        if prot.protein.residue_numbering_scheme.slug not in numbering_schemes:
            rnsn = prot.protein.residue_numbering_scheme.name
            numbering_schemes[prot.protein.residue_numbering_scheme.slug] = rnsn
    # order and convert numbering scheme dict to tuple
    return sorted(numbering_schemes.items(), key=lambda x: x[0])

def strip_html_tags(text):
    """
    Remove the html tags from a string.

    @param: text - string to clean up
    """
    return re.sub('<.*?>', '', text)

def get_format_props(freq=None, res=None):
    """
    Get the excel cell format for residues/frequencies.

    @param: freq - get the format for feature/aa frequency formatting
    @param: res - get the format for residue colouring (alignment)
    """
    residue = {
        'A': {
            'bg_color': '#E6E600',
            'font_color': '#000000',
        },
        'C': {
            'bg_color': '#B2B548',
            'font_color': '#000000',
        },
        'D': {
            'bg_color': '#E60A0A',
            'font_color': '#FDFF7B',
        },
        'E': {
            'bg_color': '#E60A0A',
            'font_color': '#FDFF7B',
        },
        'F': {
            'bg_color': '#18FF0B',
            'font_color': '#000000',
        },
        'G': {
            'bg_color': '#FF00F2',
            'font_color': '#000000',
        },
        'H': {
            'bg_color': '#0093DD',
            'font_color': '#000000',
        },
        'I': {
            'bg_color': '#E6E600',
            'font_color': '#000000',
        },
        'K': {
            'bg_color': '#145AFF',
            'font_color': '#FDFF7B',
        },
        'L': {
            'bg_color': '#E6E600',
            'font_color': '#000000',
        },
        'M': {
            'bg_color': '#E6E600',
            'font_color': '#000000',
        },
        'N': {
            'bg_color': '#A70CC6',
            'font_color': '#FDFF7B',
        },
        'P': {
            'bg_color': '#CC0099',
            'font_color': '#FDFF7B',
        },
        'Q': {
            'bg_color': '#A70CC6',
            'font_color': '#FDFF7B',
        },
        'R': {
            'bg_color': '#145AFF',
            'font_color': '#FDFF7B',
        },
        'S': {
            'bg_color': '#A70CC6',
            'font_color': '#FDFF7B',
        },
        'T': {
            'bg_color': '#A70CC6',
            'font_color': '#FDFF7B',
        },
        'V': {
            'bg_color': '#E6E600',
            'font_color': '#000000',
        },
        'W': {
            'bg_color': '#0BCF00',
            'font_color': '#000000',
        },
        'Y': {
            'bg_color': '#18FF0B',
            'font_color': '#000000',
        },
        '-': {
            'bg_color': '#FFFFFF',
            'font_color': '#000000',
        },
        '_': {
            'bg_color': '#EDEDED',
            'font_color': '#000000',
        },
    }
    properties = {
        0: {
            'bg_color': '#ff0000',
        },
        1: {
            'bg_color': '#ff3300',
        },
        2: {
            'bg_color': '#ff6600',
        },
        3: {
            'bg_color': '#ff9900',
        },
        4: {
            'bg_color': '#ffcc00',
        },
        5: {
            'bg_color': '#ffff00',
        },
        6: {
            'bg_color': '#ccff00',
        },
        7: {
            'bg_color': '#99ff00',
        },
        8: {
            'bg_color': '#66ff00',
        },
        9: {
            'bg_color': '#33ff00',
        },
        10: {
            'bg_color': '#00ff00',
        },
    }

    if freq:
        try:
            return properties[freq]
        except KeyError:
            return properties[int(freq)]
    elif res:
        return residue[res]

def prepare_excel_worksheet(workbook, worksheet_name, alignment, generic_numbers_set=None, data_block=None, data_type='features'):
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
        alignment {Alignment} -- an alignment to save the data from and to extract
                                    numbering schemes and segment names data from

    Keyword Arguments:
        generic_numbers_set {OrderedDict} -- a set of generic numbers to be used for saving.
                                Follows scheme->segment->position scheme used in Alignment
        data_block {list} -- a block of stats to be saved (features or amino acids)
                                (default: {None})
        data_type {string} -- selection of type of stats to save: "features" or "amino_acids" data
                                (default: {"features"})
    """

    if not generic_numbers_set:
        generic_numbers_set = alignment.generic_numbers

    if data_block:
        if data_type == 'features':
            props = definitions.AMINO_ACID_GROUP_NAMES.values()
            print(props)
        elif data_type == 'amino_acids':
            props = definitions.AMINO_ACIDS.keys()

    worksheet = workbook.add_worksheet(worksheet_name)

    # First column, numbering schemes
    for row, scheme in enumerate(alignment.numbering_schemes):
        worksheet.write(1 + row, 0, scheme[1])

    # First column, stats
    if data_block:
        for offset, prop in enumerate(props):
            worksheet.write(1 + len(alignment.numbering_schemes) + offset, 0, prop)

    # First column, protein list (for alignment) and line for consensus sequence
    else:
        for offset, prot in enumerate(alignment.proteins):
            worksheet.write(
                1 + len(alignment.numbering_schemes) + offset,
                0,
                prot.protein.entry_name
            )
        worksheet.write(
            1 + len(alignment.numbering_schemes) + len(alignment.proteins),
            0,
            'CONSENSUS'
            )

    # Second column and on
    # Segments
    offset = 0
    for segment in alignment.segments:
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
    if data_block:
        offset = 1 + 3 * len(alignment.numbering_schemes)

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
                cell_format = workbook.add_format(get_format_props(res[1]))
                worksheet.write(
                    row,
                    1 + col + col_offset,
                    res[0],
                    cell_format
                )
            col_offset += len(sequence.items())
