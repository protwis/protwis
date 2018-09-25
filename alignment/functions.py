"""
A set of utility functions for alignment processing.
"""
import re
from collections import OrderedDict
from common import definitions
from protein.models import Protein


def strip_html_tags(text):
    """
    Remove the html tags from a string.

    @param: text - string to clean up
    """
    return re.sub('<.*?>', '', text)

def get_format_props(freq=None, freq_gs=None, res=None, feat=None):
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
        '+': {
            'bg_color': '#FFFFFF',
            'font_color': '#000000',
        }
    }
    properties = {
        -1: {
            'bg_color': '#FFFFFF',
        },
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
    properties_gs = {
        0: {
            'bg_color': '#ffffff',
        },
        1: {
            'bg_color': '#e0e0e0',
        },
        2: {
            'bg_color': '#d0d0d0',
        },
        3: {
            'bg_color': '#c0c0c0',
        },
        4: {
            'bg_color': '#b0b0b0',
            'font_color': '#ffffff',
        },
        5: {
            'bg_color': '#a0a0a0',
            'font_color': '#ffffff',
        },
        6: {
            'bg_color': '#909090',
            'font_color': '#ffffff',
        },
        7: {
            'bg_color': '#808080',
            'font_color': '#ffffff',
        },
        8: {
            'bg_color': '#707070',
            'font_color': '#ffffff',
        },
        9: {
            'bg_color': '#606060',
            'font_color': '#ffffff',
        },
        10: {
            'bg_color': '#505050',
            'font_color': '#ffffff',
        },
    }
    property_group = {
        'HY': {
            'bg_color': '#93d050'
        },
        'HA': {
            'bg_color': '#ffff00',
        },
        'M': {
            'bg_color': '#ffff00',
        },
        'A': {
            'bg_color': '#ffff00',
        },
        'I': {
            'bg_color': '#ffff00',
        },
        'L': {
            'bg_color': '#ffff00',
        },
        'V': {
            'bg_color': '#ffff00',
        },
        'HR': {
            'bg_color': '#07b050',
        },
        'W': {
            'bg_color': '#07b050',
        },
        'Y': {
            'bg_color': '#07b050',
        },
        'F': {
            'bg_color': '#07b050',
        },
        'Hb': {
            'bg_color': '#7030a0',
            'font_color': '#ffffff',
        },
        'N': {
            'bg_color': '#7030a0',
            'font_color': '#ffffff',
        },
        'Q': {
            'bg_color': '#7030a0',
            'font_color': '#ffffff',
        },
        'S': {
            'bg_color': '#7030a0',
            'font_color': '#ffffff',
        },
        'T': {
            'bg_color': '#7030a0',
            'font_color': '#ffffff',
        },
        'Hu': {
            'bg_color': '#7030a0',
            'font_color': '#ffffff',
        },
        'Ha': {
            'bg_color': '#7030a0',
            'font_color': '#ff0000',
        },
        'Hd': {
            'bg_color': '#7030a0',
            # 'font_color': '#0070c0',
            'font_color': '#02b0f0',
        },
        '+-': {
            'bg_color': '#0070c0',
            'font_color': '#ff0000',
        },
        '+': {
            'bg_color': '#0070c0',
            'font_color': '#000000',
        },
        'H': {
            'bg_color': '#0070c0',
            'font_color': '#000000',
        },
        'K': {
            'bg_color': '#0070c0',
            'font_color': '#000000',
        },
        'R': {
            'bg_color': '#0070c0',
            'font_color': '#000000',
        },
        '-': {
            'bg_color': '#ff0000',
        },
        'D': {
            'bg_color': '#ff0000',
        },
        'E': {
            'bg_color': '#ff0000',
        },
        'Sm': {
            'bg_color': '#ffffff',
        },
        'aH': {
            'bg_color': '#d9d9d9',
        },
        'G': {
            'bg_color': '#ff02ff',
        },
        'P': {
            'bg_color': '#d603ff',
            'font_color': '#ffffff',
        },
        'C': {
            'bg_color': '#bf8f00',
        },
    }

    if freq is not None:
        try:
            return properties[freq]
        except KeyError:
            return properties[int(freq)]
    elif freq_gs is not None:
        try:
            return properties_gs[freq_gs]
        except KeyError:
            return properties_gs[int(freq_gs)]
    elif res is not None:
        return residue[res]
    elif feat is not None:
        print(feat)
        try:
            print(property_group[feat])
            return property_group[feat]
        except KeyError as msg:
            return {'bg_color': '#ffffff'}

def get_proteins_from_selection(simple_selection):

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

    return proteins

def prepare_aa_group_preference():

    pref_dict = {}
    lengths = {}
    for row, group in enumerate(definitions.AMINO_ACID_GROUPS.items()):
        tmp_len = len(group[1])
        try:
            lengths[tmp_len].append(row)
        except KeyError:
            lengths[tmp_len] = [row,]
    l_heap = sorted(lengths.keys())
    while l_heap:
        tmp = l_heap.pop()
        for feat_row in lengths[tmp]:
            pref_dict[feat_row] = []
            for pref_feat in l_heap:
                pref_dict[feat_row].extend(lengths[pref_feat])
    return pref_dict