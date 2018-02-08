"""
A set of utility functions for alignment processing.
"""
import re
from collections import OrderedDict
from common import definitions


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
        '+': {
            'bg_color': '#FFFFFF',
            'font_color': '#000000',
        }
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
