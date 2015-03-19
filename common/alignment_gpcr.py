from common.alignment import Alignment as GenericAlignment

class Alignment(GenericAlignment):
    """A class representing a protein sequence alignment, with or without a reference sequence"""

    def format_generic_number(self, generic_number):
        split_gn = generic_number.split("x")
        split_gn_helix = split_gn[0].split('.')
        formatted_gn = '<b>{:s}</b><br />.{:s}<br />x{:s}'.format(split_gn_helix[0], split_gn_helix[1], split_gn[1])
        
        return formatted_gn