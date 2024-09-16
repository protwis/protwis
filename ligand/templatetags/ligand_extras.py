from django import template

register = template.Library()

@register.filter(name='remove_species_from_entry_name')
def remove_species_from_entry_name(value):
    """
        Remove species name a UniProt entry name.
    """
    return value.split('_')[0]
