from django.conf import settings

from protein.models import Species
from protein.models import ProteinSource
from residue.models import ResidueNumberingScheme


class SimpleSelection:
    """A class representing the proteins and segments a user has selected. Can be serialized and stored in session"""
    def __init__(self):
        self.reference = []
        self.targets = []
        self.segments = []

        # species
        self.species = []

        # G proteins
        self.pref_g_proteins = []
        self.g_proteins = []

        # annotation
        ps = ProteinSource.objects.get(name='SWISSPROT') # Default protein source is SWISSPROT
        o = SelectionItem('protein_source', ps)
        self.annotation = [o]

        # numbering schemes
        gn = ResidueNumberingScheme.objects.get(slug=settings.DEFAULT_NUMBERING_SCHEME)
        o = SelectionItem('numbering_schemes', gn)
        self.numbering_schemes = [o]

        # Default values for phylogenetic tree creation
        self.tree_settings = ['0','0','0','0']

        # site residue groups (only used in site search)
        self.site_residue_groups = []
        self.active_site_residue_group = False

    def __str__(self):
        return str(self.__dict__)


class Selection(SimpleSelection):
    """A class that extends SimpleSelection, and adds methods to process the selection (these methods can not be
        serialized"""
    def importer(self, simple_selection):
        """Imports a SimpleSelection object into Selection"""
        self.reference = simple_selection.reference
        self.targets = simple_selection.targets
        self.segments = simple_selection.segments
        self.species = simple_selection.species
        self.pref_g_proteins = simple_selection.pref_g_proteins
        self.g_proteins = simple_selection.g_proteins
        self.annotation = simple_selection.annotation
        self.numbering_schemes = simple_selection.numbering_schemes
        self.tree_settings=simple_selection.tree_settings
        self.site_residue_groups = simple_selection.site_residue_groups
        self.active_site_residue_group = simple_selection.active_site_residue_group

    def exporter(self):
        """Exports the attributes of Selection to a SimpleSelection object, and returns it"""
        ss = SimpleSelection()
        ss.reference = self.reference
        ss.targets = self.targets
        ss.segments = self.segments
        ss.species = self.species
        ss.pref_g_proteins = self.pref_g_proteins
        ss.g_proteins = self.g_proteins
        ss.annotation = self.annotation
        ss.numbering_schemes = self.numbering_schemes
        ss.tree_settings=self.tree_settings
        ss.site_residue_groups = self.site_residue_groups
        ss.active_site_residue_group = self.active_site_residue_group

        return ss

    def add(self, selection_type, selection_subtype, selection_object):
        """Adds a selection item (protein, family, sequence segment etc.) to the selection"""
        # only one reference can be selected at a time, clear the selection
        if selection_type == 'reference':
            selection = []
        # for other types, add the selected item
        else:
            selection = getattr(self, selection_type)

        # make sure there is an active group
        if selection_subtype == 'site_residue':
            if not self.active_site_residue_group:
                self.active_site_residue_group = 1
        
            # update selection object with group id
            selection_object.properties['site_residue_group'] = self.active_site_residue_group
        
        # check whether the selected item is already in the selection
        if selection_object not in selection: # if not, add it
            # site residue groups
            if selection_subtype == 'site_residue':
                if not self.site_residue_groups:
                    self.site_residue_groups = [[]]
                self.site_residue_groups[self.active_site_residue_group - 1].append(1)
            
            selection.append(selection_object)
            setattr(self, selection_type, selection)

    def remove(self, selection_type, selection_subtype, selection_id):
        """Removes a selection item (protein, family, sequence segment etc.) from the selection"""
        selection = getattr(self, selection_type)
        updated_selection = []
        deleted = False

        # see if this item is part of a group
        group_id = False
        delete_group = False
        for selection_object in selection:
            if (selection_object.type == selection_subtype and selection_object.item.id == int(selection_id) and 
                'site_residue_group' in selection_object.properties and
                selection_object.properties['site_residue_group']):
                group_id = selection_object.properties['site_residue_group']
                delete_group = True

        # loop through selected objects and remove the one that matches the subtype and ID
        for selection_object in selection:
            if not (selection_object.type == selection_subtype and selection_object.item.id == int(selection_id)):
                updated_selection.append(selection_object)
                
                # check group ID
                if (group_id and 'site_residue_group' in selection_object.properties and 
                    selection_object.properties['site_residue_group'] == group_id):
                    delete_group = False
            else:
                deleted = True


        # if the deleted items group is not found anywhere else, delete the group
        if delete_group:
            del self.site_residue_groups[group_id-1]
            if self.site_residue_groups:
                self.active_site_residue_group = 1
            else:
                self.active_site_residue_group = False
        elif group_id:
            self.site_residue_groups[group_id-1].pop()
            if self.site_residue_groups[group_id - 1][0] > len(self.site_residue_groups[group_id - 1]):
                self.site_residue_groups[group_id - 1][0] -= 1

        setattr(self, selection_type, updated_selection)

        return deleted

    def clear(self, selection_type):
        """Clears a section of the selection (e.g. targets)"""
        setattr(self, selection_type, [])

        # when clearing segments, also clear residue groups
        if selection_type == 'segments':
            self.site_residue_groups = []
            self.active_site_residue_group = False

    def dict(self, selection_type):
        """Returns the selected attribute of Selection in a dictionary for rendering in templates"""
        return {
            'selection': {
                selection_type: getattr(self, selection_type),
                'site_residue_groups': getattr(self, 'site_residue_groups'),
                'active_site_residue_group': getattr(self, 'active_site_residue_group'),
            },
            'selection_type': selection_type,
        }


class SelectionItem:
    """A wrapper class for selectable objects (protein, family, sequence segment etc.) that adds a type attribute"""
    def __init__(self, selection_type, selection_object, properties={}):
        self.type = selection_type
        self.type_title = selection_type.replace('_', ' ').capitalize()
        self.item = selection_object
        self.properties = properties

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return str(self.__dict__)

    def __eq__(self, other): 
        return self.__dict__ == other.__dict__