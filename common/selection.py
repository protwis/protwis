from protein.models import Species
from protein.models import ProteinSource
from residue.models import ResidueNumberingScheme


class SimpleSelection:
    """A class representing the proteins and segments a user has selected. Can be serialized and stored in session"""
    def __init__(self):
        self.reference = []
        self.targets = []
        self.segments = []
        self.tree_settings = ['0','0','0','0'] # Default values for phylogenetic tree creation
        # species
        sp = Species.objects.get(common_name='Human') # Default species selection is human only
        o = SelectionItem('species', sp)
        self.species = [o]

        # annotation
        ps = ProteinSource.objects.get(name='SWISSPROT') # Default protein source is SWISSPROT
        o = SelectionItem('protein_source', ps)
        self.annotation = [o]

        # numbering schemes
        gn = ResidueNumberingScheme.objects.get(slug='gpcrdb')
        o = SelectionItem('numbering_schemes', gn)
        self.numbering_schemes = [o]

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
        self.annotation = simple_selection.annotation
        self.numbering_schemes = simple_selection.numbering_schemes
        self.tree_settings=simple_selection.tree_settings

    def exporter(self):
        """Exports the attributes of Selection to a SimpleSelection object, and returns it"""
        ss = SimpleSelection()
        ss.reference = self.reference
        ss.targets = self.targets
        ss.segments = self.segments
        ss.species = self.species
        ss.annotation = self.annotation
        ss.numbering_schemes = self.numbering_schemes
        ss.tree_settings=self.tree_settings

        return ss

    def add(self, selection_type, selection_subtype, selection_object):
        """Adds a selection item (protein, family, sequence segment etc.) to the selection"""
        # only one reference can be selected at a time, clear the selection
        if selection_type == 'reference':
            selection = []
        # for other types, add the selected item
        else:
            selection = getattr(self, selection_type)

        # check whether the selected item is already in the selection
        if not selection_object in selection:
            # if not, add it
            selection.append(selection_object)
            setattr(self, selection_type, selection)

    def remove(self, selection_type, selection_subtype, selection_id):
        """Removes a selection item (protein, family, sequence segment etc.) from the selection"""
        selection = getattr(self, selection_type)
        updated_selection = []
        deleted = False

        # loop through selected objects and remove the one that matches the subtype and ID
        for selection_object in selection:
            if not (selection_object.type == selection_subtype and selection_object.item.id == int(selection_id)):
                updated_selection.append(selection_object)
            else:
                deleted = True
        setattr(self, selection_type, updated_selection)

        return deleted

    def clear(self, selection_type):
        """Clears a section of the selection (e.g. targets)"""
        setattr(self, selection_type, [])

    def dict(self, selection_type):
        """Returns the selected attribute of Selection in a dictionary for rendering in templates"""
        return {
            'selection': {
                selection_type: getattr(self, selection_type),
            },
            'selection_type': selection_type,
        }


class SelectionItem:
    """A wrapper class for selectable objects (protein, family, sequence segment etc.) that adds a type attribute"""
    def __init__(self, selection_type, selection_object):
        self.type = selection_type
        self.type_title = selection_type.title()
        self.item = selection_object

    def __str__(self):
        return str(self.__dict__)

    def __eq__(self, other): 
        return self.__dict__ == other.__dict__