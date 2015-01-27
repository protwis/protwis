class SimpleSelection:
    """A class representing the proteins and segments a user has selected. Can be serialized and stored in session"""
    def __init__(self):
        self.reference = []
        self.targets = []
        self.segments = []

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

    def exporter(self):
        """Exports the attributes of Selection to a SimpleSelection object, and returns it"""
        ss = SimpleSelection()
        ss.reference = self.reference
        ss.targets = self.targets
        ss.segments = self.segments
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

        # loop through selected objects and remove the one that matches the subtype and ID
        for selection_object in selection:
            if not (selection_object.type == selection_subtype and selection_object.item.id == int(selection_id)):
                updated_selection.append(selection_object)
        setattr(self, selection_type, updated_selection)

    def expand(self):
        pass

    def render(self, selection_type):
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