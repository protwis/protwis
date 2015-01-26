class SimpleSelection:
    """A class representing the proteins and segments a user has selected. Can be serialized and stored in session"""
    def __init__(self):
        self.reference = []
        self.target = []
        self.segment = []

    def __str__(self):
        return str(self.__dict__)


class Selection(SimpleSelection):
    """A class that extends SimpleSelection, and adds methods to process the selection (these methods can not be
        serialized"""
    def importer(self, simple_selection):
        """Imports a SimpleSelection object into Selection"""
        self.reference = simple_selection.reference
        self.target = simple_selection.target
        self.segment = simple_selection.segment

    def exporter(self):
        """Exports the attributes of Selection to a SimpleSelection object, and returns it"""
        ss = SimpleSelection()
        ss.reference = self.reference
        ss.target = self.target
        ss.segment = self.segment
        return ss

    def add(self, selection_type, selection_subtype, selection_object):
        """Adds a selection item (protein, family, sequence segment etc.) to the selection"""
        selection = getattr(self, selection_type)

        # check whether the selected item is already in the selection
        if not selection_object in selection:
            # if not, add it
            selection.append(selection_object)
            setattr(self, selection_type, selection)

    def delete(self, selection_type, selection_subtype, selection_id):
        pass

    def expand(self):
        pass

    def render(self):
        """Returns the attributes of Selection in a dictionary for rendering in templates"""
        return {
            'reference': self.reference,
            'target': self.target,
            'segment': self.segment,
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