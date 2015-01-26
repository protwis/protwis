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
        self.reference = simple_selection.reference
        self.target = simple_selection.target
        self.segment = simple_selection.segment

    def exporter(self):
        ss = SimpleSelection()
        ss.reference = self.reference
        ss.target = self.target
        ss.segment = self.segment
        return ss

    def add(self, selection_type, selection_subtype, selection_object):
        selection = getattr(self, selection_type)
        if not selection_object in selection:
            selection.append(selection_object)
            setattr(self, selection_type, selection)

    def delete(self, selection_type, selection_subtype, selection_id):
        pass

    def expand(self):
        pass

    def render(self):
        return {
            'reference': self.reference,
            'target': self.target,
            'segment': self.segment,
        }


class SelectionItem:
    """Docstring"""
    def __init__(self, selection_type, selection_object):
        self.type = selection_type
        self.type_title = selection_type.title()
        self.item = selection_object

    def __str__(self):
        return str(self.__dict__)

    def __eq__(self, other): 
        return self.__dict__ == other.__dict__