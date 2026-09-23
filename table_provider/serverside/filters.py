from table_provider.serverside.querystringprocessing import DatatablesQueryStringProcessor
from django.db.models import Q

"""
 This module provides classes and functions for filtering and ordering data in a server-side processing context, specifically for use with DataTables.
 It defines a FilterSet class that processes query parameters from the request to generate filter requirements and ordering instructions based on the
 specified serializer's field mappings and data types. The module also includes specific filter classes for text and numeric range filtering,
 allowing for flexible and dynamic filtering of data based on user input.
"""

class FilterSet:
    def __init__(self, request, SerializerClass):
        """Initializes a FilterSet instance."""
        self._filters = []
        self._ordering = []
        self.request = request
        self.SerializerClass = SerializerClass
        self.field_mappings = SerializerClass.source_mapping()
        self.query_parameters = DatatablesQueryStringProcessor(self.request).query_set
        self.parse_filters()
        self.parse_ordering()

    def __iter__(self):
        return iter(self._filters)

    def has_filters(self):
        return len(self._filters) > 0

    def add_filter(self, filter_requirement):
        self._filters.append(filter_requirement)

    def parse_filters(self):
        """Convert query parameters to filter objects

        Parses query parameters from a DataTables serverside request
        to generate filter requirements based on the specified
        serializer's field mappings and data types.
        """

        if 'columns' not in self.query_parameters:
            return

        for column in self.query_parameters['columns']:
            data_field = column['data']
            field_type = self.SerializerClass.Meta.datatypes.get(data_field)
            search_value = column['search']['value']

            if not field_type:
                raise ValueError(f"No datatype mapping found for field '{data_field}' in serializer. " + \
                                 "Ensure you have declared all fields that require filtering under datatypes in meta data")

            if search_value:
                if field_type == 'text':
                    filter_req = TextFilter(self.field_mappings.get(data_field), search_value)
                elif field_type == 'numeric':
                    [filter_min, filter_max ]  = search_value.split('~')
                    if filter_min != "-Inf" or filter_max != "Inf":
                        filter_req = NumericRangeFilter(self.field_mappings.get(data_field), filter_min, filter_max)
                    else:
                        continue
                else:
                    continue

                self.add_filter(filter_req)

    def parse_ordering(self):
        """Extract columns used for row ordering from query parameters"""

        if 'order' not in self.query_parameters:
            return

        order_array = []
        for order_column in self.query_parameters['order']:
            column_idx = int(order_column['column'])
            json_name = self.query_parameters['columns'][column_idx]['data']
            db_name = self.field_mappings.get(json_name)
            order_direction = '' if order_column['dir'] == 'asc' else '-'
            order_array.append(f"{order_direction}{db_name}")
        self._ordering = order_array


    def get_filters(self, data_field=None):
        return self._filters

    def get_ordering(self):
        return self._ordering

class FilterRequirement:
    def __init__(self, db_field_name, filter_type):
        self.filter_type = filter_type
        self.db_field_name = db_field_name
        self.filter_values = None

class TextFilter(FilterRequirement):
    """A text filter requirement for a specific database field.

    Supports filtering based on exact matches or regular expressions, allowing for flexible text-based filtering.
    Interacts with Select2 "tag" system to allow partial matches and multiple values to be filtered on a single field.
    Queries are built using Django's Q objects, allowing for complex filtering logic to be constructed and applied to the database query.
    """

    def __init__(self, db_field_name, filter_value):
        super().__init__(db_field_name, 'text')
        self.tag_filters = []
        self.standard_filters = []

        filter_value = filter_value.split('|')

        #remove regex anchors added for datatables from filter values
        for fv in filter_value:
            if fv.startswith('('):
                fv = fv.strip("()")
            if fv.startswith('^'):
                fv = fv.strip("^$")
            if ".*" in fv:
                self.tag_filters.append(fv)
            else:
                fv = fv.replace("\\", "") #Remove escape characters from standard filters (added for regex compatibility with datatables)
                self.standard_filters.append(fv)

    def format_query(self):
        """Converts filter values into a Django Q object for querying the database"""
        if len(self.tag_filters) == 0 and len(self.standard_filters) == 1:
            return Q(**{f"{self.db_field_name}": self.standard_filters[0]})

        query = Q()
        if len(self.standard_filters) > 1:
            query &= Q(**{f"{self.db_field_name}__in": self.standard_filters})

        for tag_filter in self.tag_filters:
            query |= Q(**{f"{self.db_field_name}__regex": rf"{tag_filter}"})

        return query



class NumericRangeFilter(FilterRequirement):
    """Represents a numeric range filter requirement for a specific database field.

    It supports filtering based on minimum and maximum values. Values can be specified as '-Inf' or 'Inf'
    to represent unbounded ranges, allowing for filtering based on only a minimum or maximum value.
    """

    def __init__(self, db_field_name, filter_min, filter_max):
        """Initializes a NumericRangeFilter instance with the specified database field name and filter values."""
        super().__init__(db_field_name, 'numeric')

        if (filter_min == '-Inf' and filter_max == 'Inf'):
            self.filter_values = { 'min': filter_min, 'max': filter_max }
            self.filter_action = 'between'

        if (filter_min != '-Inf' and filter_max != 'Inf'):
            self.filter_values = { 'min': filter_min, 'max': filter_max }
            self.filter_action = 'between'

        elif filter_min != '-Inf':
            self.filter_values = { 'min': filter_min }
            self.filter_action = 'gt'

        elif filter_max != 'Inf':
            self.filter_values = { 'max': filter_max }
            self.filter_action = 'lt'

    def format_query(self):
        """Converts filter values into a Django Q object for querying the database"""
        if(self.filter_action == 'between'):
            return (Q(**{f"{self.db_field_name}__gt": self.filter_values.get('min')}) & Q(**{f"{self.db_field_name}__lt": self.filter_values.get('max')}))
        else:
            comp = 'min' if self.filter_action == 'gt' else 'max'
            return Q(**{f"{self.db_field_name}__{self.filter_action}": self.filter_values[comp]})
