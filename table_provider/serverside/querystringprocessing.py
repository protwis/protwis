import re

class DatatablesQueryStringProcessor:
    def __init__(self, request):
        """Initializes a DatatablesQueryStringProcessor instance and begins processing."""
        self.query_set = self._process_query_string(request.query_params)

    def _process_query_string(self, query_dict):
        """Parse the keys in the request.query_params dictionary containing the key/value pairs from DataTables AJAX query string

        Parse the data structures encoded in the query string into true python data structures, handling both dictionary and array types.
            :param query_dict: The query string parameters from the request, typically a dictionary-like object.
            :return: A nested dictionary representing the structured query parameters.
        """
        query_set = {}
        for key, value in query_dict.items():
            sub_keys = re.findall(r'([^\[\]]+)', key)
            current_level = query_set
            for (i, sub_key) in enumerate(sub_keys):
                # If it's the last sub_key, assign the value
                if i == len(sub_keys) - 1:
                    # Flatten the value if it's a list with a single element
                    if isinstance(value,list) and len(value) == 1:
                        value = value[0]
                    current_level[sub_key] = value
                else:
                    if sub_key.isdigit():
                        sub_key = int(sub_key)
                        # Ensure the current level is a list and has enough elements
                        while len(current_level) <= sub_key:
                            current_level.append(None)
                        # If the current level at sub_key is None, initialize it as a list or dict based on the next sub_key
                        if current_level[sub_key] is None:
                            current_level[sub_key] = [] if sub_keys[i+1].isdigit() else {}
                    else:
                        # Assume it's a dictionary key and initialize if it doesn't exist
                        if sub_key not in current_level:
                            current_level[sub_key] = [] if sub_keys[i+1].isdigit() else {}
                    # Move to the next level in the hierarchy
                    current_level = current_level[sub_key]
        
        return query_set