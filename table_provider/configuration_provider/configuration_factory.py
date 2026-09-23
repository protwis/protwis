import json
import pathlib

class ConfigurationFactory():
    """Factory class for providing dataTables configurations to the frontend."""

    def __init__(self, table_name, configuration_variant):
        self.table_name = table_name
        self.configuration_variant = configuration_variant

class ColumnConfigurationFactory(ConfigurationFactory):
    """Factory class for providing column configurations to the frontend."""

    def fetch(self):
        # Implement logic to fetch column configuration based on self.table_name and self.configuration_variant
        try:
            config_path = pathlib.Path(__file__).resolve().parent / "column_configurations" / f"{self.table_name}__{self.configuration_variant}.json"
            with open(config_path, 'r') as file:
                config = json.load(file)
                return config
        except FileNotFoundError:
            raise FileNotFoundError(f"Column configuration file for {self.table_name} with variant {self.configuration_variant} not found.")

class DataTablesConfigurationFactory(ConfigurationFactory):
    """Factory class for providing dataTables configurations to the frontend."""
    
    def fetch(self):
        # Implement logic to fetch dataTables configuration based on self.table_name and self.configuration_variant
        try:
            config_path = pathlib.Path(__file__).resolve().parent / "datatables_configurations" / f"{self.table_name}__{self.configuration_variant}.json"
            with open(config_path, 'r') as file:
                config = json.load(file)
                return config
        except FileNotFoundError:
            raise FileNotFoundError(f"DataTables configuration file for {self.table_name} with variant {self.configuration_variant} not found.")