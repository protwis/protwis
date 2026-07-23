import os
import logging

from django.conf import settings

from structure.model_parsers.error_handling import log_or_raise

from Bio.PDB import PDBParser

class BaseModelMetrics():
    
    """Represents the metrics associated with a structure model"""

    def __init__(self, metrics_file_path, error_handling="log"):
        self.metrics_file_path = metrics_file_path
        self.metrics = self.from_csv()
        self.logger = logging.getLogger('build')
        self.error_handling = error_handling    

    def from_csv(self):
        if os.path.exists(self.metrics_file_path):
            try:
                with open(self.metrics_file_path, 'r') as f:
                    header = f.readline().strip().split(',')
                    values = f.readline().strip().split(',')
                    if len(header) != len(values):
                        raise ValueError(f"Header and values length mismatch in metrics file {self.metrics_file_path}.")
                    return dict(zip(header, values))
                    
            except Exception as e:
                if self.error_handling == "log":
                    self.logger.error(f"Error reading metrics file {self.metrics_file_path}: {e}")
                else:
                    raise Exception(f"Error reading metrics file {self.metrics_file_path}: {e}")
        else:
            if self.error_handling == "log":
                self.logger.error(f"Metrics file {self.metrics_file_path} not found.")
            raise FileNotFoundError(f"Metrics file {self.metrics_file_path} not found.")


class BaseModel():

    """Defines a base class for a structure model, containing its PDB structure and associated metrics."""

    residue_to_one_letter = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

    def __init__(self, data_dir, error_handling="log"):
        self.data_dir = data_dir
        self.model_name = os.path.basename(self.data_dir)
        self.pdb_file_path = None  # Initialize pdb_file_path to None; subclasses should set this if needed
        self.logger = logging.getLogger('build')
        self.error_handling = error_handling

    def unpack_model_name(self):
        raise NotImplementedError("Subclasses must implement the unpack_model_name method to unpack the model name into its components.")

    def fetch_pdb_structure(self, parser_config):
        s = PDBParser(PERMISSIVE=False, get_header=True, QUIET=True) \
                                        .get_structure(self.model_name, 
                                                       self.pdb_file_path)
        if parser_config.pdb_header_override:
            s.header.update(parser_config.pdb_header_override)
        return s
    
    def fetch_pdb_content(self):
        if os.path.exists(self.pdb_file_path):
            with open(self.pdb_file_path, 'r') as f:
                return f.read()
        else:
            log_or_raise(self.logger, f"PDB file {self.pdb_file_path} not found for model {self.model_name}.", FileNotFoundError, self.error_handling)
    
    def parse_metrics(self, parser_config):
        raise NotImplementedError("Subclasses must implement the parse_metrics method to parse metrics.")


class BaseModelParser():   
    
    """Base class for parsing models organized by model set name and model name
    """

    def __init__(self, config):
        """Initialize the parser with a configuration object."""
        self.config = config
        self.logger = logging.getLogger('build')
        self.error_handling = config.error_handling

    def process_models(self):
        raise NotImplementedError("Subclasses must implement the process_models method to process models based on the configuration.")

class BaseModelParserConfig():
    """Configuration class for BaseModelParser"""
    
    def __init__(self, model_set_name, data_dir=None, pdb_header_override=None, error_handling="log"):
        self.model_set_name = model_set_name
        self.data_dir = data_dir if data_dir else os.sep.join([settings.DATA_DIR, 'structure_data', model_set_name])
        self.pdb_header_override = pdb_header_override

        self.logger = logging.getLogger('build')
        self.error_handling = error_handling