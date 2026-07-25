import json as JSON
import os
import logging

from django.conf import settings

from structure.model_parsers.error_handling import log_or_raise
from structure.model_parsers.logging import ParserVerbosity, conditional_log

from Bio.PDB import PDBParser

class BaseModelMetrics():
    
    """Represents the metrics associated with a structure model"""

    def __init__(self, metrics_file_path, error_handling="log", verbosity=ParserVerbosity.SILENT):
        self.metrics_file_path = metrics_file_path
        self.logger = logging.getLogger('build')
        self.error_handling = error_handling
        self.verbosity = verbosity
        self.metrics = self.from_csv()

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
                log_or_raise(self.logger, f"Error reading metrics file {self.metrics_file_path}: {e}", Exception, self.error_handling, parent_exception=e)
        else:
            log_or_raise(self.logger, f"Metrics file {self.metrics_file_path} not found.", FileNotFoundError, self.error_handling)

class BaseModel():

    """Defines a base class for a structure model, containing its PDB structure and associated metrics."""

    def __init__(self, data_dir, error_handling="log", verbosity=ParserVerbosity.SILENT):
        self.logger = logging.getLogger('build')
        self.error_handling = error_handling
        self.verbosity = verbosity

        conditional_log(self, f"Initializing processing of model in {data_dir}.", logging.INFO, ParserVerbosity.BASIC)

        self.data_dir = data_dir
        self.model_name = os.path.basename(self.data_dir)
        self.pdb_file_path = None  # Initialize pdb_file_path to None; subclasses should set this if needed

    def unpack_model_name(self):
        raise NotImplementedError("Subclasses must implement the unpack_model_name method to unpack the model name into its components.")

    def fetch_pdb_structure(self, parser_config):
        conditional_log(self, f"Reading PDB content as PDBStructure object for model {self.model_name}.", logging.INFO, ParserVerbosity.EVERYTHING)
        s = PDBParser(PERMISSIVE=False, get_header=True, QUIET=True) \
                                        .get_structure(self.model_name, 
                                                       self.pdb_file_path)
        if parser_config.pdb_header_override:
            s.header.update(parser_config.pdb_header_override)
        return s
    
    def fetch_pdb_content(self):
        conditional_log(self, f"Reading PDB content as text for model {self.model_name}.", logging.INFO, ParserVerbosity.EVERYTHING)
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
        self.verbosity = config.verbosity

    def process_models(self):
        raise NotImplementedError("Subclasses must implement the process_models method to process models based on the configuration.")

class BaseModelParserConfig():
    """Configuration class for BaseModelParser"""
    
    def __init__(self, model_set_name, data_dir=None, pdb_header_override=None, error_handling="log", verbosity=ParserVerbosity.SILENT):
        self.model_set_name = model_set_name
        self.data_dir = data_dir if data_dir else os.sep.join([settings.DATA_DIR, 'structure_data', model_set_name])
        self.pdb_header_override = pdb_header_override

        self.logger = logging.getLogger('build')
        self.error_handling = error_handling
        self.verbosity = verbosity

