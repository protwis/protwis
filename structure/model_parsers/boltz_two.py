
import re
import os
import csv
from django.db import transaction
import pandas as pd
from collections import OrderedDict
from django.conf import settings
from Bio.PDB import PDBParser, PDBIO, Polypeptide

import logging
from interaction.views import ligand
from structure.model_parsers.logging import ParserVerbosity, conditional_log
from structure.model_parsers.helpers import csv_to_dict
from structure.model_parsers.error_handling import log_or_raise

from structure.model_parsers.base import BaseModel, BaseModelMetrics, BaseModelParser, BaseModelParserConfig, ModelGProtienComplex, ModelLigand

from structure.models import Structure, PdbData, StructureType, StructureModelScores, StructureExtraProteins, StructureModelpLDDT
from protein.models import Protein, ProteinConformation, ProteinState, Residue

import json as JSON


from contactnetwork.cube import compute_interactions

class BoltzTwoComplexModelMetrics(BaseModelMetrics):
    
    """Represents the metrics associated with a BoltzTwoComplexModel"""

    def __init__(self, data_dir, model_name, metrics_file_prefix, model_version_number, error_handling="log", verbosity=ParserVerbosity.SILENT):
        """Build the metrics file path from the given directory, prefix, and version, then initialize the base metrics parser.

        Args:
            data_dir: Data directory of the model containing the metrics file.
            model_name: Name of the model.
            metrics_file_prefix: Prefix for the metrics file.
            model_version_number: Version number of the model (e.g. {prefix}_{version}.pdb).
            error_handling: Error handling strategy (e.g. log, raise, etc.).
            verbosity: Verbosity level for logging (ParserVerbosity.SILENT, ParserVerbosity.BASIC, ParserVerbosity.EVERYTHING).

        """
        metrics_file_path = os.sep.join([data_dir, metrics_file_prefix + '_' + model_version_number + '.csv'])
        super().__init__(metrics_file_path, error_handling=error_handling, verbosity=verbosity)

    def save(self, struct):
        """Persist the parsed metrics to the StructureModelScores record associated with the given structure, creating it if it does not already exist.

        Args:
            struct: A structure.models.Structure instance to which the metrics should be associated.
        """
        #Remove extra fields present in metrics file from metrics list before recording.
        self.metrics.pop("model_rank", None) 
        self.metrics.pop("original_model_num", None)
        self.metrics.pop("is_best_ligand_iptm", None)
        self.metrics.pop("pass_ligand_iptm_0.93", None)
        self.metrics.pop("pass_ligand_plddt_0.90", None)
        self.metrics.pop("pass_vicinity_0.70", None)

        try:
            metrics = StructureModelScores.objects.get(structure=struct)
            metrics.metrics_json = JSON.dumps(self.metrics) #Update existing metrics
        except StructureModelScores.DoesNotExist:
            metrics = StructureModelScores()
            metrics.structure = struct            
            metrics.metrics_json = JSON.dumps(self.metrics)
        metrics.save()

        return metrics


class BoltzTwoComplexModel(BaseModel):

    """Defines a BoltzTwoComplex model, containing its PDB structure and associated metrics."""

    def __init__(self, data_dir, parser_config):
        """Initialize the model with its source data directory and parser configuration.

        Args:
            data_dir: The data directory of the model.
            parser_config: An instance of a BoltzTwoComplexParserConfig object.
        """
        super().__init__(data_dir, parser_config)

    def load(self):
        """Populate the model's attributes (identifiers, ligand, PDB structure, metrics, and metadata) from the files in its data directory."""
        self.populate_identifiers_from_manifest()

        self.ligand.type = self.ligand.get_type()
        
        self.model_name = os.path.basename(self.data_dir)
        self.model_version = self.get_preferred_model_version_number()
        self.pdb_file_path = os.sep.join([self.data_dir, 
                                            self.parser_config.pdb_file_prefix + "_" + self.model_version + '.pdb'])
        
        if self.signprot:
            self.signprot_subunits = self.populate_signprot_subunits()        
        
        self.model_structure_type_name = 'Model (Boltz2)'
        self.model_structure_type_slug = f'b2{ "-signprot" if self.signprot else "" }{ "-" + self.ligand.type.replace("-", "") if self.ligand.type else "" }'

        self.pdb_raw = self.fetch_pdb_content()
        self.pdb_structure = self.fetch_pdb_structure()
        
        self.metrics = BoltzTwoComplexModelMetrics(self.data_dir, self.model_name, self.parser_config.metrics_file_prefix, self.model_version, error_handling=self.parser_config.error_handling, verbosity=self.parser_config.verbosity)
        self.metrics.load()

        self.model_date = self.parse_model_date_from_pdb_header()        

        if self.parser_config.model_receptor_state:
            self.receptor_state = self.parser_config.model_receptor_state
        else:
            log_or_raise(self.logger, "Model receptor state must be specified in the parser configuration.", ValueError, self.error_handling)

        if self.parser_config.pdb_preferred_chain:
                    self.pdb_preferred_chain = self.parser_config.pdb_preferred_chain
        else:
            log_or_raise(self.logger, "Model PDB preferred chain must be specified in the parser configuration.", ValueError, self.error_handling)    

    def populate_identifiers_from_manifest(self):
        """Read the identifiers.csv manifest file for this model and populate the receptor, ligand, and signalling protein identifiers from it."""
        try:
            manifest_file_path = os.sep.join([self.data_dir, 'identifiers.csv'])
            identifiers = csv_to_dict(manifest_file_path)

            self.receptor = identifiers.get("receptor", "") + "_" + identifiers.get("species", "")

            self.ligand = ModelLigand(error_handling=self.error_handling, verbosity=self.verbosity)                
            self.ligand.name = identifiers.get("ligand_name", "")            
            self.ligand.pubchemcid = identifiers.get("PubChemCID", "")
            self.ligand.smiles = identifiers.get("SMILES", "")
            self.ligand.inchikey = identifiers.get("InChIKey", "")
            self.ligand.pdb_chain_id = "C"  

            self.signprot = None
            has_g = identifiers.get("G_protein", "") == "with_miniG"
            if has_g:
                self.signprot = identifiers.get("mini_g_entry", None)

            self.identifiers = identifiers
            
        except Exception as e:
            log_or_raise(self.logger, f"Error reading identifiers manifest file {manifest_file_path}: {e}", Exception, self.error_handling, parent_exception=e)        

    def get_preferred_model_version_number(self):
        """Get the index of the model to use based on the parser configuration.

        If no override is specified, return the default model version number. If no default model is specified, return '1' as the default model version number.

        Returns:
            String - The preferred model version number.
        """
        if self.parser_config.default_model_version:
            if self.model_name in self.parser_config.default_model_version.get("override", {}):
                return self.parser_config.default_model_version["override"][self.model_name]       
            return self.parser_config.default_model_version.get("default_version_number", "1")
        return "1"  # Default to model version number '1' if no default model is specified
    
    def format_pdb_index(self):
        """Return the PDB index string for this model, built from the receptor, ligand, and (if present) signalling protein names.

        Returns:
            String - The PDB index string for the model.
        """
        if self.signprot:
            return f'B2M_{self.receptor.upper()}_{self.ligand.name.upper()}_{self.signprot.upper()}'
        else:
            return f'B2M_{self.receptor.upper()}_{self.ligand.name.upper()}'

    def write(self):
        """Write the loaded model, its metrics, and related records (protein, structure, ligands, extra proteins, pLDDT, contact network) to the database inside a single transaction."""
        conditional_log(self, f"Starting to write model {self.model_name} to database.", logging.INFO, ParserVerbosity.BASIC)
        try:
            with transaction.atomic():
        
                receptor_protein = self.protein_from_entry_name()

                protein_state = self.get_or_create_protein_state()

                protein_conformation = self.get_protein_conformation(receptor_protein)
                
                struct = self.get_or_initialise_structure(receptor_protein, protein_state, protein_conformation)

                signprot, signprot_conf, beta_protconf, gamma_protconf = self.get_signprot_and_conformations(struct)      

                ligands_db = self.ligand.fetch_db_entities()

                if any([ligand.ligand_type.name in ['peptide', 'protein'] for ligand in ligands_db]):
                    self.create_ligand_peptide_structure(struct, ligands_db)

                self.metrics.save(struct)

                self.create_extra_proteins(struct, signprot, signprot_conf, beta_protconf, gamma_protconf, alpha_note=self.identifiers.get("mini_gprot_construct",))

                self.store_plddt(struct, receptor_protein, signprot, beta_protconf, gamma_protconf)

                self.build_contact_network(struct, signprot)

            conditional_log(self, f"Successfully wrote model {self.model_name} to database.", logging.INFO, ParserVerbosity.BASIC)

        except Exception as e:
            log_or_raise(self.logger, f"Error writing model {self.model_name} to database: {str(e)}", Exception, self.error_handling, parent_exception=e)

class BoltzTwoComplexModelParserConfig(BaseModelParserConfig):

    """Configuration class for BoltzTwoComplexModelParser"""
    
    def __init__(self, model_set_name, data_dir=None, 
                 default_model_version={'default_version_number': '1', 'override': {}}, pdb_file_prefix="model", 
                 metrics_file_prefix="metrics", pdb_header_override=None, 
                 model_receptor_state=None, pdb_preferred_chain='A', 
                 error_handling="log", verbosity=ParserVerbosity.SILENT):
        """Initialize the parser configuration with file-naming, model-version-selection, and receptor/chain settings for BoltzTwoComplex models.

        Args:
            model_set_name: The name of the model set (also the directory name in data_dir).
            data_dir: The base data directory in which the directory named {model_set_name} is located (defaults to the structure_data directory via BaseModelParserConfig).
            default_model_version: A dictionary specifying the default model version to use (Key: default_version_number) and any specific edge cases to override (Key: override, format: { model_name : version_number }) .
            pdb_file_prefix: The prefix for the PDB files (e.g. {prefix}_{version_number}.pdb).
            metrics_file_prefix: The prefix for the metrics files (e.g. {prefix}_{version_number}.csv).
            pdb_header_override: A dictionary of fields and values to override in the PDB header.
            model_receptor_state: The state of the receptor for the model.
            pdb_preferred_chain: The preferred chain for the PDB file.
            error_handling: Error handling strategy (e.g. log, raise, etc.).
            verbosity: Verbosity level for logging (ParserVerbosity.SILENT, ParserVerbosity.BASIC, ParserVerbosity.EVERYTHING).
        """
        super().__init__(model_set_name, data_dir=data_dir, pdb_header_override=pdb_header_override, error_handling=error_handling, verbosity=verbosity)

        self.pdb_file_prefix = pdb_file_prefix
        self.model_receptor_state = model_receptor_state
        self.metrics_file_prefix = metrics_file_prefix
        self.pdb_preferred_chain = pdb_preferred_chain
        self.default_model_version = default_model_version        

class BoltzTwoComplexModelParser(BaseModelParser):   
    
    """Parses a directory of BoltzTwoComplex models organized by model-set-name and model-name
    
    The expected directory structure is as follows:
    /structure_data/{model_set_name}/{model_name}/model_{model_version_number}.pdb and metrics_{model_version_number}.csv
    where model_name is structured as {receptor}-{ligand}-{signprot(optional)}
    """

    def __init__(self, config):
        """Initialize the parser with a configuration object.

        Args:
            config: An instance of a BoltzTwoComplexModelParserConfig object.
        """
        super().__init__(config)
        self.config = config
        self.model_dirs = []
        self.models = []

    def get_model_directories(self):
        """Get a list of model directories in the model set directory."""
        try:
            self.model_dirs = list(filter(os.path.isdir, [os.path.join(self.config.data_dir, f) for f in os.listdir(self.config.data_dir)])) #/structure_data/{model_set_name}/[*]
            self.model_dirs = sorted(self.model_dirs, key=lambda x: os.path.basename(x))  # Sort the model directories by name in case OS returns them in a different order
            conditional_log(self, f"Found {len(self.model_dirs)} model directories in {self.config.data_dir}.", logging.INFO, ParserVerbosity.BASIC)
        except Exception as e:
            log_or_raise(self.logger, f"Error accessing model directories in {self.config.data_dir}: {e}", Exception, self.error_handling, parent_exception=e)

    def process_models(self, write = True, low_memory=True, offset_start=0, offset_end=None):
        """Process each model directory in the model set directory and return a list of BoltzTwoComplexModel instances.
        
        Args:
            write (bool): Whether to write the models to the database. Default is True.
            low_memory (bool): Flag indicating whether to discard models from memory after writing. Default is True. False means models will be stored in the models array (uses a lot of memory).
            offset_start (int): The starting index of the model directories to process. Default is 0.
            offset_end (int): The ending index of the model directories to process. Default is None, which means process all directories from offset_start to the end.
        """
        try:
            if low_memory and not write:
                raise ValueError("Cannot set low_memory to True when write is False. This would result in models being discarded from memory without being written to the database.")

            conditional_log(self, f"Processing models from {self.config.data_dir} with offset_start={offset_start} and offset_end={offset_end}.", logging.INFO, ParserVerbosity.BASIC)
            
            if len(self.model_dirs) == 0:
                self.get_model_directories()
            if len(self.model_dirs) == 0:
                log_or_raise(self.logger, f"No model directories found in {self.config.data_dir}.", FileNotFoundError, self.error_handling)            

            models_to_process = self.model_dirs[offset_start:offset_end] if offset_end else self.model_dirs[offset_start:]

            for model_path in models_to_process:
                model = BoltzTwoComplexModel(model_path, self.config)
                model.load()
                if not low_memory:
                    self.models.append(model)
                if write:
                    model.write()
        except Exception as e:
            log_or_raise(self.logger, f"Error processing models: {str(e)}", Exception, self.error_handling, parent_exception=e)