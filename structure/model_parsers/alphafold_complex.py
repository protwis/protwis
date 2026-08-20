import re
import os
import sys
from django.db import transaction
import pandas as pd
from django.conf import settings

import logging
from structure.model_parsers.logging import ParserVerbosity, conditional_log

from structure.model_parsers.error_handling import log_or_raise

from structure.model_parsers.base import BaseModel, BaseModelMetrics, BaseModelParser, BaseModelParserConfig, ModelLigand, LigandMultiMatchHandling

from structure.models import StructureModelScores

import json as JSON

class AlphaFoldTwoComplexModelMetrics(BaseModelMetrics):

    """Represents the metrics associated with a AlphaFoldTwoComplexModel"""

    def __init__(self, data_dir, model_name, error_handling="log", verbosity=ParserVerbosity.SILENT):
        """Locate the metrics file for the given model (trying both known naming conventions) and initialize the base metrics parser.

        Args:
            data_dir: Data directory of the model containing the metrics file.
            model_name: Name of the model.
            error_handling: Error handling strategy (e.g. log, raise, etc.).
            verbosity: Verbosity level for logging (ParserVerbosity.SILENT, ParserVerbosity.BASIC, ParserVerbosity.EVERYTHING).
        """
        path_type_1 = os.sep.join([data_dir, model_name + '_metrics.csv'])
        path_type_2 = os.sep.join([data_dir, model_name + '.csv'])

        metrics_file_path = None
        if os.path.exists(path_type_1):
            metrics_file_path = path_type_1
        elif os.path.exists(path_type_2):
            metrics_file_path = path_type_2
        else:
            self.logger = logging.getLogger('build')
            log_or_raise(self.logger, f"Metrics file not found for model {model_name}.", FileNotFoundError, error_handling)

        super().__init__(metrics_file_path, error_handling=error_handling, verbosity=verbosity)


    def save(self, struct):
        """Persist the parsed metrics to the StructureModelScores record associated with the given structure, creating it if it does not already exist.

        Args:
            struct: A structure.models.Structure instance to which the metrics should be associated.
        """
        self.metrics.pop("complex", None) #Remove complex name present in metrics file from metrics list before recording.
        try:
            metrics = StructureModelScores.objects.get(structure=struct)
            metrics.metrics_json = JSON.dumps(self.metrics) #Update existing metrics
        except StructureModelScores.DoesNotExist:
            metrics = StructureModelScores()
            metrics.structure = struct
            metrics.metrics_json = JSON.dumps(self.metrics)
        metrics.save()

        return metrics


class AlphaFoldTwoComplexModel(BaseModel):

    """Defines an object representing and AlphaFold2 model of a GPCR-ligand complex to contain its protein members, PDB structure, and associated metrics during parsing."""

    def __init__(self, data_dir, parser_config):
        """Initialize the model with its source data directory and parser configuration.

        Args:
            data_dir: The directory containing data for the model.
            parser_config: An instance of an AlphaFoldTwoComplexModelParserConfig object.
        """
        super().__init__(data_dir, parser_config)

    def load(self):
        """Populate the model's attributes (name, receptor/ligand/signprot, PDB structure, metrics, ligand sequence, and metadata) from the files in its data directory."""
        self.structure_origin = 'model'
        self.model_name = os.path.basename(self.data_dir)
        self.pdb_file_path = os.sep.join([self.data_dir,
                                            self.model_name + '.pdb'])

        self.receptor, self.ligand, self.signprot = self.unpack_model_name()
        if self.signprot:
            self.signprot_subunits = self.populate_signprot_subunits()

        self.model_structure_type_name = 'Model (AF2)'
        self.model_structure_type_slug = 'af-signprot-peptide' if self.signprot else 'af-peptide'

        self.pdb_raw = self.fetch_pdb_content()
        self.pdb_structure = self.fetch_pdb_structure()

        if not self.signprot: #ligand is on chain B in the PDB file for models without a signprot, but on chain E for models with a signprot. Remap chain B to E for consistency.
            self.pdb_raw = self.remap_chain_ids(self.pdb_raw, mapping_dict={"B":"E"})
            self.pdb_structure = self.remap_chain_ids(self.pdb_structure, mapping_dict={"B":"E"})

        self.ligand.pdb_chain_id = 'E'

        self.metrics = AlphaFoldTwoComplexModelMetrics(self.data_dir, self.model_name, error_handling=self.parser_config.error_handling, verbosity=self.parser_config.verbosity)
        self.metrics.load()

        self.annotate_ligand_sequence()
        self.model_date = self.parse_model_date_from_pdb_header()

        if self.parser_config.model_receptor_state:
            self.receptor_state = self.parser_config.model_receptor_state
        else:
            log_or_raise(self.logger, "Model receptor state must be specified in the parser configuration.", ValueError, self.error_handling)

        if self.parser_config.pdb_preferred_chain:
                    self.pdb_preferred_chain = self.parser_config.pdb_preferred_chain
        else:
            log_or_raise(self.logger, "Model PDB preferred chain must be specified in the parser configuration.", ValueError, self.error_handling)

    def annotate_ligand_sequence(self):
        """Assign the ligand's PDB chain id and populate both its original sequence and the 'cleaned' sequence with only standard amino acids"""
        self.ligand.sequence_standard_aa_only = self.ligand.get_sequence_from_pdb(self.pdb_structure)
        self.ligand.sequence = self.parser_config.old_seqs_dict.get(self.ligand.hashed_sequence, self.ligand.sequence_standard_aa_only)
        conditional_log(self, f"Annotated ligand sequences for model {self.model_name}. Original sequence: {self.ligand.sequence}. PDB sequence: {self.ligand.sequence_standard_aa_only}", logging.INFO, ParserVerbosity.EVERYTHING)

    def unpack_model_name(self):
        """Unpack the model name into receptor, ligand, and signprot components.

        Returns:
            Tuple - (receptor, ligand, signprot) parsed from the model directory name. signprot is None when not present.
        """
        parts = self.model_name.split('-')
        receptor = parts[0] if len(parts) > 0 else None
        ligand_temp = parts[1] if len(parts) > 1 else None
        signprot = parts[2] if len(parts) > 2 else None

        if ligand_temp is None:
            ligand = None
        else:
            ligand_component = re.match(r'hashedseq\[(.+)\]', ligand_temp)
            ligand = None
            if ligand_component:
                ligand = ModelLigand(self.error_handling, self.verbosity, self.parser_config.ligand_multimatch_handling)
                ligand.hashed_sequence = ligand_component.group(1)

        if not receptor or not ligand:
            log_or_raise(self.logger, f"Invalid model name format: {self.model_name}. Expected format: receptor-hashedseq[ligand_hash]-signprot (signprot is optional).", ValueError, self.error_handling)

        return receptor, ligand, signprot

    def format_pdb_index(self):
        """Return the PDB index string for this model, built from the receptor, ligand, and (if present) signalling protein names.

        Returns:
            String - The PDB index string for the model.
        """
        if 'peptide' in self.model_structure_type_slug:
            if self.signprot:
                return f'AFM_{self.receptor.upper()}_hashedseq[{self.ligand.hashed_sequence.upper()}]_{self.signprot.upper()}'
            else:
                return f'AFM_{self.receptor.upper()}_hashedseq[{self.ligand.hashed_sequence.upper()}]'
        else:
            return 'AFM_' + self.receptor.upper() + '_' + self.signprot.upper()

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

                self.create_ligand_peptide_structure(struct, ligands_db)

                self.metrics.save(struct)

                self.create_extra_proteins(struct, signprot, signprot_conf, beta_protconf, gamma_protconf)

                self.store_plddt(struct, receptor_protein, signprot, beta_protconf, gamma_protconf)

                self.build_contact_network(struct, signprot)

            conditional_log(self, f"Successfully wrote model {self.model_name} to database.", logging.INFO, ParserVerbosity.BASIC)

        except Exception as e:
            log_or_raise(self.logger, f"Error writing model {self.model_name} to database: {str(e)}", Exception, self.error_handling, parent_exception=e)

class AlphaFoldTwoComplexModelParserConfig(BaseModelParserConfig):

    """Configuration class for AlphaFoldTwoComplexModelParser"""

    def __init__(self, model_set_name, data_dir=None, cleaned_seq_csv=None, pdb_header_override=None, model_receptor_state=None, pdb_preferred_chain='A', error_handling="log", verbosity=ParserVerbosity.SILENT, ligand_multimatch_handling=LigandMultiMatchHandling.KEEP_FIRST):
        """Initialize the configuration for the AlphaFoldTwoComplex model parser.

        Args:
            model_set_name: The name of the model set (also the directory name in data_dir).
            data_dir: The base data directory in which the directory named {model_set_name} is located (defaults to the structure_data directory via BaseModelParserConfig).
            cleaned_seq_csv: Path to a CSV file mapping cleaned (hashed) sequences back to their original sequences.
            pdb_header_override: A dictionary of fields and values to override in the PDB header.
            model_receptor_state: The state of the receptor for the model.
            pdb_preferred_chain: The preferred chain for the PDB file.
            error_handling: Error handling strategy (e.g. log, raise, etc.).
            verbosity: Verbosity level for logging (ParserVerbosity.SILENT, ParserVerbosity.BASIC, ParserVerbosity.EVERYTHING).
            ligand_multimatch_handling: Strategy for handling multimatch cases for ligands (e.g. keep first, keep all, etc.).
        """
        super().__init__(model_set_name, data_dir=data_dir, pdb_header_override=pdb_header_override, 
                         error_handling=error_handling, verbosity=verbosity, 
                         ligand_multimatch_handling=ligand_multimatch_handling)

        self.cleaned_seq_csv = cleaned_seq_csv
        self.model_receptor_state = model_receptor_state
        self.pdb_preferred_chain = pdb_preferred_chain
        if cleaned_seq_csv and os.path.exists(cleaned_seq_csv):
            self.old_seqs_dict = self.generate_original_seq_lookup()
        else:
            log_or_raise(self.logger, f"Cleaned sequence CSV file not found at {cleaned_seq_csv}.", FileNotFoundError, self.error_handling)

    def generate_original_seq_lookup(self):
        """Build a lookup mapping the hashed sequence identifier found in PDB filenames to the corresponding original (pre-cleaning) sequence.

        Returns:
            Dict - Mapping of {pdb_file_hashseq: old_sequence} built from the cleaned sequence CSV.
        """
        df_cleaned_seqs = pd.read_csv(self.cleaned_seq_csv)
        df_cleaned_seqs['backwards_hex_cleaned_seq_hash_col'] = df_cleaned_seqs['cleaned_seq_hash_col'].apply(lambda x: hex(x)[2:][::-1])
        df_cleaned_seqs['pdb_file_hashseq'] = df_cleaned_seqs['cleaned_seq_hash'] + df_cleaned_seqs['backwards_hex_cleaned_seq_hash_col']
        return {k:v for k,v in zip(df_cleaned_seqs['pdb_file_hashseq'],df_cleaned_seqs['old_sequence'])}

class AlphaFoldTwoComplexModelParser(BaseModelParser):

    """Parses a directory of AlphaFoldTwoComplex models organized by model set name and model name

    The expected directory structure is as follows:
        /structure_data/{model_set_name}/{model_name}/
        where model_name is structured as {receptor}-{ligand}-{signprot(optional)}
    """

    def __init__(self, config):
        """Initialize the parser with a configuration object."""
        super().__init__(config)
        self.config = config
        self.model_dirs = []
        self.models = []

    def process_models(self, write = True, low_memory=True, offset_start=0, offset_end=None):
        """Process each model directory in the model set directory and return a list of AlphaFoldTwoComplexModel instances.

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
                model = AlphaFoldTwoComplexModel(model_path, self.config)
                model.load()
                if not low_memory:
                    self.models.append(model)
                if write:
                    model.write()
        except Exception as e:
            log_or_raise(self.logger, f"Error processing models: {str(e)}", Exception, self.error_handling, parent_exception=e)