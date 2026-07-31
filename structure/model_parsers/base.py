import json as JSON
import os
import logging
import re
from io import StringIO

from Bio.PDB import PDBParser, PDBIO, Polypeptide

from django.conf import settings

from structure.model_parsers.error_handling import log_or_raise
from structure.model_parsers.logging import ParserVerbosity, conditional_log

from protein.models import Protein, ProteinConformation, ProteinState, Residue
from structure.models import Structure, PdbData, StructureType, StructureModelScores, StructureExtraProteins, StructureModelpLDDT
from common.models import WebLink, WebResource
from ligand.models import Ligand, LigandPeptideStructure
from signprot.models import SignprotComplex

import structure.assign_generic_numbers_gpcr as generic_number_assigner

from common.definitions import G_PROTEIN_DISPLAY_NAME as g_prot_dict, ARRESTIN_DISPLAY_NAME as arr_dict
from contactnetwork.cube import compute_interactions

from structure.model_parsers.helpers import csv_to_dict

class BaseModelMetrics():
    
    """Represents the metrics associated with a structure model"""

    def __init__(self, metrics_file_path, error_handling="log", verbosity=ParserVerbosity.SILENT):
        self.metrics_file_path = metrics_file_path
        self.logger = logging.getLogger('build')
        self.error_handling = error_handling
        self.verbosity = verbosity
        self.metrics = {}

    def load(self):
        if os.path.exists(self.metrics_file_path):
            try:
               self.metrics = csv_to_dict(self.metrics_file_path)
               conditional_log(self, f"Parsed metrics from {self.metrics_file_path}: {JSON.dumps(self.metrics)}", logging.INFO, ParserVerbosity.EVERYTHING)
            except Exception as e:
                log_or_raise(self.logger, f"Error reading metrics file {self.metrics_file_path}: {e}", Exception, self.error_handling, parent_exception=e)
        else:
            log_or_raise(self.logger, f"Metrics file {self.metrics_file_path} not found.", FileNotFoundError, self.error_handling)

class BaseModel():

    """Defines a base class for a structure model, containing its PDB structure and associated metrics."""

    def __init__(self, data_dir, parser_config):
        self.parser_config = parser_config
        self.logger = logging.getLogger('build')
        self.error_handling = parser_config.error_handling
        self.verbosity = parser_config.verbosity

        conditional_log(self, f"Initializing processing of model in {data_dir}.", logging.INFO, ParserVerbosity.BASIC)

        self.data_dir = data_dir
        self.model_name = os.path.basename(self.data_dir)
        self.pdb_file_path = None  # Initialize pdb_file_path to None; subclasses should set this if needed

    def fetch_pdb_structure(self):
        conditional_log(self, f"Reading PDB content as PDBStructure object for model {self.model_name}.", logging.INFO, ParserVerbosity.EVERYTHING)
        s = PDBParser(PERMISSIVE=False, get_header=True, QUIET=True) \
                                        .get_structure(self.model_name, 
                                                       self.pdb_file_path)
        if self.parser_config.pdb_header_override:
            s.header.update(self.parser_config.pdb_header_override)
        return s
    
    def fetch_pdb_content(self):
        conditional_log(self, f"Reading PDB content as text for model {self.model_name}.", logging.INFO, ParserVerbosity.EVERYTHING)
        if os.path.exists(self.pdb_file_path):
            with open(self.pdb_file_path, 'r') as f:
                return f.read()
        else:
            log_or_raise(self.logger, f"PDB file {self.pdb_file_path} not found for model {self.model_name}.", FileNotFoundError, self.error_handling)
    
    def has_beta_gamma_complex(self):
        """Check if the model has a beta-gamma complex based on the signprot attribute and subunit list."""
        if self.signprot:
            if self.signprot_subunits.beta and self.signprot_subunits.gamma:
                return True
            else:                
                return False
        else:             
            return False

    def protein_from_entry_name(self):
        try:
            db_protein = Protein.objects.get(entry_name=self.receptor.lower())
            conditional_log(self, f"Fetched receptor object for model {self.model_name} from database with entry name {db_protein.entry_name} and id {db_protein.id}.", logging.INFO, ParserVerbosity.EVERYTHING)
            return db_protein
        except Protein.DoesNotExist:
            log_or_raise(self.logger, f"Protein object not found for entry name: {self.receptor}", ValueError, self.error_handling)

    def populate_signprot_subunits(self):
        """Populate the signprot attribute into its subunits (alpha, beta, gamma) if applicable."""
        return_subunits = ModelGProtienComplex()
        signprot_re = re.compile(r'([^_]+_[^_]+)')
        subunit_uniprot = signprot_re.findall(self.signprot)
        if subunit_uniprot:
            db_subunits = Protein.objects.filter(entry_name__in=[subunit.lower() for subunit in subunit_uniprot])
            for su in db_subunits:
                try:
                    if su.family.parent.parent.name == "Alpha":
                        return_subunits.alpha = su
                    elif su.family.parent.name == "Beta":
                        return_subunits.beta = su
                    elif su.family.parent.name == "Gamma":
                        return_subunits.gamma = su
                except Exception as e:
                    log_or_raise(self.logger, f"Error determining subunit type for {su.entry_name}: {e}", Exception, self.error_handling, parent_exception=e)
        return return_subunits

    def parse_model_date_from_pdb_header(self):
        """Extract the model date from the PDB header."""
        if self.pdb_structure.header and 'deposition_date' in self.pdb_structure.header:
            if self.pdb_structure.header['deposition_date'] != '1909-01-08':
                date = self.pdb_structure.header['deposition_date']
                conditional_log(self, f"Model date found in deposition date ({date}) field for model {self.model_name}.", logging.INFO, ParserVerbosity.EVERYTHING)
                return date
        
        if self.pdb_structure.header and 'release_date' in self.pdb_structure.header:
            if self.pdb_structure.header['release_date'] != '1909-01-08':
                date = self.pdb_structure.header['release_date']
                conditional_log(self, f"Model date found in release date ({date}) field for model {self.model_name}.", logging.INFO, ParserVerbosity.EVERYTHING)
                return date

        if self.pdb_structure.header and 'head' in self.pdb_structure.header:
            match = re.match(r'.+\s+(\d{4}-\d{2}-\d{2})', self.pdb_structure.header['head'])
            if match:
                date = match.group(1)
                conditional_log(self, f"Model date found in head ({date}) field for model {self.model_name}.", logging.INFO, ParserVerbosity.EVERYTHING)
                return date

        log_or_raise(self.logger, f"Could not parse model date from PDB header for model {self.model_name}." + 
                     "Amend the PDB header or provide a pdb_header_override in the parser configuration.", ValueError, self.error_handling)

    def get_or_initialise_structure(self, receptor_protein, protein_state, protein_conformation):
            struct = None
            try:
                struct = Structure.objects.get(protein_conformation__protein=receptor_protein, pdb_code__index=self.format_pdb_index(), structure_type__slug=self.model_structure_type_slug)
            except Structure.DoesNotExist:
                conditional_log(self, f"Structure for model {self.model_name} does not exist in the database. Creating a new Structure object.", logging.INFO, ParserVerbosity.BASIC)

            if not struct:
                try:
                    struct = Structure()
        
                    struct.structure_type = self.get_or_create_structure_type()
        
                    pdb_with_generic_numbers = self.assign_generic_numbers_to_pdb()
        
                    struct.representative = False
                    struct.state = protein_state
                    struct.author_state = protein_state
                    struct.resolution = None
                    struct.publication_date = self.model_date
                    struct.annotated = True
                    struct.refined = False
                    struct.stats_text = None
                    struct.preferred_chain = self.pdb_preferred_chain
        
                    struct.pdb_data = self.write_pdb_with_generic_numbers(pdb_with_generic_numbers)
                    struct.pdb_code = self.write_pdb_code_weblink(struct)
        
                    struct.protein_conformation = protein_conformation
        
                    struct.save()
                except Exception as e:
                    log_or_raise(self.logger, f"Failed to create Structure object for model {self.model_name}: {e}", Exception, self.error_handling, parent_exception=e)
        
                try:
                    struct.protein_conformation.generate_sites()
                except Exception as e:
                    ERROR_HANDLING_OVERRIDE = "log" #Override until we fix generate_sites.
                    log_or_raise(self.logger, f"Failed to generate sites for structure {self.model_name}: {e}", Exception, ERROR_HANDLING_OVERRIDE, parent_exception=e)
    
            return struct

    def get_or_create_protein_state(self):
        try:
            ps, created = ProteinState.objects.get_or_create(slug=self.receptor_state.lower(), defaults={'name': self.receptor_state})
            if created:
                conditional_log(self, f"Created protein state {self.receptor_state} with slug {self.receptor_state.lower()}.", logging.INFO, ParserVerbosity.BASIC)
        except Exception as e:
            log_or_raise(self.logger, f"Failed to get or create protein state {self.receptor_state}: {e}", Exception, self.error_handling, parent_exception=e)

        return ps

    def get_protein_conformation(self, receptor_protein):
        try:
            return ProteinConformation.objects.get(protein=receptor_protein)
        except ProteinConformation.DoesNotExist:
            log_or_raise(self.logger, f"Protein conformation for construct {receptor_protein.entry_name} does not exist", ValueError, self.error_handling)

    def assign_generic_numbers_to_pdb(self):
                
        try:
            pdb_struct = StringIO(self.pdb_raw)
            header = pdb_struct.readline()
            pdb_struct.seek(0)  # Reset the pointer to the beginning of the StringIO object
            assign_gn = generic_number_assigner.GenericNumbering(pdb_file=pdb_struct, blastdb=os.sep.join([settings.STATICFILES_DIRS[0], 'blast', 'protwis_human_blastdb']), sequence_parser=True)
            pdb_struct = assign_gn.assign_generic_numbers_with_sequence_parser()

            io = PDBIO()
            io.set_structure(pdb_struct)

            pdb_buffer = StringIO()
            io.save(pdb_buffer)
            return header + pdb_buffer.getvalue()
        except Exception as e:
            log_or_raise(self.logger, f"GN assignment failed for {self.model_name}.", Exception, self.error_handling, parent_exception=e)

    def write_pdb_with_generic_numbers(self, pdb_with_generic_numbers):
        try:
            pdbdata, created = PdbData.objects.get_or_create(pdb=pdb_with_generic_numbers)
            if created:
                conditional_log(self, f"Created PdbData object for model {self.model_name}.", logging.INFO, ParserVerbosity.EVERYTHING)
            return pdbdata
        except Exception as e:
            log_or_raise(self.logger, f"Failed to create PdbData object for {self.model_name}: {e}", Exception, self.error_handling, parent_exception=e )

    def write_pdb_code_weblink(self, struct):
        try:
            web_resource = WebResource.objects.get(slug='pdb')
            pdb_code = self.format_pdb_index()
            pdb_code_weblink, created = WebLink.objects.get_or_create(index=pdb_code, web_resource=web_resource)
            if created:
                conditional_log(self, f"Created WebLink object for model {self.model_name} with PDB code {pdb_code}.", logging.INFO, ParserVerbosity.EVERYTHING)                
            return pdb_code_weblink
        except Exception as e:
            log_or_raise(self.logger, f"Failed to create WebLink object for PDB code {pdb_code}: {e}", Exception, self.error_handling, parent_exception=e)

    def get_or_create_structure_type(self):
        try:
            structure_type, created = StructureType.objects.get_or_create(slug=self.model_structure_type_slug, defaults={'name': self.model_structure_type_name})
            if created:
                conditional_log(self, f"Created structure type {self.model_structure_type_name} with slug {self.model_structure_type_slug}.", logging.INFO, ParserVerbosity.BASIC)
            return structure_type
        except Exception as e:
            log_or_raise(self.logger, f"Failed to get or create structure type {self.model_structure_type_slug}: {e}", Exception, self.error_handling, parent_exception=e)

    def create_ligand_peptide_structure(self, struct, ligands_db):        
        if not self.ligand.pdb_chain_id:
            log_or_raise(self.logger, f"Ligand chain ID is not defined for in model {self.model_name}.", ValueError, self.error_handling)
        
        for ligand_db in ligands_db:
            if ligand_db.ligand_type.name not in ['peptide', 'protein']:
                continue  # Skip non-peptide ligands

            try:
                ligand_peptide_structure, created = LigandPeptideStructure.objects.get_or_create(
                    structure=struct,
                    ligand=ligand_db,
                    chain=self.ligand.pdb_chain_id,
                    defaults={'model': None} 
                )
                if created:
                    conditional_log(self, f"Created LigandPeptideStructure for ligand {ligand_db.name} in model {self.model_name}.", logging.INFO, ParserVerbosity.EVERYTHING)                
            except Exception as e:
                log_or_raise(self.logger, f"Error creating LigandPeptideStructure(s) for ligand {ligand_db.name} in model {self.model_name}: {str(e)} ", Exception, self.error_handling, parent_exception=e)

    def get_signprot_and_conformations(self, struct):
        signprot = None
        signprot_conf = None
        beta_protconf = None
        gamma_protconf = None

        if self.signprot:
            if self.signprot_subunits.alpha:
                signprot = self.signprot_subunits.alpha
            else:
                try:
                    signprot = Protein.objects.get(entry_name=self.signprot)
                except Protein.DoesNotExist as e:
                    log_or_raise(self.logger, f"Failed to fetch SignProt protein object for {self.model_name}: {str(e)}", ValueError, self.error_handling, parent_exception=e)

            try:
                signprot_conf = ProteinConformation.objects.get(protein=signprot)
            except ProteinConformation.DoesNotExist as e:
                log_or_raise(self.logger, f"Failed to fetch SignProt protein conformation for {self.model_name}: {str(e)}", ValueError, self.error_handling, parent_exception=e)

            conditional_log(self, f"Fetched Signprot protein and conformation for model {self.model_name} with entry name {signprot.entry_name} and state {signprot_conf.state}.", logging.INFO, ParserVerbosity.EVERYTHING)

            if self.has_beta_gamma_complex():
                beta_protconf = ProteinConformation.objects.get(protein__entry_name=self.signprot_subunits.beta.entry_name)
                gamma_protconf = ProteinConformation.objects.get(protein__entry_name=self.signprot_subunits.gamma.entry_name)
                sc, created = SignprotComplex.objects.get_or_create(alpha='B', protein=signprot, structure=struct,
                                                        beta_chain='C', gamma_chain='D', beta_protein=beta_protconf.protein, gamma_protein=gamma_protconf.protein)
                if created:
                    conditional_log(self, f"Created SignprotComplex for model {self.model_name} with signprot {signprot.entry_name}, beta {beta_protconf.protein.entry_name}, and gamma {gamma_protconf.protein.entry_name}.", logging.INFO, ParserVerbosity.EVERYTHING)
            else:
                sc, created = SignprotComplex.objects.get_or_create(alpha='B', protein=signprot, structure=struct,
                                                        beta_chain=None, gamma_chain=None, beta_protein=None, gamma_protein=None)
                if created:
                    conditional_log(self, f"Created SignprotComplex for model {self.model_name} with signprot {signprot.entry_name}.", logging.INFO, ParserVerbosity.EVERYTHING)

            struct.signprot_complex = sc
            struct.save()         
            conditional_log(self, f"Updated SignprotComplex on structure record for model {self.model_name} with signprot {signprot.entry_name}.", logging.INFO, ParserVerbosity.EVERYTHING)
        else:            
            signprot = None

        return signprot, signprot_conf, beta_protconf, gamma_protconf

    def create_extra_proteins(self, struct_db, signprot_db, signprot_conf_db, beta_protconf_db, gamma_protconf_db, alpha_note=None, beta_note=None, gamma_note=None):
        sep = None
        sep_beta = None
        sep_gamma = None        
        if signprot_db:        
            try:
                display_name = g_prot_dict[signprot_db.entry_name.split('_')[0].upper()]
                cat = 'G alpha'
            except:
                display_name = arr_dict[signprot_db.entry_name]
                cat = 'Arrestin'

            sep, created = StructureExtraProteins.objects.get_or_create(display_name=display_name, note=alpha_note, chain='B', category=cat, wt_coverage=100, protein_conformation=signprot_conf_db, structure=struct_db, wt_protein=signprot_db)
            if created:
                conditional_log(self, f"Created StructureExtraProteins for signprot {signprot_db.entry_name} in model {self.model_name}.", logging.INFO, ParserVerbosity.EVERYTHING)
            if self.has_beta_gamma_complex():
                sep_beta, created_beta = StructureExtraProteins.objects.get_or_create(display_name='G&beta;1', note=beta_note, chain='C', category='G beta', wt_coverage=100, protein_conformation=beta_protconf_db, structure=struct_db, wt_protein=beta_protconf_db.protein)
                sep_gamma, created_gamma = StructureExtraProteins.objects.get_or_create(display_name='G&gamma;2', note=gamma_note, chain='D', category='G gamma', wt_coverage=100, protein_conformation=gamma_protconf_db, structure=struct_db, wt_protein=gamma_protconf_db.protein)
                if created_beta:
                    conditional_log(self, f"Created StructureExtraProteins for G-beta {beta_protconf_db.protein.entry_name} in model {self.model_name}.", logging.INFO, ParserVerbosity.EVERYTHING)
                if created_gamma:
                    conditional_log(self, f"Created StructureExtraProteins for G-gamma {gamma_protconf_db.protein.entry_name} in model {self.model_name}.", logging.INFO, ParserVerbosity.EVERYTHING)

        return sep, sep_beta, sep_gamma

    def store_plddt(self, struct, receptor_protein, signprot, beta_protconf, gamma_protconf):
        #Adding plDDT for rendering
        resis = []
        for chain in self.pdb_structure.get_chains():
            for res in chain.get_residues():
                if not 'C' in res:
                    continue
                plddt = res['C'].get_bfactor()
                try:
                    if chain.get_id()=='A':
                        if not receptor_protein:                            
                            continue
                        res_obj = Residue.objects.get(protein_conformation__protein=receptor_protein, sequence_number=res.get_id()[1])
                    elif chain.get_id()=='B':
                        if not signprot:
                            continue
                        res_obj = Residue.objects.get(protein_conformation__protein=signprot, sequence_number=res.get_id()[1])
                    elif chain.get_id()=='C':
                        if not beta_protconf:
                            continue
                        res_obj = Residue.objects.get(protein_conformation__protein=beta_protconf.protein, sequence_number=res.get_id()[1])
                    elif chain.get_id()=='D':
                        if not gamma_protconf:
                            continue
                        res_obj = Residue.objects.get(protein_conformation__protein=gamma_protconf.protein, sequence_number=res.get_id()[1])
                    r = StructureModelpLDDT(structure=struct, residue=res_obj, pLDDT=plddt)
                    resis.append(r)
                except Residue.DoesNotExist:
                    continue
        try:
            StructureModelpLDDT.objects.bulk_create(resis)
            conditional_log(self, f"Stored pLDDT values for model {self.model_name}.", logging.INFO, ParserVerbosity.EVERYTHING)
        except Exception as e:
            log_or_raise(self.logger, f"Error storing pLDDT values for model {self.model_name}: {str(e)}", Exception, self.error_handling, parent_exception=e)

    def build_contact_network(self, receptor, signprot):
        if signprot:
            do_complexes = True
        else:
            do_complexes = False
        compute_interactions(self.pdb_file_path, protein=receptor, signprot=signprot, do_complexes=do_complexes, save_to_db=True, file_input=True) # add do_complexes


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

class ModelGProtienComplex():
    def __init__(self):
        self.alpha = None
        self.beta = None
        self.gamma = None

class ModelLigand():
    def __init__(self, error_handling="log", verbosity=ParserVerbosity.SILENT):
        self.name = None
        self.type = None
        self.pubchemcid = None
        self.smiles = None
        self.inchikey = None
        self.sequence_standard_aa_only = None
        self.sequence = None
        self.hashed_sequence = None
        self.pdb_chain_id = None

        self.logger = logging.getLogger('build')
        self.error_handling = error_handling
        self.verbosity = verbosity

    def get_sequence_from_pdb(self, pdb_structure):
        """Extract the ligand sequence from the PDB structure based on the ligand's chain ID."""
        if not self.pdb_chain_id:
            log_or_raise(self.logger, f"Ligand chain ID is not defined.", ValueError, self.error_handling)
        
        sequence = ""
        for model in pdb_structure:
            for chain in model:
                if chain.id == self.pdb_chain_id:
                    for residue in chain:
                        resname = residue.get_resname()
                        try:
                            one_letter = Polypeptide.protein_letters_3to1[resname]
                        except: 
                            one_letter = 'X'
                        sequence += one_letter
        return sequence

    def fetch_db_entities(self):
        ligands = None

        if self.inchikey:
            try:
                ligands = Ligand.objects.filter(inchikey=self.inchikey)
                conditional_log(self, f"Fetched ligands with InChIKey {self.inchikey} from database. Returned {len(ligands)}.", logging.INFO, ParserVerbosity.EVERYTHING)
                if ligands:
                    return ligands
            except Exception as e:
                log_or_raise(self.logger, f"Error fetching ligand from database by InChIKey: {e}", Exception, self.error_handling, parent_exception=e)

        if self.smiles:
            try:
                ligands = Ligand.objects.filter(smiles=self.smiles)
                conditional_log(self, f"Fetched ligands with SMILES {self.smiles} from database. Returned {len(ligands)}.", logging.INFO, ParserVerbosity.EVERYTHING)
                if ligands:
                    return ligands
            except Exception as e:
                log_or_raise(self.logger, f"Error fetching ligand from database by SMILES: {e}", Exception, self.error_handling, parent_exception=e)                    

        if self.sequence:
            try:
                ligands = Ligand.objects.filter(sequence=self.sequence)
                conditional_log(self, f"Fetched ligands with sequence {self.sequence} from database. Returned {len(ligands)}.", logging.INFO, ParserVerbosity.EVERYTHING)
                if ligands:
                    return ligands
            except Exception as e:
                log_or_raise(self.logger, f"Error fetching ligand from database by sequence: {e}", Exception, self.error_handling, parent_exception=e)

        if self.name:
            try:
                ligands = Ligand.objects.filter(name=self.name)
                conditional_log(self, f"Fetched ligands by name ({self.name}) from database. Returned {len(ligands)}.", logging.INFO, ParserVerbosity.EVERYTHING)
                if ligands:                    
                    return ligands
            except Exception as e:
                log_or_raise(self.logger, f"Error fetching ligand from database by name: {e}", Exception, self.error_handling, parent_exception=e)

        if not ligands:
            debug_info = f"Ligand search parameters: name={self.name}, inchikey={self.inchikey}, smiles={self.smiles}, sequence={self.sequence}"
            log_or_raise(self.logger, f"Ligand not found in database. Please ensure the ligand is present in the database or provide a valid ligand identifier or sequence. debug: {debug_info}", ValueError, self.error_handling)
            return None

    def get_type(self):
        ligands_db = self.fetch_db_entities()

        if not ligands_db:
            return None
        
        try:
            return ligands_db.first().ligand_type.name
        except Exception as e:
            log_or_raise(self.logger, f"Unable to determine ligand type.", ValueError, self.error_handling, parent_exception=e)