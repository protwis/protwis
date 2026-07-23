
from io import StringIO
import re
import os
import csv
from django.db import transaction
import pandas as pd
from collections import OrderedDict
from django.conf import settings
from Bio.PDB import PDBParser, PDBIO

from structure.model_parsers.error_handling import log_or_raise

from structure.model_parsers.base import BaseModel, BaseModelMetrics, BaseModelParser, BaseModelParserConfig

from common.models import WebLink, WebResource
from ligand.models import Ligand, LigandPeptideStructure
import structure.assign_generic_numbers_gpcr as generic_number_assigner
from structure.models import Structure, PdbData, StructureType, StructureAFScores, StructureExtraProteins, StructureModelpLDDT
from protein.models import Protein, ProteinConformation, ProteinState, Residue
from signprot.models import SignprotComplex

from common.definitions import G_PROTEIN_DISPLAY_NAME as g_prot_dict, ARRESTIN_DISPLAY_NAME as arr_dict
from contactnetwork.cube import compute_interactions

class AlphaFoldTwoComplexModelMetrics(BaseModelMetrics):
    
    """Represents the metrics associated with a AlphaFoldTwoComplexModel"""

    def save(self, struct):
        try:
            metrics = StructureAFScores.objects.get(structure=struct)
        except StructureAFScores.DoesNotExist:
            metrics = StructureAFScores()
            metrics.structure = struct
            metrics.ptm = self.metrics['ptm']
            metrics.iptm = self.metrics['iptm']
            metrics.pae_mean = self.metrics['pae_mean']
            metrics.save()


class AlphaFoldTwoComplexModel(BaseModel):

    """Defines a AlphaFoldTwoComplex model, containing its PDB structure and associated metrics."""

    def __init__(self, data_dir, parser_config):
        super().__init__(data_dir, error_handling=parser_config.error_handling)

        self.model_name = os.path.basename(self.data_dir)
        self.pdb_file_path = os.sep.join([self.data_dir, 
                                          self.model_name + '.pdb'])
        
        self.receptor, self.ligand, self.signprot = self.unpack_model_name()

        self.model_structure_type_name = 'Model (AF2)'
        self.model_structure_type_slug = 'af-signprot-peptide' if self.signprot else 'af-peptide'

        self.pdb_raw = self.fetch_pdb_content()
        self.pdb_structure = self.fetch_pdb_structure(parser_config)
        
        self.metrics = self.parse_metrics(parser_config)
        
        self.annotate_ligand_sequence(parser_config)
        self.model_date = self.parse_model_date_from_pdb_header()
        self.beta_gamma = self.has_beta_gamma_complex()

        if parser_config.model_receptor_state:
            self.receptor_state = parser_config.model_receptor_state
        else:
            log_or_raise(self.logger, "Model receptor state must be specified in the parser configuration.", ValueError, self.error_handling)

        if parser_config.pdb_preferred_chain:
                    self.pdb_preferred_chain = parser_config.pdb_preferred_chain
        else:
            log_or_raise(self.logger, "Model PDB preferred chain must be specified in the parser configuration.", ValueError, self.error_handling)


    def has_beta_gamma_complex(self):
        """Determine if the model has a beta-gamma complex based on the signprot attribute."""
        if self.signprot:
            if self.signprot.startswith('gna'):
                # Check if the model has a full heterotrimer by looking for 'gbb1_human' in the model name
                return 'gbb1_human' in self.model_name
            else: 
                log_or_raise(self.logger, f"Unexpected signprot format: {self.signprot}. " + 
                            "Expected to start with 'gna' for G-protein complexes. " + 
                            "Support for Arrestins is not implemented yet.", ValueError, self.error_handling)
        return False

    def parse_model_date_from_pdb_header(self):
        """Extract the model date from the PDB header."""
        if self.pdb_structure.header and 'deposition_date' in self.pdb_structure.header:
            if self.pdb_structure.header['deposition_date'] != '1909-01-08':
                return self.pdb_structure.header['deposition_date']
        
        if self.pdb_structure.header and 'release_date' in self.pdb_structure.header:
            if self.pdb_structure.header['release_date'] != '1909-01-08':
                return self.pdb_structure.header['release_date']

        if self.pdb_structure.header['head']:
            match = re.match(r'.+\s+(\d{4}-\d{2}-\d{2})', self.pdb_structure.header['head'])
            if match:
                return match.group(1)

        log_or_raise(self.logger, f"Could not parse model date from PDB header for model {self.model_name}." + 
                     "Amend the PDB header or provide a pdb_header_override in the parser configuration.", ValueError, self.error_handling)

    def annotate_ligand_sequence(self, parser_config):
        # Default to 'E' for complexes with a signprot, otherwise default to 'B' for complexes without a signprot
        self.ligand["pdb_chain_id"] = 'E' if self.signprot else 'B'  
        self.ligand["sequence_pdb"] = self.get_ligand_sequence_from_pdb()
        self.ligand["sequence_original"] = parser_config.old_seqs_dict.get(self.ligand.get("hash"), self.ligand["sequence_pdb"])

    def get_ligand_sequence_from_pdb(self):
        """Extract the ligand sequence from the PDB structure based on the ligand's chain ID."""
        chain_id = self.ligand.get("pdb_chain_id")
        if not chain_id:
            log_or_raise(self.logger, f"Ligand chain ID is not defined for model {self.model_name}.", ValueError, self.error_handling)
        
        sequence = ""
        for model in self.pdb_structure:
            for chain in model:
                if chain.id == chain_id:
                    for residue in chain:
                        resname = residue.get_resname()
                        one_letter = self.residue_to_one_letter.get(resname, 'X')
                        sequence += one_letter
        return sequence

    def unpack_model_name(self):
        """Unpack the model name into receptor, ligand, and signprot components."""
        parts = self.model_name.split('-')
        receptor = parts[0] if len(parts) > 0 else None
        ligand_temp = parts[1] if len(parts) > 1 else None
        signprot = parts[2] if len(parts) > 2 else None

        ligand_component = re.match(r'hashedseq\[(.+)\]', ligand_temp)
        ligand = None
        if ligand_component:
            ligand = {'hash': ligand_component.group(1)}

        if not receptor or not ligand:
            log_or_raise(self.logger, f"Invalid model name format: {self.model_name}. Expected format: receptor-hashedseq[ligand_hash]-signprot (signprot is optional).", ValueError, self.error_handling)

        return receptor, ligand, signprot

    def get_preferred_model_version_number(self, model_name, parser_config):
        
        """Get the index of the model to use based on the parser configuration
        
        If no override is specified, return the default model version number. If no default model is specified, return '1' as the default model version number.
        """

        if parser_config.default_model_version:
            if model_name in parser_config.default_model_version.get("override", {}):
                return parser_config.default_model_version["override"][model_name]       
            return parser_config.default_model_version.get("default_version_number", "1")
        return "1"  # Default to model version number '1' if no default model is specified

    def parse_metrics(self, parser_config):
        path_type_1 = os.sep.join([self.data_dir, self.model_name + '_metrics.csv'])
        path_type_2 = os.sep.join([self.data_dir, self.model_name + '.csv'])
        
        if os.path.exists(path_type_1):
            metrics_file_path = path_type_1
        elif os.path.exists(path_type_2):
            metrics_file_path = path_type_2
        else:   
            log_or_raise(self.logger, f"Metrics file not found for model {self.model_name}.", FileNotFoundError, self.error_handling)
        return AlphaFoldTwoComplexModelMetrics(metrics_file_path, error_handling=self.error_handling)

    def protein_from_entry_name(self):
        db_protein = Protein.objects.get(entry_name=self.receptor.lower())
        if not db_protein:
            log_or_raise(self.logger, f"Protein object not found for entry name: {self.receptor}", ValueError, self.error_handling)
        return db_protein

    def format_pdb_index(self):
        if 'peptide' in self.model_structure_type_slug:
            if self.signprot:
                return f'AFM_{self.receptor.upper()}_hashedseq[{self.ligand["hash"].upper()}]_{self.signprot.upper()}'
            else:
                return f'AFM_{self.receptor.upper()}_hashedseq[{self.ligand["hash"].upper()}]'
        else:
            return 'AFM_' + self.receptor.upper() + '_' + self.signprot.upper()

    def get_or_initialise_structure(self, receptor_protein, protein_state, protein_conformation):
        try:
            struct = Structure.objects.get(protein_conformation__protein=receptor_protein, pdb_code__index=self.format_pdb_index(), structure_type__slug=self.model_structure_type_slug)
        except Structure.DoesNotExist:
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

            try:
                struct.protein_conformation.generate_sites()
            except:
                pass

        return struct

    def get_or_create_protein_state(self):
        try:
            ps, created = ProteinState.objects.get_or_create(slug=self.receptor_state.lower(), defaults={'name': self.receptor_state})
            if created:
                self.logger.info('Created protein state {}'.format(ps.name))
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
            return pdbdata
        except Exception as e:
            log_or_raise(self.logger, f"Failed to create PdbData object for {self.model_name}: {e}", Exception, self.error_handling, parent_exception=e )

    def write_pdb_code_weblink(self, struct):
        try:
            web_resource = WebResource.objects.get(slug='pdb')
            pdb_code_weblink, created = WebLink.objects.get_or_create(index=self.format_pdb_index(), web_resource=web_resource)
            return pdb_code_weblink
        except Exception as e:
            log_or_raise(self.logger, f"Failed to create WebLink object for PDB code {self.format_pdb_index()}: {e}", Exception, self.error_handling, parent_exception=e)

    def get_or_create_structure_type(self):
        try:
            structure_type, created = StructureType.objects.get_or_create(slug=self.model_structure_type_slug, defaults={'name': self.model_structure_type_name})
            if created:
                self.logger.info('Created structure type {}'.format(structure_type))
            return structure_type
        except Exception as e:
            log_or_raise(self.logger, f"Failed to get or create structure type {self.model_structure_type_slug}: {e}", Exception, self.error_handling, parent_exception=e)

    def fetch_ligands_from_db(self):
        ligands = None
        if self.ligand.get("name", None):
            try:
                ligands = Ligand.objects.filter(name=self.ligand["name"]).first()
                if ligands:                    
                    return ligands
            except Exception as e:
                log_or_raise(self.logger, f"Error fetching ligand from database by name for model {self.model_name}: {e}", Exception, self.error_handling, parent_exception=e)

        if self.ligand.get("sequence_original", None):
            try:
                ligands = Ligand.objects.filter(sequence=self.ligand["sequence_original"])
                if ligands:
                    return ligands
            except Exception as e:
                log_or_raise(self.logger, f"Error fetching ligand from database by sequence for model {self.model_name}: {e}", Exception, self.error_handling, parent_exception=e)

        if not ligands:
            log_or_raise(self.logger, f"Ligand not found in database for model {self.model_name}. Please ensure the ligand is present in the database or provide a valid ligand name or sequence.", ValueError, self.error_handling)

    def create_ligand_peptide_structure(self, struct, ligands):
        for ligand in ligands:        
            try:
                ligand_peptide_structure, created = LigandPeptideStructure.objects.get_or_create(
                    structure=struct,
                    ligand=ligand,
                    chain=self.ligand.get("pdb_chain_id", None),
                    defaults={'model': None} 
                )
            except Exception as e:
                log_or_raise(self.logger, f"Error creating LigandPeptideStructure(s) for ligand {ligand.name} in model {self.model_name}: {str(e)} ", Exception, self.error_handling, parent_exception=e)

    def get_signprot_and_conformations(self, struct):
        signprot = None
        signprot_conf = None
        beta_protconf = None
        gamma_protconf = None

        if self.signprot:
            signprot = Protein.objects.get(entry_name=self.signprot.lower())
            signprot_conf = ProteinConformation.objects.get(protein=signprot)
            if self.beta_gamma:
                beta_protconf = ProteinConformation.objects.get(protein__entry_name='gbb1_human')
                gamma_protconf = ProteinConformation.objects.get(protein__entry_name='gbg2_human')
                sc = SignprotComplex.objects.get_or_create(alpha='B', protein=signprot, structure=struct,
                                                        beta_chain='C', gamma_chain='D', beta_protein=beta_protconf.protein, gamma_protein=gamma_protconf.protein).first()
            else:
                sc = SignprotComplex.objects.get_or_create(alpha='B', protein=signprot, structure=struct,
                                                        beta_chain=None, gamma_chain=None, beta_protein=None, gamma_protein=None).first()

            struct.signprot_complex = sc
            struct.save()         
        else:            
            signprot = None

        return signprot, signprot_conf, beta_protconf, gamma_protconf

    def create_extra_proteins(self, struct, signprot, signprot_conf, beta_protconf, gamma_protconf):
        sep = None
        sep_beta = None
        if signprot:        
            try:
                display_name = g_prot_dict[signprot.entry_name.split('_')[0].upper()]
                cat = 'G alpha'
            except:
                display_name = arr_dict[signprot.entry_name]
                cat = 'Arrestin'

            sep = StructureExtraProteins.objects.get_or_create(display_name=display_name, note=None, chain='B', category=cat, wt_coverage=100, protein_conformation=signprot_conf, structure=struct, wt_protein=signprot)
            if self.beta_gamma:
                sep_beta = StructureExtraProteins.objects.get_or_create(display_name='G&beta;1', note=None, chain='C', category='G beta', wt_coverage=100, protein_conformation=beta_protconf, structure=struct, wt_protein=beta_protconf.protein)
                sep_beta = StructureExtraProteins.objects.get_or_create(display_name='G&gamma;2', note=None, chain='D', category='G gamma', wt_coverage=100, protein_conformation=gamma_protconf, structure=struct, wt_protein=gamma_protconf.protein)
            # g beta - TO BE ADDED
            # g gamma - TO BE ADDED
        return sep, sep_beta

    def store_plddt(self, struct, receptor_protein, signprot, beta_protconf, gamma_protconf):
        #Adding plDDT for rendering
        resis = []
        for chain in self.pdb_structure.get_chains():
            for res in chain.get_residues():
                plddt = res['C'].get_bfactor()
                try:
                    if chain.get_id()=='A':
                        res_obj = Residue.objects.get(protein_conformation__protein=receptor_protein, sequence_number=res.get_id()[1])
                    elif chain.get_id()=='B' and signprot:
                        res_obj = Residue.objects.get(protein_conformation__protein=signprot, sequence_number=res.get_id()[1])
                    elif chain.get_id()=='C':
                        res_obj = Residue.objects.get(protein_conformation__protein=beta_protconf.protein, sequence_number=res.get_id()[1])
                    elif chain.get_id()=='D':
                        res_obj = Residue.objects.get(protein_conformation__protein=gamma_protconf.protein, sequence_number=res.get_id()[1])
                    r = StructureModelpLDDT(structure=struct, residue=res_obj, pLDDT=plddt)
                    resis.append(r)
                except Residue.DoesNotExist:
                    continue
        try:
            StructureModelpLDDT.objects.bulk_create(resis)
        except Exception as e:
            log_or_raise(self.logger, f"Error storing pLDDT values for model {self.model_name}: {str(e)}", Exception, self.error_handling, parent_exception=e)

    def build_contact_network(self, receptor, signprot):
        if signprot:
            do_complexes = True
        else:
            do_complexes = False
        compute_interactions(self.pdb_file_path, protein=receptor, signprot=signprot, do_complexes=do_complexes, save_to_db=True, file_input=True) # add do_complexes

    def write(self):

        try:
            with transaction.atomic():
        
                receptor_protein = self.protein_from_entry_name()

                protein_state = self.get_or_create_protein_state()

                protein_conformation = self.get_protein_conformation(receptor_protein)
                
                struct = self.get_or_initialise_structure(receptor_protein, protein_state, protein_conformation)

                signprot, signprot_conf, beta_protconf, gamma_protconf = self.get_signprot_and_conformations(struct)      

                ligands = self.fetch_ligands_from_db()

                self.create_ligand_peptide_structure(struct, ligands)

                self.metrics.save(struct)

                self.create_extra_proteins(struct, signprot, signprot_conf, beta_protconf, gamma_protconf)

                self.store_plddt(struct, receptor_protein, signprot, beta_protconf, gamma_protconf)

                self.build_contact_network(struct, signprot)

        except Exception as e:
            log_or_raise(self.logger, f"Error writing model {self.model_name} to database: {str(e)}", Exception, self.error_handling, parent_exception=e)
   
class AlphaFoldTwoComplexModelParserConfig(BaseModelParserConfig):

    """Configuration class for AlphaFoldTwoComplexModelParser"""
    
    def __init__(self, model_set_name, data_dir=None, cleaned_seq_csv=None, pdb_header_override=None, model_receptor_state=None, pdb_preferred_chain='A', error_handling="log"):
        super().__init__(model_set_name, data_dir=data_dir, pdb_header_override=pdb_header_override, error_handling=error_handling)

        self.cleaned_seq_csv = cleaned_seq_csv
        self.model_receptor_state = model_receptor_state
        self.pdb_preferred_chain = pdb_preferred_chain
        if cleaned_seq_csv and os.path.exists(cleaned_seq_csv):
            self.old_seqs_dict = self.generate_original_seq_lookup()
        else:
            log_or_raise(self.logger, f"Cleaned sequence CSV file not found at {cleaned_seq_csv}.", FileNotFoundError, self.error_handling)

    def generate_original_seq_lookup(self):
        df_cleaned_seqs = pd.read_csv(self.cleaned_seq_csv)
        df_cleaned_seqs['backwards_hex_cleaned_seq_hash_col'] = df_cleaned_seqs['cleaned_seq_hash_col'].apply(lambda x: hex(x)[2:][::-1])
        df_cleaned_seqs['pdb_file_hashseq'] = df_cleaned_seqs['cleaned_seq_hash'] + df_cleaned_seqs['backwards_hex_cleaned_seq_hash_col']
        return {k:v for k,v in zip(df_cleaned_seqs['pdb_file_hashseq'],df_cleaned_seqs['old_sequence'])}

class AlphaFoldTwoComplexModelParser(BaseModelParser):   
    
    """Parses a directory of AlphaFoldTwoComplex models organized by model set name and model name
    
    The expected directory structure is as follows:
    /structure_data/{model_set_name}/{receptor}-{ligand}-{signprot}(optional)/model_{model_version_number}.pdb and metrics_{model_version_number}.csv
    """

    def __init__(self, config):
        """Initialize the parser with a configuration object."""
        super().__init__(config)
        self.config = config
        self.models = []
        
    def load_models(self):
        """Process each model directory in the model set directory and return a list of AlphaFoldTwoComplexModel instances."""
        try:
            model_dirs = list(filter(os.path.isdir, [os.path.join(self.config.data_dir, f) for f in os.listdir(self.config.data_dir)])) #/structure_data/{model_set_name}/[*]
            models = []
            for model_path in model_dirs:
                models.append(AlphaFoldTwoComplexModel(model_path, self.config))
            self.models = models
        except Exception as e:
            log_or_raise(self.logger, f"Error processing models: {str(e)}", Exception, self.error_handling, parent_exception=e)

    def write_models(self):
        """Write each model to the database."""
        for model in self.models:
            model.write()
    

#DEBUG
# sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
# class settings:
#     DATA_DIR = "/home/rnd457/gpcr"


def tester():

    config = AlphaFoldTwoComplexModelParserConfig(model_set_name="AlphaFold_multimer",
                                                  cleaned_seq_csv="/home/rnd457/gpcr/structure_data/AlphaFold_multimer/cleaned_seqs.csv",
                                                  model_receptor_state="Active",
                                                  pdb_preferred_chain="A",
                                                  error_handling="raise")
    af_parser = AlphaFoldTwoComplexModelParser(config)
    af_parser.load_models()
    af_parser.write_models()


tester()
#DEBUG