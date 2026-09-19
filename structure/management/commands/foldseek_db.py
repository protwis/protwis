from build.management.commands.base_build import Command as BaseBuild
from structure.models import Structure, StructureModel
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Model import Model
from structure.assign_generic_numbers_gpcr import GenericNumbering as as_gn
from io import StringIO
import os
import shutil
from protwis import settings

# Number of DB rows to keep in memory at once while streaming each queryset
# (pdb_data__pdb is a large TextField, so this bounds peak memory).
QUERYSET_CHUNK_SIZE = 200

class ResidueBFactorSelect(Select):
    """
    A selection class for filtering residues based on the B-factor of their CA atoms.
    Only residues with CA atom B-factors within the specified range are selected.
    """
    def __init__(self, bfactor_range=(-8.0, 8.0)):
        """
        Initializes the selection object with the specified B-factor range.

        Parameters:
        -----------
        bfactor_range : tuple, optional
            The inclusive range of B-factors to select residues. Default is (-8.0, 8.0).
            Negative range is to catch bulges.
        """
        self.bfactor_range = bfactor_range

    def accept_residue(self, residue):
        """
        Checks if a residue should be accepted based on the B-factor of its CA atom.

        Parameters:
        -----------
        residue : Bio.PDB.Residue
            The residue to check.

        Returns:
        --------
        bool
            True if the residue should be accepted (CA atom's B-factor is within the range), False otherwise.
        """
        for atom in residue:
            if atom.get_id() == 'CA':
                bfactor = atom.get_bfactor()
                if self.bfactor_range[0] < bfactor < self.bfactor_range[1]:
                    return True
        return False

def create_db_dir(structure_type, output_dir):
    """
    Creates an output directory for the specified structure type.

    Parameters:
    -----------
    structure_type : str
        The structure type ('raw', 'ref', or 'af').
    output_dir : str
        The base output directory.

    Returns:
    --------
    str
        The path to the created output directory.
    """
    db_dir = f'{structure_type}_foldseek_db_trim'
    output_db_dir = os.path.join(output_dir, db_dir)
    if os.path.exists(output_db_dir):
        shutil.rmtree(output_db_dir)
    os.mkdir(output_db_dir)
    return output_db_dir

def process_structure(pdb_code, pdb_data, output_filename, preferred_chain=None, residue_cache=None):
    """
    Processes a single structure: parses PDB data, assigns generic numbers, saves annotated structure.

    Parameters:
    -----------
    pdb_code : str
        The PDB code or identifier of the structure.
    pdb_data : str
        The PDB data as a string.
    output_filename : str
        The filename to save the annotated structure.
    preferred_chain : str, optional
        The identifier of the preferred chain to extract. If None, the entire structure is used.
    residue_cache : dict, optional
        A dict (keyed by reference protein id) reused across multiple calls within the same
        worker process, so repeated structures of the same receptor don't re-fetch identical
        reference-residue data from the DB. Safe to reuse because it only ever lives for the
        duration of this one command invocation. See GenericNumbering.__init__.
    """
    parser = PDBParser(QUIET=True)
    pdb_io = StringIO(pdb_data)
    structure = parser.get_structure(pdb_code, pdb_io)

    if preferred_chain:
        # Extract the preferred chain
        model = structure[0]  # Assuming only one model
        if preferred_chain in model:
            chain = model[preferred_chain]
            # Create a new model with the chain
            new_model = Model(model.id)
            new_model.add(chain)
            structure_to_use = new_model
        else:
            print(f"Chain {preferred_chain} not found in structure {pdb_code}")
            return
    else:
        # Use the entire structure
        structure_to_use = structure[0]  # Assuming first model

    # Assign generic numbers
    gn = as_gn(structure=structure_to_use, pdb_code=pdb_code, residue_cache=residue_cache)
    gn.assign_generic_numbers()
    annotated_structure = gn.get_annotated_structure()

    # Save annotated structure
    io = PDBIO()
    io.set_structure(annotated_structure)
    io.save(output_filename, ResidueBFactorSelect())
    print(f"Saved selected residues to {output_filename}")

def build_items(structure_type, output_dir):
    """
    Creates the output directory for the given structure_type and builds the list of
    work items (one per structure) to be processed, without doing any of the (expensive,
    per-structure) generic-numbering work itself. Kept separate from processing so that
    the work items can be handed out to worker processes.

    Parameters:
    -----------
    structure_type : str
        The structure type to process ('raw', 'ref', or 'af').
    output_dir : str
        The base output directory.

    Returns:
    --------
    list of tuple
        Each tuple is (pdb_code, pdb_data, output_filename, preferred_chain).
    """
    output_db_dir = create_db_dir(structure_type, output_dir)
    items = []

    if structure_type == 'raw':
        # Query experimental structures from the database
        exp_structures = Structure.objects.filter(
            structure_type__slug__in=['x-ray-diffraction', 'electron-microscopy', 'electron-crystallography']
        ).values_list(
            'pdb_code__index',
            'pdb_data__pdb',
            'preferred_chain'
        ).iterator(chunk_size=QUERYSET_CHUNK_SIZE)

        for pdb_code, pdb_data, preferred_chain in exp_structures:
            output_filename = f"{output_db_dir}/{pdb_code}_raw_info.pdb"
            items.append((pdb_code, pdb_data, output_filename, preferred_chain))

    elif structure_type == 'af':
        af_structures = StructureModel.objects.filter(main_template_id__isnull=True).values_list(
            'protein__entry_name',
            'pdb_data__pdb',
            'state__slug'
        ).iterator(chunk_size=QUERYSET_CHUNK_SIZE)

        for pdb_code, pdb_data, state in af_structures:
            output_filename = f"{output_db_dir}/{pdb_code}_{state}_af_info.pdb"
            items.append((pdb_code, pdb_data, output_filename, None))

    elif structure_type == 'ref':
        inactive_structures = StructureModel.objects.filter(main_template_id__isnull=False).values_list(
            'protein__entry_name',
            'pdb_data__pdb'
        ).iterator(chunk_size=QUERYSET_CHUNK_SIZE)

        active_structures = Structure.objects.filter(
            structure_type__slug__in=['af-signprot-refined-cem', 'af-signprot-refined-xray']
        ).values_list(
            'pdb_code__index',
            'pdb_data__pdb',
            'preferred_chain'
        ).iterator(chunk_size=QUERYSET_CHUNK_SIZE)

        # Active structures
        for pdb_code, pdb_data, preferred_chain in active_structures:
            output_filename = f"{output_db_dir}/{pdb_code}_ref_info.pdb"
            items.append((pdb_code, pdb_data, output_filename, preferred_chain))

        # Inactive structures
        for pdb_code, pdb_data in inactive_structures:
            output_filename = f"{output_db_dir}/{pdb_code}_ref_info.pdb"
            items.append((pdb_code, pdb_data, output_filename, None))

    return items

class Command(BaseBuild):
    help = "Assigns generic numbers to structures and extracts residues based on B-factors, for Foldseek databases"

    def add_arguments(self, parser):
        super().add_arguments(parser=parser)
        parser.add_argument(
            '--structure_type',
            nargs='*',
            choices=['raw', 'ref', 'af'],
            help='Specify which structure types to process. If omitted, all types are processed.',
        )

    def handle(self, *args, **options):
        # Define the output directory
        output_dir = os.path.join(settings.DATA_DIR, 'structure_data')

        # Get the structure types to process
        structure_types = options['structure_type'] or ['raw', 'ref', 'af']

        # Build the full list of work items up front (this also creates the output
        # directories), then hand them out to worker processes via prepare_input/main_func.
        self.items = []
        for structure_type in structure_types:
            print(f"Collecting structures of type: {structure_type}")
            self.items.extend(build_items(structure_type, output_dir))

        print(f"Processing {len(self.items)} structures with {options['proc']} process(es)")
        self.prepare_input(options['proc'], self.items)

    def main_func(self, positions, iteration, count, lock):
        # One residue_cache per worker process, reused across every item this worker
        # handles (see process_structure's docstring / GenericNumbering.__init__).
        residue_cache = {}

        # Work-stealing loop: each worker grabs the next unclaimed item from self.items
        # (shared via fork's copy-on-write) using the shared counter/lock, so that workers
        # finishing early (structures vary a lot in size/chain count) pick up more work
        # rather than sitting idle.
        while count.value < len(self.items):
            with lock:
                if count.value >= len(self.items):
                    break
                index = count.value
                count.value += 1

            pdb_code, pdb_data, output_filename, preferred_chain = self.items[index]
            if preferred_chain:
                print(f"Processing structure {pdb_code} with preferred chain {preferred_chain}")
            else:
                print(f"Processing structure {pdb_code}")
            process_structure(pdb_code, pdb_data, output_filename, preferred_chain, residue_cache=residue_cache)
