from protein.models import Protein, ProteinSegment, Residue, ProteinConformation, ProteinCouplings
from structure.models import Structure, StructureModel
from structure.assign_gn_foldseek import *
import csv
from django.core.management.base import *
import pandas as pd 
from Bio.PDB import PDBParser, PDBIO, Select
from io import StringIO

import re
import time
from django.db.models import Q, OuterRef, Prefetch
import shutil
from pathlib import Path
from django.conf import settings

class ResidueBFactorSelect(Select):
    """A selection class for filtering residues based on the B-factor of their CA atoms."""

    def __init__(self, bfactor_range=(1.0, 9.0)):
        self.bfactor_range = bfactor_range

    def accept_residue(self, residue):
        for atom in residue:
            if atom.get_id() == 'CA' and self.bfactor_range[0] <= atom.get_bfactor() < self.bfactor_range[1]:
                return True
        return False


class Command(BaseCommand):
    def handle(self, *args, **options):

        output_dir = Path(os.path.join(settings.DATA_DIR, 'structure_data'))
        db_id = 'foldseek_db'
        raw_types = ['1', '2', '3']
        ref_types = ['4', '5']
        af_types = ['7']

        # Clear and recreate output directories
        shutil.rmtree(output_dir, ignore_errors=True)
        dirs = [f'{prefix}_{db_id}{suffix}' for prefix in ['raw', 'ref', 'af'] for suffix in ['', '_trim']]
        for dir in dirs:
            (output_dir / dir).mkdir(parents=True, exist_ok=True)

        # Process raw and reference structures
        structures = Structure.objects.filter(structure_type__in=raw_types + ref_types).select_related("pdb_data", "pdb_code")
        for structure in structures:
            process_structure(structure, raw_types, ref_types, af_types, db_id, output_dir)

        # Process AF multistate structures
        af_multistate = StructureModel.objects.filter(main_template__isnull=True).select_related('pdb_data', 'state', 'protein').prefetch_related(
            Prefetch('protein__proteinconformation_set', queryset=ProteinConformation.objects.filter(protein__isnull=False), to_attr='protein_conformation')
        )
        for structure in af_multistate:
            process_af_structure(structure, db_id, output_dir)

def process_structure(structure, raw_types, ref_types, af_types, db_id, output_dir):
    pdb_data = structure.pdb_data.pdb
    pdb_code = structure.pdb_code.index
    s_type = str(structure.structure_type.id)
    structure_type = "raw" if s_type in raw_types else "ref" if s_type in ref_types else "af" if s_type in af_types else 'empty'
    assign_grn = GenericNumberingFromDB(structure, pdb_data)
    save_structure(assign_grn.assign_generic_numbers(), structure_type, db_id, pdb_code, output_dir)
    

def process_af_structure(structure, db_id, output_dir):
    pdb_data = structure.pdb_data.pdb
    pdb_code = structure.protein.entry_name
    state = structure.state.slug
    assign_grn = GenericNumberingFromDB1(structure, pdb_data)
    save_structure(assign_grn.assign_generic_numbers(), 'af', db_id, f'{pdb_code}_{state}', output_dir)


def save_structure(structure, structure_type, db_id, pdb_code, output_dir):
    io = PDBIO()
    io.set_structure(structure)
    output_pdb = f'{output_dir}/{structure_type}_{db_id}/{pdb_code}_{structure_type}_info.pdb'
    io.save(output_pdb)
    trimmed_output_pdb = f'{output_dir}/{structure_type}_{db_id}_trim/{pdb_code}_{structure_type}_info.pdb'
    process_pdb(output_pdb, trimmed_output_pdb)


def process_pdb(input_pdb, output_pdb, bfactor_range=(1.0, 9.0)):
    parser = PDBParser()
    structure = parser.get_structure('structure', input_pdb)
    selector = ResidueBFactorSelect(bfactor_range)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, selector)