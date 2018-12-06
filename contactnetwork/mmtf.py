import urllib.request
import os

from Bio.PDB import *
from Bio.PDB.mmtf import *
from Bio.PDB.PDBExceptions import PDBException

from io import StringIO

# Download the specified MMTF file from RCSB.
def mmtf_get_structure(pdb_name):
    # Download the PDB
    response = urllib.request.urlopen('http://mmtf.rcsb.org/v1.0/full/%s' % pdb_name)

    return MMTFParser.get_structure_from_url(pdb_name)
