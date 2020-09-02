import urllib.request
import os

from Bio.PDB import *
from Bio.PDB.PDBExceptions import PDBException

from io import StringIO

# Download the specified PDB file from RCSB.
def pdb_get_structure(pdb_name):
    # Download the PDB
    response = urllib.request.urlopen('https://www.rcsb.org/pdb/files/%s.pdb' % pdb_name)

    # Check that the requested PDB name is valid
    if response.getcode() == 404:
        raise PDBException('No PDB with name %s in RCSB.' % pdb_name)
    # And that it could be retrieved
    elif response.getcode() != 200:
        raise PDBException('Could not retrieve PDB file from RCSB.')

    # Decode the response
    contents = response.read().decode("utf-8")

    # Create filehandle
    f = StringIO(contents)

    # Return the structure
    p = PDBParser(QUIET=True)
    return p.get_structure(pdb_name, f)
