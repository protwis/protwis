from interaction import *
import numpy

# Amino acid names.
AA = {'ALA', 'ARG', 'ASN', 'ASP',
      'CYS', 'GLN', 'GLU', 'GLY',
      'HIS', 'ILE', 'LEU', 'LYS',
      'MET', 'PHE', 'PRO', 'SER',
      'THR', 'TRP', 'TYR', 'VAL'}

# Side-chain H-bond donors
HBD = {'ARG', 'ASN', 'GLN', 'HIS',
       'LYS', 'SER', 'THR', 'TRP',
       'TYR'}

# Side-chain H-bond acceptors
HBA = {'ASN', 'ASP', 'GLN', 'GLU',
       'HIS', 'SER', 'THR', 'TYR'}

# Aromatic amino acid names.
AROMATIC_AA = {'TYR', 'TRP', 'PHE', 'HIS'}

# Positively charged amino acid names.
POS_CHARGED_AA = {'ARG', 'HIS', 'LYS'}

# Negatively charged amino acid names.
NEG_CHARGED_AA = {'ASP', 'GLU'}

# Hydrophobic amino acid shortnames.
HYDROPHOBIC_AA = {'ALA', 'CYS', 'PHE',
                  'ILE', 'LEU', 'MET',
                  'PRO', 'VAL', 'TRP',
                  'TYR'}

# Hydrogen placement properties taken from CHARMM36
ANGLE_REFERENCES = {'ARG':
                        {'NH1': [['CZ', 120.61, 0.9903, 'NE'],      # HH11
                                ['CZ', 116.29, 1.0023, 'NH2']],     # HH12
                        'NH2': [['CZ', 119.91, 0.9899, 'NH1'],      # HH21
                                ['CZ', 116.88, 0.9914, 'NE']],      # HH22
                        'NE': [['CD', 113.14, 1.0065, 'NH1']]},      # HE
                    'ASN':
                        {'ND2': [['CG', 117.35, 0.9963, 'CB'],      # HD21
                                ['CG', 120.05, 0.9951, 'OD1']]},     # HD22
                    'CYS':
                        {'SG': [['CB', 97.15, 1.3341]]},             # HG1
                    'GLN':
                        {'NE2': [['CD', 116.86, 0.9959, 'CG'],      # HE21
                                ['CD', 119.83, 0.9943, 'OE1']]},     # HE22
                    'HIS':
                        {'ND1': [['CG', 126.09, 1.0020, 'CD2']],    # HD1
                        'NE2': [['CD2', 125.52, 1.0020, 'CG']]},     # HE2
                    'LYS':
                        {'NZ': [['CE', 110.020, 1.0404]]},           # HZ1-3
                    'SER':
                        {'OG': [['CB', 107.08, 0.9655]]},            # HG1
                    'THR':
                        {'OG1': [['CB', 105.45, 0.9633]]},           # HG1
                    'TRP':
                        {'NE1': [['CD1', 124.68, 0.9767, 'CG']]},    # HE1
                    'TYR':
                        {'OH': [['CZ', 107.47, 0.9594]]},			# HH
                    'HOH':
                        {'O': []},                                   # H1/H2
                    'BB':
                        {'N': [['CA', 116.67, 0.9973, 'C']]},		# NH - C-atom from previous connected residue
                }

ACCEPTING_REFERENCES = {'ASN':
                            {'OD1': 'CG'},
                    'ASP':
                        {'OD1': 'CG',
                    	'OD2': 'CG'},
                    'CYS':
                        {'SG': 'CB'},
                    'GLN':
                        {'OE1': 'CD'},
                    'GLU':
                        {'OE1': 'CD',
                    	'OE2': 'CD' },
                    'HIS':
                        {'ND1': ['CG', 'CE1'],
                    	'NE2': ['CD2', 'CE1']},
                    'SER':
                        {'OG': 'CB' },
                    'THR':
                        {'OG1': 'CB' },
                    'TYR':
                        {'OH': 'CZ' },
                    'HOH':
                        {'O': []},
                    'BB':
                        {'O': 'C' }}

# Van-der-Waals radii of the elements.
VDW_RADII = {'H': 1.2, 'C': 1.7, 'N': 1.55, 'O': 1.52, 'P': 1.80, 'S': 1.80, 'F': 1.47, 'CL': 1.75}

# Factor multiplied with VDW radii when determining distance between elements.
VDW_TRESHOLD_FACTOR = 1.1


def is_aromatic_aa(res):
    return res.get_resname() in AROMATIC_AA


def is_aa(res):
    return res.get_resname() in AA

def is_charged(res):
    return is_pos_charged(res) or is_neg_charged(res)

def is_hba(res):
    return res.get_resname() in HBA or is_water(res)

def is_hbd(res):
    return res.get_resname() in HBD or is_water(res)

def is_pos_charged(res):
    return res.get_resname() in POS_CHARGED_AA

def is_neg_charged(res):
    return res.get_resname() in NEG_CHARGED_AA

# TODO: add is solvent

def is_water(res):
    return res.get_resname() == "HOH"

# Get lists of atoms of all rings in residue
def get_ring_atom_name_lists(res):
    atom_name_lists = []

    res_name = res.get_resname()

    if res_name == 'TYR':
        atom_name_lists.append(['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'])
    elif res_name == 'TRP':
        # Tryptophane has two rings
        atom_name_lists.append(['CG', 'CD1', 'CD2', 'CE2', 'NE1'])
        atom_name_lists.append(['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'])
    elif res_name == 'PHE':
        atom_name_lists.append(['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'])
    elif res_name == 'HIS':
        atom_name_lists.append(['CG', 'ND1', 'CE1', 'NE2', 'CD2'])
    else:
        return []

    # Atom objects
    res_atoms = res.child_dict

    # Filter lists where all atoms are not present
    atom_name_lists = [ l for l in atom_name_lists if all(k in res_atoms for k in l) ]

    # Get lists of atoms in for each ring; one list containing one list for TYR, PHE and HIS, two for TRP.
    return [[a for a in res_atoms.values() if a.name in a_l] for a_l in atom_name_lists]


# Returns a list of pairs, [(center_coord, normal)], representing
# the coordinates of the center of the rings in a residue
# as well as a normal vector of the plane in which the ring lies.
def get_ring_descriptors(res):
    try:
        # Get the coords of each atom in the rings; a list containing one list for TYR, PHE and HIS, two for TRP.
        ring_atom_coords = [[a.coord for a in a_l] for a_l in get_ring_atom_name_lists(res)]

        # Atoms are equidistantly places on the circumference on a circle, so the average will yield the center.
        ring_centers = [numpy.divide(numpy.sum(a_c_l, axis=0), len(a_c_l)) for a_c_l in ring_atom_coords]

        # Get normals by taking the cross product of two vectors in the ring plane. All atom coordinates are co-planar.
        ring_normals = [numpy.cross(numpy.subtract(a_c_l[0], a_c_l[1]), numpy.subtract(a_c_l[0], a_c_l[2])) for a_c_l in ring_atom_coords]

        return list(zip(ring_centers, ring_normals))
    except:
        return []

# FROM: https://stackoverflow.com/questions/38987
# Given two dicts, merge them into a new dict as a shallow copy.
def merge_two_dicts(x, y):
    z = x.copy()
    z.update(y)
    return z

 # Returns hydrogen placement angles/distances (CHARMM36) + reference atoms
 # for the calculation of hydrogen bonds
def get_hbond_donor_references(res):
    resname = res.get_resname()
    if resname == "HOH":
        return ANGLE_REFERENCES[resname]
    elif resname in ANGLE_REFERENCES:
        # return backbone + residue specific
        # Python >= 3.5
        #return {**ANGLE_REFERENCES[resname], **ANGLE_REFERENCES['BB']}
        # Python < 3.5
        return merge_two_dicts(ANGLE_REFERENCES[resname], ANGLE_REFERENCES['BB'])
    else:
        # Only return backbone
        return ANGLE_REFERENCES['BB']

# Returns reference atoms for hydrogen bond acceptors
# for the calculation of hydrogen bonds
def get_hbond_acceptors(res):
    resname = res.get_resname()
    if resname == "HOH":
        return ACCEPTING_REFERENCES[resname]
    elif resname in ACCEPTING_REFERENCES:
        # return backbone + residue specific
        # Python >= 3.5
        #return {**ACCEPTING_REFERENCES[resname], **ACCEPTING_REFERENCES['BB']}
        # Python < 3.5
        return merge_two_dicts(ACCEPTING_REFERENCES[resname], ACCEPTING_REFERENCES['BB'])
    else:
        # Only return backbone
        return ACCEPTING_REFERENCES['BB']

# redefine the unit_vector function to replace the internal function
def get_unit_vector(vector):
    return vector / numpy.linalg.norm(vector)

# Returns a list of positively charges atoms in a residue
def get_pos_charged_atom_names(res):
    resname = res.get_resname()

    # For now: simple assumption that they are always charged
    # in vicinity of acidic residues
    atomnames = []
    if resname == 'ARG':
        atomnames = ['CZ', 'NE', 'NH1', 'NH2']
    elif resname == 'LYS':
        atomnames = ['NZ']
    elif resname == 'HIS':
        # TODO: Implement using e.g. ProPka
        atomnames = ['ND1', 'NE2']
    else:
        return []

    return match_atomselection_residue(res, atomnames)


# Returns a list of positively charges atoms in a residue
def get_neg_charged_atom_names(res):
    resname = res.get_resname()

    # For now: simple assumption that they are always charged
    # in vicinity of basic residues
    atomnames = []
    if resname == 'ASP':
        atomnames = ['OD1', 'OD2']
    elif resname == 'GLU':
        atomnames = ['OE1', 'OE2']
    else:
        return atomnames

    return match_atomselection_residue(res, atomnames)

def get_charged_atom_names(res):
    if is_pos_charged(res):
        return get_pos_charged_atom_names(res)
    elif is_neg_charged(res):
        return get_neg_charged_atom_names(res)
    else:
        return []

# Returns the list of atom IDs that are actually present in for the residue
def match_atomselection_residue(res, atomnames):
    res_atoms = res.child_dict
    # print(res_atoms)
    return [ name for name in atomnames if name in res_atoms ]
