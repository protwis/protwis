from interaction import *
import numpy

AA = {'ALA', 'ARG', 'ASN', 'ASP',
      'CYS', 'GLN', 'GLU', 'GLY',
      'HIS', 'ILE', 'LEU', 'LYS',
      'MET', 'PHE', 'PRO', 'SER',
      'THR', 'TRP', 'TYR', 'VAL'}

AROMATIC_AA = {'TYR', 'TRP', 'PHE', 'HIS'}

POS_CHARGED_AA = {'ARG', 'LYS'}  # skip ,'HIS'
NEG_CHARGED_AA = {'ASP', 'GLU'}

HYDROPHOBIC_AA = {'ALA', 'CYS', 'PHE',
                  'ILE', 'LEU', 'MET',
                  'PRO', 'VAL', 'TRP',
                  'TYR'}

VDW_RADII = {'H': 1.2, 'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.80}


def is_aromatic_aa(res):
    return res.get_resname() in AROMATIC_AA


def is_aa(res):
    return res.get_resname() in AA


def is_charged(res):
    return is_pos_charged(res) or is_neg_charged(res)


def is_pos_charged(res):
    return res.get_resname() in POS_CHARGED_AA


def is_neg_charged(res):
    return res.get_resname() in NEG_CHARGED_AA


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

        return zip(ring_centers, ring_normals)
    except:
        return []


# Returns a list of positively charges atoms in a residue
def get_pos_charged_atom_names(res):
    resname = res.get_resname()

    # Only ARG, LYS and HIS are positively charged
    if resname == 'ARG':
        return ['NH1', 'NH2', 'NE']
    elif resname == 'LYS':
        return ['N']
    elif resname == 'HIS':
        # TODO: Implement using e.g. ProPka
        return []


# Returns a list of positively charges atoms in a residue
def get_neg_charged_atom_names(res):
    resname = res.get_resname()

    # Only ASP and GLU are negatively charged
    if resname == 'ASP':
        return ['OD1', 'OD2']
    elif resname == 'GLU':
        return ['OE1', 'OE2']
