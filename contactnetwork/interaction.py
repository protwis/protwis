from contactnetwork.residue import *
from Bio.PDB.Polypeptide import *

import math

class InteractingPair:
    'Common base class for all interactions'
    def __init__(self, res1, res2, interactions):
        self.res1 = res1
        self.res2 = res2
        self.interactions = interactions

    def add_interaction(self, interaction):
        self.interactions.append(interaction)

    def get_interactions(self):
        return self.interactions

    def get_residue_1(self):
        return self.res1

    def get_residue_2(self):
        return self.res2

    def get_pymol_selection_code(self):
        # TODO fix static chain selection
        return "select chain A and (resi {0}+{1}); zoom sele; show sticks, sele".format(str(self.res1.id[1]), str(self.res2.id[1]))

    def get_interaction_text(self):
        text = 'Interaction (\'{0}\',\'{1}\') at ({2},{3}) ('.format(self.res1.resname, self.res2.resname, str(self.res1.id[1]), str(self.res2.id[1]))
        first = True

        if self.interactions:
            for i in self.interactions:
                if not first:
                    text += ', '
                first = False
                text += i.get_name()

        text += ')'

        return text

    def get_interaction_json(self, generic):
        text = '[\'{0}\',\'{1}\',{2},\'{3}\',\'{4}\',\'{5}\',{6}, ["'.format(self.res1.parent.id, three_to_one(self.res1.get_resname()), self.res1.id[1], generic, self.res2.parent.id, three_to_one(self.res2.get_resname()), self.res2.id[1])
        first = True

        if self.interactions:
            for i in self.interactions:
                if not first:
                    text += ', "'
                first = False
                text += i.get_name() + '"'

        text += ']]'

        return text


class Interaction(object):
    def __init__(self):
        pass


class VanDerWaalsInteraction(Interaction):
    def get_name(self):
        return "van-der-waals"


class HydrophobicInteraction(Interaction):
    def get_name(self):
        return "hydrophobic"

class HydrogenBondInteraction(Interaction):
    def get_name(self):
        return "h-bond"

class IonicInteraction(Interaction):
    def get_name(self):
        return "ionic"

class PolarInteraction(Interaction):
    def __init__(self, is_charged_res1, is_charged_res2):
        Interaction.__init__(self)
        self.is_charged_res1 = is_charged_res1
        self.is_charged_res2 = is_charged_res2

    def get_name(self):
        return "polar"


class PolarSidechainSidechainInteraction(PolarInteraction):
    def get_name(self):
        return "polar-sidechain-sidechain"


class PolarBackboneSidechainInteraction(PolarInteraction):
    def get_name(self):
        return "polar-backbone-sidechain"


class PolarSideChainBackboneInteraction(PolarInteraction):
    def get_name(self):
        return "polar-sidechain-backbone"

class PolarWaterInteraction(Interaction):
    def get_name(self):
        return "polar-water"

class AromaticInteraction(Interaction):
    def get_name(self):
        return "aromatic"


class FaceToFaceInteraction(AromaticInteraction):
    def get_name(self):
        return "face-to-face"


class EdgeToFaceInteraction(AromaticInteraction):
    def get_name(self):
        return "edge-to-face"


class FaceToEdgeInteraction(AromaticInteraction):
    def get_name(self):
        return "face-to-edge"


class PiCationInteraction(AromaticInteraction):
    def get_name(self):
        return "pi-cation"


class CationPiInteraction(AromaticInteraction):
    def get_name(self):
        return "cation-pi"

class WaterMediated(Interaction):
    def get_name(self):
        return "water-mediated"

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / numpy.linalg.norm(vector)


def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return numpy.arccos(numpy.clip(numpy.dot(v1_u, v2_u), -1.0, 1.0))


# Return the acute angles between two vectors
def angle_between_plane_normals(v1, v2):
    return min([angle_between(v1, v2), angle_between(v1, numpy.multiply(v2, -1.0))])


def distance_between(v1, v2):
    return numpy.linalg.norm(numpy.subtract(v1, v2))


# Checks if two residues have a face to face interaction
def has_face_to_face_interaction(res1, res2):
    rings_res1 = get_ring_descriptors(res1)
    rings_res2 = get_ring_descriptors(res2)

    '''
    if any([(angle_between_plane_normals(r1[1], r2[1]) < 0.34906585
            and numpy.linalg.norm(numpy.subtract(r1[0], r2[0])) < 5.0)
            for r1 in rings_res1 for r2 in rings_res2]):
        print "Distance between ring centers: %f" % numpy.linalg.norm(numpy.subtract(rings_res1[0][0], rings_res2[0][0]))
        print "Center of first ring is %f,%f,%f, second is %f,%f,%f!".format()
    '''

    # Make sure that the acute angle between the planes are less than (or eq.) 20
    # degrees and that the distance between centers is less than  (or eq.) 5 Angstrom.
    #return any([(angle_between_plane_normals(r1[1], r2[1]) <= 0.34906585)
    return any([(math.degrees(angle_between_plane_normals(r1[1], r2[1])) <= 20)
                and (numpy.linalg.norm(numpy.subtract(r1[0], r2[0])) <= 5.0)
                for r1 in rings_res1 for r2 in rings_res2])


# Checks if two residues have an edge to face interaction
def has_edge_to_face_interaction(res1, res2):
    res1_desc = get_ring_descriptors(res1)
    res2_desc = get_ring_descriptors(res2)

    # Make sure the edge atom closest to the center is closer than
    # 4.5 angstroms and that the perpendicular angles is with
    # +/- 30 degrees, i.e. the acute angle is greater than 60
    # degrees = 1.04719755 radians.
    #return any([(abs(angle_between_plane_normals(r1[1], r2[1]) - 1.5707963267) < 0.523598776)
    return any([(math.degrees(abs(angle_between_plane_normals(r1[1], r2[1]) - 1.5707963267)) <= 30)
                and (numpy.linalg.norm(numpy.subtract(r1[0], r2[0])) <= 5.2)
                for r1 in res1_desc for r2 in res2_desc])


def has_pi_cation_interaction(res1, res2):
    # Only aromatic residues will have ring descriptors
    res1_descs = get_ring_descriptors(res1)

    # Only positively charged atoms will have charged atoms
    res2_pos_atom_names = get_pos_charged_atom_names(res2)

    # Check if any charged residue atom is closer to any ring centroid than 6 angstroms.
    # TODO: also adjust this to add an angle check of 30 degrees
    try:
        close_enough = any(
            [any(
                [distance_between(res2.child_dict[atom_name].coord, desc[0]) <= 6 for desc in res1_descs]
            ) for atom_name in res2_pos_atom_names])
    except KeyError:
        return False

    # print "select chain A and (resi {0}+{1}); zoom sele; show sticks, sele".format(str(res1.id[1]), str(res2.id[1]))

    # The epsilon carbon of LYS is 2.4x as likely to be closer
    # to a center of a ring than the protonated nitrogen, see:
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC22230/
    # consider using this information!

    return close_enough

def get_polar_hbonds_interactions(res1, res2):
    if has_hbond_interaction(res1, res2):
        return [HydrogenBondInteraction()]
    else:
        return []

# TODO: cleanup + optimize this function
def has_hbond_interaction(res1, res2, switch = 1):
    # Initially focus on side-chain and water H-bonds
    if is_hbd(res1) and is_hba(res2):
        hbd = res1
        hba = res2

        # Get acceptors
        acceptors = get_hbond_acceptors(hba)

        # Get H-bond donor information
        donors = get_hbond_donor_references(hbd)
        for donor in donors:

            # Pairs within 3.5A?
            pairs = [ [donor, acceptor] for acceptor in acceptors if donor in hbd.child_dict and acceptor in hba.child_dict and distance_between(hbd.child_dict[donor].coord, hba.child_dict[acceptor].coord) <= 3.5 ]

            # Check angles
            if len(pairs) > 0:
                if is_water(hbd):
                    for pair in pairs:
                        return True
                else:
                    for pair in pairs:
                        donor = pair[0]
                        acceptor = pair[1]

                        for set in donors[donor]:
                            p1 = hbd.child_dict[set[0]].coord
                            p2 = hbd.child_dict[donor].coord

                            p3 = hba.child_dict[acceptor].coord
                            if len(set) == 4: # secondary
                                p3 = hbd.child_dict[set[3]].coord

                            # calculate optimal H-bonding vector to acceptor
                            d=get_unit_vector(p2-p1)
                            v=p3-p1
                            t=numpy.dot(v, d)
                            p4=p1+t*d
                            best_vector=get_unit_vector(p3-p4)
                            if len(set) == 4: # secondary
                                best_vector=-1 * best_vector

                            angle=math.radians(set[1]-90)
                            x=abs(math.cos(angle)*set[2])
                            y=abs(math.sin(angle)*set[2])
                            hydrogen=p2+y*d+x*best_vector

                            # check angle
                            if 180 - math.degrees(angle_between(hydrogen - p2, hba.child_dict[acceptor].coord - hydrogen)) >= 120:
                                return True

#                            print("POTENTIAL pair " + res1.get_resname() + str(res1.id[1]) + " - " + res2.get_resname() + str(res2.id[1]))
#                            print("ANGLE3 is " + str(180-math.degrees(angle_between(hydrogen-p2,hba.child_dict[acceptor].coord-hydrogen))))
#                            print("HETATM{:5d} C    HYD A{:4d}    {:8.3f}{:8.3f}{:8.3f}   1.00  1.00           C  ".format(1, 1, hydrogen[0], hydrogen[1], hydrogen[2]))

    if switch == 1:
        return has_hbond_interaction(res2, res1, 0)
    else:
        return False

def has_cation_pi_interaction(res1, res2):
    return has_pi_cation_interaction(res2, res1)


# Checks if two residues have an edge to face interaction
def has_face_to_edge_interaction(res1, res2):
    return has_edge_to_face_interaction(res2, res1)


def get_aromatic_interactions(res1, res2):
    interactions = []

    if is_aromatic_aa(res1) and is_aromatic_aa(res2):
        if has_face_to_face_interaction(res1, res2):
            interactions.append(FaceToFaceInteraction())

        if has_edge_to_face_interaction(res1, res2):
            interactions.append(EdgeToFaceInteraction())

        if has_face_to_edge_interaction(res1, res2):
            interactions.append(FaceToEdgeInteraction())

    if is_aromatic_aa(res1) and is_pos_charged(res2):
        if has_pi_cation_interaction(res1, res2):
            interactions.append(PiCationInteraction())

    if is_pos_charged(res1) and is_aromatic_aa(res2):
        if has_cation_pi_interaction(res1, res2):
            interactions.append(CationPiInteraction())

    return interactions


# Check if residues have any hydrophobic, i.e. C-C interactions
def get_hydrophobic_interactions(res1, res2):
    res1_carbons = [atom for atom in res1.child_list if atom.element == 'C']
    res2_carbons = [atom for atom in res2.child_list if atom.element == 'C']

    if any([any([distance_between(a1.coord, a2.coord) <= 4.5 for a1 in res1_carbons]) for a2 in res2_carbons]):
        return [HydrophobicInteraction()]
    else:
        return []


# Check if the sidechain atoms of res1 interacts with the backbone atoms of res2
def get_polar_sidechain_backbone_interactions(res1, res2):
    # Sidechain nitrogens and oxygens of residue 1
    res1_sidechain_nitrogens = [atom for atom in res1.child_list if atom.element == 'N' and atom.name != 'N']
    res1_sidechain_oxygens = [atom for atom in res1.child_list if atom.element == 'O' and atom.name != 'O']
    res1_sidechain_sulfurs = [atom for atom in res1.child_list if atom.element == 'S']

    # Backbone nitrogen and oxygen of residue 2
    try:
        res2_backbone_nitrogen = [res2.child_dict['N']]
        res2_backbone_oxygen = [res2.child_dict['O']]
    except KeyError:
        return []

    # Compute all polar interacting pairs
    sidechain_atoms = res1_sidechain_nitrogens + res1_sidechain_oxygens + res1_sidechain_sulfurs
    backbone_atoms = res2_backbone_nitrogen + res2_backbone_oxygen

    polarInteraction = any([distance_between(sca.coord, bba.coord) <= 4.5 for sca in sidechain_atoms for bba in backbone_atoms])

    if polarInteraction:
        return [PolarSideChainBackboneInteraction(is_charged(res1), is_charged(res2))]
    else:
        return []


# Check if the sidechain atoms of res2 interacts with the backbone atoms of res1
def get_polar_backbone_sidechain_interactions(res1, res2):
    interaction = get_polar_sidechain_backbone_interactions(res2, res1)

    if (interaction):
        return [PolarBackboneSidechainInteraction(is_charged(res1), is_charged(res2))]
    else:
        return []


# Check if the backbone atoms of res1 interacts with backbone atoms of res2
def get_polar_sidechain_sidechain_interactions(res1, res2):
    # Sidechain nitrogens and oxygens of residue 1
    res1_sidechain_nitrogens = [atom for atom in res1.child_list if atom.element == 'N' and atom.name != 'N']
    res1_sidechain_oxygens = [atom for atom in res1.child_list if atom.element == 'O' and atom.name != 'O']
    res1_sidechain_sulfurs = [atom for atom in res1.child_list if atom.element == 'S']

    # Sidechain nitrogens and oxygens of residue 2
    res2_sidechain_nitrogens = [atom for atom in res2.child_list if atom.element == 'N' and atom.name != 'N']
    res2_sidechain_oxygens = [atom for atom in res2.child_list if atom.element == 'O' and atom.name != 'O']
    res2_sidechain_sulfurs = [atom for atom in res2.child_list if atom.element == 'S']

    # Compute all polar interacting pairs
    res1_sidechain_atoms = res1_sidechain_nitrogens + res1_sidechain_oxygens + res1_sidechain_sulfurs
    res2_sidechain_atoms = res2_sidechain_nitrogens + res2_sidechain_oxygens + res2_sidechain_sulfurs

    polarInteraction = any([
       distance_between(sca1.coord, sca2.coord) <= 4.5 for sca1 in res1_sidechain_atoms for
       sca2 in res2_sidechain_atoms
    ])

    if polarInteraction:
        return [PolarSidechainSidechainInteraction(is_charged(res1), is_charged(res2))]
    else:
        return []

# TODO: cleanup
# Check if a water atom interacts with polar atoms of res2
def get_polar_water_interactions(res1, res2):
    if res1.resname == "HOH" or res2.resname == "HOH":
        # nitrogens and oxygens of residue 1
        res1_nitrogens = [atom for atom in res1.child_list if atom.element == 'N']
        res1_oxygens = [atom for atom in res1.child_list if atom.element == 'O']
        res1_sulfurs = [atom for atom in res1.child_list if atom.element == 'S']

        # nitrogens and oxygens of residue 2
        res2_nitrogens = [atom for atom in res2.child_list if atom.element == 'N']
        res2_oxygens = [atom for atom in res2.child_list if atom.element == 'O']
        res2_sulfurs = [atom for atom in res2.child_list if atom.element == 'S']

        # Compute all polar interacting pairs
        res1_atoms = res1_nitrogens + res1_oxygens + res1_sulfurs
        res2_atoms = res2_nitrogens + res2_oxygens + res2_sulfurs

        polarInteraction = any([
           distance_between(sca1.coord, sca2.coord) <= 4.5 for sca1 in res1_atoms for
           sca2 in res2_atoms
        ])

        if polarInteraction:
            return [PolarWaterInteraction()]
        else:
            return []
    else:
        return []

# Get polar contacts and interactions between 2 residues
def get_polar_interactions(res1, res2):
    polar_interactions = []
    polar_interactions += get_polar_backbone_sidechain_interactions(res1, res2)
    polar_interactions += get_polar_sidechain_backbone_interactions(res1, res2)
    polar_interactions += get_polar_sidechain_sidechain_interactions(res1, res2)
    polar_interactions += get_polar_water_interactions(res1, res2)
    polar_interactions += get_polar_hbonds_interactions(res1, res2)
    return polar_interactions


# Get Van der Waals inteactions between 2 residues
def get_van_der_waals_interactions(res1, res2):
    res1_atoms = res1.child_list
    res2_atoms = res2.child_list
    # Consider just labeling all unlabeled interacting pairs as VDW
    if any([distance_between(a1.coord, a2.coord) <= ((VDW_RADII[a1.element] + VDW_RADII[a2.element]) * VDW_TRESHOLD_FACTOR) for a1 in res1_atoms for a2 in res2_atoms]):
        return [VanDerWaalsInteraction()]
    else:
        return []


# Returns a list of interactions between two pairs of residues
def get_interactions(res1, res2):
    # Found interactions
    interactions = []

    # Aromatic interactions
    interactions += get_aromatic_interactions(res1, res2)

    # Hydrophobic interactions
    interactions += get_hydrophobic_interactions(res1, res2)

    # Polar interactions
    interactions += get_polar_interactions(res1, res2)

    # TODO: extend with C or N-terms and backbones?
    if (is_hba(res1) and is_hbd(res2)) or (is_hbd(res1) and is_hba(res2)):
        has_hbond_interaction(res1, res2)

    # Van der Waals interactions
    interactions += get_van_der_waals_interactions(res1, res2)

    return interactions
