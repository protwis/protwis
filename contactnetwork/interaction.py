from contactnetwork.residue import *
from Bio.PDB.Polypeptide import *
from contactnetwork.models import *

from residue.models import Residue

import math

class InteractingPair:

    'Common base class for all interactions'
    def __init__(self, res1, res2, interactions,dbres1,dbres2,s):
        self.res1 = res1
        self.res2 = res2
        self.interactions = interactions
        self.dbres1 = dbres1
        self.dbres2 = dbres2
        self.structure = s

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

    def get_interaction_json(self, generic, pdb):
        # Temporary mapping G-proteins using static variable mapping
        # G-prot numbering_schemes
        gprot = self.mapping[pdb][self.res2.id[1]];
        selected = Residue.objects.filter(sequence_number=gprot, protein_conformation__protein__entry_name="gnas2_human")
        gprot_gn = selected[0].display_generic_number
        text = '[\'{0}\',\'{1}\',{2},\'{3}\',\'{4}\',\'{5}\',{6},\'{7}\', ["'.format(self.res1.parent.id, three_to_one(self.res1.get_resname()), self.res1.id[1], generic, self.res2.parent.id, three_to_one(self.res2.get_resname()), self.res2.id[1], gprot_gn)
        first = True

        if self.interactions:
            for i in self.interactions:
                if not first:
                    text += ', "'
                first = False
                text += i.get_name() + '"'

        text += ']]'

        return text

    def save_into_database(self):
        # Save the pair
        pair,created = InteractingResiduePair.objects.get_or_create(res1=self.dbres1, res2=self.dbres2, referenced_structure=self.structure)

        # Add the interactions to the pair'
        bulk = []
        for i in self.get_interactions():
            # if type(i) is VanDerWaalsInteraction:
            #     ni = Interaction(interaction_type='VanDerWaals', interacting_pair=pair)
            # elif type(i) is HydrophobicInteraction:
            #     ni = Interaction(interaction_type='Hydrophobic', interacting_pair=pair)
            # elif type(i) is PolarSidechainSidechainInteraction:
            #     ni = Interaction(interaction_type='Polar',specific_type='PolarSidechainSidechain', interacting_pair=pair)
            # elif type(i) is PolarBackboneSidechainInteraction:
            #     ni = Interaction(interaction_type='Polar',specific_type='PolarBackboneSidechain', interacting_pair=pair)
            # elif type(i) is PolarSideChainBackboneInteraction:
            #     ni = Interaction(interaction_type='Polar',specific_type='PolarSideChainBackbone', interacting_pair=pair)
            # elif type(i) is PosNegIonicInteraction:
            #     ni = Interaction(interaction_type='Polar',specific_type='IonicPositiveNegative', interacting_pair=pair)
            # elif type(i) is NegPosIonicInteraction:
            #     ni = Interaction(interaction_type='Polar',specific_type='IonicNegativePositive', interacting_pair=pair)
            # elif type(i) is HydrogenBondDAInteraction:
            #     ni = Interaction(interaction_type='Polar',specific_type='HBondDonorAcceptor', interacting_pair=pair)
            # elif type(i) is HydrogenBondADInteraction:
            #     ni = Interaction(interaction_type='Polar',specific_type='HBondAcceptorDonor', interacting_pair=pair)
            # elif type(i) is WaterMediated:
            #     ni = Interaction(interaction_type='Polar',specific_type='HBondWaterMediated', interacting_pair=pair)
            # elif type(i) is FaceToFaceInteraction:
            #     ni = Interaction(interaction_type='Aromatic',specific_type='FaceToFace', interacting_pair=pair)
            # elif type(i) is FaceToEdgeInteraction:
            #     ni = Interaction(interaction_type='Aromatic',specific_type='FaceToEdge', interacting_pair=pair)
            # elif type(i) is EdgeToFaceInteraction:
            #     ni = Interaction(interaction_type='Aromatic',specific_type='EdgeToFace', interacting_pair=pair)
            # elif type(i) is PiCationInteraction:
            #     ni = Interaction(interaction_type='Aromatic',specific_type='PiCation', interacting_pair=pair)
            # elif type(i) is CationPiInteraction:
            #     ni = Interaction(interaction_type='Aromatic',specific_type='CationPi', interacting_pair=pair)
            #
            # else:
            #     ni = Interaction(interaction_type=i.get_name() , interacting_pair=pair)
                #ni.res1_has_pi = False

            ni = Interaction(interaction_type=i.get_type(),specific_type=i.get_details(), interacting_pair=pair)
            bulk.append(ni)
        Interaction.objects.bulk_create(bulk)


# Make type and detail variables and default functions
class CI(object):
    def __init__(self):
        self.type = "interaction"
        self.detail = ""
        pass

    def get_type(self):
        return self.type

    def get_details(self):
        return self.detail

class VanDerWaalsInteraction(CI):
    def __init__(self):
        self.type = "van-der-waals"
        self.detail = ""
        pass

class HydrophobicInteraction(CI):
    def __init__(self):
        self.type = "hydrophobic"
        self.detail = ""
        pass

class IonicInteraction(CI):
    def __init__(self):
        super().__init__()
        self.type = "ionic"
        self.detail = ""

class PolarInteraction(CI):
    def __init__(self):
        self.type = "polar"
        self.detail = ""
        pass

class AromaticInteraction(CI):
    def __init__(self):
        self.type = "aromatic"
        self.detail = ""
        pass

class NegPosIonicInteraction(IonicInteraction):
    def __init__(self):
        super().__init__()
        self.detail = "negative-positive"

class PosNegIonicInteraction(IonicInteraction):
    def __init__(self):
        super().__init__()
        self.detail = "positive-negative"

class HydrogenBondInteraction(PolarInteraction):
    def __init__(self):
        super().__init__()
        self.detail = "H-bond"

class HydrogenBondDAInteraction(PolarInteraction):
    def __init__(self):
        super().__init__()
        self.detail = "donor-acceptor"

class HydrogenBondADInteraction(PolarInteraction):
    def __init__(self):
        super().__init__()
        self.detail = "acceptor-donor"

class WaterMediated(PolarInteraction):
    def __init__(self):
        super().__init__()
        self.detail = "water-mediated"

class PolarSidechainSidechainInteraction(PolarInteraction):
    def __init__(self):
        super().__init__()
        self.detail = "sidechain-sidechain"

class PolarBackboneSidechainInteraction(PolarInteraction):
    def __init__(self):
        super().__init__()
        self.detail = "backbone-sidechain"

class PolarSideChainBackboneInteraction(PolarInteraction):
    def __init__(self):
        super().__init__()
        self.detail = "sidechain-backbone"

# AK: I vote to REMOVE
class PolarWaterInteraction(PolarInteraction):
    def __init__(self):
        super().__init__()
        self.detail = "water"

class FaceToFaceInteraction(AromaticInteraction):
    def __init__(self):
        super().__init__()
        self.detail = "face-to-face"

class EdgeToFaceInteraction(AromaticInteraction):
    def __init__(self):
        super().__init__()
        self.detail = "edge-to-face"

class FaceToEdgeInteraction(AromaticInteraction):
    def __init__(self):
        super().__init__()
        self.detail = "edge-to-edge"

class PiCationInteraction(AromaticInteraction):
    def __init__(self):
        super().__init__()
        self.detail = "pi-cation"

class CationPiInteraction(AromaticInteraction):
    def __init__(self):
        super().__init__()
        self.detail = "cation-pi"

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
    res1_desc = get_ring_descriptors(res1)
    res2_desc = get_ring_descriptors(res2)

    # Make sure that the acute angle between the planes are less than (or eq.) 30
    # degrees and that the distance between centers is less than  (or eq.) 4.4 Angstrom.
    return any([(distance_between(r1[0], r2[0]) <= 4.4)
                and (math.degrees(angle_between_plane_normals(r1[1], r2[1])) <= 30)
                for r1 in res1_desc for r2 in res2_desc])


# Checks if two residues have an edge to face interaction
def has_edge_to_face_interaction(res1, res2):
    res1_desc = get_ring_descriptors(res1)
    res2_desc = get_ring_descriptors(res2)

    # Make sure the ring centers are closer than 5.5 angstroms
    # and that the perpendicular angle is within +/- 30 degrees
    # i.e. the acute angle is between 60-120 degrees
    return any([(distance_between(r1[0], r2[0]) <= 5.5)
                and (abs(math.degrees(angle_between_plane_normals(r1[1], r2[1])) - 90) <= 30)
                for r1 in res1_desc for r2 in res2_desc])


def has_pi_cation_interaction(res1, res2):
    # Only aromatic residues will have ring descriptors
    res1_desc = get_ring_descriptors(res1)

    # Only positively charged atoms will have charged atoms
    res2_pos_atom_names = get_pos_charged_atom_names(res2)

    # Check if any charged atom is within 6.6 angstroms of any ring centroid
    # and that the angle is within +/- 30 degrees
    return any([distance_between(res2.child_dict[atom_name].coord, r1[0]) <= 6.6
        and abs(math.degrees(angle_between_plane_normals(r1[1], res2.child_dict[atom_name].coord-r1[0]))) <= 30
        for r1 in res1_desc for atom_name in res2_pos_atom_names])

def get_polar_hbonds_interactions(res1, res2):
    # TODO: extend with loose H-bond definitions as new polar contacts
    if has_hbond_interaction(res1, res2):
        return [HydrogenBondDAInteraction()]
    elif has_hbond_interaction(res2, res1):
        return [HydrogenBondADInteraction()]
    else:
        return []

# TODO: cleanup + optimize this function
def has_hbond_interaction(res1, res2):
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
                            if set[0] in hbd.child_dict and (len(set) < 4 or set[3] in hbd.child_dict):
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

    return False

def has_cation_pi_interaction(res1, res2):
    return has_pi_cation_interaction(res2, res1)


# Checks if two residues have an edge to face interaction
def has_face_to_edge_interaction(res1, res2):
    return has_edge_to_face_interaction(res2, res1)


def get_aromatic_interactions(res1, res2):
    interactions = []

    if is_aromatic_aa(res1) and is_pos_charged(res2):
        if has_pi_cation_interaction(res1, res2):
            interactions.append(PiCationInteraction())

    if is_pos_charged(res1) and is_aromatic_aa(res2):
        if has_cation_pi_interaction(res1, res2):
            interactions.append(CationPiInteraction())

    if is_aromatic_aa(res1) and is_aromatic_aa(res2):
        if has_face_to_face_interaction(res1, res2):
            interactions.append(FaceToFaceInteraction())

        elif has_edge_to_face_interaction(res1, res2):
            interactions.append(EdgeToFaceInteraction())

        # AK 28-11-2018 - now redundant given the calculation method
        #if has_face_to_edge_interaction(res1, res2):
        #    interactions.append(FaceToEdgeInteraction())

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
        return [PolarSideChainBackboneInteraction()]
    else:
        return []


# Check if the sidechain atoms of res2 interacts with the backbone atoms of res1
def get_polar_backbone_sidechain_interactions(res1, res2):
    interaction = get_polar_sidechain_backbone_interactions(res2, res1)

    if (interaction):
        return [PolarBackboneSidechainInteraction()]
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
        return [PolarSidechainSidechainInteraction()]
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
    polar_interactions += get_polar_hbonds_interactions(res1, res2)
    # AK: as discussed disabled for now
    #polar_interactions += get_polar_water_interactions(res1, res2)
    #polar_interactions += get_polar_sidechain_sidechain_interactions(res1, res2)
    #polar_interactions += get_polar_sidechain_backbone_interactions(res1, res2)
    #polar_interactions += get_polar_backbone_sidechain_interactions(res1, res2)

    return polar_interactions


def get_ionic_interactions(res1, res2):
    # Only oppositely charged residues
    if (is_pos_charged(res1) and is_pos_charged(res2)) or (is_neg_charged(res1) and is_neg_charged(res2)):
        return []

    res1_atom_names = get_charged_atom_names(res1)

    res2_atom_names = get_charged_atom_names(res2)

    # Check if any charged atoms are within 4.5 angstroms
    if any([distance_between(res1.child_dict[name1].coord, res2.child_dict[name2].coord) <= 4.5 for name1 in res1_atom_names for name2 in res2_atom_names]):
        if is_pos_charged(res1):
            return [NegPosIonicInteraction()]
        else:
            return [PosNegIonicInteraction()]
    else:
        return []

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

    # TODO: extend with C or N-terms and backbones
    #if (is_hba(res1) and is_hbd(res2)) or (is_hbd(res1) and is_hba(res2)):
    #    has_hbond_interaction(res1, res2)

    # Ionic interactions
    if (is_charged(res1) and is_charged(res2)):
        interactions += get_ionic_interactions(res1, res2)

    # Polar interactions
    interactions += get_polar_interactions(res1, res2)

    # Aromatic interactions
    interactions += get_aromatic_interactions(res1, res2)

    # Van der Waals interactions
    interactions += get_van_der_waals_interactions(res1, res2)

    # Hydrophobic interactions
    interactions += get_hydrophobic_interactions(res1, res2)

    return interactions
