from contactnetwork.residue import *
from Bio.PDB.Polypeptide import *
from contactnetwork.models import *

from residue.models import Residue

import math

class InteractingPair:

    # Distance for encoding solely side-chain interactions
    NUM_SKIP_BB_INTERACTIONS = 4

    'Common base class for all interactions'
    def __init__(self, res1, res2, dbres1, dbres2, structure):
        self.res1 = res1
        self.res2 = res2
        self.dbres1 = dbres1
        self.dbres2 = dbres2
        self.structure = structure
        self.interactions = []
        self.compute_interactions()

    def __repr__(self):
        return '<InteractingPair: {}-{}-{}>'.format(self.dbres1, self.dbres2, self.structure)

    def add_interactions(self, interaction):
        self.interactions.append(interaction)

    # Creates a list of interactions between two pairs of residues
    def compute_interactions(self):
        # Ionic interactions
        self.ionic_interactions()

        # Polar interactions
        self.hbond_interactions()

        # Aromatic interactions
        self.aromatic_interactions()

        # Hydrophobic interactions
        self.hydrophobic_interactions()

        # Van der Waals interactions
        self.van_der_waals_interactions()

        # DISABLED for now - is a toggle in the browser
        # Filter interactions between backbone atoms if residue distance within cut-off set in cube.py
        # backbone_atoms = ["C", "O", "N", "CA"]
        # if self.interactions and abs(self.res1.id[1] - self.res2.id[1]) <= self.NUM_SKIP_BB_INTERACTIONS:
        #     for i in self.interactions:
        #         backbone_res1 = i.atomname_residue1 in backbone_atoms
        #         backbone_res2 = i.atomname_residue2 in backbone_atoms
        #
        #         # Remove interaction if both are backbone atoms
        #         if backbone_res1 and backbone_res2:
        #             #print("REMOVING this interaction", self.res1.id[1], self.res2.id[1], i.atomname_residue1, i.atomname_residue2)
        #             self.interactions.remove(i)


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
            ni = Interaction(interaction_type=i.get_type(),specific_type=i.get_details(), interacting_pair=pair, atomname_residue1=i.atomname_residue1, atomname_residue2=i.atomname_residue2, interaction_level=i.get_level())
            bulk.append(ni)
        Interaction.objects.bulk_create(bulk)

    def save_peptide_interactions(self, peptide):
        pair, _ = InteractingPeptideResiduePair.objects.get_or_create(peptide_amino_acid_three_letter=self.dbres2.three_letter, peptide_amino_acid=self.dbres2.amino_acid, peptide_sequence_number=self.dbres2.sequence_number,
                                                                            peptide=peptide, receptor_residue=self.dbres1)
        bulk = []
        for i in self.get_interactions():
            ip = InteractionPeptide(interacting_peptide_pair=pair, peptide_atom=i.atomname_residue2, receptor_atom=i.atomname_residue1, interaction_type=i.get_type(), specific_type=i.get_details(), interaction_level=i.get_level())
            bulk.append(ip)
        InteractionPeptide.objects.bulk_create(bulk)

    def ionic_interactions(self):
        # Only oppositely charged residues
        if (is_pos_charged(self.res1) and is_neg_charged(self.res2)) or (is_neg_charged(self.res1) and is_pos_charged(self.res2)):

            res1_atom_names = get_charged_atom_names(self.res1)
            res2_atom_names = get_charged_atom_names(self.res2)

            # Check if any charged atoms are within 4.5 angstroms
            for match1, match2 in [[name1, name2] for name1 in res1_atom_names for name2 in res2_atom_names if distance_between(self.res1.child_dict[name1].coord, self.res2.child_dict[name2].coord) <= 4.5 ]:
                if is_pos_charged(self.res1):
                    self.add_interactions(PosNegIonicInteraction(match1, match2))
                else:
                    self.add_interactions(NegPosIonicInteraction(match1, match2))

    def strict_hbond_interactions(self, switch = False):
        hbd = self.res1
        hba = self.res2
        if switch:
            hbd = self.res2
            hba = self.res1

        found = False

        # Get acceptors
        acceptors = get_hbond_acceptors(hba)

        # Get H-bond donor information
        donors = get_hbond_donor_references(hbd)
        for donor in donors:
            # Pairs within 3.5A?
            pairs = [ [donor, acceptor] for acceptor in acceptors if donor in hbd.child_dict and acceptor in hba.child_dict and distance_between(hbd.child_dict[donor].coord, hba.child_dict[acceptor].coord) <= 3.5 ]

            # Check angles
            if len(pairs) > 0:
                #if is_water(hbd):
                #    for pair in pairs:
                #        return True
                #else:
                for pair in pairs:
                    donor = pair[0]
                    acceptor = pair[1]

                    if InteractingPair.verify_hbond_angle(hbd, donor, hba, acceptor):
                        if switch:
                            self.add_interactions(HydrogenBondADInteraction(acceptor, donor))
                            found = True
                        else:
                            self.add_interactions(HydrogenBondDAInteraction(donor, acceptor))
                            found = True
        if not switch:
            found = self.strict_hbond_interactions(switch = True) or found
        return found

    @staticmethod
    def verify_water_hbond(residue, water):
        atom_names = []

        # Water as donor
        if is_hba(residue):
            acceptors = get_hbond_acceptors(residue)
            matches = [ acceptor for acceptor in acceptors if acceptor in residue.child_dict and distance_between(water.coord, residue.child_dict[acceptor].coord) <= 3.5 ]
            if len(matches) > 0:
                atom_names = matches

        # Water as acceptor
        if is_hbd(residue):
            donors = get_hbond_donor_references(residue)
            matches = [ donor for donor in donors if donor in residue.child_dict and distance_between(water.coord, residue.child_dict[donor].coord) <= 3.5 ]

            for donor in matches:
                if InteractingPair.verify_hbond_angle(residue, donor, water.get_parent(), water.get_name()):
                    atom_names.append(donor)

        # return unique atom names (e.g. when atom is both acceptor + donor)
        return list(set(atom_names))

    @staticmethod
    def verify_hbond_angle(hbd_residue, hbd_atomname, hba_residue, hba_atomname):
        donors = get_hbond_donor_references(hbd_residue)
        for donor_set in donors[hbd_atomname]:
            if donor_set[0] in hbd_residue.child_dict and (len(donor_set) < 4 or donor_set[3] in hbd_residue.child_dict):
                p1 = hbd_residue.child_dict[donor_set[0]].coord
                p2 = hbd_residue.child_dict[hbd_atomname].coord

                # Option 1: freely rotating donor, take acceptor as reference
                p3 = hba_residue[hba_atomname].coord

                # Option 2: fixed donor orientation (dihedral), internal reference
                if len(donor_set) == 4: # secondary
                    # backbone nitrogen - second reference atom comes from the previous residue
                    if hbd_atomname == 'N':
                        # 1. grab previous connected residue (if numbering switches, this might give an error)
                        # 2. grab coordinate atom previous residue
                        hbd_chain = hbd_residue.get_parent()
                        if hbd_chain.has_id(hbd_residue.id[1]-1):
                            p3 = hbd_chain[hbd_residue.id[1]-1].child_dict[donor_set[3]].coord
                        else:
                            continue
                    else:
                        p3 = hbd_residue.child_dict[donor_set[3]].coord

                # calculate optimal H-bonding vector to acceptor
                d=get_unit_vector(p2-p1)
                v=p3-p1
                t=numpy.dot(v, d)
                p4=p1+t*d
                best_vector=get_unit_vector(p3-p4)
                if len(donor_set) == 4: # secondary
                    best_vector=-1 * best_vector

                angle=math.radians(donor_set[1]-90)
                x=abs(math.cos(angle)*donor_set[2])
                y=abs(math.sin(angle)*donor_set[2])
                hydrogen=p2+y*d+x*best_vector

                # check angle
                if 180 - math.degrees(angle_between(hydrogen - p2, hba_residue[hba_atomname].coord - hydrogen)) >= 120:
                    return True

        return False


    def loose_hbond_interactions(self, switch = False):
        hbd = self.res1
        hba = self.res2
        if switch:
            hbd = self.res2
            hba = self.res1

        # Get acceptors
        acceptors = get_hbond_acceptors(hba)

        # Get H-bond donor information
        donors = get_hbond_donor_references(hbd)
        for donor in donors:
            # Matching donor-acceptor pairs within 4A - no angle requirements
            pairs = [ [donor, acceptor] for acceptor in acceptors if donor in hbd.child_dict and acceptor in hba.child_dict and distance_between(hbd.child_dict[donor].coord, hba.child_dict[acceptor].coord) <= 4 ]
            for pair in pairs:
                donor = pair[0]
                acceptor = pair[1]

                if switch:
                    self.add_interactions(LooseHydrogenBondADInteraction(acceptor, donor))
                else:
                    self.add_interactions(LooseHydrogenBondDAInteraction(donor, acceptor))

        if not switch:
            self.loose_hbond_interactions(switch = True)

    def hbond_interactions(self):
        # calculate interactions for pair
        found = self.strict_hbond_interactions()
        if not found:
            self.loose_hbond_interactions()


    def aromatic_interactions(self):
        aromatic_count = is_aromatic_aa(self.res1) + is_aromatic_aa(self.res2)

        if aromatic_count == 1:
            self.pi_cation_interactions()
        elif aromatic_count == 2:
            found = self.face_to_face_interactions()
            found = self.edge_to_face_interactions() or found
            if not found:
                self.loose_aromatic_interactions()

    def pi_cation_interactions(self, switch = False):
        aromatic = self.res1
        cation = self.res2
        if switch:
            aromatic = self.res2
            cation = self.res1

        if is_aromatic_aa(aromatic) and is_pos_charged(cation):
            res1_desc = get_ring_descriptors(aromatic)
            res2_pos_atom_names = get_pos_charged_atom_names(cation)

            # Check if any charged atom is within 6.6 angstroms of any ring centroid
            # and that the angle is within +/- 30 degrees
            for match1, match2 in [[index1, atom_name] for index1, r1 in enumerate(res1_desc)
                        for atom_name in res2_pos_atom_names
                        if distance_between(cation.child_dict[atom_name].coord, r1[0]) <= 6.6
                        and abs(math.degrees(angle_between_plane_normals(r1[1], cation.child_dict[atom_name].coord-r1[0]))) <= 30]:
                if switch:
                    self.add_interactions(CationPiInteraction("RN"+str(match1+1), match2))
                else:
                    self.add_interactions(PiCationInteraction(match2, "RN"+str(match1+1)))

        if not switch:
            self.pi_cation_interactions(switch = True)

    # Checks if two residues have a face to face interaction
    def face_to_face_interactions(self):
        res1_desc = get_ring_descriptors(self.res1)
        res2_desc = get_ring_descriptors(self.res2)

        # Make sure that the acute angle between the planes are less than (or eq.) 30
        # degrees and that the distance between centers is less than  (or eq.) 4.4 Angstrom.
        for match1, match2 in [[index1, index2] for index1, r1 in enumerate(res1_desc)
                for index2, r2 in enumerate(res2_desc)
                    if (distance_between(r1[0], r2[0]) <= 4.4)
                    and (math.degrees(angle_between_plane_normals(r1[1], r2[0])) <= 30)]:
            self.add_interactions(FaceToFaceInteraction("RN"+str(match1 + 1), "RN"+str(match2 + 1)))

    # Checks if two residues have an edge to face interaction
    def edge_to_face_interactions(self, switch = False):
        res1_desc = get_ring_descriptors(self.res1)
        res2_desc = get_ring_descriptors(self.res2)
        if switch:
            res1_desc = res2_desc
            res2_desc = get_ring_descriptors(self.res1)


        # Make sure the ring centers are closer than 5.5 angstroms
        # and that the perpendicular angle is within +/- 30 degrees
        for match1, match2 in [[index1, index2] for index1, r1 in enumerate(res1_desc)
                for index2, r2 in enumerate(res2_desc)
                    if (distance_between(r1[0], r2[0]) <= 5.5)
                    and (math.degrees(angle_between_plane_normals(r1[1], r2[1])) > 30)]:

            # angle = math.degrees(angle_between_plane_normals(res1_desc[match1][1], res2_desc[match2][1]))
            angle = abs(90 - math.degrees(angle_between_plane_normals(res1_desc[match1][1], res2_desc[match2][0]-res1_desc[match1][0])))
            if angle <= 30:
                if switch:
                    self.add_interactions(FaceToEdgeInteraction("RN"+str(match1 + 1), "RN"+str(match2 + 1)))
                else:
                    self.add_interactions(EdgeToFaceInteraction("RN"+str(match1 + 1), "RN"+str(match2 + 1)))

        if not switch:
            self.edge_to_face_interactions(switch = True)

    def loose_aromatic_interactions(self):
        res1_desc = get_ring_descriptors(self.res1)
        res2_desc = get_ring_descriptors(self.res2)

        # Make sure the ring centers are closer than 5.5 angstroms, not additional requirements
        for match1, match2 in [[index1, index2] for index1, r1 in enumerate(res1_desc)
                for index2, r2 in enumerate(res2_desc)
                    if (distance_between(r1[0], r2[0]) <= 5.5)]:
                        self.add_interactions(LooseAromaticInteraction("RN"+str(match1 + 1), "RN"+str(match2 + 1)))

    def hydrophobic_interactions(self):
        res1_carbons = [atom for atom in self.res1.child_list if atom.element == 'C' or atom.element == 'S']
        res2_carbons = [atom for atom in self.res2.child_list if atom.element == 'C' or atom.element == 'S']

        for match1, match2 in [[a1, a2] for a1 in res1_carbons for a2 in res2_carbons if distance_between(a1.coord, a2.coord) <= 4.5]:
            self.add_interactions(HydrophobicInteraction(match1.name, match2.name))

    # Get Van der Waals inteactions between 2 residues
    def van_der_waals_interactions(self):
        res1_atoms = self.res1.child_list
        res2_atoms = self.res2.child_list

        for match1, match2 in [[a1, a2] for a1 in res1_atoms for a2 in res2_atoms if distance_between(a1.coord, a2.coord) <= ((VDW_RADII[a1.element] + VDW_RADII[a2.element]) * VDW_TRESHOLD_FACTOR)]:
            self.add_interactions(VanDerWaalsInteraction(match1.name, match2.name))

# Make type and detail variables and default functions
class CI(object):
    def __init__(self, name1, name2):
        self.type = "interaction"
        self.detail = ""
        self.atomname_residue1 = name1
        self.atomname_residue2 = name2
        self.interaction_level = 0

    def get_type(self):
        return self.type

    def get_details(self):
        return self.detail

    def get_level(self):
        return self.interaction_level

    def set_atomname_residue1(self, name):
        self.atomname_residue1 = name

    def set_atomname_residue2(self, name):
        self.atomname_residue2 = name

# PRIMARY TYPES
class VanDerWaalsInteraction(CI):
    def __init__(self, name1, name2):
        self.type = "van-der-waals"
        self.detail = ""
        self.atomname_residue1 = name1
        self.atomname_residue2 = name2
        self.interaction_level = 0

class HydrophobicInteraction(CI):
    def __init__(self, name1, name2):
        self.type = "hydrophobic"
        self.detail = ""
        self.atomname_residue1 = name1
        self.atomname_residue2 = name2
        self.interaction_level = 0

class IonicInteraction(CI):
    def __init__(self, name1, name2):
        self.type = "ionic"
        self.detail = ""
        self.atomname_residue1 = name1
        self.atomname_residue2 = name2
        self.interaction_level = 0

class PolarInteraction(CI):
    def __init__(self, name1, name2):
        self.type = "polar"
        self.detail = ""
        self.atomname_residue1 = name1
        self.atomname_residue2 = name2
        self.interaction_level = 0

class AromaticInteraction(CI):
    def __init__(self, name1, name2):
        self.type = "aromatic"
        self.detail = ""
        self.atomname_residue1 = name1
        self.atomname_residue2 = name2
        self.interaction_level = 0

# SECONDARY TYPES
class NegPosIonicInteraction(IonicInteraction):
    def __init__(self, name1, name2):
        super().__init__(name1, name2)
        self.detail = "negative-positive"

class PosNegIonicInteraction(IonicInteraction):
    def __init__(self, name1, name2):
        super().__init__(name1, name2)
        self.detail = "positive-negative"

class HydrogenBondDAInteraction(PolarInteraction):
    def __init__(self, name1, name2):
        super().__init__(name1, name2)
        self.detail = "donor-acceptor"

class HydrogenBondADInteraction(PolarInteraction):
    def __init__(self, name1, name2):
        super().__init__(name1, name2)
        self.detail = "acceptor-donor"

class WaterMediated(PolarInteraction):
    def __init__(self, name1, name2):
        super().__init__(name1, name2)
        self.detail = "water-mediated"

class LooseHydrogenBondDAInteraction(PolarInteraction):
    def __init__(self, name1, name2):
        super().__init__(name1, name2)
        self.detail = "donor-acceptor"
        self.interaction_level = 1

class LooseHydrogenBondADInteraction(PolarInteraction):
    def __init__(self, name1, name2):
        super().__init__(name1, name2)
        self.detail = "acceptor-donor"
        self.interaction_level = 1

class FaceToFaceInteraction(AromaticInteraction):
    def __init__(self, name1, name2):
        super().__init__(name1, name2)
        self.detail = "face-to-face"

class EdgeToFaceInteraction(AromaticInteraction):
    def __init__(self, name1, name2):
        super().__init__(name1, name2)
        self.detail = "edge-to-face"

class FaceToEdgeInteraction(AromaticInteraction):
    def __init__(self, name1, name2):
        super().__init__(name1, name2)
        self.detail = "face-to-edge"

class LooseAromaticInteraction(AromaticInteraction):
    def __init__(self, name1, name2):
        super().__init__(name1, name2)
        self.detail = "loose-aromatic"
        self.interaction_level = 1

class PiCationInteraction(AromaticInteraction):
    def __init__(self, name1, name2):
        super().__init__(name1, name2)
        self.detail = "pi-cation"

class CationPiInteraction(AromaticInteraction):
    def __init__(self, name1, name2):
        super().__init__(name1, name2)
        self.detail = "cation-pi"

# class PolarSidechainSidechainInteraction(PolarInteraction):
#     def __init__(self):
#         super().__init__()
#         self.detail = "sidechain-sidechain"
#
# class PolarBackboneSidechainInteraction(PolarInteraction):
#     def __init__(self):
#         super().__init__()
#         self.detail = "backbone-sidechain"
#
# class PolarSideChainBackboneInteraction(PolarInteraction):
#     def __init__(self):
#         super().__init__()
#         self.detail = "sidechain-backbone"
#
# class PolarWaterInteraction(PolarInteraction):
#     def __init__(self):
#         super().__init__()
#         self.detail = "water"


# Return the unit vector for a given vector
def unit_vector(vector):
    return vector / numpy.linalg.norm(vector)


# Return the angle between two given vectors
def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return numpy.arccos(numpy.clip(numpy.dot(v1_u, v2_u), -1.0, 1.0))


# Return the acute angles between two given vectors
def angle_between_plane_normals(v1, v2):
    return min([angle_between(v1, v2), angle_between(v1, numpy.multiply(v2, -1.0))])


# Return the distance between two given points
def distance_between(v1, v2):
    return numpy.linalg.norm(numpy.subtract(v1, v2))

### DEPRECATED ###
# def get_polar_hbonds_interactions(res1, res2):
#     # TODO: extend with loose H-bond definitions as new polar contacts
#     if has_hbond_interactions(res1, res2):
#         return [HydrogenBondDAInteraction()]
#     elif has_hbond_interactions(res2, res1):
#         return [HydrogenBondADInteraction()]
#     else:
#         return []
#
# # Check if the sidechain atoms of res1 interacts with the backbone atoms of res2
# def get_polar_sidechain_backbone_interactions(res1, res2):
#     # Sidechain nitrogens and oxygens of residue 1
#     res1_sidechain_nitrogens = [atom for atom in res1.child_list if atom.element == 'N' and atom.name != 'N']
#     res1_sidechain_oxygens = [atom for atom in res1.child_list if atom.element == 'O' and atom.name != 'O']
#     res1_sidechain_sulfurs = [atom for atom in res1.child_list if atom.element == 'S']
#
#     # Backbone nitrogen and oxygen of residue 2
#     try:
#         res2_backbone_nitrogen = [res2.child_dict['N']]
#         res2_backbone_oxygen = [res2.child_dict['O']]
#     except KeyError:
#         return []
#
#     # Compute all polar interacting pairs
#     sidechain_atoms = res1_sidechain_nitrogens + res1_sidechain_oxygens + res1_sidechain_sulfurs
#     backbone_atoms = res2_backbone_nitrogen + res2_backbone_oxygen
#
#     polarInteraction = any([distance_between(sca.coord, bba.coord) <= 4.5 for sca in sidechain_atoms for bba in backbone_atoms])
#
#     if polarInteraction:
#         return [PolarSideChainBackboneInteraction()]
#     else:
#         return []
#
#
# # Check if the sidechain atoms of res2 interacts with the backbone atoms of res1
# def get_polar_backbone_sidechain_interactions(res1, res2):
#     interaction = get_polar_sidechain_backbone_interactions(res2, res1)
#
#     if (interaction):
#         return [PolarBackboneSidechainInteraction()]
#     else:
#         return []
#
#
# # Check if the backbone atoms of res1 interacts with backbone atoms of res2
# def get_polar_sidechain_sidechain_interactions(res1, res2):
#     # Sidechain nitrogens and oxygens of residue 1
#     res1_sidechain_nitrogens = [atom for atom in res1.child_list if atom.element == 'N' and atom.name != 'N']
#     res1_sidechain_oxygens = [atom for atom in res1.child_list if atom.element == 'O' and atom.name != 'O']
#     res1_sidechain_sulfurs = [atom for atom in res1.child_list if atom.element == 'S']
#
#     # Sidechain nitrogens and oxygens of residue 2
#     res2_sidechain_nitrogens = [atom for atom in res2.child_list if atom.element == 'N' and atom.name != 'N']
#     res2_sidechain_oxygens = [atom for atom in res2.child_list if atom.element == 'O' and atom.name != 'O']
#     res2_sidechain_sulfurs = [atom for atom in res2.child_list if atom.element == 'S']
#
#     # Compute all polar interacting pairs
#     res1_sidechain_atoms = res1_sidechain_nitrogens + res1_sidechain_oxygens + res1_sidechain_sulfurs
#     res2_sidechain_atoms = res2_sidechain_nitrogens + res2_sidechain_oxygens + res2_sidechain_sulfurs
#
#     polarInteraction = any([
#        distance_between(sca1.coord, sca2.coord) <= 4.5 for sca1 in res1_sidechain_atoms for
#        sca2 in res2_sidechain_atoms
#     ])
#
#     if polarInteraction:
#         return [PolarSidechainSidechainInteraction()]
#     else:
#         return []
#
# # Check if a water atom interacts with polar atoms of res2
# def get_polar_water_interactions(res1, res2):
#     if res1.resname == "HOH" or res2.resname == "HOH":
#         # nitrogens and oxygens of residue 1
#         res1_nitrogens = [atom for atom in res1.child_list if atom.element == 'N']
#         res1_oxygens = [atom for atom in res1.child_list if atom.element == 'O']
#         res1_sulfurs = [atom for atom in res1.child_list if atom.element == 'S']
#
#         # nitrogens and oxygens of residue 2
#         res2_nitrogens = [atom for atom in res2.child_list if atom.element == 'N']
#         res2_oxygens = [atom for atom in res2.child_list if atom.element == 'O']
#         res2_sulfurs = [atom for atom in res2.child_list if atom.element == 'S']
#
#         # Compute all polar interacting pairs
#         res1_atoms = res1_nitrogens + res1_oxygens + res1_sulfurs
#         res2_atoms = res2_nitrogens + res2_oxygens + res2_sulfurs
#
#         polarInteraction = any([
#            distance_between(sca1.coord, sca2.coord) <= 4.5 for sca1 in res1_atoms for
#            sca2 in res2_atoms
#         ])
#
#         if polarInteraction:
#             return [PolarWaterInteraction()]
#         else:
#             return []
#     else:
#         return []
#
# # Get polar contacts and interactions between 2 residues
# def get_polar_interactions(res1, res2):
#     polar_interactions = []
#     polar_interactions += get_hbond_interactions(res1, res2)
#     # AK: as discussed disabled for now
#     #polar_interactions += get_polar_water_interactions(res1, res2)
#     #polar_interactions += get_polar_sidechain_sidechain_interactions(res1, res2)
#     #polar_interactions += get_polar_sidechain_backbone_interactions(res1, res2)
#     #polar_interactions += get_polar_backbone_sidechain_interactions(res1, res2)
#
#     return polar_interactions
