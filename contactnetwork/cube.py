from django.conf import settings

from Bio.PDB import Selection, PDBParser, Polypeptide
from Bio.PDB.NeighborSearch import NeighborSearch

from contactnetwork.interaction import *
from contactnetwork.pdb import *
from contactnetwork.models import *
from contactnetwork.residue import is_aa
from ligand.models import LigandPeptideStructure
from residue.functions import DummyResidue
from io import StringIO

from protein.models import ProteinConformation

from structure.models import Structure, StructureExtraProteins

from signprot.models import SignprotComplex

import copy
import yaml
import os

# Distance between residues in peptide
NUM_SKIP_RESIDUES = 1

def compute_interactions(pdb_name, protein=None, signprot=None, lig=None, do_interactions=False, do_complexes=False, do_peptide_ligand=False, save_to_db=False, file_input=False):
    classified = []
    classified_complex = []
    with open(os.sep.join([settings.DATA_DIR, 'residue_data', 'unnatural_amino_acids.yaml']), 'r') as f_yaml:
        unnatural_amino_acids = yaml.safe_load(f_yaml)
        unnatural_amino_acids = {str(x):unnatural_amino_acids[x] for x in unnatural_amino_acids}

    if file_input:
        # Get the preferred chain
        if protein:
            preferred_chain = protein.preferred_chain
        else:
            preferred_chain = 'A' #I guess?
        # Get the Biopython structure for the PDB
        parser = PDBParser()
        s = parser.get_structure("ComplexModel", pdb_name)[0]
        #s = pdb_get_structure(pdb_name)[0]
        chain = s[preferred_chain]

        # fetching residues for wild type structure of protein
        if do_complexes:
            residues = protein.protein_conformation.residue_set.all().prefetch_related('generic_number')
        else:
            residues = protein.protein_conformation.residue_set.exclude(generic_number=None).all().prefetch_related('generic_number')
        struc = protein
        # residues = ProteinConformation.objects.get(protein__entry_name=protein.protein.entry_name).residue_set.exclude(generic_number=None).all().prefetch_related('generic_number')
        # CREATING THE DUMMY RESIDUES FOR ALL THE RESIDUES OF CHAIN A (aka protein)
        # for res in chain:
        #     try:
        #         one_letter = Polypeptide.three_to_one(res.get_resname())
        #     except KeyError:
        #         if res.get_resname() in unnatural_amino_acids:
        #             one_letter = unnatural_amino_acids[res.get_resname()]
        #         else:
        #             print('WARNING: {} residue in structure {} is missing from unnatural amino acid definitions (data/protwis/gpcr/residue_data/unnatural_amino_acids.yaml)'.format(res, struc))
        #             continue
        #     res_obj = DummyResidue(res.get_id()[1], one_letter, res.get_resname())
        #     dbres[res_obj.sequence_number] = res_obj
    else:
        # Ensure that the PDB name is lowercase
        pdb_name = pdb_name.lower()
        struc = Structure.objects.get(protein_conformation__protein__entry_name=pdb_name)
        # Get the pdb structure
        pdb_io = StringIO(struc.pdb_data.pdb)
        # Get the preferred chain
        preferred_chain = struc.preferred_chain.split(',')[0]
        # Get the Biopython structure for the PDB
        s = PDBParser(PERMISSIVE=True, QUIET=True).get_structure('ref', pdb_io)[0]
        #s = pdb_get_structure(pdb_name)[0]
        chain = s[preferred_chain]
        # remove residues without GN and only those matching receptor.
        if do_complexes:
            residues = struc.protein_conformation.residue_set.all().prefetch_related('generic_number')
        else:
            residues = struc.protein_conformation.residue_set.exclude(generic_number=None).all().prefetch_related('generic_number')

    dbres = {}
    dblabel = {}
    for r in residues:
        dbres[r.sequence_number] = r
        if r.generic_number:
            dblabel[r.sequence_number] = r.generic_number.label
        else:
            dblabel[r.sequence_number] = '-'
    ids_to_remove = []
    for res in chain:
        if not res.id[1] in dbres.keys() and res.get_resname() != "HOH":
            ids_to_remove.append(res.id)
    for i in ids_to_remove:
        chain.detach_child(i)
    if do_interactions:
        atom_list = Selection.unfold_entities(s[preferred_chain], 'A')
        # Search for all neighbouring residues
        ns = NeighborSearch(atom_list)
        all_neighbors = ns.search_all(6.6, "R")
        # Filter all pairs containing non AA residues
        all_aa_neighbors = [pair for pair in all_neighbors if is_aa(pair[0]) and is_aa(pair[1])]
        # Only include contacts between residues more than NUM_SKIP_RESIDUES sequence steps apart
        all_aa_neighbors = [pair for pair in all_aa_neighbors if abs(pair[0].id[1] - pair[1].id[1]) > NUM_SKIP_RESIDUES]
        # For each pair of interacting residues, determine the type of interaction
        interactions = [InteractingPair(res_pair[0], res_pair[1], dbres[res_pair[0].id[1]], dbres[res_pair[1].id[1]], struc) for res_pair in all_aa_neighbors if not is_water(res_pair[0]) and not is_water(res_pair[1])]
        # Split unto classified and unclassified.
        classified = [interaction for interaction in interactions if len(interaction.get_interactions()) > 0]

    if do_complexes:
        try:
            if file_input:
                if protein:
                    signprot_chain = protein.signprot_complex.alpha
                else:
                    signprot_chain = 'B' ### Hardcoded for now
                extension = ''
                pdb_name = signprot.entry_name
            else:
                # check if structure in signprot_complex
                signprot_chain = ""
                extension = ""
                if StructureExtraProteins.objects.filter(structure=struc, category="Arrestin"):
                    signprot_chain = StructureExtraProteins.objects.get(structure=struc, category="Arrestin").chain
                    extension = "_arrestin"
                else:
                    signprot_chain = SignprotComplex.objects.get(structure=struc).alpha
                    extension = "_a"

                # Workaround for fused receptor - signaling proteins constructs
                if signprot_chain == preferred_chain:
                    pdb_io = StringIO(struc.pdb_data.pdb)
                    s = PDBParser(PERMISSIVE=True, QUIET=True).get_structure('ref', pdb_io)[0]

            # Get all GPCR residue atoms based on preferred chain
            gpcr_atom_list = [ atom for residue in Selection.unfold_entities(s[preferred_chain], 'R') if is_aa(residue) and residue.get_id()[1] in dbres \
                            for atom in residue.get_atoms()]

            # TOFIX: Current workaround is forcing _a to pdb for indicating alpha-subunit
            residues_sign = ProteinConformation.objects.get(protein__entry_name=pdb_name+extension).residue_set.all().prefetch_related('generic_number')

            # grab labels from sign protein
            dbres_sign = {}
            dblabel_sign = {}
            for r in residues_sign:
                dbres_sign[r.sequence_number] = r
                dblabel_sign[r.sequence_number] = r.generic_number.label

            # Get all residue atoms from the coupled protein (e.g. G-protein)
            sign_atom_list = [ atom for residue in Selection.unfold_entities(s[signprot_chain], 'R') if is_aa(residue) and residue.get_id()[1] in dbres_sign \
                                for atom in residue.get_atoms()]

            ns_gpcr = NeighborSearch(gpcr_atom_list)
            ns_sign = NeighborSearch(sign_atom_list)

            # For each GPCR atom perform the neighbor search on the signaling protein
            all_neighbors = {(gpcr_atom.parent, match_res) for gpcr_atom in gpcr_atom_list
                            for match_res in ns_sign.search(gpcr_atom.coord, 4.5, "R")}

            # Find interactions
            interactions = [InteractingPair(res_pair[0], res_pair[1], dbres[res_pair[0].id[1]], dbres_sign[res_pair[1].id[1]], struc) for res_pair in all_neighbors if res_pair[0].id[1] in dbres and res_pair[1].id[1] in dbres_sign ]

            # Filter unclassified interactions
            classified_complex = [interaction for interaction in interactions if len(interaction.get_interactions()) > 0]

            # Convert to dictionary for water calculations
            interaction_pairs = {}
            for pair in classified_complex:
                res_1 = pair.get_residue_1()
                res_2 = pair.get_residue_2()
                key =  res_1.get_parent().get_id()+str(res_1.get_id()[1]) + "_" + res_2.get_parent().get_id()+str(res_2.get_id()[1])
                interaction_pairs[key] = pair

                # TODO: DEBUG AND VERIFY this code as water-mediated interactions were present at this time
                # 1. UPDATE complexes to include also mini Gs and peptides (e.g. 4X1H/6FUF/5G53)
                # 2. Run and verify water-mediated do_interactions
                # 3. Improve the intersection between the two hit lists

                ## TODO: cleaner intersection between hits from the two Lists
                # see new code below
#                for water_pair_one in water_neighbors_gpcr:
#                    for water_pair_two in water_neighbors_sign:
#                        if water_pair_one[0]==water_pair_two[0]:
#                            res_1 = water_pair_one[1]
#                            res_2 = water_pair_two[1]
#                            key =  res_1.get_parent().get_id()+str(res_1.get_id()[1]) + "_" + res_2.get_parent().get_id()+str(res_2.get_id()[1])

                            # Check if interaction is polar
#                            if any(get_polar_interactions(water_pair_one[0].get_parent(), water_pair_one[1])) and any(get_polar_interactions(water_pair_two[0].get_parent(), water_pair_two[1])):
                                # TODO Check if water interaction is already present (e.g. multiple waters)
                                # TODO Is splitting of sidechain and backbone-mediated interactions desired?
#                                if not key in interaction_pairs:
#                                    interaction_pairs[key] = InteractingPair(res_1, res_2, dbres[res_1.id[1]], dbres_sign[res_2.id[1]], struc)

                                # TODO: fix assignment of interacting atom labels (now seems limited to residues)
#                                interaction_pairs[key].interactions.append(WaterMediated(a + "|" + str(water_pair_one[0].get_parent().get_id()[1]), b))

        except SignprotComplex.DoesNotExist:
#            print("No complex definition found for", pdb_name)
            log = "No complex definition found for " + pdb_name
        except ProteinConformation.DoesNotExist:
            print("No protein conformation definition found for signaling protein of ", pdb_name)
#            log = "No protein conformation definition found for signaling protein of " + pdb_name

    # TODO for peptides:
    # - support deviating atom names in unnatural AAs for detecting HB-donors and acceptors
    # - support for terminal carboxyl (OXT) + amidation etc for ionic+hbond interactions
    if do_peptide_ligand:
        if file_input:
            #Fetching from the correct model
            ligands = LigandPeptideStructure.objects.filter(structure=protein, ligand=lig)
        else:
            ligands = LigandPeptideStructure.objects.filter(structure=struc)
        if len(ligands)>0:
            # Get all GPCR residue atoms based on preferred chain
            gpcr_atom_list = [atom for residue in Selection.unfold_entities(s[preferred_chain], 'R') if is_aa(residue) and residue.get_id()[1] in dbres \
                            for atom in residue.get_atoms()]
            ns_gpcr = NeighborSearch(gpcr_atom_list)

            classified_ligand_complex = {}
            for ligand in ligands:
                pep_chain = ligand.chain

                # Peptide chain not resolved/modelled
                if pep_chain.strip() == "":
                    classified_ligand_complex[ligand] = []
                    continue

                dbres_pep = {}
                for res in s[pep_chain]:
                    try:
                        one_letter = Polypeptide.three_to_one(res.get_resname())
                    except KeyError:
                        if res.get_resname() in unnatural_amino_acids:
                            one_letter = unnatural_amino_acids[res.get_resname()]
                        else:
                            print('WARNING: {} residue in structure {} is missing from unnatural amino acid definitions (data/protwis/gpcr/residue_data/unnatural_amino_acids.yaml)'.format(res, struc))
                            continue

                    res_obj = DummyResidue(res.get_id()[1], one_letter, res.get_resname())
                    dbres_pep[res.get_id()[1]] = res_obj

                pep_atom_list = [ atom for residue in Selection.unfold_entities(s[pep_chain], 'R') if residue.get_id()[1] in dbres_pep \
                                    for atom in residue.get_atoms()]

                ns_pep = NeighborSearch(pep_atom_list)

                # For each GPCR atom perform the neighbor search on the signaling protein
                all_neighbors = {(gpcr_atom.parent, match_res) for gpcr_atom in gpcr_atom_list
                                for match_res in ns_pep.search(gpcr_atom.coord, 4.5, "R")}

                # Find interactions
                interactions = [InteractingPair(res_pair[0], res_pair[1], dbres[res_pair[0].id[1]], dbres_pep[res_pair[1].id[1]], protein) for res_pair in all_neighbors if res_pair[0].id[1] in dbres and res_pair[1].id[1] in dbres_pep ]
                # Filter unclassified interactions
                classified_ligand_complex[ligand] = [interaction for interaction in interactions if len(interaction.get_interactions()) > 0]
    if save_to_db:

        if do_interactions:
            # Delete previous for faster load in
            InteractingResiduePair.objects.filter(referenced_structure=struc).all().delete()

            # Create interaction dictionary
            interaction_pairs = {}
            for pair in classified:
                res_1 = pair.get_residue_1()
                res_2 = pair.get_residue_2()
                key =  res_1.get_parent().get_id()+str(res_1.get_id()[1]) + "_" + res_2.get_parent().get_id()+str(res_2.get_id()[1])
                interaction_pairs[key] = pair

            # POSSIBLE ADDON: support for multiple water-mediated bonds
            ## Obtain list of water molecules
            water_list = { water for residue in s[preferred_chain] if residue.get_resname() == "HOH" for water in residue.get_atoms() }
            if len(water_list) > 0:
                ## Iterate water molecules over residue atom list
                water_neighbors = [(water, match_res) for water in water_list
                                for match_res in ns.search(water.coord, 3.5, "R") if not is_water(match_res) and (is_hba(match_res) or is_hbd(match_res))]

                # intersect between residues sharing the same interacting water
                for index_one in range(len(water_neighbors)):
                    water_pair_one = water_neighbors[index_one]

                    for index_two in [ index for index in range(index_one+1, len(water_neighbors)) if water_pair_one[0]==water_neighbors[index][0] ]:
                        water_pair_two = water_neighbors[index_two]
                        res_1 = water_pair_one[1]
                        res_2 = water_pair_two[1]

                        # TODO: order residues + check minimum spacing between residues
                        key =  res_1.get_parent().get_id()+str(res_1.get_id()[1]) + "_" + res_2.get_parent().get_id()+str(res_2.get_id()[1])

                        # Verify h-bonds between water and both residues
                        matches_one = InteractingPair.verify_water_hbond(water_pair_one[1], water_pair_one[0])
                        matches_two = InteractingPair.verify_water_hbond(water_pair_two[1], water_pair_two[0])
                        if len(matches_one) > 0 and len(matches_two) > 0:
                            # if not exists, create residue pair without interactions
                            if not key in interaction_pairs:
                                interaction_pairs[key] = InteractingPair(res_1, res_2, dbres[res_1.id[1]], dbres[res_2.id[1]], struc)

                            for a,b in zip(matches_one, matches_two):
                                # HACK: store water ID as part of first atom name
                                interaction_pairs[key].interactions.append(WaterMediated(a + "|" + str(water_pair_one[0].get_parent().get_id()[1]), b))

            for p in classified:
                p.save_into_database()

        if do_complexes:
            for pair in classified_complex:
                pair.save_into_database()

        if do_peptide_ligand and len(ligands)>0:
            for ligand in ligands:
                for pair in classified_ligand_complex[ligand]:
                    # SAME HERE, SOME CHANGES WERE NEEDED
                    pair.save_peptide_interactions(ligand)
    return classified

# def compute_interactions(pdb_name,save_to_db = False):
#
#     do_distances = False ## Distance calculation moved to build_structure_angles
#     do_interactions = True
#     do_complexes = True
#     do_peptide_ligand = True
#     distances = []
#     classified = []
#     classified_complex = []
#     with open(os.sep.join([settings.DATA_DIR, 'residue_data', 'unnatural_amino_acids.yaml']), 'r') as f_yaml:
#         unnatural_amino_acids = yaml.safe_load(f_yaml)
#         unnatural_amino_acids = {str(x):unnatural_amino_acids[x] for x in unnatural_amino_acids}
#
#     # Ensure that the PDB name is lowercase
#     pdb_name = pdb_name.lower()
#
#     # Get the pdb structure
#     struc = Structure.objects.get(protein_conformation__protein__entry_name=pdb_name)
#     pdb_io = StringIO(struc.pdb_data.pdb)
#     # Get the preferred chain
#     preferred_chain = struc.preferred_chain.split(',')[0]
#
#     # Get the Biopython structure for the PDB
#     s = PDBParser(PERMISSIVE=True, QUIET=True).get_structure('ref', pdb_io)[0]
#     #s = pdb_get_structure(pdb_name)[0]
#     chain = s[preferred_chain]
#     #return classified, distances
#
#     # remove residues without GN and only those matching receptor.
#     residues = struc.protein_conformation.residue_set.exclude(generic_number=None).all().prefetch_related('generic_number')
#     dbres = {}
#     dblabel = {}
#     for r in residues:
#         dbres[r.sequence_number] = r
#         dblabel[r.sequence_number] = r.generic_number.label
#     ids_to_remove = []
#     for res in chain:
#         if not res.id[1] in dbres.keys() and res.get_resname() != "HOH":
#             ids_to_remove.append(res.id)
#     for i in ids_to_remove:
#         chain.detach_child(i)
#
#     # if do_distances:
#     #     for i1,res1 in enumerate(chain,1):
#     #         if not is_water(res1):
#     #             for i2,res2 in enumerate(chain,1):
#     #                 if i2>i1 and not is_water(res2):
#     #                     # Do not calculate twice.
#     #                     distance = res1['CA']-res2['CA']
#     #                     distances.append((dbres[res1.id[1]],dbres[res2.id[1]],distance,dblabel[res1.id[1]],dblabel[res2.id[1]]))
#
#     if do_interactions:
#         atom_list = Selection.unfold_entities(s[preferred_chain], 'A')
#
#         # Search for all neighbouring residues
#         ns = NeighborSearch(atom_list)
#         all_neighbors = ns.search_all(6.6, "R")
#
#         # Filter all pairs containing non AA residues
#         all_aa_neighbors = [pair for pair in all_neighbors if is_aa(pair[0]) and is_aa(pair[1])]
#
#         # Only include contacts between residues more than NUM_SKIP_RESIDUES sequence steps apart
#         all_aa_neighbors = [pair for pair in all_aa_neighbors if abs(pair[0].id[1] - pair[1].id[1]) > NUM_SKIP_RESIDUES]
#
#         # For each pair of interacting residues, determine the type of interaction
#         interactions = [InteractingPair(res_pair[0], res_pair[1], dbres[res_pair[0].id[1]], dbres[res_pair[1].id[1]], struc) for res_pair in all_aa_neighbors if not is_water(res_pair[0]) and not is_water(res_pair[1]) ]
#
#         # Split unto classified and unclassified.
#         classified = [interaction for interaction in interactions if len(interaction.get_interactions()) > 0]
#
#     if do_complexes:
#         try:
#             # check if structure in signprot_complex
#             signprot_chain = ""
#             extension = ""
#             if StructureExtraProteins.objects.filter(structure=struc, category="Arrestin"):
#                 signprot_chain = StructureExtraProteins.objects.get(structure=struc, category="Arrestin").chain
#                 extension = "_arrestin"
#             else:
#                 signprot_chain = SignprotComplex.objects.get(structure=struc).alpha
#                 extension = "_a"
#
#             # Workaround for fused receptor - signaling proteins constructs
#             if signprot_chain == preferred_chain:
#                 pdb_io = StringIO(struc.pdb_data.pdb)
#                 s = PDBParser(PERMISSIVE=True, QUIET=True).get_structure('ref', pdb_io)[0]
#
#             # Get all GPCR residue atoms based on preferred chain
#             gpcr_atom_list = [ atom for residue in Selection.unfold_entities(s[preferred_chain], 'R') if is_aa(residue) and residue.get_id()[1] in dbres \
#                             for atom in residue.get_atoms()]
#
#             # For each pair of interacting residues, determine the type of interaction
#             #residues_sign = ProteinConformation.objects.get(protein__entry_name=pdb_name+"_"+complex.alpha.lower()).residue_set.exclude(generic_number=None).all().prefetch_related('generic_number')
#             # TOFIX: Current workaround is forcing _a to pdb for indicating alpha-subunit
#             residues_sign = ProteinConformation.objects.get(protein__entry_name=pdb_name+extension).residue_set.exclude(generic_number=None).all().prefetch_related('generic_number')
#
#             # grab labels from sign protein
#             dbres_sign = {}
#             dblabel_sign = {}
#             for r in residues_sign:
#                 dbres_sign[r.sequence_number] = r
#                 dblabel_sign[r.sequence_number] = r.generic_number.label
#
#             # Get all residue atoms from the coupled protein (e.g. G-protein)
#             sign_atom_list = [ atom for residue in Selection.unfold_entities(s[signprot_chain], 'R') if is_aa(residue) and residue.get_id()[1] in dbres_sign \
#                                 for atom in residue.get_atoms()]
#
#             ns_gpcr = NeighborSearch(gpcr_atom_list)
#             ns_sign = NeighborSearch(sign_atom_list)
#
#             # For each GPCR atom perform the neighbor search on the signaling protein
#             all_neighbors = {(gpcr_atom.parent, match_res) for gpcr_atom in gpcr_atom_list
#                             for match_res in ns_sign.search(gpcr_atom.coord, 4.5, "R")}
#
#             # Find interactions
#             interactions = [InteractingPair(res_pair[0], res_pair[1], dbres[res_pair[0].id[1]], dbres_sign[res_pair[1].id[1]], struc) for res_pair in all_neighbors if res_pair[0].id[1] in dbres and res_pair[1].id[1] in dbres_sign ]
#
#             # Filter unclassified interactions
#             classified_complex = [interaction for interaction in interactions if len(interaction.get_interactions()) > 0]
#
#             # Convert to dictionary for water calculations
#             interaction_pairs = {}
#             for pair in classified_complex:
#                 res_1 = pair.get_residue_1()
#                 res_2 = pair.get_residue_2()
#                 key =  res_1.get_parent().get_id()+str(res_1.get_id()[1]) + "_" + res_2.get_parent().get_id()+str(res_2.get_id()[1])
#                 interaction_pairs[key] = pair
#
#             # # Obtain list of all water molecules in the structure
#             # water_list = { water for chain in s for residue in chain
#             #                 if residue.get_resname() == "HOH" for water in residue.get_atoms() }
#             #
#             # # If waters are present calculate water-mediated interactions
#             # if len(water_list) > 0:
#             #     ## Iterate water molecules over coupled and gpcr atom list
#             #     water_neighbors_gpcr = {(water, match_res) for water in water_list
#             #                     for match_res in ns_gpcr.search(water.coord, 3.5, "R")}
#             #
#             #     water_neighbors_sign = {(water, match_res) for water in water_list
#             #                     for match_res in ns_sign.search(water.coord, 3.5, "R")}
#
#
#                 # TODO: DEBUG AND VERIFY this code as water-mediated interactions were present at this time
#                 # 1. UPDATE complexes to include also mini Gs and peptides (e.g. 4X1H/6FUF/5G53)
#                 # 2. Run and verify water-mediated do_interactions
#                 # 3. Improve the intersection between the two hit lists
#
#                 ## TODO: cleaner intersection between hits from the two Lists
#                 # see new code below
# #                for water_pair_one in water_neighbors_gpcr:
# #                    for water_pair_two in water_neighbors_sign:
# #                        if water_pair_one[0]==water_pair_two[0]:
# #                            res_1 = water_pair_one[1]
# #                            res_2 = water_pair_two[1]
# #                            key =  res_1.get_parent().get_id()+str(res_1.get_id()[1]) + "_" + res_2.get_parent().get_id()+str(res_2.get_id()[1])
#
#                             # Check if interaction is polar
# #                            if any(get_polar_interactions(water_pair_one[0].get_parent(), water_pair_one[1])) and any(get_polar_interactions(water_pair_two[0].get_parent(), water_pair_two[1])):
#                                 # TODO Check if water interaction is already present (e.g. multiple waters)
#                                 # TODO Is splitting of sidechain and backbone-mediated interactions desired?
# #                                if not key in interaction_pairs:
# #                                    interaction_pairs[key] = InteractingPair(res_1, res_2, dbres[res_1.id[1]], dbres_sign[res_2.id[1]], struc)
#
#                                 # TODO: fix assignment of interacting atom labels (now seems limited to residues)
# #                                interaction_pairs[key].interactions.append(WaterMediated(a + "|" + str(water_pair_one[0].get_parent().get_id()[1]), b))
#
#         except SignprotComplex.DoesNotExist:
# #            print("No complex definition found for", pdb_name)
#             log = "No complex definition found for " + pdb_name
#         except ProteinConformation.DoesNotExist:
#             print("No protein conformation definition found for signaling protein of ", pdb_name)
# #            log = "No protein conformation definition found for signaling protein of " + pdb_name
#
#     # TODO for peptides:
#     # - support deviating atom names in unnatural AAs for detecting HB-donors and acceptors
#     # - support for terminal carboxyl (OXT) + amidation etc for ionic+hbond interactions
#     if do_peptide_ligand:
#         ligands = LigandPeptideStructure.objects.filter(structure=struc)
#         if len(ligands)>0:
#             # Get all GPCR residue atoms based on preferred chain
#             gpcr_atom_list = [ atom for residue in Selection.unfold_entities(s[preferred_chain], 'R') if is_aa(residue) and residue.get_id()[1] in dbres \
#                             for atom in residue.get_atoms()]
#             ns_gpcr = NeighborSearch(gpcr_atom_list)
#
#             classified_ligand_complex = {}
#             for ligand in ligands:
#                 pep_chain = ligand.chain
#
#                 # Peptide chain not resolved/modelled
#                 if pep_chain.strip() == "":
#                     classified_ligand_complex[ligand] = []
#                     continue
#
#                 dbres_pep = {}
#                 for res in s[pep_chain]:
#                     try:
#                         one_letter = Polypeptide.three_to_one(res.get_resname())
#                     except KeyError:
#                         if res.get_resname() in unnatural_amino_acids:
#                             one_letter = unnatural_amino_acids[res.get_resname()]
#                         else:
#                             print('WARNING: {} residue in structure {} is missing from unnatural amino acid definitions (data/protwis/gpcr/residue_data/unnatural_amino_acids.yaml)'.format(res, struc))
#                             continue
#
#                     res_obj = DummyResidue(res.get_id()[1], one_letter, res.get_resname())
#                     dbres_pep[res.get_id()[1]] = res_obj
#
#                 pep_atom_list = [ atom for residue in Selection.unfold_entities(s[pep_chain], 'R') if residue.get_id()[1] in dbres_pep \
#                                     for atom in residue.get_atoms()]
#
#                 ns_pep = NeighborSearch(pep_atom_list)
#
#                 # For each GPCR atom perform the neighbor search on the signaling protein
#                 all_neighbors = {(gpcr_atom.parent, match_res) for gpcr_atom in gpcr_atom_list
#                                 for match_res in ns_pep.search(gpcr_atom.coord, 4.5, "R")}
#
#                 # Find interactions
#                 interactions = [InteractingPair(res_pair[0], res_pair[1], dbres[res_pair[0].id[1]], dbres_pep[res_pair[1].id[1]], struc) for res_pair in all_neighbors if res_pair[0].id[1] in dbres and res_pair[1].id[1] in dbres_pep ]
#
#                 # Filter unclassified interactions
#                 classified_ligand_complex[ligand] = [interaction for interaction in interactions if len(interaction.get_interactions()) > 0]
#
#     if save_to_db:
#
#         if do_interactions:
#             # Delete previous for faster load in
#             InteractingResiduePair.objects.filter(referenced_structure=struc).all().delete()
#
#             # bulk_pair = []
#             # for d in distances:
#             #     pair = InteractingResiduePair(res1=d[0], res2=d[1], referenced_structure=struc)
#             #     bulk_pair.append(pair)
#
#             # Create interaction dictionary
#             interaction_pairs = {}
#             for pair in classified:
#                 res_1 = pair.get_residue_1()
#                 res_2 = pair.get_residue_2()
#                 key =  res_1.get_parent().get_id()+str(res_1.get_id()[1]) + "_" + res_2.get_parent().get_id()+str(res_2.get_id()[1])
#                 interaction_pairs[key] = pair
#
#             # POSSIBLE ADDON: support for multiple water-mediated bonds
#             ## Obtain list of water molecules
#             water_list = { water for residue in s[preferred_chain] if residue.get_resname() == "HOH" for water in residue.get_atoms() }
#             if len(water_list) > 0:
#                 ## Iterate water molecules over residue atom list
#                 water_neighbors = [(water, match_res) for water in water_list
#                                 for match_res in ns.search(water.coord, 3.5, "R") if not is_water(match_res) and (is_hba(match_res) or is_hbd(match_res))]
#
#                 # intersect between residues sharing the same interacting water
#                 for index_one in range(len(water_neighbors)):
#                     water_pair_one = water_neighbors[index_one]
#
#                     for index_two in [ index for index in range(index_one+1, len(water_neighbors)) if water_pair_one[0]==water_neighbors[index][0] ]:
#                         water_pair_two = water_neighbors[index_two]
#                         res_1 = water_pair_one[1]
#                         res_2 = water_pair_two[1]
#
#                         # TODO: order residues + check minimum spacing between residues
#                         key =  res_1.get_parent().get_id()+str(res_1.get_id()[1]) + "_" + res_2.get_parent().get_id()+str(res_2.get_id()[1])
#
#                         # Verify h-bonds between water and both residues
#                         matches_one = InteractingPair.verify_water_hbond(water_pair_one[1], water_pair_one[0])
#                         matches_two = InteractingPair.verify_water_hbond(water_pair_two[1], water_pair_two[0])
#                         if len(matches_one) > 0 and len(matches_two) > 0:
#                             # if not exists, create residue pair without interactions
#                             if not key in interaction_pairs:
#                                 interaction_pairs[key] = InteractingPair(res_1, res_2, dbres[res_1.id[1]], dbres[res_2.id[1]], struc)
#
#                             for a,b in zip(matches_one, matches_two):
#                                 # HACK: store water ID as part of first atom name
#                                 interaction_pairs[key].interactions.append(WaterMediated(a + "|" + str(water_pair_one[0].get_parent().get_id()[1]), b))
#
#             for p in classified:
#                 p.save_into_database()
#
#         if do_complexes:
#             for pair in classified_complex:
#                 pair.save_into_database()
#
#         if do_peptide_ligand and len(ligands)>0:
#             for ligand in ligands:
#                 for pair in classified_ligand_complex[ligand]:
#                     pair.save_peptide_interactions(ligand)
#
#         # if do_distances:
#         #     # Distance.objects.filter(structure=struc).all().delete()
#         #     bulk_distances = []
#         #     for i,d in enumerate(distances):
#         #         distance = Distance(distance=int(10000*d[2]),res1=d[0], res2=d[1],gn1=d[3], gn2=d[4], gns_pair='_'.join([d[3],d[4]]), structure=struc)
#         #         bulk_distances.append(distance)
#         #         if len(bulk_distances)>1000:
#         #             pairs = Distance.objects.bulk_create(bulk_distances)
#         #             bulk_distances = []
#         #
#         #     pairs = Distance.objects.bulk_create(bulk_distances)
#     return classified, distances
