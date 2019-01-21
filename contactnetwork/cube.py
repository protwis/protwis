from Bio.PDB import Selection, PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch

from contactnetwork.interaction import *
from contactnetwork.pdb import *
from contactnetwork.models import *

from structure.models import Structure

# Distance between residues in peptide
NUM_SKIP_RESIDUES = 4


def compute_interactions(pdb_name,save_to_db = False):
    # Ensure that the PDB name is lowercase
    pdb_name = pdb_name.lower()

    # Get the pdb structure
    struc = Structure.objects.get(protein_conformation__protein__entry_name=pdb_name)

    # Get the preferred chain
    preferred_chain = struc.preferred_chain.split(',')[0]
    # Get the Biopython structure for the PDB
    s = pdb_get_structure(pdb_name)[0]
    chain = s[preferred_chain]


     # remove residues without GN and only those matching receptor.
    residues = struc.protein_conformation.residue_set.exclude(generic_number=None).all()
    dbres = {}
    for r in residues:
        dbres[r.sequence_number] = r
    ids_to_remove = []
    for res in chain:
        if not res.id[1] in dbres.keys():
            ids_to_remove.append(res.id)
    for i in ids_to_remove:
        chain.detach_child(i)


    distances = []
    for i1,res1 in enumerate(chain,1):
        for i2,res2 in enumerate(chain,1):
            if i2>i1:
                # Do not calculate twice.
                distance = res1['CA']-res2['CA']
                distances.append((dbres[res1.id[1]],dbres[res2.id[1]],distance))


    atom_list = Selection.unfold_entities(s[preferred_chain], 'A')
    # Search for all neighbouring residues
    ns = NeighborSearch(atom_list)
    all_neighbors = ns.search_all(4.5, "R")

    # Filter all pairs containing non AA residues
    all_aa_neighbors = [pair for pair in all_neighbors if is_aa(pair[0]) and is_aa(pair[1])]

    # Only include contacts between residues less that NUM_SKIP_RESIDUES sequence steps apart
    all_aa_neighbors = [pair for pair in all_aa_neighbors if abs(pair[0].id[1] - pair[1].id[1]) > NUM_SKIP_RESIDUES]

    # For each pair of interacting residues, determine the type of interaction
    interactions = [InteractingPair(res_pair[0], res_pair[1], get_interactions(res_pair[0], res_pair[1]),dbres[res_pair[0].id[1]], dbres[res_pair[1].id[1]],struc) for res_pair in all_aa_neighbors]

    # Split unto classified and unclassified.
    classified = [interaction for interaction in interactions if len(interaction.get_interactions()) > 0]

    if save_to_db: 

        # Delete previous for faster load in
        InteractingResiduePair.objects.filter(referenced_structure=struc).all().delete()

        bulk_pair = []
        for d in distances:
            pair = InteractingResiduePair(res1=d[0], res2=d[1], referenced_structure=struc)
            bulk_pair.append(pair)

        pairs = InteractingResiduePair.objects.bulk_create(bulk_pair)

        bulk_distances = []
        for i,d in enumerate(distances):
            distance = Distance(distance=d[2],interacting_pair=pairs[i])
            bulk_distances.append(distance)

        pairs = Distance.objects.bulk_create(bulk_distances)

        for p in classified:
            p.save_into_database()


    return classified, distances
