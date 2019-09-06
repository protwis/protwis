from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from django.db.models import Max, Min

from contactnetwork.distances import *
from contactnetwork.models import *
from protein.models import *
from signprot.models import *
from structure.models import *

import logging, json, os

import pandas as pd
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd
from statistics import mean

class Command(BaseCommand):

    help = "Evaluate activation states and level of IC opening."

    logger = logging.getLogger(__name__)

    def handle(self, *args, **options):
        # Loop over classes
        class_slugs = list(ProteinFamily.objects.filter(parent__slug="000") \
                            .filter(slug__startswith="00").values_list("slug"))

        for slug in class_slugs:
            print("Processing class {}".format(slug[0]))

            # grab all PDB-codes for this class
            structure_ids = list(Structure.objects.filter(protein_conformation__protein__family__slug__startswith=slug[0]) \
                                .values_list("pdb_code__index"))
            structure_ids = [x[0] for x in structure_ids]

#            print("Identified the following PDBs")
#            print(structure_ids)

            # Get all PDB-codes for G-protein coupled structures in this class
            # extra: filter 6CMO (unfit reference, see hardcoded exceptions)
            active_ids = list(SignprotComplex.objects.filter(structure__pdb_code__index__in=structure_ids) \
                                .exclude(structure__pdb_code__index="6CMO") \
                                .values_list("structure__pdb_code__index"))
            active_ids = [x[0] for x in active_ids] # flatten
            #                print("The following PDBs are G-prot complex structures:")
            #                print(active_ids)

            # Grab most inactive PDB per ligandType -> 2x46 - 6x37 distance present and <13Å (all classes)
            inactive_ids = list(Distance.objects.filter(distance__lt=1300) \
                                .filter(gn1="2x46").filter(gn2="6x37") \
                                .filter(structure__pdb_code__index__in=structure_ids) \
                                .order_by("structure__protein_conformation__protein__family__parent__parent__name", "distance") \
                                .distinct("structure__protein_conformation__protein__family__parent__parent__name") \
                                .values_list("structure__pdb_code__index"))
            inactive_ids = [x[0] for x in inactive_ids]

            if len(structure_ids) > 0 and len(active_ids) > 0 and len(inactive_ids) > 0:

#                print("The following PDBs are inactive state structures:")
#                print(inactive_ids)

                # create distance matrix for given structures on lower half TM + G-Prot only residues
                dis = Distances()
                dis.lower_only = True # only lower half TM + G-prot only helices
                dis.load_pdbs(structure_ids)
                distance_matrix = dis.get_distance_matrix(True, False) # normalize, but don't use the cache
                distance_matrix = pd.DataFrame(distance_matrix, columns=dis.pdbs, index=dis.pdbs)

                # # Calculate score per pdbs directly based on distance matrix
                # scoring_results = {}
                # for pdb in dis.pdbs:
                #     print("Processing {}".format(pdb))
                #     min_active_distance = min(distance_matrix.loc[pdb, active_ids])
                #     min_inactive_distance = min(distance_matrix.loc[pdb, inactive_ids])
                #
                #     scoring_results[pdb] = min_active_distance-min_inactive_distance


                # hierarchical clustering -> create distance matrix from tree
                hclust = sch.linkage(ssd.squareform(distance_matrix), method='average')
                tree = sch.to_tree(hclust, False)
                tree_distance = getDistanceMatrix(tree, dis.pdbs)
                finalMap = {}
                for d in tree_distance:
                    finalMap.update(d)

                # Calculate score per pdbs
                scoring_results = {}
                for pdb in dis.pdbs:
#                    print("Processing {}".format(pdb))
                    min_active_distance = mean([ finalMap[pdb+"_"+x] for x in active_ids ])
                    min_inactive_distance = mean([ finalMap[pdb+"_"+x] for x in inactive_ids ])
                    scoring_results[pdb] = min_active_distance-min_inactive_distance

#                    print("{}|{}|{}|{}".format(pdb, scoring_results[pdb], min_active_distance, min_inactive_distance))

                # Hardcoded annotations
                hardcoded = {
                    "6CMO" : "active", # Complex with G prot - irregular conformation
                    "5ZKP" : "unknown" # Unknown state (auto-inhibited with H8?)
                }

                distances = list(Distance.objects.filter(gn1="2x46").filter(gn2="6x37") \
                                    .filter(structure__pdb_code__index__in=structure_ids) \
                                    .distinct("gns_pair", "structure") \
                                    .values_list("structure__pdb_code__index", "distance"))

                range_distance = Distance.objects.filter(gn1="2x46").filter(gn2="6x37") \
                                    .filter(structure__pdb_code__index__in=structure_ids) \
                                    .aggregate(Max('distance'), Min('distance'))

                min_open = range_distance['distance__min']
                max_open = range_distance['distance__max']

                # find smallest distance between any active structure and any inactive structure
                lowest_inactive_distance = min([ finalMap[y+"_"+x] for y in inactive_ids for x in active_ids ])
                for entry in distances:
                    # Percentage score
                    percentage = int(round((entry[1]-min_open)/(max_open-min_open)*100))
                    if percentage < 0:
                        percentage = 0
                    elif percentage > 100:
                        percentage = 100

                    # Classification
                    score = scoring_results[entry[0]]
                    structure_state = "inactive"
                    if entry[1] >= 13: # below this distance always inactive structure
                        if score < -0.95*lowest_inactive_distance:
                            structure_state = "active"
                        elif score < 0:
                            structure_state = "intermediate"

                    # UGLY: knowledge-based hardcoded corrections
                    if entry[0] in hardcoded:
                        structure_state = hardcoded[entry[0]]

                    # Store for structure
                    struct = Structure.objects.get(pdb_code__index=entry[0])
                    struct.state, created = ProteinState.objects.get_or_create(slug=structure_state, defaults={'name': structure_state.capitalize()})
                    struct.tm6_angle = percentage
                    struct.save()

                    #print("Class {}: structure {} to state {} and opening is {}%".format(slug, entry[0], structure_state, percentage))
            elif len(structure_ids) > 0:
                distances = list(Distance.objects.filter(gn1="2x46").filter(gn2="6x37") \
                                    .filter(structure__pdb_code__index__in=structure_ids) \
                                    .distinct("gns_pair", "structure") \
                                    .values_list("structure__pdb_code__index", "distance"))

                range_distance = Distance.objects.filter(gn1="2x46").filter(gn2="6x37") \
                                    .aggregate(Max('distance'), Min('distance'))

                min_open = range_distance['distance__min']
                max_open = range_distance['distance__max']
                for entry in distances:
                    # Percentage score
                    percentage = int(round((entry[1]-min_open)/(max_open-min_open)*100))
                    if percentage < 0:
                        percentage = 0
                    elif percentage > 100:
                        percentage = 100

                    # Store for structure
                    struct = Structure.objects.get(pdb_code__index=entry[0])
                    struct.tm6_angle = percentage

                    # Definitely an inactive state structure When distance is smaller than 13Å
                    if entry[1] < 13:
                        struct.state, created = ProteinState.objects.get_or_create(slug="inactive", defaults={'name': "Inactive"})

                    # Save changes
                    struct.save()

                    #print("Class {}: structure {} to state ? and opening is {}%".format(slug, entry[0], percentage))

def getDistanceMatrix(node, pdbs):
    if node.is_leaf():
#        print("Found {} with {}".format(node.id, node.dist))
        return [{pdbs[node.id]: node.dist}, {pdbs[node.id]+"_"+pdbs[node.id]: node.dist}]
    else:
        # Get the left half
        left = getDistanceMatrix(node.get_left(), pdbs)

        # Get the right half
        right = getDistanceMatrix(node.get_right() , pdbs)

        # Combined distances
        combined_results = {}
        own_distance = {}
        for i in left[0]:
            own_distance[i] = left[0][i] + node.dist
            for j in right[0]:
                combined_results[i+"_"+j] = left[0][i] + right[0][j]
                combined_results[j+"_"+i] = combined_results[i+"_"+j]
        for j in right[0]:
            own_distance[j] = right[0][j] + node.dist

        combined_results.update(left[1])
        combined_results.update(right[1])

        return [own_distance, combined_results]
