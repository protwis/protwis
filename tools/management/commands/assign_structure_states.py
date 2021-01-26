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
        self.logger.info('ASSIGNING the "Degree Active" levels and activation states')
        # Loop over classes
        class_slugs = list(ProteinFamily.objects.filter(parent__slug="000") \
                            .filter(slug__startswith="00").values_list("slug"))

        for slug in class_slugs:
            print("Processing class {}".format(slug[0]))

            # grab all PDB-codes for this class
            structure_ids = list(Structure.objects.exclude(refined=True).filter(protein_conformation__protein__family__slug__startswith=slug[0]) \
                                .values_list("pdb_code__index"))

            structure_ids = [x[0] for x in structure_ids]

            if len(structure_ids) > 0:

                # Get all PDB-codes for G-protein coupled structures in this class
                # extra: filter 6CMO (unfit reference, see hardcoded exceptions)
                active_ids = list(SignprotComplex.objects.filter(structure__pdb_code__index__in=structure_ids) \
                                    .exclude(structure__pdb_code__index="6CMO") \
                                    .values_list("structure__pdb_code__index"))
                active_ids = [x[0] for x in active_ids if x[0][0].isnumeric()] # flatten

                # Hardcoded active structures
                if slug[0] == "001":
                    active_ids.extend(["6LI3"])
                elif slug[0] == "004":
                    active_ids = ["7C7Q"]

                # print("The following PDBs are G-prot complex structures:")
                # print(slug[0], active_ids)


                # V1: Grab most inactive PDB per ligandType -> 2x46 - 6x37 distance should be present and < 13Å (all classes)
                # V2: Grab most inactive PDB per Receptor family -> 2x46 - 6x37 distance should be present and < 13Å (cut-off valid for all classes)
                # V3: Just grab most inactive PDBs -> 2x46 - 6x37 distance should be present and < 13Å (cut-off valid for all classes at this point 10-01-2020)
                class_pair_inactives = {}
                class_pair_inactives['001'] = ["2x46_6x37", 11.9] #A
                class_pair_inactives['002'] = ["2x46_6x37", 13] #B1
                class_pair_inactives['003'] = ["2x47_6x37", 1000] #B2 PLACEHOLDER
                class_pair_inactives['004'] = ["2x47_6x37", 14.5] #C
                class_pair_inactives['005'] = ["2x47_6x37", 1000] #D PLACEHOLDER
                class_pair_inactives['006'] = ["2x44_6x31", 13] #F

                inactive_ids = list(Distance.objects.filter(distance__lt=class_pair_inactives[slug[0]][1]*distance_scaling_factor) \
                                    .filter(gns_pair=class_pair_inactives[slug[0]][0]) \
                                    .filter(structure__pdb_code__index__in=structure_ids) \
                                    .exclude(structure__pdb_code__index__in=active_ids) \
                                    .values_list("structure__pdb_code__index"))

                inactive_ids = [x[0] for x in inactive_ids if x[0][0].isnumeric()]

                # HARDCODED INACTIVE STRUCTURES
                if slug[0] == "004":
                    inactive_ids.extend(["4OR2", "6FFI"])
                elif slug[0] == "006":
                    inactive_ids.extend(["4JKV", "6D35", "6D32", "5V57", "5V56", "4N4W"])

                if "6FJ3" in inactive_ids:
                    inactive_ids.remove("6FJ3")

                #print(slug[0], len(inactive_ids),len(active_ids))
                if len(active_ids) > 0 and len(inactive_ids) > 0:
                    # create distance matrix for given structures on lower half TM + G-Prot only residues
                    dis = Distances()
                    dis.filtered_gns = True # only lower half TM + G-prot only helices
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
                    all_scoring_results = {}
                    for pdb in dis.pdbs:
                        min_active_distance = mean([ abs(finalMap[pdb+"_"+x]) for x in active_ids ])
                        min_inactive_distance = mean([ abs(finalMap[pdb+"_"+x]) for x in inactive_ids ])
                        scoring_results[pdb] = min_active_distance-min_inactive_distance
                        all_scoring_results[pdb] = [min_active_distance-min_inactive_distance, min_active_distance, min_inactive_distance]
                        # print("{}|{}|{}|{}".format(pdb, scoring_results[pdb], min_active_distance, min_inactive_distance))

                    min_score = min(scoring_results.items(), key=lambda x: x[1])[1]
                    max_score = max(scoring_results.items(), key=lambda x: x[1])[1]

                    # Hardcoded annotations
                    hardcoded = {
                        "6CMO" : "active", # Complex with G prot - irregular conformation
                        "5ZKP" : "other", # Unknown/other activation states (in this case auto-inhibited with H8?)
                        "5LWE" : "inactive", # Cannot be determined using this method because of missing TM in annotation
                        "6KUX" : "inactive", #
                        "5NX2" : "intermediate", # Closer to active + groups together but internally more inactive
                        "6N51" : "intermediate", # Holds middle between active and inactive
                        "7CA5" : "intermediate"  # Apo state holds middle between active and inactive
                    }

                    # Percentage score for TM2-TM6 opening
                    #range_distance = Distance.objects.filter(gn1="2x46").filter(gn2="6x37") \
                    range_distance = Distance.objects.filter(gns_pair=class_pair_inactives[slug[0]][0]) \
                                        .filter(structure__pdb_code__index__in=structure_ids) \
                                        .aggregate(Max('distance'), Min('distance'))

                    min_open = range_distance['distance__min']
                    max_open = range_distance['distance__max']

                    #distances = list(Distance.objects.filter(gn1="2x46").filter(gn2="6x37") \
                    distances = list(Distance.objects.filter(gns_pair=class_pair_inactives[slug[0]][0]) \
                                        .filter(structure__pdb_code__index__in=structure_ids) \
                                        .distinct("gns_pair", "structure") \
                                        .values_list("structure__pdb_code__index", "distance"))

                    opening_percentage = {}
                    for entry in distances:
                        percentage = int(round((entry[1]-min_open)/(max_open-min_open)*100))
                        if percentage < 0:
                            percentage = 0
                        elif percentage > 100:
                            percentage = 100
                        opening_percentage[entry[0]] = percentage

                    # find smallest distance between any active structure and any inactive structure
                    # lowest_inactive_distance = min([ finalMap[y+"_"+x] for y in inactive_ids for x in active_ids ])
                    for pdb in structure_ids:
                        #if entry[1] >= 13: # below this distance always inactive structure
                        #    if score < -0.95*lowest_inactive_distance:
                        #        structure_state = "active"
                        #    elif score < 0:
                        #        structure_state = "intermediate"

                        # Classification
                        score = scoring_results[pdb]
                        structure_state = "inactive"
                        if score < 40 and slug[0] == "001": # above this score always inactive structure
                            structure_state = "active"
                            if slug[0] == "001" and score > -70:
                                structure_state = "intermediate"
                        elif score < 0 and slug[0] == "004": # above this score always inactive structure
                            structure_state = "active"
                        elif score < 0 and slug[0] == "006": # above this score always inactive structure
                            structure_state = "active"
                        elif score < 0 and slug[0] not in ["001", "004", "006"]: # above this score always inactive structure
                            structure_state = "active"

                                #print(slug[0], entry[0], structure_state, score)
                            #if slug=="002" and score > -20:
                            #    structure_state = "intermediate"
                            #elif slug=="004" and score > 0:
                            #    structure_state = "intermediate"
                            #elif slug=="006" and score > 0 :
                            #    structure_state = "intermediate"

                        # UGLY: knowledge-based hardcoded corrections
                        if pdb in hardcoded:
                            structure_state = hardcoded[pdb]

                        # Percentage score for TM2-TM6 opening
                        percentage = None
                        if pdb in opening_percentage:
                            percentage = opening_percentage[pdb]

                        # Percentage Gprot-bound likeness
                        gprot_likeness = 100 - int(round((score-min_score)/(max_score-min_score)*100))

                        if pdb in active_ids:
                            gprot_likeness = 100
                            structure_state = "active"
                        elif structure_state == "other":
                            gprot_likeness = None
                            percentage = None

                        #print(slug[0], pdb, score, min_score, max_score, gprot_likeness, structure_state)

                        # Store for structure
                        struct = Structure.objects.get(pdb_code__index=pdb)
                        struct.state, created = ProteinState.objects.get_or_create(slug=structure_state, defaults={'name': structure_state.capitalize()})
                        struct.tm6_angle = percentage
                        struct.gprot_bound_likeness = gprot_likeness
                        struct.save()
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
                        if entry[1] < 13*distance_scaling_factor:
                            struct.state, created = ProteinState.objects.get_or_create(slug="inactive", defaults={'name': "Inactive"})

                        # UGLY: knowledge-based hardcoded corrections
                        if entry[0] in hardcoded:
                            structure_state = hardcoded[entry[0]]
                            struct.state, created = ProteinState.objects.get_or_create(slug=structure_state, defaults={'name': structure_state.capitalize()})

                        # Save changes
                        struct.save()

        self.logger.info('DONE assiging the "Degree Active" levels and activation states')

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
