from django.core.management.base import BaseCommand
from django.contrib.postgres.aggregates import ArrayAgg
from django.db.models import Max, Min

from contactnetwork.distances import Distance, Distances
from contactnetwork.models import distance_scaling_factor
from protein.models import ProteinFamily, ProteinState
from residue.models import Residue
from signprot.models import SignprotComplex
from structure.models import Structure

import logging
import pandas as pd
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd
from statistics import mean

class Command(BaseCommand):
    help = "Evaluate activation states and level of IC \"opening\"."
    logger = logging.getLogger(__name__)

    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--states_only', help='Only (re)determine the activation states of GPCR structures.', default=False, action='store_true')
        parser.add_argument('--representatives_only', help='Only (re)assign the representative tags to GPCR structures.', default=False, action='store_true')

    def handle(self, *args, **options):
        self.logger.info("ASSIGNING the \"Degree Active\" levels and activation states")

        # Loop over classes
        class_slugs = []
        if not options["representatives_only"]:
            class_slugs = list(ProteinFamily.objects.filter(parent__slug="000") \
                                .filter(slug__startswith="00").values_list("slug"))

        for slug in class_slugs:
            print("Processing class {}".format(slug[0]))

            # grab all PDB-codes for this class
            structure_ids = list(Structure.objects.filter(protein_conformation__protein__family__slug__startswith=slug[0]) \
                                .exclude(structure_type__slug__startswith='af-').values_list("pdb_code__index"))

            structure_ids = [x[0] for x in structure_ids]

            ### Skipping class C as author based state annotation is used here instead
            if slug[0].startswith('004'):
                continue

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
                    active_ids.extend(["7C7Q"])

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
                class_pair_inactives['007'] = ["2x46_6x37", 1000] #T PLACEHOLDER

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
                if "7XBX" in active_ids:
                    active_ids.remove("7XBX")
                if "7XBW" in active_ids:
                    active_ids.remove("7XBW")

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
                        "7CA5" : "intermediate",  # Apo state holds middle between active and inactive
                        "7CA3" : "active",  # GABA with PAM groups to active state structures within family, but not class
                        "7M3E" : "inactive", #
                        "7M3J" : "inactive", #
                        "7DD5" : "inactive", #
                        "4IAR" : "inactive", #
                        "4IAQ" : "inactive", #
                        "5V54" : "inactive", #
                        "6IQL" : "inactive", #
                        "7C61" : "inactive", #
                        "7XBX" : "active", # CX3R1 binds to G protein differently with no outwards movement of TM6
                        "7XBW" : "active", # CX3R1 binds to G protein differently with no outwards movement of TM6
                        "6Z66" : "intermediate",
                        "6Z4V" : "intermediate",
                        "6Z8N" : "intermediate",
                        "6ZA8" : "intermediate"
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
                        # Classification
                        score = scoring_results[pdb]
                        structure_state = "inactive"
                        if score < 55 and slug[0] == "001": # above this score always inactive structure
                            structure_state = "active"
                            if slug[0] == "001" and score > -15:
                                structure_state = "intermediate"
                        elif score < -10 and slug[0] == "004": # above this score always inactive structure
                            structure_state = "active"
                        elif score < 0 and slug[0] == "006": # above this score always inactive structure
                            structure_state = "active"
                        elif score < 0 and slug[0] not in ["001", "004", "006"]: # above this score always inactive structure
                            structure_state = "active"

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
                        struct.state = ProteinState.objects.get_or_create(slug=structure_state, defaults={'name': structure_state.capitalize()})[0]
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
                            struct.state = ProteinState.objects.get_or_create(slug="inactive", defaults={'name': "Inactive"})[0]

                        # UGLY: knowledge-based hardcoded corrections
                        if entry[0] in hardcoded:
                            structure_state = hardcoded[entry[0]]
                            struct.state = ProteinState.objects.get_or_create(slug=structure_state, defaults={'name': structure_state.capitalize()})[0]

                        # Save changes
                        struct.save()

        self.logger.info("DONE assiging the \"Degree Active\" levels and activation states")

        if not options["states_only"]:
            # Assigning the representative tag for the "best" structure with the same receptor-state combinations
            self.logger.info("ASSIGNING the \"representative\" tag for unique structure-state complexes")

            # Set the representative state of all GPCR structure to False
            Structure.objects.filter(protein_conformation__protein__family__slug__startswith="00").exclude(structure_type__slug__startswith='af-').update(representative=False)


            # Select all GPCR structures and get unique slug-state combinations
            struct_combs = list(Structure.objects.filter(protein_conformation__protein__family__slug__startswith="00") \
                                .exclude(structure_type__slug__startswith='af-').values_list("protein_conformation__protein__family__slug", "state", "pk", "resolution", "protein_conformation__pk", "protein_conformation__protein__parent__pk"))

            # Grab protein conformations IDs and receptor slugs
            struct_conf_pks = [struct[4] for struct in struct_combs]
            struct_parent_pks = set([struct[5] for struct in struct_combs])

            # Collect unique wildtype GNs per parent protein ID
            wildtype_gns = Residue.objects.filter(protein_conformation__protein_id__in=struct_parent_pks,
                protein_conformation__protein__sequence_type__slug="wt") \
                .exclude(generic_number=None) \
                .values_list("protein_conformation__protein_id") \
                .order_by("protein_conformation__protein_id") \
                .annotate(gns=ArrayAgg("generic_number_id"))

            wildtype_gns_dict = {entry[0]:set(entry[1]) for entry in wildtype_gns}

            # Collect bundle GNs + H8 per structure
            gns_bundle = Residue.objects.filter(protein_conformation_id__in=struct_conf_pks) \
                .filter(protein_segment__slug__in=['TM1','TM2','TM3','TM4','TM5','TM6','TM7','H8']) \
                .exclude(generic_number=None) \
                .values_list("protein_conformation_id") \
                .order_by("protein_conformation_id") \
                .annotate(gns=ArrayAgg("generic_number_id"))

            gns_bundle_dict = {entry[0]:set(entry[1]) for entry in gns_bundle}

            # Collect remaining (mainly loop) GNs per structure
            gns_other = Residue.objects.filter(protein_conformation_id__in=struct_conf_pks) \
                .exclude(protein_segment__slug__in=['TM1','TM2','TM3','TM4','TM5','TM6','TM7','H8']) \
                .exclude(generic_number=None) \
                .values_list("protein_conformation_id") \
                .order_by("protein_conformation_id") \
                .annotate(gns=ArrayAgg("generic_number_id"))

            gns_other_dict = {entry[0]:set(entry[1]) for entry in gns_other}

            # Find and annotate unique combinations with decision criteria
            unique_combs = {}
            for struct in struct_combs:
                combo_id = struct[0] + "_" + str(struct[1])
                if combo_id not in unique_combs:
                    unique_combs[combo_id] = []

                gns_bundle_count = len(wildtype_gns_dict[struct[5]].intersection(gns_bundle_dict[struct[4]]))

                # A couple of structures have no GNs outside of the main bundle
                gns_other_count = 0
                if struct[4] in gns_other_dict:
                    gns_other_count = len(wildtype_gns_dict[struct[5]].intersection(gns_other_dict[struct[4]]))

                unique_combs[combo_id].append([struct[2], gns_bundle_count, gns_other_count, float(struct[3])])

            # Rank and assign representative tag
            for combo in unique_combs:
                if len(unique_combs[combo]) > 1:
                    # Rank by 1. wildtype GNs bundle, 2. # wildtype GNs other, 3. structure resolution
                    unique_combs[combo].sort(
                      key = lambda l: (-1*l[1], -1*l[2], l[3])
                    )

                # Assign "representative" tag to #1
                ref_struct = Structure.objects.get(pk=unique_combs[combo][0][0])
                ref_struct.representative = True
                ref_struct.save()

            self.logger.info("DONE assigning the \"representative\" tag for unique structure-state complexes")

def getDistanceMatrix(node, pdbs):
    if node.is_leaf():
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
