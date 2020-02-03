"""
A set of utility functions for processing distances
"""
from contactnetwork.distances import *

import math
import scipy
import numpy as np
from itertools import combinations

from Bio.SVDSuperimposer import SVDSuperimposer
from numpy.linalg import norm

# convert 3D points to 2D points on best fitting plane
def calculatePlane(tm_points, intracellular):
    # Find best fit plane through set of points (singular value decomposition)
    # Based on https://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points
    centroid = np.sum(tm_points,0) / len(tm_points)
    points = np.transpose(tm_points - centroid) # subtract out the centroid
    svd = np.transpose(np.linalg.svd(points)) # singular value decomposition
    # the corresponding left singular vector is the normal vector of the best-fitting plane
    normal = np.transpose(svd[0])[2]

    return normal, centroid

# convert 3D points to 2D points on best fitting plane
def convert3D_to_2D_plane(tm_points, intracellular, normal, centroid):
    # 3D project points onto plane -> 2D points
    points_plane = [[0]] * 7
    points_z = [0] * 7
    for i in range(0,len(tm_points)):
        # point relative to centroid
        v = tm_points[i] - centroid
        # distance from point to plane along normal vector
        dist = sum(v*normal)
        # point onto plane
        points_plane[i] = tm_points[i] - dist * normal
        points_z[i] = dist

    # Create X and Y axes on plane
    locx = (points_plane[0] - centroid) # TM1 right side
    locy = np.cross(normal, locx)

    locx = locx/norm(locx)
    locy = locy/norm(locy)
    points_2d = [[np.dot(p - centroid, locx), np.dot(p - centroid, locy)] for p in points_plane]

    # Needs flipping?
    if (intracellular and (points_2d[3][1]-points_2d[0][1]) < 0) or (not intracellular and (points_2d[3][1]-points_2d[0][1]) > 0):
        points_2d = np.array([[x[0], -1*x[1]] for x in points_2d])

    return np.array(points_2d), points_z

# Find best fit plane through set of points (least squares)
# Based on https://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points
def convert2D_PLS(tm_points):
    tmp_A = []
    tmp_b = []
    for i in range(len(tm_points)):
        tmp_A.append([tm_points[i][0], tm_points[i][1], 1])
        tmp_b.append(tm_points[i][2])

    b = np.matrix(tmp_b).T
    A = np.matrix(tmp_A)
    fit = (A.T * A).I * A.T * b
    #errors = b - A * fit
    #residual = norm(errors)
    return fit

def determineOrder(distance_matrix):
    ref_distances = np.triu(np.array([np.array(x) for x in distance_matrix]))
    ranking = []
    while len(ranking) < 7:
        max_ref = np.unravel_index(ref_distances.argmax(), ref_distances.shape)
        # highest = ref_distances[max_ref[0]][max_ref[1]]

        ref_distances[max_ref[0]][max_ref[1]] = 0
        filtered_ref = [ x for x in max_ref if x not in ranking ]
        ranking.extend(filtered_ref)
        # print(ranking, max_ref, highest)

    return ranking


# 3D reconstruction based on average distances
def recreate3Dorder(distance_matrix):
    to = determineOrder(distance_matrix)

    tms = [[0]] * 7
    tms[0] = np.array([0, 0, 0])
    tms[1] = np.array([distance_matrix[to[0]][to[1]], 0, 0])

    # place TM3 - relative to TM1 and TM2
    d = tms[1][0]
    a = (math.pow(distance_matrix[to[0]][to[2]],2) - math.pow(distance_matrix[to[1]][to[2]],2) + math.pow(d,2)) / (2*d)
    h = math.sqrt(math.pow(distance_matrix[to[0]][to[2]],2) - math.pow(a,2))
    x3 = a+h*tms[1][1]/d
    y3 = abs(0-h*tms[1][0]/d)

    tms[2] = np.array([x3, y3, 0])
    for i in range(3,7):
        #print("calculating for TM", str(i+1))
        if i == 3:
            # just take the first solution (mirrored solution)
            tms[i] = trilaterate(tms[0], tms[1], tms[2], distance_matrix[to[0]][to[i]], distance_matrix[to[1]][to[i]], distance_matrix[to[2]][to[i]])[0]
        else:
            if False:
                # Select most distant references to this point (for more accurate placement)
                ref_distances = [distance_matrix[to[i]][x] for x in to[:i]]
                distance_order = [sorted(ref_distances, reverse=True).index(x) for x in ref_distances]
                selected_references = [distance_order.index(0), distance_order.index(1), distance_order.index(2)]
                validation_point = distance_order.index(3)

                # Place point using the identified references
                sr = selected_references
                ref_dist = distance_matrix[to[validation_point]][to[i]]
                tms[i] = trilaterate(tms[sr[0]], tms[sr[1]], tms[sr[2]], distance_matrix[to[sr[0]]][to[i]], distance_matrix[to[sr[1]]][to[i]], distance_matrix[to[sr[2]]][to[i]])
                changes = [abs(np.linalg.norm(tms[i][j]-tms[validation_point]) - ref_dist) for j in range(0,2)]
                tms[i] = tms[i][changes.index(min(changes))]
            else:
                # Use all reference to optimize point placement
                reference_tms = range(0,i)
                reference_points = [ comb for comb in combinations(reference_tms, 3)]

                ref_worked = 0
                tms[i] = np.array([0,0,0], dtype='f')
                for x in reference_points:
                    # Find TM not in list
                    check_tm = [j for j in reference_tms if j not in x][0]

                    # Find best coordinate matching the fourth reference point
                    ref_dist = distance_matrix[to[check_tm]][to[i]]
                    scaler = 1
                    while True:
                        try:
                            combi_pos = trilaterate(tms[x[0]], tms[x[1]], tms[x[2]], distance_matrix[to[x[0]]][to[i]]*scaler, distance_matrix[to[x[1]]][to[i]]*scaler, distance_matrix[to[x[2]]][to[i]]*scaler)
                            break
                        except:
                            scaler += 0.01
                            # print("Scaling up", scaler)
                            if scaler > 1.1: # Way off, only if DB data is wrong
                                break
                            pass

                    changes = [abs(np.linalg.norm(combi_pos[j]-tms[check_tm])-ref_dist) for j in range(0,2)]
                    lowest = changes.index(min(changes))
                    tms[i] += combi_pos[lowest]/len(reference_points)

    # Rearrange to correct order
    tms = [tms[to.index(i)] for i in range(0,7)]

    # DEBUG distances
    # for i in range(0,6):
    #     for j in range(i+1, 7):
    #         ref_dist = distance_matrix[i][j]
    #         print (i+1,j+1, round(np.linalg.norm(tms[i] - tms[j]),3), round(ref_dist,3), round(np.linalg.norm(tms[i] - tms[j]) - ref_dist,3))

    return np.array(tms)

# Find the intersection of three spheres
# Based on https://stackoverflow.com/questions/1406375/finding-intersection-points-between-3-spheres
# P1,P2,P3 are the centers, r1,r2,r3 are the radii
# Implementaton based on Wikipedia Trilateration article.
def trilaterate(P1,P2,P3,r1,r2,r3):
    temp1 = P2-P1
    e_x = temp1/norm(temp1)
    temp2 = P3-P1
    i = np.dot(e_x,temp2)
    temp3 = temp2 - i*e_x
    e_y = temp3/norm(temp3)
    e_z = np.cross(e_x,e_y)
    d = norm(P2-P1)
    j = np.dot(e_y,temp2)
    x = (r1*r1 - r2*r2 + d*d) / (2*d)
    y = (r1*r1 - r3*r3 -2*i*x + i*i + j*j) / (2*j)
    temp4 = r1*r1 - x*x - y*y
    if temp4<0:
        raise Exception("The three spheres do not intersect!");
    z = math.sqrt(temp4)
    p_12_a = P1 + x*e_x + y*e_y + z*e_z
    p_12_b = P1 + x*e_x + y*e_y - z*e_z
    return p_12_a,p_12_b

def tm_movement_2D(pdbs1, pdbs2, intracellular, data, gn_dictionary):
    print("COMPARISON")
    print(pdbs1)
    print("VS")
    print(pdbs2)

    distances_set1 = Distances()
    distances_set1.load_pdbs(pdbs1)
    distances_set1.filtered_gns = True

    distances_set2 = Distances()
    distances_set2.load_pdbs(pdbs2)
    distances_set2.filtered_gns = True

    conserved_set1 = distances_set1.fetch_conserved_gns_tm()
    conserved_set2 = distances_set2.fetch_conserved_gns_tm()

    conserved = [x for x in conserved_set2 if x in conserved_set1]

    gns = [[]] * 7
    for i in range(0,7):
        tm_only = [x for x in conserved if x[0]==str(i+1)]
        if intracellular and i % 2 == 0: #all uneven TMs (as # = i+1)
            tm_only.reverse()
        elif not intracellular and i % 2 == 1: # all even TMs (as # i+1)
            tm_only.reverse()
        if len(tm_only) < 3:
            print("too few residues")
            return []
        gns[i] = tm_only[0:3]
        # gns[i] = tm_only[0:2]

    # INCLUDING References points from membrane middle of GPCR
    ref_membrane_mid = {}
    ref_membrane_mid["001"] = [['1x43', '1x44','1x45'], ['2x51', '2x52','2x53'], ['3x35', '3x36', '3x37'], ['4x53', '4x54', '4x55'], ['5x45', '5x46', '5x47'], ['6x47', '6x48', '6x49'], ['7x42', '7x43', '7x44']]
    ref_membrane_mid["002"] = [['1x50', '1x51', '1x52'], ['2x57', '2x58', '2x59'], ['3x40','3x41','3x42'], ['4x53', '4x54', '4x55'], ['5x44', '5x45', '5x46'], ['6x48', '6x49', '6x50'], ['7x49', '7x50', '7x51']]
    ref_membrane_mid["003"] = ref_membrane_mid["002"]
    ref_membrane_mid["004"] = [['1x48', '1x49', '1x50'], ['2x47', '2x48', '2x49'], ['3x39', '3x40', '3x41'], ['4x40', '4x41', '4x42'], ['5x47', '5x48', '5x49'], ['6x47', '6x48', '6x49'], ['7x39', '7x40', '7x41']]
    ref_membrane_mid["005"] = [['1x42', '1x43', '1x44'], ['2x52', '2x53', '2x54'], ['3x37', '3x38', '3x39'], ['4x52', '4x53', '4x54'], ['5x52', '5x53', '5x54'], ['6x42', '6x43', '6x44'], ['7x46', '7x47', '7x48']]

    # FILTER not conserved GNs
    middle_gpcr = [list(filter(lambda x: x in conserved, tm_list)) for tm_list in ref_membrane_mid[data['gpcr_class']]]
    ends_and_middle = gns[:]
    ends_and_middle.extend(middle_gpcr)


    distances_set1.filter_gns.extend([y for x in ends_and_middle for y in x])
    distances_set2.filter_gns = distances_set1.filter_gns
    distances_set1.fetch_distances_tm()
    distances_set2.fetch_distances_tm()

    membrane_data1 = [x[:] for x in [[0] * 14] * 14]
    membrane_data2 = [x[:] for x in [[0] * 14] * 14]
    for i in range(0,13):
        for j in range(i+1, 14):
            filter_keys = [x+"_"+y if right_gn_order(x,y) else y+"_"+x for x in ends_and_middle[i] for y in ends_and_middle[j] if x != y]
            if len(filter_keys) == 0:
                print("no filter keys")
                return []
            membrane_data1[i][j] = sum([k for x in filter_keys for k in distances_set1.data[x]])/(len(filter_keys)*len(pdbs1))
            membrane_data1[j][i] = membrane_data1[i][j]
            membrane_data2[i][j] = sum([k for x in filter_keys for k in distances_set2.data[x]])/(len(filter_keys)*len(pdbs2))
            membrane_data2[j][i] = membrane_data2[i][j]

    # Identify most stable TMs by ranking the variations to all other helices
    diff_distances = [x[:] for x in [[0] * len(ends_and_middle)] * len(ends_and_middle)]
    # real_differences = [x[:] for x in [[0] * 7] * 7]
    for i in range(0,13):
        for j in range(i+1, 14):
            # Calculate movements for each TM relative to their "normal" distance
            diff_distances[i][j] = abs(membrane_data1[i][j] - membrane_data2[i][j])/((membrane_data1[i][j]+membrane_data2[i][j])/2)*100
            diff_distances[j][i] = diff_distances[i][j]
            # real_differences[i][j] = membrane_data2[i][j] - membrane_data1[i][j]
            # real_differences[j][i] = real_differences[i][j]

    # Ranking for each TM
    sum_differences = [sum(x) for x in diff_distances]
    normalized_differences = [((sum_differences[i]-min(sum_differences[0:7]))/(max(sum_differences[0:7])-min(sum_differences[0:7])))**2 for i in range(0,7)]
#    print("Measured diff")
    for i in range(0,7):
#        print(diff_distances[i])
        # Real variations
        diff_distances[i] = [sorted(diff_distances[i]).index(x) for x in diff_distances[i]]
        #diff_distances[i] = [sorted(diff_distances[i]).index(x) for j in range(0,7) for x in diff_distances[i]]
        #print(i,diff_distances[i], real_differences[i])

    # TODO Simplify
    final_rank = [0] * 7
    for i in range(0,7):
        final_rank[i] = sum([diff_distances[j][i] for j in range(0,7)])

    # Grab stable TMs
    tm_ranking = [0] * 7
    sorted_rank = sorted(final_rank)
    for i in range(0,7):
        tm_ranking[i] = final_rank.index(sorted_rank[i])
        final_rank[tm_ranking[i]] = 100 # make sure this TM isn't repeated

    stable_one, stable_two = tm_ranking[0], tm_ranking[1]
    print("TM stable ranking")
    print(tm_ranking)


    # Possible workaround for two consecutive helices if combined with rotation/translation
    # if abs(stable_two-stable_one) == 1 or abs(stable_two-stable_one) == 6:
    #     if abs(tm_ranking[2]-stable_one) > 1:
    #         stable_two = tm_ranking[2]
    #     else:
    #         stable_one = tm_ranking[2]
    # if stable_one > stable_two:
    #     stable_one, stable_two = stable_two, stable_one

    # Recalculate 3D network and populate 2D views
    # Note: rebuilding in order stable helices resulted in lower quality reconstruction, therefore back to TM order (1-7)
    #tms_set1 = recreate3Dorder([x[:7] for x in membrane_data1[:7]], range(0,7))
    tms_set1 = recreate3Dorder([x[:7] for x in membrane_data1[:7]])
    # tms_set1 = recreate3Dorder([x[:7] for x in membrane_data1[:7]], tm_ranking)
    #tms_set2 = recreate3Dorder([x[:7] for x in membrane_data2[:7]], range(0,7))
    tms_set2 = recreate3Dorder([x[:7] for x in membrane_data2[:7]])
    # tms_set2 = recreate3Dorder([x[:7] for x in membrane_data2[:7]], tm_ranking)

    #normal_set1, mid_set1 = calculatePlane(tms_set1, intracellular)
    normal_set1, mid_set1 = calculatePlane(np.concatenate((tms_set1, tms_set2)), intracellular)
    #plane_set1, z_set1, normal_set1, mid_set1 = convert2D_SVD(tms_set1, intracellular)
    plane_set1, z_set1 = convert3D_to_2D_plane(tms_set1, intracellular, normal_set1, mid_set1)

    # Align 3D points of set2 with 3D points of set1
    imposer = SVDSuperimposer()
    imposer.set(tms_set1[tm_ranking[0:3]], tms_set2[tm_ranking[0:3]])
    imposer.run()
    print("RMSD")
    print(imposer.get_rms())
    rot, trans = imposer.get_rotran()

    # Based on Biopython SVDSuperimposer
    # coords = tms_set2
    # reference_coords = tms_set1

    # OLD centroid calcalation
    # av1 = sum(coords) / len(coords)
    # av2 = sum(reference_coords) / len(reference_coords)

    # NEW weighted centroid calculation
    # print(normalized_differences)
    # av1, av2 = 0, 0
    # totalweight = 0
    # for i in range(0,7):
    #     # print("Round",i)
    #     #weight = 1+(7-tm_ranking.index(i))/7
    #     weight = (1-normalized_differences[i]+0.1)/1.1
    #     totalweight += weight
    #     print("TM", str(i+1), "weight",weight)
    #     av1 += coords[i]*weight
    #     av2 += reference_coords[i]*weight
    #
    # av1 = av1/totalweight
    # av2 = av2/totalweight
    #
    # coords = coords - av1
    # reference_coords = reference_coords - av2
    #
    # # correlation matrix
    # a = np.dot(np.transpose(coords), reference_coords)
    # u, d, vt = np.linalg.svd(a)
    # rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
    # # check if we have found a reflection
    # if np.linalg.det(rot) < 0:
    #     vt[2] = -vt[2]
    #     rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
    # trans = av2 - np.dot(av1, rot)

    # Apply
    tms_set2 = np.dot(tms_set2, rot) + trans

    # Convert the 3D points of set2 to 2D for the plane of set1
    plane_set2, z_set2 = convert3D_to_2D_plane(tms_set2, intracellular, normal_set1, mid_set1)

    # Add angles
    rotations = [0] * 7
    for i in range(0,7):
        # Ca-angle to axis core

        rotations[i] = [data['tab4'][gn_dictionary[x]]['angles_set1'][1]-data['tab4'][gn_dictionary[x]]['angles_set2'][1] if abs(data['tab4'][gn_dictionary[x]]['angles_set1'][1]-data['tab4'][gn_dictionary[x]]['angles_set2'][1]) < 180 else -1*data['tab4'][gn_dictionary[x]]['angles_set2'][1]-data['tab4'][gn_dictionary[x]]['angles_set1'][1] for x in gns[i]]
        if intracellular:
            rotations[i] = -1*sum(rotations[i])/3
        else:
            rotations[i] = sum(rotations[i])/3

    # TODO: check z-coordinates orientation
    # Step 1: collect movement relative to membrane mid
    # Step 2: find min and max TM
    # Step 3: check if orientation of min/max TM matches the z-scales + intra/extra - if not invert z-coordinates
    labeled_set1 = [{"label": "TM"+str(i+1), "x": float(plane_set1[i][0]), "y": float(plane_set1[i][1]), "z": float(z_set1[i]), "rotation" : 0} for i in range(0,7)]
    labeled_set2 = [{"label": "TM"+str(i+1), "x": float(plane_set2[i][0]), "y": float(plane_set2[i][1]), "z": float(z_set2[i]), "rotation" : rotations[i]} for i in range(0,7)]

    return {"coordinates_set1" : labeled_set1, "coordinates_set2": labeled_set2, "gns_used": gns}


def right_gn_order(pos1, pos2):
    res1 = pos1.split("x")
    res2 = pos2.split("x")
    return int(res1[0]) < int(res2[0]) or (res1[0] == res2[0] and int(res1[1]) < int(res2[1]))
