"""
A set of utility functions for processing distances
"""
from contactnetwork.distances import *

import math
import random
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
    points_plane = [[0]] * len(tm_points)
    points_z = [0] * len(tm_points)
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
    points_2d = [np.array([np.dot(p - centroid, locx), np.dot(p - centroid, locy)]) for p in points_plane]

    # Needs flipping?
    #if (intracellular and (points_2d[3][1]-points_2d[0][1]) < 0) or (not intracellular and (points_2d[3][1]-points_2d[0][1]) > 0):
    #    points_2d = np.array([np.array([x[0], -1*x[1]]) for x in points_2d])
    # UPDATE 21-02-2020 - no mirroring for intracellular side, top-down slicing through GPCR
    if (intracellular and (points_2d[3][1]-points_2d[0][1]) > 0) or (not intracellular and (points_2d[3][1]-points_2d[0][1]) > 0):
        points_2d = np.array([np.array([x[0], -1*x[1]]) for x in points_2d])

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

def determine_order(distance_matrix):
    ref_distances = np.triu(np.array([np.array(x) for x in distance_matrix]))
    ranking = []
    while len(ranking) < 7:
        max_ref = np.unravel_index(ref_distances.argmax(), ref_distances.shape)

        ref_distances[max_ref[0]][max_ref[1]] = 0
        filtered_ref = [ x for x in max_ref if x not in ranking ]
        ranking.extend(filtered_ref)

    return ranking

def determine_order_group(distance_matrix, gn_grouping):
    first_index_group = [gn_grouping.index(x) for x in range(0,max(gn_grouping)+1)]
    ref_distances = np.triu(np.array([np.array(x)[first_index_group] for x in distance_matrix])[first_index_group])

    ranking = []
    while len(ranking) <= max(gn_grouping):
        max_ref = np.unravel_index(ref_distances.argmax(), ref_distances.shape)

        ref_distances[max_ref[0]][max_ref[1]] = 0
        filtered_ref = [ x for x in max_ref if x not in ranking ]
        ranking.extend(filtered_ref)

    ranking = [first_index_group[i] for i in ranking]
    ranking_groups = ranking[:]
    for looping in range(0,2):
        for i in ranking:
            add = [x for x in range(0,len(gn_grouping)) if gn_grouping[x] == gn_grouping[i] and x not in ranking_groups]
            if len(add) > 0:
                ranking_groups.append(add[0])

    return ranking_groups

def consecutive_group_order(gn_grouping):
    first_index_group = [gn_grouping.index(x) for x in range(0,max(gn_grouping)+1)]
    ranking = [first_index_group[i] for i in range(0,max(gn_grouping)+1)]
    ranking_groups = ranking[:]
    for looping in range(0,2):
        for i in ranking:
            add = [x for x in range(0,len(gn_grouping)) if gn_grouping[x] == gn_grouping[i] and x not in ranking_groups]
            if len(add) > 0:
                ranking_groups.append(add[0])

    return ranking_groups

def reconstruction_error(distance_matrix, tm_points):
    error = sum([abs(np.linalg.norm(tm_points[i] - tm_points[j]) - distance_matrix[i][j]) for i in range(0,len(tm_points)-1) for j in range(i+1, len(tm_points))])
    return error, error/(len(tm_points)*(len(tm_points)-1))

# 3D reconstruction based on average distances
def recreate3Dorder(distance_matrix, gn_grouping, mirror = False):
    # Reorder with respect to distances
    best_set = []
    best_error = [100000, 100000]

    #to = determine_order_group(distance_matrix, gn_grouping) # based on distance
    to_ref = consecutive_group_order(gn_grouping) # based on groups - same initial four points, same plane
    for repeat in range(0,10):
        to = to_ref[:]
        skip_this_loop = False

        if repeat > 0:
            random.seed(repeat)
            random.shuffle(to)
            #print("Randomized reconstruction order", to)

        reorder_dist = np.array([np.array(x)[to] for x in distance_matrix])[to]

        tms = [[0]] * len(gn_grouping)
        tms[0] = np.array([0, 0, 0])
        tms[1] = np.array([reorder_dist[0][1], 0, 0])

        # place TM3 - relative to TM1 and TM2
        d = tms[1][0]
        a = (math.pow(reorder_dist[0][2],2) - math.pow(reorder_dist[1][2],2) + math.pow(d,2)) / (2*d)
        h = math.sqrt(math.pow(reorder_dist[0][2],2) - math.pow(a,2))
        x3 = a+h*tms[1][1]/d
        y3 = abs(0-h*tms[1][0]/d)

        tms[2] = np.array([x3, y3, 0])
        for i in range(3,len(gn_grouping)):
            #print("calculating for TM", str(i+1))
            if i == 3:
                # just take the first solution (mirrored solution)
                try:
                    if mirror:
                        tms[i] = trilaterate(tms[0], tms[1], tms[2], reorder_dist[0][i], reorder_dist[1][i], reorder_dist[2][i])[1]
                    else:
                        tms[i] = trilaterate(tms[0], tms[1], tms[2], reorder_dist[0][i], reorder_dist[1][i], reorder_dist[2][i])[0]
                except:
                    skip_this_loop = True
                    break
            else:
                # Alternative 1: Place point using just the previous three points
                # sr = [i-3,i-2,i-1]
                # ref_dist = reorder_dist[i-4][i]
                # scaler = 1
                # while True:
                #     try:
                #         tms[i] = trilaterate(tms[sr[0]], tms[sr[1]], tms[sr[2]], reorder_dist[sr[0]][i]*scaler, reorder_dist[sr[1]][i]*scaler, reorder_dist[sr[2]][i]*scaler)
                #         changes = [abs(np.linalg.norm(tms[i][j]-tms[i-4]) - ref_dist) for j in range(0,2)]
                #         print("CHANGES", i, changes)
                #         tms[i] = tms[i][changes.index(min(changes))]
                #         break
                #     except:
                #         scaler += 0.01
                #         print("Scaling up", scaler)
                #         if scaler > 1.1: # Way off, only if DB data is wrong (TODO: throw and handle error )
                #             break
                #         pass

                # Alternative 2: Place point using just 1-3
                sr = [0,1,2]
                ref_dist = reorder_dist[3][i]
                scaler = 1
                while True:
                    try:
                        tms[i] = trilaterate(tms[sr[0]], tms[sr[1]], tms[sr[2]], reorder_dist[sr[0]][i]*scaler, reorder_dist[sr[1]][i]*scaler, reorder_dist[sr[2]][i]*scaler)
                        changes = [abs(np.linalg.norm(tms[i][j]-tms[3]) - ref_dist) for j in range(0,2)]
                        # print("CHANGES", i, changes)
                        tms[i] = tms[i][changes.index(min(changes))]
                        break
                    except:
                        scaler += 0.01
                        # print("Scaling up", scaler)
                        if scaler > 1.1: # Way off, only if DB data is wrong (TODO: throw and handle error )
                            skip_this_loop = True
                            break
                        pass


                # Alternative 3: place point using the most distant references (should minimize placement rounding error)
                # ref_distances = [reorder_dist[i][x] for x in range(0,min(i,max(gn_grouping)+1))]
                # distance_order = [sorted(ref_distances, reverse=True).index(x) for x in ref_distances[:(max(gn_grouping)+1)]]
                #
                # sr = [distance_order.index(0), distance_order.index(1), distance_order.index(2)]
                # ref_point = distance_order.index(3)
                # try:
                #     # tms[i] = trilaterate(tms[sr[0]], tms[sr[1]], tms[sr[2]], reorder_dist[sr[0]][i]*scaler, reorder_dist[sr[1]][i]*scaler, reorder_dist[sr[2]][i]*scaler)
                #     tms[i] = trilaterate(tms[sr[0]], tms[sr[1]], tms[sr[2]], reorder_dist[sr[0]][i], reorder_dist[sr[1]][i], reorder_dist[sr[2]][i])
                # except:
                #     print("Using backup scenario") # Scaler could also be used, increased error
                #     sr = [0,1,2]
                #     ref_point = 3
                #     tms[i] = trilaterate(tms[sr[0]], tms[sr[1]], tms[sr[2]], reorder_dist[sr[0]][i], reorder_dist[sr[1]][i], reorder_dist[sr[2]][i])
                #
                # ref_dist = reorder_dist[ref_point][i]
                # changes = [abs(np.linalg.norm(tms[i][j]-tms[ref_point]) - ref_dist) for j in range(0,2)]
                # print("CHANGES", i, changes)
                # tms[i] = tms[i][changes.index(min(changes))]

        # Verify if this loop could reconstruct all points
        if skip_this_loop:
            continue

        # Rearrange to correct order
        tms = [tms[to.index(i)] for i in range(0,len(gn_grouping))]

        # # DEBUG residue-pair distances
        # for i in range(0,len(gn_grouping)-1):
        #     for j in range(i+1, len(gn_grouping)):
        #         ref_dist = distance_matrix[i][j]
        #         print (i+1,j+1, round(np.linalg.norm(tms[i] - tms[j]),3), round(ref_dist,3), round(np.linalg.norm(tms[i] - tms[j]) - ref_dist,3))

        total_error, point_error = reconstruction_error(distance_matrix,tms)
        #print(repeat, "Total error for", len(tms), "points:", total_error, "averaging", round(point_error, 4), "Å per distance")
        if total_error < best_error[0]:
            best_set = np.array(tms)
            best_error = [total_error, point_error]
            if point_error < 0.1:
                break

    print("Total error for", len(tms), "points:", best_error[0], "averaging", round(best_error[1], 4), "Å per distance")
    print("------")

    tms = best_set

    # Create centroids for points in the same group
    tm_centroids = [[0.0,0.0,0.0]] * (max(gn_grouping)+1)
    for group in range(0, max(gn_grouping)+1):
        # Grab points from same group
        ref_points = [x for x in range(0,len(gn_grouping)) if gn_grouping[x] == group]
        for ref_point in ref_points:
            tm_centroids[group] += np.array(tms[ref_point])/np.float(len(ref_points))

    # DEBUG centroid distances
    # for i in range(0,len(tm_centroids)-1):
    #     for j in range(i+1, len(tm_centroids)):
    #         print (i+1,j+1, round(np.linalg.norm(tm_centroids[i] - tm_centroids[j]),3))

    return np.array(tm_centroids), np.array(tms)

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


# mode => 0 - extracellular, 1 intracellular, 2 major pocket (class A), 3 middle
def tm_movement_2D(pdbs1, pdbs2, mode, data, gn_dictionary):
    string_mode = ["extracellular", "intracellular", "pocket", "middle"]
    intracellular = (mode == 1)
    print("COMPARISON", string_mode[mode])
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
    middle_gpcr = [[]] * 7
    if mode <= 1: # Intracellular or Extracellular
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

            for upwards in range(12, 6, -1):
                if len(tm_only) >= upwards:
                    middle_gpcr[i] = tm_only[(upwards-3):upwards]
                    break

        # INCLUDING References points from membrane middle of GPCR
        # ref_membrane_mid = {}
        # ref_membrane_mid["001"] = [['1x43', '1x44','1x45'], ['2x51', '2x52','2x53'], ['3x35', '3x36', '3x37'], ['4x53', '4x54', '4x55'], ['5x45', '5x46', '5x47'], ['6x47', '6x48', '6x49'], ['7x42', '7x43', '7x44']] # A
        # #ref_membrane_mid["002"] = [['1x50', '1x51', '1x52'], ['2x57', '2x58', '2x59'], ['3x40','3x41','3x42'], ['4x53', '4x54', '4x55'], ['5x44', '5x45', '5x46'], ['6x48', '6x49', '6x50'], ['7x49', '7x50', '7x51']] # B1
        # ref_membrane_mid["002"] = [['1x50', '1x51', '1x52'], ['2x57', '2x58', '2x59'], ['3x40','3x41','3x42'], ['4x55', '4x56'], ['5x42', '5x43', '5x44'], ['7x47', '7x49']] # B1
        # ref_membrane_mid["003"] = ref_membrane_mid["002"] # B2
        # ref_membrane_mid["004"] = [['1x48', '1x49', '1x50'], ['2x47', '2x48', '2x49'], ['3x39', '3x40', '3x41'], ['4x40', '4x41', '4x42'], ['5x47', '5x48', '5x49'], ['6x47', '6x48', '6x49'], ['7x39', '7x40', '7x41']] # C
        # ref_membrane_mid["006"] = [['1x42', '1x43', '1x44'], ['2x52', '2x53', '2x54'], ['3x37', '3x38', '3x39'], ['4x52', '4x53', '4x54'], ['5x52', '5x53', '5x54'], ['6x42', '6x43', '6x44'], ['7x46', '7x47', '7x48']] # F
        #
        # middle_gpcr = ref_membrane_mid[data['gpcr_class']]
    elif mode == 2: # Major pocket (class A)
        ligand_references = [['1x39', '1x40','1x41'], ['2x56', '2x57','2x58'], ['3x31', '3x32', '3x33'], ['4x56', '4x57', '4x58'], ['5x43', '5x44', '5x45'], ['6x51', '6x52', '6x53'], ['7x39', '7x40', '7x41']]
        for i in range(0,7):
            gns[i] = [x for x in ligand_references[i] if x in conserved]
            tm_only = [x for x in conserved if x[0]==str(i+1)]
            if i % 2 == 1: #all uneven TMs (as # = i+1)
                tm_only.reverse()
            if len(gns[i]) > 0:
                if i % 2 == 1: #all uneven TMs (as # = i+1)
                    start_pos = tm_only.index(gns[i][-1])
                else:
                    start_pos = tm_only.index(gns[i][0])

                gns[i] = tm_only[start_pos:(start_pos+3)]

                # Stay close for this as references
                #middle_gpcr[i] = tm_only[(start_pos+6):(start_pos+9)]
                for upwards in range(9, 6, -1):
                   if len(tm_only) >= (start_pos+upwards):
                       middle_gpcr[i] = tm_only[(start_pos+upwards-3):(start_pos+upwards)]
                       continue
            else:
                if len(tm_only) < 9:
                    print("too few residues")
                    return []
                else:
                    #print("Refind",i, gns[i])
                    gns[i] = tm_only[0:3]
                    middle_gpcr[i] = tm_only[6:9]

                    # for upwards in range(15, 6, -1):
                    #     if len(tm_only) >= upwards:
                    #         middle_gpcr[i] = tm_only[(upwards-3):upwards]

        # # FILTER not conserved GNs
        # middle_gpcr = [[]] * 7
        # for i in range(0,7):
        #     tm_only = [x for x in conserved if x[0]==str(i+1)]
        #     if i % 2 == 0: #all uneven TMs (as # = i+1)
        #         tm_only.reverse()
        #
        #     if len(tm_only) < 3:
        #         print("too few residues")
        #         return []
        #
        #     middle_gpcr[i] = tm_only[0:3]
        #print(middle_gpcr)

    elif mode == 3: # Middle
        # References points from membrane middle of GPCR
        ref_membrane_mid = {}
        ref_membrane_mid["001"] = [['1x43', '1x44','1x45'], ['2x51', '2x52','2x53'], ['3x35', '3x36', '3x37'], ['4x53', '4x54', '4x55'], ['5x45', '5x46', '5x47'], ['6x47', '6x48', '6x49'], ['7x42', '7x43', '7x44']] # A
        #ref_membrane_mid["002"] = [['1x50', '1x51', '1x52'], ['2x57', '2x58', '2x59'], ['3x40','3x41','3x42'], ['4x53', '4x54', '4x55'], ['5x44', '5x45', '5x46'], ['6x48', '6x49', '6x50'], ['7x49', '7x50', '7x51']] # B1
        #ref_membrane_mid["002"] = [['1x50', '1x51', '1x52'], ['2x57', '2x58', '2x59'], ['3x40','3x41','3x42'], ['4x55', '4x56'], ['5x42', '5x43', '5x44'], ['7x47', '7x49']] # B1
        ref_membrane_mid["002"] = [['1x50', '1x51', '1x52'], ['2x57', '2x58', '2x59'], ['3x40','3x41','3x42'], ['4x55', '4x56'], ['5x42', '5x43', '5x44'], ['6x48', '6x49', '6x50'], ['7x47', '7x49']] # B1
        ref_membrane_mid["003"] = ref_membrane_mid["002"] # B2
        ref_membrane_mid["004"] = [['1x48', '1x49', '1x50'], ['2x47', '2x48', '2x49'], ['3x39', '3x40', '3x41'], ['4x40', '4x41', '4x42'], ['5x47', '5x48', '5x49'], ['6x47', '6x48', '6x49'], ['7x39', '7x40', '7x41']] # C
        ref_membrane_mid["006"] = [['1x42', '1x43', '1x44'], ['2x52', '2x53', '2x54'], ['3x37', '3x38', '3x39'], ['4x52', '4x53', '4x54'], ['5x52', '5x53', '5x54'], ['6x42', '6x43', '6x44'], ['7x46', '7x47', '7x48']] # F

        membrane_mid = ref_membrane_mid[data['gpcr_class']]

        if data['gpcr_class'] != "001":
            inv_gn_dictionary = {v: k for k, v in gn_dictionary.items()}
            for index in range(len(membrane_mid)):
                membrane_mid[index] = [inv_gn_dictionary[res] for res in membrane_mid[index]]

        for i in range(0,7):
            gns[i] = [x for x in membrane_mid[i] if x in conserved]
            tm_only = [x for x in conserved if x[0]==str(i+1)]
            if i % 2 == 1: #all uneven TMs (as # = i+1)
                tm_only.reverse()
            if len(gns[i]) > 0:
                if i % 2 == 1: #all uneven TMs (as # = i+1)
                    start_pos = tm_only.index(gns[i][-1])
                else:
                    start_pos = tm_only.index(gns[i][0])

                gns[i] = tm_only[start_pos:(start_pos+3)]

                # Stay close for this as references
                #middle_gpcr[i] = tm_only[(start_pos+6):(start_pos+9)]
                for upwards in range(6, 3, -1):
                   if len(tm_only) >= (start_pos+upwards):
                       middle_gpcr[i] = tm_only[(start_pos+upwards-3):(start_pos+upwards)]
                       continue
            else:
                if len(tm_only) < 6:
                    print("too few residues")
                    return []
                else:
                    #print("Refind",i, gns[i])
                    gns[i] = tm_only[0:3]
                    middle_gpcr[i] = tm_only[3:6]

                    # for upwards in range(15, 6, -1):
                    #     if len(tm_only) >= upwards:
                    #         middle_gpcr[i] = tm_only[(upwards-3):upwards]

    # Merge the reference and the helper points
    gns_flat = [y for x in gns for y in x]
    middle_gpcr = [list(filter(lambda x: x in conserved and x not in gns_flat, tm_list)) for tm_list in middle_gpcr]
    # print(gns)
    # print(middle_gpcr)

    ends_and_middle = gns[:]
    ends_and_middle.extend(middle_gpcr)
    ends_and_middle_flat = [y for x in ends_and_middle for y in x]
    ends_and_middle_grouping = [x for x in range(0, len(ends_and_middle)) for y in ends_and_middle[x]]
    segment_order = [int(ends_and_middle[x][0][0])-1 for x in range(0, len(ends_and_middle))]

    distances_set1.filter_gns.extend([y for x in ends_and_middle for y in x])
    distances_set2.filter_gns = distances_set1.filter_gns
    distances_set1.fetch_distances_tm(distance_type = "HC")
    distances_set2.fetch_distances_tm(distance_type = "HC")


    membrane_data1 = [x[:] for x in [[0] * len(ends_and_middle_flat)] * len(ends_and_middle_flat)]
    membrane_data2 = [x[:] for x in [[0] * len(ends_and_middle_flat)] * len(ends_and_middle_flat)]
    for i in range(0,len(ends_and_middle_flat)-1):
        for j in range(i+1, len(ends_and_middle_flat)):
            if right_gn_order(ends_and_middle_flat[i], ends_and_middle_flat[j]):
                filter_key = ends_and_middle_flat[i] + "_" + ends_and_middle_flat[j]
            else:
                filter_key = ends_and_middle_flat[j] + "_" + ends_and_middle_flat[i]

            if ends_and_middle_flat[i] != ends_and_middle_flat[j]:
                membrane_data1[i][j] = sum(distances_set1.data[filter_key])/len(pdbs1)
                membrane_data1[j][i] = membrane_data1[i][j]
                membrane_data2[i][j] = sum(distances_set2.data[filter_key])/len(pdbs2)
                membrane_data2[j][i] = membrane_data2[i][j]

    # Identify most stable TMs by ranking the variations to all other helices
    membrane_data1 = np.array([np.array(x) for x in membrane_data1])
    membrane_data2 = np.array([np.array(x) for x in membrane_data2])
    diff_distances = [x[:] for x in [[0] * len(ends_and_middle)] * len(ends_and_middle)]
    for i in range(0,max(ends_and_middle_grouping)):
        for j in range(i+1, max(ends_and_middle_grouping)+1):
            # Calculate movements for each TM relative to their "normal" distance
            # selected residues for group 1 and 2
            group_1 = [x for x in range(0,len(ends_and_middle_grouping)) if ends_and_middle_grouping[x] == i]
            group_2 = [x for x in range(0,len(ends_and_middle_grouping)) if ends_and_middle_grouping[x] == j]

            diff_distances[i][j] = np.sum(abs(membrane_data1[group_1][:, group_2] - membrane_data2[group_1][:, group_2]))/(np.sum(membrane_data1[group_1][:, group_2]+membrane_data2[group_1][:, group_2])/2)*100
            diff_distances[j][i] = diff_distances[i][j]

    # Ranking for each TM
    sum_differences = [sum(x) for x in diff_distances]
    # normalized_differences = [((sum_differences[i]-min(sum_differences[0:7]))/(max(sum_differences[0:7])-min(sum_differences[0:7])))**2 for i in range(0,7)]
    for i in range(0,7):
        diff_distances[i] = [sorted(diff_distances[i]).index(x) for x in diff_distances[i]]
    final_rank = [sum([diff_distances[j][i] for j in range(0,7)]) for i in range(0,7)]

    # Grab stable TMs
    tm_ranking = [0] * 7
    sorted_rank = sorted(final_rank)
    for i in range(0,7):
        tm_ranking[i] = final_rank.index(sorted_rank[i])
        final_rank[tm_ranking[i]] = 100 # make sure this TM isn't repeated

    # Calculate 3D coordinates from distance matrix
    tms_centroids_set1, tms_set1 = recreate3Dorder(membrane_data1, ends_and_middle_grouping)
    tms_centroids_set2, tms_set2 = recreate3Dorder(membrane_data2, ends_and_middle_grouping)

    # Align 3D points of set2 with 3D points of set1 using the most stable reference points
    best_rmsd = 1000
    best_set = []
    # Disabled the testing RMSD for now
    for comb in combinations(tm_ranking[:3], 3):
    #for comb in combinations(tm_ranking[:4], 3):
        sel_refs = [x for x in range(0,len(segment_order)) if segment_order[x] in comb]
        #print(sel_refs)

        tms_reference_set1 = np.array(tms_centroids_set1[sel_refs], copy = True)
        tms_reference_set2 = np.array(tms_centroids_set2[sel_refs], copy = True)

        imposer = SVDSuperimposer()
        imposer.set(tms_reference_set1, tms_reference_set2)
        imposer.run()
        rot, trans = imposer.get_rotran()
        rmsd = imposer.get_rms()

        print("RMSD", round(rmsd,2), tm_ranking)
        if rmsd < best_rmsd:
            best_set = comb
            best_rmsd = rmsd

    # Check for possible mirroring error
    test_set2 = np.dot(tms_centroids_set2, rot) + trans
    error = 0
    for i in tm_ranking[3:7]:
        if np.linalg.norm(test_set2[i] - tms_centroids_set1[i]) > 5:
            error += 1

    #if rmsd > 2:
    #if error >= 3 or rmsd > 2:
    if True:
        for i in range(0,len(tms_centroids_set2)):
            tms_centroids_set2[i][2] = tms_centroids_set2[i][2]*-1

        # Align 3D points of set2 with 3D points of set1 using the most stable reference points
        tms_reference_set1 = tms_centroids_set1[[x for x in range(0,len(segment_order)) if segment_order[x] in tm_ranking[0:3]]]
        tms_reference_set2 = tms_centroids_set2[[x for x in range(0,len(segment_order)) if segment_order[x] in tm_ranking[0:3]]]

        imposer = SVDSuperimposer()
        imposer.set(tms_reference_set1, tms_reference_set2)
        imposer.run()
        new_rot, new_trans = imposer.get_rotran()
        new_rmsd = imposer.get_rms()
        print("RMSD2", round(new_rmsd,2))

        if new_rmsd < rmsd:
            rot = new_rot
            trans = new_trans
            rmsd = new_rmsd
        else:
            for i in range(0,len(tms_centroids_set2)):
                tms_centroids_set2[i][2] = tms_centroids_set2[i][2]*-1

    # test_set2 = np.dot(tms_reference_set2, rot) + trans
    # for i in range(0,len(test_set2)):
    #     print("pseudoatom s1_tm" + str(i+1), ", pos=[", ','.join([str(x) for x in tms_reference_set1[i]]), "]")
    # for i in range(0,len(test_set2)):
    #     print("pseudoatom s2_tm" + str(i+1), ", pos=[", ','.join([str(x) for x in test_set2[i]]), "]")
    #
    # print("############")
    # #test_set2 = np.dot(tms_centroids_set2, rot) + trans
    # test_set2 = np.array(tms_centroids_set2, copy = True)
    # for i in range(0,len(tms_centroids_set1)):
    #     print("pseudoatom s1_tm" + str(i+1), ", pos=[", ','.join([str(x) for x in tms_centroids_set1[i]]), "]")
    # for i in range(0,len(tms_centroids_set2)):
    #     print("pseudoatom s2_tm" + str(i+1), ", pos=[", ','.join([str(x) for x in tms_centroids_set2[i]]), "]")

    # if rmsd > 2:
    #     for i in range(0,len(tms_centroids_set2)):
    #         tms_centroids_set2[i][2] = tms_centroids_set2[i][2]*-1
    #     # Huge error during alignment of "stable" helices, just use the references not the helper points
    #     tms_reference_set1 = tms_centroids_set1[[x for x in range(0,7) if segment_order[x] in tm_ranking[0:4]]]
    #     tms_reference_set2 = tms_centroids_set2[[x for x in range(0,7) if segment_order[x] in tm_ranking[0:4]]]
    #     imposer = SVDSuperimposer()
    #     imposer.set(tms_reference_set1, tms_reference_set2)
    #     imposer.run()
    #     rot, trans = imposer.get_rotran()
    #     rmsd = imposer.get_rms()
    #     print("RMSD3", round(rmsd,2))
    #

    tms_centroids_set2 = np.dot(tms_centroids_set2, rot) + trans
    tms_set2 = np.dot(tms_set2, rot) + trans

    # Calculate optimal plane through points in both sets and convert to 2D
    # Try normal based on TM7
    # tm7_centroids = tms_centroids_set1[[x for x in range(0,len(segment_order)) if segment_order[x] == 6]]
    # if len(tm7_centroids) == 2:
    #     normal = (tm7_centroids[1] - tm7_centroids[0])/np.linalg.norm(tm7_centroids[1] - tm7_centroids[0])
    # else:
    #     # Using TM mid as reference plane
    #     normal, midpoint = calculatePlane(np.concatenate((tms_centroids_set1[7:], tms_centroids_set2[7:])), intracellular)

    # Alternative: use center of helical ends and center of helical middle
    #    normal = tms_centroids_set1[:7].mean(axis=0)  - tms_centroids_set1[7:].mean(axis=0)
    #    normal = normal/np.linalg.norm(normal)

    # 7TM references
    tm_centroids = {y:[] for y in range(0,7)}
    [tm_centroids[y].append(tms_centroids_set1[x]) for y in range(0,7) for x in range(0,len(segment_order)) if segment_order[x] == y]
    count = 0
    normal = np.array([0.0,0.0,0.0])
    for y in range(0,7):
        #if len(tm_centroids[y]) == 2 and (mode != 1 or y != 5):
        if len(tm_centroids[y]) == 2:
            normal += np.array((tm_centroids[y][1] - tm_centroids[y][0])/np.linalg.norm(tm_centroids[y][1] - tm_centroids[y][0]))
            count += 1
    normal = normal/count

    midpoint = tms_centroids_set1[:7].mean(axis=0)

    #plane_set1, z_set1 = convert3D_to_2D_plane(tms_centroids_set1[:7], intracellular, normal, midpoint)
    #plane_set2, z_set2 = convert3D_to_2D_plane(tms_centroids_set2[:7], intracellular, normal, midpoint)
    plane_set, z_set = convert3D_to_2D_plane(np.concatenate((tms_centroids_set1[:7], tms_centroids_set2[:7]), axis = 0), intracellular, normal, midpoint)
    plane_set1 = plane_set[:7]
    plane_set2 = plane_set[7:]
    z_set1 = z_set[:7]
    z_set2 = z_set[7:]

    # DO NOT REMOVE: possibly we want to upgrade to weighted superposing
    # Based on Biopython SVDSuperimposer
    # coords = tms_centroids_set2
    # reference_coords = tms_centroids_set1

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
    # rot, trans = imposer.get_rotran()
    # tms_set2 = np.dot(tms_set2, rot) + trans

    # CURRENT: Ca-angle to axis core
    rotations = [0] * 7
    for i in range(0,7):
        try:
            # rotations[i] = [data['tab4'][gn_dictionary[x]]['angles_set1'][1]-data['tab4'][gn_dictionary[x]]['angles_set2'][1] if abs(data['tab4'][gn_dictionary[x]]['angles_set1'][1]-data['tab4'][gn_dictionary[x]]['angles_set2'][1]) < 180 else -1*data['tab4'][gn_dictionary[x]]['angles_set2'][1]-data['tab4'][gn_dictionary[x]]['angles_set1'][1] for x in gns[i]]
            angles1 = [data['tab4'][gn_dictionary[x]]['angles_set1'][11] for x in gns[i]]
            angles1 = [angle if angle > 0 else angle + 360 for angle in angles1 ]
            angles2 = [data['tab4'][gn_dictionary[x]]['angles_set2'][11] for x in gns[i]]
            angles2 = [angle if angle > 0 else angle + 360 for angle in angles2 ]

            rotations[i] = [angles1[x] - angles2[x] for x in range(3)]
            rotations[i] = [value if abs(value) <= 180 else value-360 if value > 0 else value+360 for value in rotations[i]]

            # count=0
            # for x in gns[i]:
            #     print(i, x, data['tab4'][gn_dictionary[x]]['angles_set1'][11], data['tab4'][gn_dictionary[x]]['angles_set2'][11], rotations[i][count])
            #     count += 1

        except:
            rotations[i] = [0.0, 0.0, 0.0]  # TODO: verify other class B errors

        # UPDATE 20-02-2020 No mirroring but top-down through GPCR
        rotations[i] = sum(rotations[i])/3
        # if intracellular:
        #     rotations[i] = -1*sum(rotations[i])/3
        # else:
        #     rotations[i] = sum(rotations[i])/3


    # ALTERNATIVE: utilize TM tip alignment (needs debugging as some angles seem off, e.g. GLP-1 active vs inactive TM2)
    # Add rotation angle based on TM point placement
    # tms_2d_set1, junk = convert3D_to_2D_plane(tms_set1, intracellular, normal, midpoint)
    # tms_2d_set2, junk = convert3D_to_2D_plane(tms_set2, intracellular, normal, midpoint)

    # rotations = [0] * 7
    # for i in range(0,7):
    #     positions = [x for x in range(0, len(ends_and_middle_grouping)) if ends_and_middle_grouping[x] == i]
    #     turn_set1 = tms_2d_set1[positions]
    #     turn_set2 = tms_2d_set2[positions]
    #
    #     # set to middle
    #     turn_set1 = turn_set1 - turn_set1.mean(axis=0)
    #     turn_set2 = turn_set2 - turn_set2.mean(axis=0)
    #
    #     # Calculate shift per residue and take average for this TM
    #     for j in range(0,len(turn_set1)):
    #         v1 = turn_set1[j]/np.linalg.norm(turn_set1[j])
    #         v2 = turn_set2[j]/np.linalg.norm(turn_set2[j])
    #         angle = np.degrees(np.arctan2(v2[1], v2[0]) - np.arctan2(v1[1],v1[0]))
    #
    #         if abs(angle) > 180:
    #             angle = 360 - abs(angle)
    #
    #         rotations[i] += angle/len(turn_set1)

    # TODO: check z-coordinates orientation
    # Step 1: collect movement relative to membrane mid
    # Step 2: find min and max TM
    # Step 3: check if orientation of min/max TM matches the z-scales + intra/extra - if not invert z-coordinates
    labeled_set1 = [{"label": "TM"+str(i+1), "x": float(plane_set1[i][0]), "y": float(plane_set1[i][1]), "z": float(z_set1[i]), "rotation" : 0} for i in range(0,7)]
    labeled_set2 = [{"label": "TM"+str(i+1), "x": float(plane_set2[i][0]), "y": float(plane_set2[i][1]), "z": float(z_set2[i]), "rotation" : rotations[i]} for i in range(0,7)]

    # Convert used GNs to right numbering
    gns_used = gns[:]
    for i in range(0,len(gns)):
        for j in range(0,len(gns[i])):
            gns_used[i][j] = gn_dictionary[gns[i][j]]
    return {"coordinates_set1" : labeled_set1, "coordinates_set2": labeled_set2, "gns_used": gns_used}


def right_gn_order(pos1, pos2):
    res1 = pos1.split("x")
    res2 = pos2.split("x")

    # Handling bulges
    pos1 = int(res1[1])
    if pos1 < 100:
        pos1 = pos1 * 10
    pos2 = int(res2[1])
    if pos2 < 100:
        pos2 = pos2 * 10
    return int(res1[0]) < int(res2[0]) or (res1[0] == res2[0] and pos1 < pos2)
