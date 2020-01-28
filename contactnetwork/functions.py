"""
A set of utility functions for processing distances
"""
from contactnetwork.distances import *

import math
import numpy as np
from numpy.linalg import norm

# convert 3D points to 2D points on best fitting plane
def convert2D_SVD(tm_points, intracellular):
    # Find best fit plane through set of points (singular value decomposition)
    # Based on https://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points
    centroid = np.sum(tm_points,0) / len(tm_points)
    points = np.transpose(tm_points - centroid) # subtract out the centroid
    svd = np.transpose(np.linalg.svd(points)) # singular value decomposition
    # the corresponding left singular vector is the normal vector of the best-fitting plane
    normal = np.transpose(svd[0])[2]

    # 3D project points onto plane -> 2D points
    points_plane = [[0]] * 7
    for i in range(0,len(tm_points)):
        # point relative to centroid
        v = tm_points[i] - centroid
        # distance from point to plane along normal vector
        dist = sum(v*normal)
        # point onto plane
        points_plane[i] = tm_points[i] - dist * normal

    # Create X and Y axes on plane
    locx = (points_plane[0] - centroid) # TM1 right side
    locy = np.cross(normal, locx)
    if not intracellular:
        locy = -1 * locy
    locx = locx/norm(locx)
    locy = locy/norm(locy)
    points_2d = [[np.dot(p - centroid, locx), np.dot(p - centroid, locy)] for p in points_plane]

    return np.array(points_2d)

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

# 3D reconstruction based on average distances
def recreate3D(distance_matrix, intracellular = True):
    tms = [[0]] * 7
    tms[0] = np.array([0, 0, 0])
    tms[1] = np.array([distance_matrix[0][1], 0, 0])

    # place TM3 - relative to TM1 and TM2
    d = tms[1][0]
    a = (math.pow(distance_matrix[0][2],2) - math.pow(distance_matrix[1][2],2) + math.pow(d,2)) / (2*d)
    h = math.sqrt(math.pow(distance_matrix[0][2],2) - math.pow(a,2))
    x2 = a*tms[1][0]/d
    y2 = a*tms[1][1]/d
    x3 = x2+h*tms[1][1]/d
    y3 = y2-h*tms[1][0]/d

    # Either positive or negative (switch for IC/EC)
    if intracellular:
        y3 = abs(y3)
    else:
        y3 = -1*abs(y3)

    tms[2] = np.array([x3, y3, 0])

    for i in range(3,7):
        #print("calculating for TM", str(i+1))
        tms[i] = trilaterate(tms[0], tms[1], tms[2], distance_matrix[0][i], distance_matrix[1][i], distance_matrix[2][i])
        if i == 4:
            # Determine placement TM4 and TM5 using each other
            ref_dist = distance_matrix[(i-1)][i]
            changes = [abs(math.sqrt(sum([math.pow(x,2) for x in tms[4][j]-tms[3][k]]))-ref_dist) for j in range(0,2) for k in range(0,2)]
            #print("Differences for TM", str(i+1))
            #print(changes)
            lowest = changes.index(min(changes))
            tms[3] = tms[3][math.floor(lowest/2)]
            tms[4] = tms[4][lowest%2]
        elif i > 4:
            ref_dist = distance_matrix[3][i]
            changes = [abs(math.sqrt(sum([math.pow(x,2) for x in tms[i][j]-tms[3]]))-ref_dist) for j in range(0,2)]
            #print("Differences for TM", str(i+1))
            #print(changes)
            lowest = changes.index(min(changes))
            tms[i] = tms[i][lowest]

    return tms

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

def tm_movement_2D(pdbs1, pdbs2, intracellular):
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

    distances_set1.filter_gns = [y for x in gns for y in x]
    distances_set2.filter_gns = distances_set1.filter_gns
    distances_set1.fetch_distances_tm()
    distances_set2.fetch_distances_tm()

    distance_data1 = [x[:] for x in [[0] * 7] * 7]
    distance_data2 = [x[:] for x in [[0] * 7] * 7]
    for i in range(0,6):
        for j in range(i+1, 7):
            filter_keys = [x+"_"+y for x in gns[i] for y in gns[j]]
            if len(filter_keys) == 0:
                print("no filter keys")
                return []
            distance_data1[i][j] = sum([k for x in filter_keys for k in distances_set1.data[x]])/(len(filter_keys)*len(pdbs1))
            distance_data2[i][j] = sum([k for x in filter_keys for k in distances_set2.data[x]])/(len(filter_keys)*len(pdbs2))

    tms_set1 = recreate3D(distance_data1, 1)
    plane_set1 = convert2D_SVD(tms_set1, intracellular)

    tms_set2 = recreate3D(distance_data2, 1)
    plane_set2 = convert2D_SVD(tms_set2, intracellular)

    # Identify most stable TMs by ranking the variations to all other helices
    diff_distances = [x[:] for x in [[0] * 7] * 7]
    real_differences = [x[:] for x in [[0] * 7] * 7]
    for i in range(0,6):
        for j in range(i+1, 7):
            # Calculate movements for each TM relative to their "normal" distance
            diff_distances[i][j] = abs(distance_data2[i][j] - distance_data1[i][j])/distance_data1[i][j]*100
            diff_distances[j][i] = diff_distances[i][j]
            real_differences[i][j] = distance_data2[i][j] - distance_data1[i][j]
            real_differences[j][i] = real_differences[i][j]

    # Ranking for each TM
    for i in range(0,7):
        diff_distances[i] = [sorted(diff_distances[i]).index(x) for i in range(0,7) for x in diff_distances[i]]

    # TODO Make oneliner from this
    final_rank = [0] * 7
    for i in range(0,7):
        final_rank[i] = sum([diff_distances[j][i] for j in range(0,7)])

    # Grab stable TMs
    final_rank = [sorted(final_rank).index(x) for x in final_rank]
    stable_one, stable_two = final_rank.index(0), final_rank.index(1)
    if stable_one > stable_two:
        stable_one, stable_two = stable_two, stable_one


    # Rescale positions to match true length
    scale = math.sqrt(math.pow(plane_set1[stable_two][0]-plane_set1[stable_one][0],2) + math.pow(plane_set1[stable_two][1]-plane_set1[stable_one][1],2))/distance_data1[stable_one][stable_two]
    plane_set1 = plane_set1/scale
    plane_set2 = plane_set2/scale

    vector = plane_set1[stable_two] - plane_set1[stable_one]
    vector = vector/np.linalg.norm(vector)
    length_stable2 = math.sqrt(math.pow(plane_set2[stable_two][0]-plane_set2[stable_one][0],2) + math.pow(plane_set2[stable_two][1]-plane_set2[stable_one][1],2))

    diff_correction = (length_stable2 - distance_data1[stable_one][stable_two])/2

    # Apply the translation
    new_one = plane_set1[stable_one] - vector * diff_correction
    new_two = new_one + vector * length_stable2
    plane_set2 = plane_set2 - plane_set2[stable_one] + new_one

    # Apply rotation
    old_vector = plane_set2[stable_two] - plane_set2[stable_one]
    old_vector = old_vector/np.linalg.norm(old_vector)
    theta = np.arctan2(np.linalg.norm(np.cross(old_vector, vector)), np.dot(old_vector, vector))
    r = np.array(( (np.cos(theta), -np.sin(theta)),
                   (np.sin(theta),  np.cos(theta)) ))
    plane_set2 = plane_set2 - new_one
    plane_set2 = [r.dot(x) for x in plane_set2]
    plane_set2 = plane_set2 + new_one

    labeled_set1 = [{"label": "TM"+str(i+1), "x": plane_set1[i][0], "y": plane_set1[i][1], "rotation" : 0} for i in range(0,7)]
    labeled_set2 = [{"label": "TM"+str(i+1), "x": plane_set2[i][0], "y": plane_set2[i][1], "rotation" : 0} for i in range(0,7)]

    return {"coordinates_set1" : labeled_set1, "coordinates_set2": labeled_set2, "rotation": [], "gns_used": gns}
