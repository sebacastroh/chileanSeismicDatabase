#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 11:34:42 2019

@author: srcastro
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy.spatial as spsp
import scipy.io as scpio
import multiprocessing as mp
from functools import partial
import time
import gc
import shapely.geometry

#%% Functions
def LatLonDepth2XYZNum(lat, lon, depth):
    a = 6378.137 # km
    f = 1./298.257223563
    b = a*(1. - f)
    d = math.sqrt(a**2*math.cos(lat*math.pi/180.)**2 + b**2*math.sin(lat*math.pi/180.)**2)
    
    xo = (a**2/d + depth)*math.cos(lat*math.pi/180.)*math.cos(lon*math.pi/180.)
    yo = (a**2/d + depth)*math.cos(lat*math.pi/180.)*math.sin(lon*math.pi/180.)
    zo = (b**2/d + depth)*math.sin(lat*math.pi/180.)
    
    xyz = np.array([xo, yo, zo])
    
    return xyz

def LatLonDepth2XYZNumPy(lat, lon, depth):
    a = 6378.137 # km
    f = 1./298.257223563
    b = a*(1. - f)
    d = np.sqrt(a**2*np.cos(lat*np.pi/180.)**2 + b**2*np.sin(lat*np.pi/180.)**2)
    
    xo = (a**2/d + depth)*np.cos(lat*np.pi/180.)*np.cos(lon*np.pi/180.)
    yo = (a**2/d + depth)*np.cos(lat*np.pi/180.)*np.sin(lon*np.pi/180.)
    zo = (b**2/d + depth)*np.sin(lat*np.pi/180.)
    
    xyz = np.c_[xo, yo, zo]
    
    return xyz

def perp(a):
    b = np.empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b

def seg_intersect(a1,a2, b1,b2) :
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = perp(da)
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    return (num / denom.astype(float))*db + b1

def distance_point2segment(A, B, P):
    """ segment line AB, point P, where each one is an array([x, y]) """
    
    distance = np.zeros(len(A))
    
    val1 = (P - A) / np.linalg.norm(P - A, axis=1)[:,np.newaxis]
    val2 = (B - A) / np.linalg.norm(B - A, axis=1)[:,np.newaxis]
    angles = np.arccos(np.sum(val1*val2, axis=1))
    posA = np.where(angles > np.pi/2.)[0]
    
    distance[posA] = np.linalg.norm(P-A, axis=1)[posA]
    
    val1 = (P - B) / np.linalg.norm(P - B, axis=1)[:,np.newaxis]
    val2 = (A - B) / np.linalg.norm(A - B, axis=1)[:,np.newaxis]
    angles = np.arccos(np.sum(val1*val2, axis=1))
    posB = np.where(angles > np.pi/2.)[0]
    
    distance[posB] = np.linalg.norm(P-B, axis=1)[posB]
    
    pos = np.ones(len(distance), dtype=bool)
    pos[posA] = False
    pos[posB] = False
    
    distance[pos] = ((B[pos,1] - A[pos,1])*P[0] - (B[pos,0] - A[pos,0])*P[1]\
            + B[pos,0]*A[pos,1] - B[pos,1]*A[pos,0])/\
            np.sqrt((B[pos,1] - A[pos,1])**2 + (B[pos,0] - A[pos,0])**2)
    
    return distance

#%% Parameters and data
plt.close('all')
plt.ioff()
#@profile
def SurfaceFault(inputs, f, all_depths, contour_lines, lon, lat, dep, extraOutputs=False):
    plt.ioff()
    lonO, latO, depthO, S, W = inputs
    L = 0.
    
    #%% Hypocenter's contour line    
    if depthO not in all_depths:
        depths = np.sort(np.hstack((all_depths, depthO)))[::-1]
        if depthO > -10.:
            line10 = contour_lines[0]
            distance = distance_point2segment(line10[:-1], line10[1:], np.array([lonO, latO]))
            pos = distance.argmin()
            
            vector = line10[pos+1] - line10[pos]
            vector = vector/math.sqrt(vector[0]**2 + vector[1]**2)
            
            vectorP = perp(vector)
            
            line = line10 - vectorP*distance[pos]
            
        else:
            contour = plt.contour(lon, lat, dep, levels=[depthO])
            line = np.vstack([path.vertices for path in contour.collections[0].get_paths()])#contour.collections[0].get_paths()[0].vertices
        index = np.where(depths==depthO)[0][0]
        lines = contour_lines[:index] + [line] + contour_lines[index:]
    else:
        index = np.where(depths==depthO)[0][0]
        lines = contour_lines
        line = contour_lines[index]
        
    n = len(depths)
    
    pos = np.argmin(np.sqrt((line[:,0]-lonO)**2 + (line[:,1]-latO)**2))
    
    if latO >= line[pos,1]:
        line = np.vstack((line[:pos], [lonO, latO], line[pos:]))
    else:
        line = np.vstack((line[:pos+1], [lonO, latO], line[pos+1:]))
        pos += 1
    
    #%% Left contour line at mid point
    d = 0
    this_W = 0.
    
    pointO = np.array([lonO, latO])
    
    normalO1 = line[pos+1]-line[pos]
    normalO1 = normalO1/np.sqrt(np.sum(normalO1**2))
    
    normalO2 = line[pos]-line[pos-1]
    normalO2 = normalO2/np.sqrt(np.sum(normalO2**2))
    
    normalO = 0.5*(normalO1 + normalO2)
    normalO = normalO/np.sqrt(np.sum(normalO**2))
    
    borderO1 = pointO + perp(normalO)*1000.
    borderO2 = pointO - perp(normalO)*1000.
    
    point1 = LatLonDepth2XYZNum(pointO[1], pointO[0], depthO)
    
    valid_left = True
    length_to_add = 0.
    
    k = 1
    mid_line_W = [[pointO[0], pointO[1], depthO]]
    while d < W/2 and index > 0:
        l2 = lines[index-k]
        new_depth = depths[index-k]
        
        side = np.sign(np.dot((l2-pointO), normalO))
        p = ((side[1:] - side[:-1]) != 0).astype(int)
        p1 = np.where(p==1)[0][0]
        p2 = p1+1
        
        intersection = seg_intersect(l2[p1], l2[p2], borderO1, borderO2)
        point2 = LatLonDepth2XYZNum(intersection[1], intersection[0], new_depth)
        
        d += np.sqrt(np.sum((point1 - point2)**2))
        
        point1 = point2.copy()
        k += 1
        mid_line_W.append([intersection[0], intersection[1], new_depth])
        if index - k <= -1:
            length_to_add = max(W/2. - d, 0.)
            valid_left = False
            break
    
    if index == 0:
        l2 = lines[index]
        length_to_add = W/2. - d
        valid_left = False
    
    this_W += d
    
    leftBor = l2.copy()
    
    left_k = index - k + 1
    left_dep = depths[left_k]
    
    aux = point1.copy()
    auxk = k - 1
    
    #%% Right contour line at mid point
    d = 0.
    point1 = LatLonDepth2XYZNum(pointO[1], pointO[0], depthO)
    add_right = False
    
    k = 1
    while d < (W/2. + length_to_add) and index < (n-1):
        l2 = lines[index+k]
        new_depth = depths[index+k]
        
        side = np.sign(np.dot((l2-pointO), normalO))
        p = ((side[1:] - side[:-1]) != 0).astype(int)
        p1 = np.where(p==1)[0][0]
        p2 = p1+1
        
        intersection = seg_intersect(l2[p1], l2[p2], borderO1, borderO2)
        point2 = LatLonDepth2XYZNum(intersection[1], intersection[0], new_depth)
        
        d += np.sqrt(np.sum((point1 - point2)**2))
        
        point1 = point2.copy()
        k += 1
        mid_line_W.append([intersection[0], intersection[1], new_depth])
        if index + k >= n:
            length_to_add = max(W/2. + length_to_add - d, 0.)
            if length_to_add != 0.:
                add_right = True
            break
    
    if index == n-1:
        l2 = lines[index]
        length_to_add = W/2. + length_to_add
        add_right = True
    
    this_W += d
    
    rightBor = l2.copy()
    right_k = index + k -1
    #%% Add to left contour line at mid point
    if valid_left and add_right:
        d = 0
        new_depth = left_dep
        point1 = aux.copy()
        k = auxk + 1
        
        while d < length_to_add:
            l2 = lines[index-k]
            new_depth = depths[index-k]
            
            side = np.sign(np.dot((l2-pointO), normalO))
            p = ((side[1:] - side[:-1]) != 0).astype(int)
            p1 = np.where(p==1)[0][0]
            p2 = p1+1
            
            intersection = seg_intersect(l2[p1], l2[p2], borderO1, borderO2)
            point2 = LatLonDepth2XYZNum(intersection[1], intersection[0], new_depth)
            
            d += np.sqrt(np.sum((point1 - point2)**2))
            
            point1 = point2.copy()
            k += 1
            mid_line_W.append([intersection[0], intersection[1], new_depth])
            if index - k <= -1:
                break
    
        leftBor = l2.copy()
        
        left_k = index - k + 1
        left_dep = depths[left_k]
        
        this_W += d
    
    #%% Superior point and normal
    i = 1
    pointI = LatLonDepth2XYZNum(latO, lonO, depthO)
    pointJ = LatLonDepth2XYZNum(line[pos-i,1], line[pos-i,0], depthO)
    dL = np.sqrt(np.sum((pointI-pointJ)**2)) # kilometers
    s = this_W*dL
    this_S = 0.
    L += dL
    mid_line_L = [[lonO, latO, depthO]]
    while s < S/2.:
        this_S += this_W*dL
        point1 = pointJ.copy()
        auxpoint1 = point1.copy()
        pointI = line[pos-i,:2]
        pointJ = line[pos-i-1,:2]
        
        normalI = pointJ - pointI
        normalI = normalI/np.sqrt(np.sum(normalI**2))
        
        borderI1 = pointI + perp(normalI)*1000.
        borderI2 = pointI - perp(normalI)*1000.
        
        mid_line_L.append([line[pos-i,0], line[pos-i,1], depthO])
        
        k = 1
        this_W = 0.
        
        while index - k >= left_k:
            l2 = lines[index-k]
            new_depth = depths[index-k]
            
            side = np.sign(np.dot((l2-pointI), normalI))
            p = ((side[1:] - side[:-1]) != 0).astype(int)
            p1 = np.where(p==1)[0][0]
            p2 = p1+1
            
            intersection = seg_intersect(l2[p1], l2[p2], borderI1, borderI2)
            point2 = LatLonDepth2XYZNum(intersection[1], intersection[0], new_depth)
            
            this_W += np.sqrt(np.sum((point1 - point2)**2))
            
            point1 = point2.copy()
            k += 1
        
        point1 = auxpoint1.copy()
        k = 1
        while index + k <= right_k:
            l2 = lines[index+k]
            new_depth = depths[index+k]
            
            side = np.sign(np.dot((l2-pointI), normalI))
            p = ((side[1:] - side[:-1]) != 0).astype(int)
            p1 = np.where(p==1)[0][0]
            p2 = p1+1
            
            intersection = seg_intersect(l2[p1], l2[p2], borderI1, borderI2)
            point2 = LatLonDepth2XYZNum(intersection[1], intersection[0], new_depth)
            
            this_W += np.sqrt(np.sum((point1 - point2)**2))
            
            point1 = point2.copy()
            k += 1
        
        # Calculate dL
        pointI = LatLonDepth2XYZNum(line[pos-i,1], line[pos-i,0], depthO)
        pointJ = LatLonDepth2XYZNum(line[pos-i-1,1], line[pos-i-1,0], depthO)
        
        dL = np.sqrt(np.sum((pointI-pointJ)**2)) # kilometers
        s += this_W*dL
        L += dL
        i += 1
        
    posA = pos - i + 1
    normalA = line[posA] - line[posA+1]
    pointA = line[posA]
    borderA1 = pointA + perp(normalA)*1000.
    borderA2 = pointA - perp(normalA)*1000.
    
    #%% Inferior point and normal
    j = 1
    pointI = LatLonDepth2XYZNum(latO, lonO, depthO)
    pointJ = LatLonDepth2XYZNum(line[pos+j,1], line[pos+j,0], depthO)
    dL = np.sqrt(np.sum((pointI-pointJ)**2)) # kilometers
    s = this_W*dL
    L += dL
    
    mid_line_L = mid_line_L[::-1]
    
    while s < S/2.:
        this_S += this_W*dL
        point1 = pointJ.copy()
        auxpoint1 = point1.copy()
        pointI = line[pos+j,:2]
        pointJ = line[pos+j+1,:2]
        
        normalI = pointJ - pointI
        normalI = normalI/np.sqrt(np.sum(normalI**2))
        
        borderI1 = pointI + perp(normalI)*1000.
        borderI2 = pointI - perp(normalI)*1000.
        
        mid_line_L.append([line[pos+j,0], line[pos+j,1], depthO])
        
        k = 1
        this_W = 0.
        while index - k >= left_k:
            l2 = lines[index-k]
            new_depth = depths[index-k]
            
            side = np.sign(np.dot((l2-pointI), normalI))
            p = ((side[1:] - side[:-1]) != 0).astype(int)
            p1 = np.where(p==1)[0][0]
            p2 = p1+1
            
            intersection = seg_intersect(l2[p1], l2[p2], borderI1, borderI2)
            point2 = LatLonDepth2XYZNum(intersection[1], intersection[0], new_depth)
            
            this_W += np.sqrt(np.sum((point1 - point2)**2))
            
            point1 = point2.copy()
            k += 1
        
        point1 = auxpoint1.copy()
        k = 1
        while index + k <= right_k:
            l2 = lines[index+k]
            new_depth = depths[index+k]
            
            side = np.sign(np.dot((l2-pointI), normalI))
            p = ((side[1:] - side[:-1]) != 0).astype(int)
            p1 = np.where(p==1)[0][0]
            p2 = p1+1
            
            intersection = seg_intersect(l2[p1], l2[p2], borderI1, borderI2)
            point2 = LatLonDepth2XYZNum(intersection[1], intersection[0], new_depth)
            
            this_W += np.sqrt(np.sum((point1 - point2)**2))
            
            point1 = point2.copy()
            k += 1
        
        # Calculate dL
        pointI = LatLonDepth2XYZNum(line[pos+j,1], line[pos+j,0], depthO)
        pointJ = LatLonDepth2XYZNum(line[pos+j+1,1], line[pos+j+1,0], depthO)
        
        dL = np.sqrt(np.sum((pointI-pointJ)**2)) # kilometers
        s += this_W*dL
        L += dL
        j += 1
        
    posB = pos + j - 1
    normalB = line[posB] - line[posB-1]
    pointB = line[posB]
    borderB1 = pointB + perp(normalB)*1000.
    borderB2 = pointB - perp(normalB)*1000.
    #%% Left superior and inferior border
    k = 1
    
    lineA = [[pointA[0], pointA[1], depthO]]
    lineB = [[pointB[0], pointB[1], depthO]]
    
    while index - k > left_k - 1:
        l2 = lines[index-k]
        new_depth = depths[index-k]
        
        side = np.sign(np.dot((l2-pointA), normalA))
        p = ((side[1:] - side[:-1]) != 0).astype(int)
        p1 = np.where(p==1)[0][0]
        p2 = p1+1
        
        intersection = seg_intersect(l2[p1], l2[p2], borderA1, borderA2)
        lineA.append([intersection[0], intersection[1], new_depth])
        
        side = np.sign(np.dot((l2-pointB), normalB))
        p = ((side[1:] - side[:-1]) != 0).astype(int)
        p1 = np.where(p==1)[0][0]
        p2 = p1+1
        
        intersection = seg_intersect(l2[p1], l2[p2], borderB1, borderB2)
        lineB.append([intersection[0], intersection[1], new_depth])
        
        k += 1
    
    lineA = np.array(lineA)[::-1].tolist()
    lineB = np.array(lineB)[::-1].tolist()
    
    #%% Right superior and inferior border
    k = 1
    
    while index + k <= right_k:
        l2 = lines[index+k]
        new_depth = depths[index+k]
        
        side = np.sign(np.dot((l2-pointA), normalA))
        p = ((side[1:] - side[:-1]) != 0).astype(int)
        p1 = np.where(p==1)[0][0]
        p2 = p1+1
        
        intersection = seg_intersect(l2[p1], l2[p2], borderA1, borderA2)
        lineA.append([intersection[0], intersection[1], new_depth])
        
        side = np.sign(np.dot((l2-pointB), normalB))
        p = ((side[1:] - side[:-1]) != 0).astype(int)
        p1 = np.where(p==1)[0][0]
        p2 = p1+1
        
        intersection = seg_intersect(l2[p1], l2[p2], borderB1, borderB2)
        
        lineB.append([intersection[0], intersection[1], new_depth])
        
        k += 1
    
    lineA = np.array(lineA)
    lineB = np.array(lineB)
    
    p1 = np.argmin(np.sqrt(np.sum((lineA[0,:2]-leftBor)**2, axis=1)))
    alpha = (leftBor[p1,0]-lineA[0,0])*(lineA[1,1]-lineA[0,1]) - (leftBor[p1,1]-lineA[0,1])*(lineA[1,0]-lineA[0,0])
    if alpha < 0.:
        p1 += 1
        
    p2 = np.argmin(np.sqrt(np.sum((lineB[0,:2]-leftBor)**2, axis=1)))
    alpha = (leftBor[p2,0]-lineB[0,0])*(lineB[1,1]-lineB[0,1]) - (leftBor[p2,1]-lineB[0,1])*(lineB[1,0]-lineB[0,0])
    if alpha > 0.:
        p2 -= 1
    
    leftBorder = leftBor[p1:p2+1]
    
    p1 = np.argmin(np.sqrt(np.sum((lineA[-1,:2]-rightBor)**2, axis=1)))
    alpha = (rightBor[p1,0]-lineA[-1,0])*(lineA[-1,1]-lineA[-2,1]) - (rightBor[p1,1]-lineA[-1,1])*(lineA[-1,0]-lineA[-2,0])
    if alpha < 0.:
        p1 += 1
        
    p2 = np.argmin(np.sqrt(np.sum((lineB[-1,:2]-rightBor)**2, axis=1)))
    alpha = (rightBor[p2,0]-lineB[-1,0])*(lineB[-1,1]-lineB[-2,1]) - (rightBor[p2,1]-lineB[-1,1])*(lineB[-1,0]-lineB[-2,0])
    if alpha > 0.:
        p2 -= 1
    
    rightBorder = rightBor[p1:p2+1]
    
    border = np.vstack((lineB[:,:2], rightBorder[::-1], lineA[:,:2][::-1], leftBorder))
    depths_surface = f(border[:,[1,0]])
    border = np.hstack((border, depths_surface[:,np.newaxis]))
    
    plt.close('all')
    gc.collect()
    
    if extraOutputs:
        mid_line_W = np.array(mid_line_W)[np.argsort(np.array(mid_line_W)[:,2])]
        mid_line_L = np.array(mid_line_L)
        
        return border, mid_line_W, mid_line_L
    else:
        return border

def border2matrix(border, points, shape, indices):
    polygon = shapely.geometry.Polygon(border)
    subpoints = np.where((points[:,0]>= border[:,0].min()) &\
                         (points[:,0]<= border[:,0].max()) &\
                         (points[:,1]>= border[:,1].min()) &\
                         (points[:,1]<= border[:,1].max()))[0]
    inside = [polygon.contains(shapely.geometry.Point(point)) for point in points[subpoints]]
    
    surface_matrix = np.zeros(shape, dtype=bool)
    surface_matrix[np.unravel_index(subpoints, shape)] = inside
    surface_matrix[indices] = False
    
    gc.collect()
    return surface_matrix
    
##%% Main
#if __name__ == '__main__':
#    
#    compute_borders = False
#    compute_matrices = False
#    read_matrices = False
#    compute_distances = False
#    plotting = True
#    
#    if compute_distances and not compute_matrices:
#        read_matrices = True
#    
#    data = np.load('sam_slab1.npz')
#    
#    lat = data['lat']
#    lon = data['lon']
#    dep = data['dep']
#    
#    data = scpio.loadmat('10000_EarthquakeScenarios_TE', squeeze_me=True)
#    
#    m = 10000
#    
#    nx = len(np.unique(lon))
#    ny = len(np.unique(lat))
#    
#    dep = np.reshape(dep, (ny,nx))
#    
#    lonO = data['longitude'][:m]
#    latO = data['latitude'][:m]
#    magO = data['magnitude'][:m]
#    eleO = data['elevation'][:m]
#    
##    aW = np.zeros_like(magO)
##    bW = np.zeros_like(magO)
#    
#    Mwlim = 8.67
#    pos = np.zeros_like(magO, dtype=bool)
#    pos[magO>Mwlim] = True
#    
##    aW[pos] = 2.29
##    aW[~pos] = -1.91
##    bW[~pos] = 0.48
#    aW = -0.86
#    bW = 0.35
#    
##    aS = np.zeros_like(magO)
##    bS = np.zeros_like(magO)
#    
##    Mwlim = 8.63#(math.log10(74000.)+5.62)/1.22
##    pos = np.zeros_like(magO, dtype=bool)
##    pos[magO>Mwlim] = True
#    
##    aS[pos] = 2.23
##    bS[pos] = 0.31
##    aS[~pos] = -5.62
##    bS[~pos] = 1.22
#    
#    aS = -3.63
#    bS = 0.96
#    
#    SO = 10**(aS + bS*magO)
#    WO = 10**(aW + bW*magO)
#    
#    lat = np.reshape(lat, (ny,nx))
#    lon = np.reshape(lon, (ny,nx))
#    
#    f = interpolate.RegularGridInterpolator((lat[:,0][::-1], lon[0]), values=dep[::-1], method='linear')
#    depthO = f((latO, lonO))
#    
#    pos = np.isnan(depthO)
#    depthO[pos] = eleO[pos]
#    
#    n = 499 # n = 29, 99, 499, 999
#    new_contour_lines = False
#    
#    inputs = np.hstack((lonO[:,np.newaxis], latO[:,np.newaxis],\
#                        depthO[:,np.newaxis], SO[:,np.newaxis], WO[:,np.newaxis]))
#    
#    if new_contour_lines:
#        min_depth = 61.
#        max_depth = 10.
#        all_depths = -np.logspace(math.log10(min_depth), math.log10(max_depth), n)
#        contour_lines = plt.contour(lon, lat, dep, levels=all_depths)
#        contour_lines = [np.vstack([path.vertices for path in collection.get_paths()])\
#                         for collection in reversed(contour_lines.collections)]
#        np.savez_compressed('contour_lines_n%i.npz' %n, contour_lines=contour_lines, all_depths=all_depths[::-1])
#    else:
#        data = np.load('contour_lines_n%i.npz' %n, allow_pickle=True)
#        all_depths = data['all_depths']
#        contour_lines = data['contour_lines'].tolist()
#    
#    if compute_borders:
#        pool = mp.Pool(processes=48)
#        t = time.time()
#        borders = pool.map(partial(SurfaceFault, f=f, all_depths=all_depths,\
#                                   contour_lines=contour_lines, lon=lon,\
#                                   lat=lat, dep=dep), inputs.tolist())
#        pool.close()
#        print(time.time() - t)
#        np.savez_compressed('borders_W1S1_n%i.npz' %n, borders=borders)
#    else:
#        borders = np.load('borders_W1S1_n%i.npz' %n, allow_pickle=True)
#        borders = borders['borders']
#    
#    if compute_matrices:
#        points = np.hstack((lon.ravel()[:,np.newaxis], lat.ravel()[:,np.newaxis]))
#        shape = lon.shape
#        indices = np.isnan(dep)
#        
#        pool = mp.Pool(processes=48)
#        t = time.time()
#        matrices = pool.map(partial(border2matrix, points=points, shape=shape, indices=indices), borders)
#        pool.close()
#        print(time.time() - t)
#        
#        np.savez_compressed('matrices_W1S1.npz', matrices=matrices)
#    else:
#        if read_matrices:
#            matrices = np.load('matrices_W1S1.npz')
#            matrices = matrices['matrices']
#        
#    if compute_distances:
#        iquiquePoint = LatLonDepth2XYZNumPy(-20.221011, -70.152011, 0.) # lat, lon, depth
#        distances = []
#        
#        for i in range(m):
#            lat_fault = lat[matrices[i]]
#            lon_fault = lon[matrices[i]]
#            dep_fault = dep[matrices[i]]
#            
#            points = LatLonDepth2XYZNumPy(lat_fault, lon_fault, dep_fault)
#            
#            distances.append(np.min(spsp.distance.cdist(points, iquiquePoint)))
#            
#        distances = np.array(distances)
#        np.save('distanceFaults2Iquique.npy', distances)
#    else:
#        distances = np.load('distanceFaults2Iquique.npy')
#    
#    if plotting:
#        k = 7301#1530#788
#        border, mid_line_W, mid_line_L = SurfaceFault(inputs[k], f, all_depths,\
#                                                      contour_lines, lon, lat,\
#                                                      dep, extraOutputs=True)
#        plt.ion()
#        fig = plt.figure(figsize=(14.4,4.8))
#        
#        # Left figure
#        from mpl_toolkits.mplot3d import Axes3D
#        import matplotlib.cm as cm
#        from matplotlib.gridspec import GridSpec
#        
#        gs = GridSpec(1, 3, figure=fig)
#        
#        colors = cm.viridis((all_depths-all_depths.min())/(all_depths.max() - all_depths.min()))
#        
#        a1 = fig.add_subplot(gs[0,:-1], projection='3d')
#        
#        p1 = np.argmin(np.abs(all_depths-np.nanmin(border[:,2])))
#        p2 = np.argmin(np.abs(all_depths-np.nanmax(border[:,2])))
#        
#        p1 = np.where(np.abs(border[:,2]-all_depths[p1])<1e-5)[0]
#        p1 = p1[np.where(p1[1:]-p1[:-1]==1)]
#        p2 = np.where(np.abs(border[:,2]-all_depths[p2])<1e-5)[0]
#        p2 = p2[np.where(p2[1:]-p2[:-1]==1)]
#        
#        lineA = border[p1[-1]:p2[0]+1]
#        lineB = border[:p1[0]+1]
#        
#        surface = np.empty((0,3))
#        
#        for i in range(n):
#            line_displaced = lineA + np.array([-0.2, 0.2, 0.])
#            p1 = line_displaced[0,[0,1]]
#            p2 = line_displaced[-1,[0,1]]
#            
#            signs_up = (p2[1]-p1[1])*contour_lines[i][:,0] - (p2[0]-p1[0])*contour_lines[i][:,1] + p2[0]*p1[1] - p2[1]*p1[0]
#            
#            line_displaced = lineB + np.array([-0.7, -0.7, 0.])
#            p1 = line_displaced[0,[0,1]]
#            p2 = line_displaced[-1,[0,1]]
#            
#            signs_down = (p2[1]-p1[1])*contour_lines[i][:,0] - (p2[0]-p1[0])*contour_lines[i][:,1] + p2[0]*p1[1] - p2[1]*p1[0]
#            
#            pos = np.where((signs_down < 0) & (signs_up < 0))[0]
#            this_line = contour_lines[i][pos]
#            line_xyz = np.vstack([LatLonDepth2XYZNum(this_line[j,1],\
#                                                     this_line[j,0],\
#                                                     all_depths[i])\
#            for j in range(len(this_line))])
#            a1.plot(line_xyz[:,0], line_xyz[:,1], line_xyz[:,2], c=colors[i])
#            surface = np.vstack((surface, line_xyz))
#        
#        points = np.hstack((lon.ravel()[:,np.newaxis], lat.ravel()[:,np.newaxis]))
#        shape = lon.shape
#        indices = np.isnan(dep)
#        
#        matrix = border2matrix(border, points, shape, indices)
#        
#        for i in range(1101):
#            this_lon = lon[matrix[:,i],i]
#            this_lat = lat[matrix[:,i],i]
#            this_dep = dep[matrix[:,i],i]
#            slab_xyz = LatLonDepth2XYZNumPy(this_lat, this_lon, this_dep)
#            if len(slab_xyz) > 0:
#                a1.plot(slab_xyz[:,0], slab_xyz[:,1], slab_xyz[:,2], 'gray', lw=0.5)
#                
#        for i in range(2501):
#            this_lon = lon[i,matrix[i,:]]
#            this_lat = lat[i,matrix[i,:]]
#            this_dep = dep[i,matrix[i,:]]
#            slab_xyz = LatLonDepth2XYZNumPy(this_lat, this_lon, this_dep)
#            if len(slab_xyz) > 0:
#                a1.plot(slab_xyz[:,0], slab_xyz[:,1], slab_xyz[:,2], 'gray', lw=0.5)
#        
#        mid_line_W_xyz = LatLonDepth2XYZNumPy(mid_line_W[:,1], mid_line_W[:,0], mid_line_W[:,2])
#        mid_line_L_xyz = LatLonDepth2XYZNumPy(mid_line_L[:,1], mid_line_L[:,0], mid_line_L[:,2])
#        
#        a1.plot(mid_line_W_xyz[:,0], mid_line_W_xyz[:,1], mid_line_W_xyz[:,2], 'k')
#        a1.plot(mid_line_L_xyz[:,0], mid_line_L_xyz[:,1], mid_line_L_xyz[:,2], 'k')
#        
#        border_xyz = np.vstack([LatLonDepth2XYZNum(border[i,1], border[i,0], border[i,2]) for i in range(len(border))])
#        a1.plot(border_xyz[:,0], border_xyz[:,1], border_xyz[:,2], lw=2)
#        
#        textW_xyz = LatLonDepth2XYZNum(-19.52, -71.12, -20.)
#        a1.text3D(textW_xyz[0], textW_xyz[1], textW_xyz[2], 'W', fontsize=20)
#        
#        textL_xyz = LatLonDepth2XYZNum(-18.52, -71.2, -20.)
#        a1.text3D(textL_xyz[0], textL_xyz[1], textL_xyz[2], 'L', fontsize=20)
#        
#        a1.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
#        a1.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
#        a1.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
#        a1.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
#        a1.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
#        a1.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
#        a1.set_xticks([])
#        a1.set_yticks([])
#        a1.set_zticks([])
#        a1._axis3don = False
#        
#        a1.view_init(-65, -80)
#        a1.dist = 6.4 # menor valor, mayor zoom -- duh, obvio, es distancia
#        
#        cbaxes = fig.add_axes([0.1, 0.1, 0.03, 0.8])
#        scalarMappable = cm.ScalarMappable()
#        scalarMappable.set_clim(vmin=-60., vmax=-10.)
#        cbar = fig.colorbar(scalarMappable, cax=cbaxes, fraction=0.05)
#        cbar.set_label('Depth [km]')
#        cbar.set_ticks(np.linspace(-60, -10, 6))
#        cbar.set_ticklabels(['%i' %abs(tick) for tick in cbar.get_ticks()])
#        cbar.ax.yaxis.set_label_position('left')
#        
#        # Right figure
#        import cartopy.crs as ccrs
#        from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
#        
#        resolution = 'large'#4096px'
#        
#        a2 = fig.add_subplot(gs[0,-1], projection=ccrs.PlateCarree())
#        a2.set_extent([-74., -68., -23., -14.], crs=ccrs.PlateCarree())
#        a2.background_img(name="natural-earth-1", resolution=resolution)
#        a2.plot(border[:,0], border[:,1])
#        
#        a2.set_xticks(np.arange(-74, -67, 3), crs=ccrs.PlateCarree())
#        lon_formatter = LongitudeFormatter(zero_direction_label=True)
#        a2.xaxis.set_major_formatter(lon_formatter)
#        
#        a2.set_yticks(np.arange(-23, -13, 2), crs=ccrs.PlateCarree())
#        lat_formatter = LatitudeFormatter()
#        a2.yaxis.set_major_formatter(lat_formatter)
#        
#        a2.scatter(-70.152011, -20.221011, s=200, marker='*', c='#FF9900')
#
#        plt.show()
#        
#        fig.savefig('surfaceFault.pdf', fmt='pdf', bbox_inches='tight', pad_inches=0, dpi=300)