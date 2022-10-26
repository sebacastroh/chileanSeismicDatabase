#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 09:51:43 2020

@author: srcastro
"""
import os
import math
import pickle
import pyproj
import geojson
import numpy as np
import scipy.io as spio


def vinc_pt( phi1, lembda1, alpha12, s ) :
    """

    Returns the lat and long of projected point and reverse azimuth
    given a reference point and a distance and azimuth to project.
    lats, longs and azimuths are passed in decimal degrees

    Returns ( phi2,  lambda2,  alpha21 ) as a tuple 

    """
    f = 1.0 / 298.257223563		# WGS84
    a = 6378137.0 			# metres
    piD4 = math.atan( 1.0 )
    two_pi = piD4 * 8.0

    phi1    = phi1    * piD4 / 45.0
    lembda1 = lembda1 * piD4 / 45.0
    alpha12 = alpha12 * piD4 / 45.0
    if ( alpha12 < 0.0 ) : 
            alpha12 = alpha12 + two_pi
    if ( alpha12 > two_pi ) : 
            alpha12 = alpha12 - two_pi

    b = a * (1.0 - f)

    TanU1 = (1-f) * math.tan(phi1)
    U1 = math.atan( TanU1 )
    sigma1 = math.atan2( TanU1, math.cos(alpha12) )
    Sinalpha = math.cos(U1) * math.sin(alpha12)
    cosalpha_sq = 1.0 - Sinalpha * Sinalpha

    u2 = cosalpha_sq * (a * a - b * b ) / (b * b)
    A = 1.0 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * \
            (320 - 175 * u2) ) )
    B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2) ) )

    # Starting with the approximation
    sigma = (s / (b * A))

    last_sigma = 2.0 * sigma + 2.0	# something impossible

    # Iterate the following three equations 
    #  until there is no significant change in sigma 

    # two_sigma_m , delta_sigma
    while ( abs( (last_sigma - sigma) / sigma) > 1.0e-9 ) :
            two_sigma_m = 2 * sigma1 + sigma

            delta_sigma = B * math.sin(sigma) * ( math.cos(two_sigma_m) \
                    + (B/4) * (math.cos(sigma) * \
                    (-1 + 2 * pow( math.cos(two_sigma_m), 2 ) -  \
                    (B/6) * math.cos(two_sigma_m) * \
                    (-3 + 4 * pow(math.sin(sigma), 2 )) *  \
                    (-3 + 4 * pow( math.cos (two_sigma_m), 2 ))))) \

            last_sigma = sigma
            sigma = (s / (b * A)) + delta_sigma

    phi2 = math.atan2 ( (math.sin(U1) * math.cos(sigma) + math.cos(U1) * math.sin(sigma) * math.cos(alpha12) ), \
            ((1-f) * math.sqrt( pow(Sinalpha, 2) +  \
            pow(math.sin(U1) * math.sin(sigma) - math.cos(U1) * math.cos(sigma) * math.cos(alpha12), 2))))

    lembda = math.atan2( (math.sin(sigma) * math.sin(alpha12 )), (math.cos(U1) * math.cos(sigma) -  \
            math.sin(U1) *  math.sin(sigma) * math.cos(alpha12)))

    C = (f/16) * cosalpha_sq * (4 + f * (4 - 3 * cosalpha_sq ))

    omega = lembda - (1-C) * f * Sinalpha *  \
            (sigma + C * math.sin(sigma) * (math.cos(two_sigma_m) + \
            C * math.cos(sigma) * (-1 + 2 * pow(math.cos(two_sigma_m),2) )))

    lembda2 = lembda1 + omega

    alpha21 = math.atan2 ( Sinalpha, (-math.sin(U1) * math.sin(sigma) +  \
            math.cos(U1) * math.cos(sigma) * math.cos(alpha12)))

    alpha21 = alpha21 + two_pi / 2.0
    if ( alpha21 < 0.0 ) :
            alpha21 = alpha21 + two_pi
    if ( alpha21 > two_pi ) :
            alpha21 = alpha21 - two_pi

    phi2       = phi2       * 45.0 / piD4
    lembda2    = lembda2    * 45.0 / piD4
    alpha21    = alpha21    * 45.0 / piD4

    return phi2,  lembda2,  alpha21 

def generateSubFaults(lat, lon, depth, strike, dip, length, width):
    nx = int(length) #  ~ subfaults 1km
    ny = int(width)
    
    L = length * 1000
    W = width * 1000
    
    # Procedemos a calcular las coordenadas geograficas de los 4 extremos que daran lugar a nuestro plano de falla
    latt1, lont1, alpha21 = vinc_pt(lat, lon, (strike + 90), (W/2)*math.cos(math.radians(dip)))
    
    # Los 4 puntos que definen la falla grande correponden a los siguientes:
    lat1, lon1, alpha21 = vinc_pt(latt1,lont1, strike+0., L/2. )
    
    # Pasamos a calcular las coordenadas de cada una de las subfallas
    l = L/(nx-1)
    w = W/(ny-1)
    Hmax = depth + (width/2)*math.sin(math.radians(dip))
    longitudes = []
    latitudes = []
    depths = []
    
    for i in range(nx):
        if i == 0:
            this_lon = lon1
            this_lat = lat1
        else:
            this_lat, this_lon, alpha21 = vinc_pt(lat1, lon1, (strike - 180), l*i)
            
        longitudes.append(this_lon)
        latitudes.append(this_lat)
        depths.append(Hmax)
        for j in range(1,ny):
            lla2, llo2, alpha21 = vinc_pt(this_lat, this_lon, (strike - 90), w*j)
            if lla2 < -35:
                print(i,j)
            longitudes.append(llo2)
            latitudes.append(lla2)
            depths.append(Hmax - w*math.sin(math.radians(dip))*j/1000.)
            
    return np.c_[longitudes, latitudes, depths]

def inPolygon(point, polygon):
    nvert = len(polygon)
    inside = 0
    j = nvert - 1
    for i in range(nvert):
        if ((polygon[i,1] > point[1]) != (polygon[j,1] > point[1])):
            if polygon[i,1] == polygon[j,1]:
                inside = not inside
            else:
                check = (polygon[j,0] - polygon[i,0]) * (point[1] - polygon[i,1]) / (polygon[j,1] - polygon[i,1]) + polygon[i,0]
                if point[0] < check:
                    inside = not inside
        j = i
    
    return bool(inside)


def LatLonDepth2XYZNumPy(lat, lon, depth):
    """Transform latitude, longitude and depth into xyz system. Uses NumPy.

    Parameters
    ----------
    lat : ndarray
        Latitude points
    lon : ndarray
        Longitude points
    depth : ndarray
        Depth points. In kilometers

    Returns
    -------
    xyz : ndarray
        Points in xyz system. In kilometers
    """
    a = 6378.137 # km
    f = 1./298.257223563
    b = a*(1. - f)
    d = np.sqrt(a**2*np.cos(lat*np.pi/180.)**2 + b**2*np.sin(lat*np.pi/180.)**2)
    
    xo = (a**2/d - depth)*np.cos(lat*np.pi/180.)*np.cos(lon*np.pi/180.)
    yo = (a**2/d - depth)*np.cos(lat*np.pi/180.)*np.sin(lon*np.pi/180.)
    zo = (b**2/d - depth)*np.sin(lat*np.pi/180.)
    
    xyz = np.c_[xo, yo, zo]
    
    return xyz

def LatLonDepth2XYZNum(lat, lon, depth):
    """Transform latitude, longitude and depth into xyz system.

    Parameters
    ----------
    lat : float
        Latitude point
    lon : float
        Longitude point
    depth : float
        Depth point. In kilometers

    Returns
    -------
    xyz : ndarray
        Point in xyz system in kilometers
    """
    a = 6378.137 # km
    f = 1./298.257223563
    b = a*(1. - f)
    d = math.sqrt(a**2*math.cos(lat*math.pi/180.)**2 + b**2*math.sin(lat*math.pi/180.)**2)
    
    xo = (a**2/d - depth)*math.cos(lat*math.pi/180.)*math.cos(lon*math.pi/180.)
    yo = (a**2/d - depth)*math.cos(lat*math.pi/180.)*math.sin(lon*math.pi/180.)
    zo = (b**2/d - depth)*math.sin(lat*math.pi/180.)
    
    xyz = np.array([xo, yo, zo])
    
    return xyz

def distancesPoint(filename, seismic_path):
    event = spio.loadmat(os.path.join(seismic_path, filename + '.mat'),
                         struct_as_record=False, squeeze_me=True)
    event.pop('__header__')
    event.pop('__globals__')
    event.pop('__version__')
    
    geod = pyproj.Geod(ellps='WGS84')
    
    hypocenter = LatLonDepth2XYZNum(event['st00'].hypocenter_lat,\
                                    event['st00'].hypocenter_lon,\
                                    event['st00'].depth)
    distances = {}
    for key, station in event.items():
        lon = station.lon
        lat = station.lat
        station_xyz = LatLonDepth2XYZNum(lat, lon, 0.)
        
        Rhypo = np.sqrt(np.sum((hypocenter - station_xyz)**2))
        azimuth1, azimuth2, Repi = geod.inv(lon, lat,\
                                            station.hypocenter_lon, station.hypocenter_lat)
        Repi /= 1000.
        
        distances[key] = [Rhypo, Repi, Rhypo, Repi]
        
    return distances
    

def distancesWithFFM(filename, seismic_path, ffm_path):
    event = spio.loadmat(os.path.join(seismic_path, filename + '.mat'),
                         struct_as_record=False, squeeze_me=True)
    event.pop('__header__')
    event.pop('__globals__')
    event.pop('__version__')
    
    geod = pyproj.Geod(ellps='WGS84')
    
    if filename.startswith('1985'):
        fault = np.loadtxt(os.path.join(ffm_path, filename + '.fsp'), comments='%')
        ffm = {}
        ffm['features'] = []
        
        fault = np.reshape(fault, (11,17,6))
        for i in range(10):
            for j in range(16):
                p1 = fault[i,j]
                p2 = fault[i,j+1]
                p3 = fault[i+1,j+1]
                p4 = fault[i+1,j]
                data = np.vstack((p1,p2,p3,p4))
                data[:,4] *= 1000.
                slip = data[:,5].mean()
                points = [data[:,[1,0,4]].tolist()]
                ffm['features'].append({'geometry': {'coordinates':points}, 'properties': {'slip': slip}})
    else:
        with open(os.path.join(ffm_path, filename + '.geojson'), 'r') as f:
            ffm = geojson.load(f)
        
    max_slip = -np.inf
    min_slip =  np.inf
    for feature in ffm['features']:
        slip = feature['properties']['slip']
        max_slip = max(max_slip, slip)
        min_slip = min(min_slip, slip)
    
    dslip = 0.15*(max_slip - min_slip) + min_slip
    
    hypocenter = LatLonDepth2XYZNum(event['st00'].hypocenter_lat,\
                                    event['st00'].hypocenter_lon,\
                                    event['st00'].depth)
    
    distances = {}
        
    for key, station in event.items():
        lon = station.lon
        lat = station.lat
        station_xyz = LatLonDepth2XYZNum(lat, lon, 0.)
        
        Rhypo = np.sqrt(np.sum((hypocenter - station_xyz)**2))
        azimuth1, azimuth2, Repi = geod.inv(lon, lat,\
                                            station.hypocenter_lon, station.hypocenter_lat)
        Repi /= 1000.
        
        Rrup = np.inf
        Rjb = np.inf
        for feature in ffm['features']:
            slip = feature['properties']['slip']
            if slip < dslip:
                continue
            
            points = np.array(feature['geometry']['coordinates'][0])[:-1]
            sub_fault = LatLonDepth2XYZNumPy(points[:,1], points[:,0], points[:,2]/1000.)
            dist = np.sqrt(np.sum((sub_fault - station_xyz)**2, axis=1)).min()
                
            Rrup = min(Rrup, dist)
            
            if inPolygon(station_xyz[:2], sub_fault[:,:2]):
                Rjb = 0.
            elif Rjb != 0.:
                azimuth1, azimuth2, dist = geod.inv(lon*np.ones(len(points)),\
                                                    lat*np.ones(len(points)),\
                                                    points[:,0], points[:,1])
                Rjb = min(Rjb, dist.min()/1000.)
        
        distances[key] = [Rhypo, Repi, Rrup, Rjb]
        
    return distances

def distancesWithoutFFM(filename, seismic_path):
    event = spio.loadmat(os.path.join(seismic_path, filename + '.mat'),
                         struct_as_record=False, squeeze_me=True)
    event.pop('__header__')
    event.pop('__globals__')
    event.pop('__version__')
    
    geod = pyproj.Geod(ellps='WGS84')
    
    with open('fault_plane_properties.pkl', 'rb') as f:
        data = pickle.load(f)
        
    if filename in data.keys():
        strike, dip, rake = data[filename]
        
        Mw = event['st00'].magnitude
        L = 10**(-2.9 + 0.63*Mw) # Allen et al 2017, Table 2
        W = 10**(-0.86 + 0.35*Mw) # Allen et al 2017, Table 2
        
        hypocenter_lat = event['st00'].hypocenter_lat
        hypocenter_lon = event['st00'].hypocenter_lon
        depth = event['st00'].depth
        subfaults = generateSubFaults(hypocenter_lat, hypocenter_lon, depth, strike, dip, L, W)
        
        valid = True
    else:
        Rrup = 'Finite fault model required'
        Rjb = 'Finite fault model required'
        valid = False
    
    hypocenter = LatLonDepth2XYZNum(event['st00'].hypocenter_lat,\
                                    event['st00'].hypocenter_lon,\
                                    event['st00'].depth)
    
    distances = {}
        
    for key, station in event.items():
        lon = station.lon
        lat = station.lat
        station_xyz = LatLonDepth2XYZNum(lat, lon, 0.)
        
        Rhypo = np.sqrt(np.sum((hypocenter - station_xyz)**2))
        azimuth1, azimuth2, Repi = geod.inv(lon, lat,\
                                            station.hypocenter_lon, station.hypocenter_lat)
        Repi /= 1000.
        
        if valid:           
            faultXYZ = LatLonDepth2XYZNumPy(subfaults[:,1],\
                                            subfaults[:,0],\
                                            subfaults[:,2])
                
            Rrup = np.sqrt(np.sum((faultXYZ - station_xyz)**2, axis=1)).min()
            
            nx = int(L)
            ny = int(W)
            polygon = subfaults[[0,ny-1,nx*ny-1,nx*ny-ny]]
            
            if inPolygon([lon,lat], polygon):
                Rjb = 0.
            else:
                azimuth1, azimuth2, dist = geod.inv(lon*np.ones(nx*ny),\
                                                    lat*np.ones(nx*ny),\
                                                    subfaults[:,0], subfaults[:,1])
                    
                Rjb = dist.min()/1000.
        
        distances[key] = [Rhypo, Repi, Rrup, Rjb]
        
    return distances

seismic_path = os.path.join(os.getcwd(), 'events_mat_uncorrected')
# seismic_path = os.path.join(os.getcwd(), 'events_mat_corrected_v2')
ffm_path = os.path.join(os.getcwd(), 'ffm')

filenames = os.listdir(seismic_path)
ffm_files = os.listdir(ffm_path)
ffm_files = [ffm.split('M.')[0] + 'M' for ffm in ffm_files]

slab = np.load('sam_slab2.npz')

for filename in filenames:
    if not filename.startswith('202107'):
    # if not filename.startswith('20210527_4.3'):
        continue
    print(filename)
    date, magnitude, lon, lat, depth = filename.split('_')
    magnitude = float(magnitude[:-2])
    
    if magnitude < 7.1:
        distances = distancesPoint(filename[:-4], seismic_path)
    else:
        if filename[:-4] in ffm_files:
            distances = distancesWithFFM(filename[:-4], seismic_path, ffm_path)
        else:
            distances = distancesWithoutFFM(filename[:-4], seismic_path)
    
    event = spio.loadmat(os.path.join(seismic_path, filename),\
                          struct_as_record=False, squeeze_me=True)
    event.pop('__globals__')
    event.pop('__version__')
    event.pop('__header__')
    new_event = {}
    
    lon = event['st00'].hypocenter_lon
    lat = event['st00'].hypocenter_lat
    dep_pos = np.nanargmin(np.sqrt((slab['lon'] - lon)**2 + (slab['lat'] - lat)**2))
    depth = -slab['dep'][dep_pos]
    this_depth = event['st00'].depth
    difference = depth - this_depth
    
    if difference > 20:
        event_type = 'crustal'
    elif difference > -20:
        event_type = 'interface'
    else:
        event_type = 'intraslab'
    
    for key, station in event.items():
        new_station = {}
        for attribute in sorted(dir(station)):
            if attribute.startswith('_') or attribute.startswith('R') or attribute.startswith('mecha'):
                continue
            new_station[attribute] = getattr(station, attribute)
        
        new_station['event_type'] = event_type
        new_station['Rhypo'] = distances[key][0]
        new_station['Repi'] = distances[key][1]
        new_station['Rrup'] = distances[key][2]
        new_station['Rjb'] = distances[key][3]
        new_station['units'] = 'Acceleration: m/s/s; Magnitude: Mw; dt: s; Depth: km; Vs30: m/s; Rhypo: km; Repi: km; Rrup: km; Rjb: km'
        new_event[key] = new_station
        
    spio.savemat(os.path.join(seismic_path, filename), new_event)
    # break