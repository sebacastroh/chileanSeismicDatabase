# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 16:13:26 2023

@author: sebacastroh
"""

import os
import json
import pyproj
import datetime
import numpy as np
import pandas as pd
import scipy.io as spio
import lib.computeDistances as computeDistances

def raw2Uncorrected():

    licensing = 'This SIBER-RISK Strong Motion Database is made available '
    licensing += 'under the Creative Commons Attribution-NonCommercial-'
    licensing += 'ShareAlike 4.0 International Public License: '
    licensing += 'https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode. '
    licensing += 'Any rights in individual contents of the database are '
    licensing += 'licensed under the Creative Commons Attribution-'
    licensing += 'NonCommercial-ShareAlike 4.0 International Public License: '
    licensing += 'https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode'
    
    cite = 'Sebastián Castro, Roberto Benavente, Jorge G. F. Crempien, Gabriel Candia, Juan Carlos de la Llera; '
    cite += 'A Consistently Processed Strong‐Motion Database for Chilean Earthquakes. '
    cite += 'Seismological Research Letters 2022;; 93 (5): 2700–2718. doi: https://doi.org/10.1785/0220200336'
    database_doi = 'https://doi.org/10.7764/datasetUC/ING-UC.1170836_1'
    
    df = pd.read_csv('siberrisk.csv', dtype={'Identificador': str})
    
    geod = pyproj.Geod(ellps='WGS84')
    
    ffm_files = os.listdir('ffm')
    
    with open('stationsInfo.json') as f:
        sinfo = json.load(f)
    
    with open('fault_plane_properties.json') as f:
        fpp = json.load(f)
    
    xChannels = ['LONGITUDINAL', '350', 'NORTH-SOUTH', '340', '290 DEG', '100',
                 '140', '170', '280', '160', '290', '100 DEGREES', 'L', 'NS', 'N-S',
                 'HNN', 'HN1', 'HLN']
    
    yChannels = ['TRANSVERSE', '80', 'EAST-WEST', '70', '200 DEG', '190', '10', '080',
                 '50', '200', '010 DEGREES', 'T', 'EW', 'E-W', 'HLE', 'HN2', 'HNE']
    
    zChannels = ['VERTICAL', 'UP', 'V', 'Z', 'NZ', 'HLZ', 'HNZ']
    
    omitChannel = ['INTC']
    
    event_ids = df['ID'].unique()
    
    if not os.path.exists('seismicDatabase'):
        os.mkdir('seismicDatabase')
    
    if not os.path.exists(os.path.join('seismicDatabase', 'npz')):
        os.mkdir(os.path.join('seismicDatabase', 'npz'))
    
    if not os.path.exists(os.path.join('seismicDatabase', 'mat')):
        os.mkdir(os.path.join('seismicDatabase', 'mat'))
    
    slab = np.load('sam_slab2.npz')
    for event_id in event_ids.tolist():
        
        if os.path.exists(os.path.join('seismicDatabase', 'npz', event_id + '.npz')):
            continue
        
        info = df[df['ID'] == event_id]
        st = 0
        event = {
            'event_id': event_id,
            'licensing': licensing,
            'cite': cite,
            'databaseURL': database_doi
        }
        
        for r, row in info.iterrows():
            filename = row['Identificador'] + '.npz'
            
            with np.load(os.path.join('rawData', filename), allow_pickle=True) as f:
                data = {}
                for key, value in f.items():
                    data[key] = value.item()
            
            event_mag = row['Magnitud [*]']
            event_lon = row['Longitud']
            event_lat = row['Latitud']
            event_dep = row['Profundidad [km]']
            
            if np.any(np.isnan([event_mag, event_lon, event_lat, event_dep])):
                event_type = 'Undetermined'
                Rhypo = np.nan
                Repi  = np.nan
                Rrup  = np.nan
                Rjb   = np.nan
                
                hypocenter = np.nan
            else:
                # Event type
                dep_pos = np.nanargmin(np.sqrt((slab['lon'] - event_lon)**2 + (slab['lat'] - event_lat)**2))
                slab_depth = -slab['dep'][dep_pos]
                difference = slab_depth - event_dep
                
                if difference > 10:
                    event_type = 'crustal'
                elif difference >= -10:
                    event_type = 'interface'
                else:
                    event_type = 'intraslab'
                    
                # Distance
                hypocenter = computeDistances.LatLonDepth2XYZNum(event_lat, event_lon, event_dep)
            
            for i, (stationCode, station) in enumerate(data.items()):
                acc_1 = np.empty(0)
                acc_2 = np.empty(0)
                acc_3 = np.empty(0)
                
                x1 = np.inf
                x2 = np.inf
                x3 = np.inf
                
                for channelCode, channel in station.items():
                    if channelCode.strip() in omitChannel:
                        continue
                    
                    location = channel.get('location')
                    if location is None:
                        location = channel.get('loc')
                    stationLon = location.get('lon')
                    stationLat = location.get('lat')
                    
                    metadata = channel.get('metadata')
                    if metadata is None:
                        metadata = channel.get('m')
                    
                    stationStarttime = metadata.get('starttime')
                    if stationStarttime is None:
                        stationStarttime = row['Fecha (UTC)']
                        
                    stationDt = metadata.get('delta')
                    if stationDt is None:
                        x = channel.get('x')
                        stationDt = np.mean(x[1:] - x[:-1])
                    
                    if channelCode.strip() in xChannels:
                        x1 = channel.get('x')[0]
                        acc_1 = channel.get('y').copy()
                        component_1 = channelCode.strip()
                    elif channelCode.strip() in yChannels:
                        x2 = channel.get('x')[0]
                        acc_2 = channel.get('y').copy()
                        component_2 = channelCode.strip()
                    elif channelCode.strip() in zChannels:
                        x3 = channel.get('x')[0]
                        acc_3 = channel.get('y').copy()
                        component_3 = channelCode.strip()
                
                xini = np.min([x1, x2, x3])
                if xini > x1:
                    n = int((xini - x1)/stationDt)
                    if n > 0:
                        acc_1 = np.hstack((np.zeros(n), acc_1))
                
                if xini > x2:
                    n = int((xini - x2)/stationDt)
                    if n > 0:
                        acc_2 = np.hstack((np.zeros(n), acc_2))
                
                if xini > x3:
                    n = int((xini - x3)/stationDt)
                    if n > 0:
                        acc_3 = np.hstack((np.zeros(n), acc_3))
                
                # Distances
                if not np.all(np.isnan(hypocenter)):
                    Rhypo = None
                    Repi  = None
                    Rrup  = None
                    Rjb   = None
                    station_xyz = computeDistances.LatLonDepth2XYZNum(stationLat, stationLon, 0.)
                    if event_mag < 7.1:
                        Rhypo = np.sqrt(np.sum((hypocenter - station_xyz)**2))
                        azimuth1, azimuth2, Repi = geod.inv(stationLon, stationLat, event_lon, event_lat)
                        Repi /= 1000.
                        Rrup = Rhypo
                        Rjb = Repi
                    else:
                        if event_id + '.geojson' in ffm_files:
                            with open(os.path.join('ffm', event_id + '.geojson')) as f:
                                ffm = json.load(f)
                                
                            max_slip = -np.inf
                            min_slip =  np.inf
                            for feature in ffm['features']:
                                slip = feature['properties']['slip']
                                max_slip = max(max_slip, slip)
                                min_slip = min(min_slip, slip)
                            
                            dslip = 0.15*(max_slip - min_slip) + min_slip
                            
                            Rhypo = np.sqrt(np.sum((hypocenter - station_xyz)**2))
                            azimuth1, azimuth2, Repi = geod.inv(stationLon, stationLat, event_lon, event_lat)
                            Repi /= 1000.
                            
                            Rrup = np.inf
                            Rjb = np.inf
                            for feature in ffm['features']:
                                slip = feature['properties']['slip']
                                if slip < dslip:
                                    continue
                                
                                points = np.array(feature['geometry']['coordinates'][0])[:-1]
                                sub_fault = computeDistances.LatLonDepth2XYZNumPy(points[:,1], points[:,0], points[:,2]/1000.)
                                dist = np.sqrt(np.sum((sub_fault - station_xyz)**2, axis=1)).min()
                                    
                                Rrup = min(Rrup, dist)
                                
                                if computeDistances.inPolygon(station_xyz[:2], sub_fault[:,:2]):
                                    Rjb = 0.
                                elif Rjb != 0.:
                                    azimuth1, azimuth2, dist = geod.inv(stationLon*np.ones(len(points)),\
                                                                        stationLat*np.ones(len(points)),\
                                                                        points[:,0], points[:,1])
                                    Rjb = min(Rjb, dist.min()/1000.)
                                    
                        else:
                            Rhypo = np.sqrt(np.sum((hypocenter - station_xyz)**2))
                            azimuth1, azimuth2, Repi = geod.inv(stationLon, stationLat, event_lon, event_lat)
                            Repi /= 1000.
                            
                            properties = fpp.get(event_id)
                            if properties is None:
                                Rrup = 'Finite fault model required'
                                Rjb = 'Finite fault model required'
                            else:
                                strike, dip, rake = properties
                                
                                L = 10**(-2.9  + 0.63*event_mag) # Allen et al 2017, Table 2
                                W = 10**(-0.86 + 0.35*event_mag) # Allen et al 2017, Table 2
                                subfaults = computeDistances.generateSubFaults(event_lat, event_lon, event_dep, strike, dip, L, W)
                                
                                faultXYZ = computeDistances.LatLonDepth2XYZNumPy(subfaults[:,1],\
                                                                                 subfaults[:,0],\
                                                                                 subfaults[:,2])
                                    
                                Rrup = np.sqrt(np.sum((faultXYZ - station_xyz)**2, axis=1)).min()
                                
                                nx = int(L)
                                ny = int(W)
                                polygon = subfaults[[0,ny-1,nx*ny-1,nx*ny-ny]]
                                
                                if computeDistances.inPolygon([stationLon, stationLat], polygon):
                                    Rjb = 0.
                                else:
                                    azimuth1, azimuth2, dist = geod.inv(stationLon*np.ones(nx*ny),\
                                                                        stationLat*np.ones(nx*ny),\
                                                                        subfaults[:,0], subfaults[:,1])
                                        
                                    Rjb = dist.min()/1000.
                
                # Station properties
                properties = sinfo.get(stationCode)
                if properties is None:
                    vs30 = np.nan
                    azimuth = np.nan
                    hvsr = 'Undetermined'
                    station_name = 'Unknown'
                else:
                    vs30 = properties[2]
                    azimuth = properties[4]
                    hvsr = 'Undetermined'
                    station_name = 'Unknown'
                
                # Save results
                station_dict = {
                    'starttime': stationStarttime,
                    'magnitude': event_mag,
                    'hypocenter_lon': event_lon,
                    'hypocenter_lat': event_lat,
                    'depth': event_dep,
                    'event_type': event_type,
                    'station_name': station_name,
                    'station_code': stationCode,
                    'station_lon': stationLon,
                    'station_lat': stationLat,
                    'acc_uncorrected_1': acc_1,
                    'acc_corrected_1': np.empty(0),
                    'acc_filtered_1': np.empty(0),
                    'component_1': component_1,
                    'acc_uncorrected_2': acc_2,
                    'acc_corrected_2': np.empty(0),
                    'acc_filtered_2': np.empty(0),
                    'component_2': component_2,
                    'acc_uncorrected_3': acc_3,
                    'acc_corrected_3': np.empty(0),
                    'acc_filtered_3': np.empty(0),
                    'component_3': component_3,
                    'dt': stationDt,
                    'p_wave': -1,
                    'units': 'Acceleration: m/s/s; Magnitude: M; dt: s; Depth: km; Vs30: m/s; Rhypo: km; Repi: km; Rrup: km; Rjb: km',
                    'Rhypo': Rhypo,
                    'Repi': Repi,
                    'Rrup': Rrup,
                    'Rjb': Rjb,
                    'vs30': vs30,
                    'hvsr': hvsr,
                    'corner_freqs_1': np.array([np.nan, np.nan]),
                    'corner_freqs_2': np.array([np.nan, np.nan]),
                    'corner_freqs_3': np.array([np.nan, np.nan]),
                    'azimuth': azimuth,
                    'last_update': datetime.datetime.now().isoformat()
                }
                event['st%0.2i' %st] = station_dict
                st += 1
            
            np.savez_compressed(os.path.join('seismicDatabase', 'npz', event_id), **event)
            spio.savemat(os.path.join('seismicDatabase', 'mat', event_id + '.mat'), event, do_compression=True)
    slab.close()
    