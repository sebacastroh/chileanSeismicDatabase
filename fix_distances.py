#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 17:52:42 2019

@author: srcastro
"""

import os
import pickle
import scipy.io as spio
import distanceSurfaceFault

with open('fault_plane_properties.pkl', 'rb') as f:
    data = pickle.load(f)

mat_uncorrected = os.path.join(os.getcwd(), 'events_mat_uncorrected')
mat_corrected = os.path.join(os.getcwd(), 'events_mat_corrected_v2')

filenames = sorted(os.listdir(os.path.join(os.getcwd(), 'events')))
for filename in filenames:
    if filename in data.keys():
        save = True
        event_uncorrected = spio.loadmat(os.path.join(mat_uncorrected, filename + '.mat'),
                                         squeeze_me=True, struct_as_record=False)
        
        event_uncorrected.pop('__header__')
        event_uncorrected.pop('__version__')
        event_uncorrected.pop('__globals__')
        
        strike, dip, rake = data[filename]
        
        new_event_uncorrected = {}
        
        for key, station in event_uncorrected.iteritems():
            if station.Rrup == 'To determine':
                distance = distanceSurfaceFault.MinDistToFault(station.hypocenter_lat,
                                                               station.hypocenter_lon,
                                                               station.depth,
                                                               strike,
                                                               dip,
                                                               station.magnitude,
                                                               station.lon,
                                                               station.lat)
                station.Rrup = distance
            else:
                save = False
                break
            
            new_event_uncorrected[key] = station
        
        if save:
            spio.savemat(os.path.join(mat_uncorrected, filename + '.mat'), new_event_uncorrected)
        
        if os.path.exists(os.path.join(mat_corrected, filename + '.mat')):
            save = True
            event_corrected = spio.loadmat(os.path.join(mat_corrected, filename + '.mat'),
                                           squeeze_me=True, struct_as_record=False)
            
            event_corrected.pop('__header__')
            event_corrected.pop('__version__')
            event_corrected.pop('__globals__')
        
            strike, dip, rake = data[filename]
            
            new_event_corrected = {}
        
            for key,station in event_corrected.iteritems():
                if station.Rrup == 'To determine':
                    distance = distanceSurfaceFault.MinDistToFault(station.hypocenter_lat,
                                                                   station.hypocenter_lon,
                                                                   station.depth,
                                                                   strike,
                                                                   dip,
                                                                   station.magnitude,
                                                                   station.lon,
                                                                   station.lat)
                    station.Rrup = distance
                else:
                    save = False
                    break
                
                new_event_corrected[key] = station
            
            if save:
                spio.savemat(os.path.join(mat_corrected, filename + '.mat'), new_event_corrected)
