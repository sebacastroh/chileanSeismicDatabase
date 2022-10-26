#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 12:14:42 2020

@author: srcastro
"""
import os
import numpy as np
import scipy.io as spio

path = 'events_mat_uncorrected/'
files = os.listdir('/home/srcastro/projects/correctSeismicDatabase/newMethodology/' + path)

slab = np.load('sam_slab2.npz')

for filename in files:
    event = spio.loadmat('./' + path + filename, struct_as_record=False, squeeze_me=True)
    event.pop('__globals__')
    event.pop('__version__')
    event.pop('__header__')
    new_event = {}
    
    lon = event['st00'].hypocenter_lon
    lat = event['st00'].hypocenter_lat
    
    if isinstance(lon, str):
        event_type = 'N/A'
    else:
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
        for attribute in dir(station):
            if attribute.startswith('_'):
                continue
            if attribute == 'mechanism':
                continue
            new_station[attribute] = getattr(station, attribute)
        
        new_station['event_type'] = event_type
        new_event[key] = new_station
    
    spio.savemat('./new_mat_uncorrected/' + filename, new_event)
    # break

path = 'events_mat_corrected_v2/'
files = os.listdir('/home/srcastro/projects/correctSeismicDatabase/newMethodology/' + path)

slab = np.load('sam_slab2.npz')

for filename in files:
    event = spio.loadmat('./' + path + filename, struct_as_record=False, squeeze_me=True)
    event.pop('__globals__')
    event.pop('__version__')
    event.pop('__header__')
    new_event = {}
    
    lon = event['st00'].hypocenter_lon
    lat = event['st00'].hypocenter_lat
    
    if isinstance(lon, str):
        event_type = 'N/A'
    else:
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
        for attribute in dir(station):
            if attribute.startswith('_'):
                continue
            if attribute == 'mechanism':
                continue
            new_station[attribute] = getattr(station, attribute)
        
        new_station['event_type'] = event_type
        new_event[key] = new_station
    
    spio.savemat('./new_mat_corrected_v2/' + filename, new_event)
    # break