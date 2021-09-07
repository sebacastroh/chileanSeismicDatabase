#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 17:03:03 2020

@author: srcastro
"""
import os
import pickle
import scipy.io as spio

filenames = os.listdir('./events_mat_corrected_v2/')

if __name__ == '__main__':
    for filename in sorted(filenames):
        
        if int(filename[:8]) < 202107:
        # if not filename.startswith('20210527_4.3'):
            continue
        
        event = spio.loadmat('./events_mat_corrected_v2/' + filename, struct_as_record=False, squeeze_me=True)
        
        new_event = {}
        for key, station in event.items():
            if key.startswith('_'):
                continue
            
            new_station = {}
            
            for attribute in dir(station):
                if attribute.startswith('_'):
                    continue
                
                new_station[attribute] = getattr(station, attribute)
            
            new_event[key] = new_station
            
        with open('./events_pkl_corrected_v2/' + filename[:-4] + '.pkl', 'wb') as f:    
            pickle.dump(new_event, f)