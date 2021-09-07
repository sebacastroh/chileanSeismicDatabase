#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 19:47:54 2020

@author: srcastro
"""
import scipy.io as spio

folder = './events_mat_corrected_v2/'

filenames = ['20180121_6.2Mw_18.88S_69.61W_129KM.mat',
             '20181101_6.3Mw_19.65S_69.41W_101KM.mat',
             '20200717_5.7Mw_20.24S_70.08W_76KM.mat']

for filename in filenames:    
    event = spio.loadmat(folder + filename, struct_as_record=False, squeeze_me=True)
    event.pop('__header__')
    event.pop('__globals__')
    event.pop('__version__')
    
    new_event = {}
    this = False
    for key, station in event.items():
        new_station = {}
        name = station.name
        if name == 'PX02':
            this = True
            print('Yes')
        else:
            this = False
        
        for attribute in dir(station):
            if attribute.startswith('_'):
                continue
            if attribute == 'acc_1' and this:
                continue
            if attribute == 'acc_3' and this:
                continue
            
            new_station[attribute] = getattr(station, attribute)
        
        if this:
            new_station['acc_1'] = station.acc_3
            new_station['acc_3'] = station.acc_1
            
        new_event[key] = new_station
     
    spio.savemat(folder + filename, new_event)