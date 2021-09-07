#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 17:53:27 2019

@author: srcastro
"""
import os
import pickle
import scipy.io as spio
import SeismicCorrectionLibrary

path = os.path.join(os.getcwd(), 'events_mat_uncorrected')

with open('p_waves.pkl', 'rb') as f:
    data = pickle.load(f)
    
for i,row in enumerate(data):
    if row[-2] == 'Valid':
        event = spio.loadmat(row[0], struct_as_record=False, squeeze_me=False)
        
        station = event[row[1]]
        pos = SeismicCorrectionLibrary.PWaveDetection(station.acc_3, 1./station.dt)
        
        data[i][4] = pos
