#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 14:55:24 2020

@author: srcastro
"""
import os
import pickle
import numpy as np
import correctRecords
import scipy.io as spio
import scipy.integrate as spin
import matplotlib.pyplot as plt

plt.ioff()

with open('p_waves.pkl', 'rb') as f:
    data = pickle.load(f)
    
choices = np.random.choice(np.arange(len(data)), 100, replace=True)

path = '/home/srcastro/projects/correctSeismicDatabase/newMethodology/events_mat_uncorrected'
path_c = '/home/srcastro/projects/correctSeismicDatabase/newMethodology/events_mat_corrected_v2'
# filename = '19850303_7.9Mw_33.95S_71.71W_40.7KM.mat'
# filename = '20200430_5.3Mw_22.67S_68.22W_126KM.mat'
filename = '20161225_7.6Mw_43.52S_74.39W_30KM.mat'

for choice in choices:
    filename, key, name, p_wave, status, pType, value = data[choice]
    if status != 'Valid':
        continue
    full_path = os.path.join(path, filename)
    full_path_c = os.path.join(path_c, filename)
    
    event = spio.loadmat(full_path, struct_as_record=False, squeeze_me=True)
    event_c = spio.loadmat(full_path_c, struct_as_record=False, squeeze_me=True)
    
    station = event[key]
    station_c = event_c[key]
    dt = station.dt
    n = len(station.acc_1)
    t = np.linspace(0., (n-1)*dt, n)
    for i in range(3):
        acc = getattr(station, 'acc_%i' %(i+1))
        
        if np.allclose(acc, np.zeros_like(acc)):
            continue
        acc_corr = correctRecords.correct_seismogram(acc, p_wave, dt)
        
        acc_c = getattr(station_c, 'acc_%i' %(i+1))
        
        vel_corr = spin.cumtrapz(acc_corr, t, initial=0.)
        vel_c = spin.cumtrapz(acc_c, t, initial=0.)
        
        fig1 = plt.figure()
        a11 = fig1.add_subplot(211)
        a11.plot(t, acc_corr, lw=1)
        a12 = fig1.add_subplot(212)
        a12.plot(t, acc_c, lw=1)
        fig1.suptitle(filename[:-4] + '\nStation: %s - Component %i' %(station.name, i+1))

        fig2 = plt.figure()
        a21 = fig2.add_subplot(211)
        a21.plot(t, vel_corr, lw=1)
        a22 = fig2.add_subplot(212)
        a22.plot(t, vel_c, lw=1)
        fig2.suptitle(filename[:-4] + '\nStation: %s - Component %i' %(station.name, i+1))
        
        plt.show()
        
        a11.cla()
        a12.cla()
        fig1.clf()
        
        a21.cla()
        a22.cla()
        fig2.clf()
        
        plt.close('all')
        
        
with open('p_waves.pkl', 'rb') as f:
    data = pickle.load(f)

path = '/home/srcastro/projects/correctSeismicDatabase/newMethodology/events_mat_uncorrected'
path_c = '/home/srcastro/projects/correctSeismicDatabase/newMethodology/events_mat_corrected_v2'

for row in data:
    filename, key, name, p_wave, status, pType, value = row
    
    if status != 'Valid':
        continue
    
    full_path = os.path.join(path, filename)
    full_path_c = os.path.join(path_c, filename)
    
    if not os.path.exists(full_path):
        continue
    
    event = spio.loadmat(full_path, struct_as_record=False, squeeze_me=True)
    # event_c = spio.loadmat(full_path_c, struct_as_record=False, squeeze_me=True)
    
    for thisK, station in event.items():
        if thisK.startswith('__'):
            continue
        if station.name == name:
            if thisK != key:
                print(row)
            
            
new_data = []
for row in data:
    if row[0] == '20200421_4Mw_22.06S_70.44W_31KM.mat':
        continue
    new_data.append(row)
        
new_data.append(['20200421_4Mw_22.06S_70.44W_31KM.mat', 'st00', 'A07F', 12129, 'Valid', 'Manual', -1.0])
new_data.append(['20200421_4Mw_22.06S_70.44W_31KM.mat', 'st01', 'A14F', 11239, 'Valid', 'Manual', -1.0])
with open('p_waves.pkl','wb') as f:
    pickle.dump(new_data, f)