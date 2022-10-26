#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 18:09:24 2020

@author: srcastro
"""
import os
import pickle
import numpy as np
import multiprocessing
import scipy.io as spio
import matplotlib.pyplot as plt
import SeismicCorrectionLibrary

with open('p_waves.pkl', 'rb') as f:
    data = pickle.load(f)
    

def computeError(filename, key, name, pos, status, method, num):
    if status == 'Valid':
        if filename.startswith('1985'):
            return -999999
        path = './events_mat_uncorrected/'
        if not os.path.exists(path + filename):
            return -999999
        
        event = spio.loadmat(path + filename, struct_as_record=False, squeeze_me=True)
                
        station = event[key]
        
        this_pos = SeismicCorrectionLibrary.PWaveDetection(station.acc_3, 1./station.dt)
        
        return (pos - this_pos)*station.dt
    else:
        return -999999

if __name__ == '__main__':
    
    pool = multiprocessing.Pool(processes=56)
    errors = pool.starmap(computeError, data)
    pool.close()
    
    filtered_error = np.array([this_error for this_error in errors if this_error > -999999])
    
    plt.figure(figsize=(10.62, 9.82*0.725))
    plt.subplot(311)
    plt.hist(filtered_error, ec='k', bins= np.hstack([0., 10**np.linspace(np.log10(0.005), np.log10(100), 49)]))
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.xlabel('Absolute error  [s]')
    plt.ylabel('Number of records')
    # plt.xticks(np.hstack([10**np.linspace(np.log10(0.005), np.log10(100), 49)]))