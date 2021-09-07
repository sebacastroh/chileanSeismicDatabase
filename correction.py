#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 13:07:11 2019

@author: srcastro
"""
import os
import numpy as np
import scipy.io as spio
import scipy.signal as spsig
import scipy.integrate as spin
import SeismicCorrectionLibrary
from obspy.signal import filter as flt
import rjmcmc
import matplotlib.pyplot as plt

plt.close('all')

events_path = '/home/srcastro/projects/downloadCSN/events_uncorrected_mat/'
filenames = os.listdir(events_path)

import pickle

name_pkl = 'p_waves.pkl'
with open(name_pkl, 'rb') as f:
    p_waves = pickle.load(f)

#for filename in filenames:
nsim = 2
k = 0
for value in p_waves:
    filename = value[0]
    key = value[1]
    stations = spio.loadmat(os.path.join(events_path, filename),
                            squeeze_me=True,
                            struct_as_record=False)
    station = stations[key]
    
    p_wave = value[3]
    status = value[4]
#    stations.pop('__header__')
#    stations.pop('__version__')
#    stations.pop('__globals__')
#    
#    for station in stations.itervalues():
    if status == 'Valid':
        k += 1
        acc = station.acc_1.copy()
        dt = station.dt
        n = len(acc)
        t = np.linspace(0., (n-1)*dt, n)
        vel = spin.cumtrapz(acc, dx=dt, initial=0.)
        dis = spin.cumtrapz(vel, dx=dt, initial=0.)
        
        plt.figure()
        plt.subplot(311)
        plt.plot(t, acc)
        plt.ylabel(r'Acceleration [m/$s^2$]')
        plt.title('Original record')
        plt.grid()
        
        plt.subplot(312)
        plt.plot(t, vel)
        plt.ylabel('Velocity [m/s]')
        plt.grid()
        
        plt.subplot(313)
        plt.plot(t, dis)
        plt.xlabel('Time [s]')
        plt.ylabel('Displacement [m]')
        plt.grid()
        
        # P wave detection
#        p_wave = SeismicCorrectionLibrary.PWaveDetection(acc, 1./dt)
        
        plt.figure()
        plt.plot(t, acc)
        plt.plot(t[p_wave], acc[p_wave], 'o')
        plt.xlabel('Time [s]')
        plt.ylabel(r'Acceleration [m/$s^2$]')
        plt.title('P wave detection')
        plt.grid()
        
        # Zero completation and band-pass filtering
        Tmax = 100.
        nn = len(acc) - p_wave
        mm = int(Tmax/dt)
        window = spsig.tukey(nn, alpha = 0.005)
        
        new_acc = np.hstack((np.zeros(mm), acc[p_wave:]*window))
        
        freq_min = 0.01 # Hz
        freq_max = 20. # Hz
        fsamp = 1./dt # Hz
        order = 4
        zerophase = False

        fil_acc = flt.bandpass(new_acc, freq_min, freq_max, fsamp, corners=order, zerophase=zerophase)[(mm-p_wave):]
        fil_vel = spin.cumtrapz(fil_acc, dx=dt, initial=0.)
        fil_dis = spin.cumtrapz(fil_vel, dx=dt, initial=0.)
        
        plt.figure()
        plt.subplot(311)
        plt.plot(t, fil_acc)
        plt.ylabel(r'Acceleration [m/$s^2$]')
        if zerophase:
            plt.title('Filtered record\nButterworth, order %i, acausal, freq_min = %0.2f Hz, freq_max = %0.2f Hz' %(order, freq_min, freq_max))
        else:
            plt.title('Filtered record\nButterworth, order %i, causal, freq_min = %0.2f Hz, freq_max = %0.2f Hz' %(order, freq_min, freq_max))
        plt.grid()
        
        plt.subplot(312)
        plt.plot(t, fil_vel)
        plt.ylabel('Velocity [m/s]')
        plt.grid()
        
        plt.subplot(313)
        plt.plot(t, fil_dis)
        plt.xlabel('Time [s]')
        plt.ylabel('Displacement [m]')
        plt.grid()
        
        # Base line correction on velocity
        vel_mean = np.convolve(fil_vel, np.ones(int(fsamp/2.)+1)/(int(fsamp/2.)+1), 'same')
        vel_mean2 = np.convolve(fil_vel**2, np.ones(int(fsamp/2.)+1)/(int(fsamp/2.)+1), 'same')
        std = np.sqrt(vel_mean2 - vel_mean**2)
        peaks = spsig.find_peaks(std, distance=int(2./dt))[0]
        smooth_std = np.interp(t, t[peaks], std[peaks])
        smooth_std[:p_wave] = 0.
        
        dataset = rjmcmc.dataset1d(*np.c_[t[:nn], fil_vel[p_wave:], smooth_std[p_wave:]].T.tolist())
        
        i = 0
        j = 0
        while i < 1:
            print(j)
            j += 1
            sample_x = []
            sample_y = []
            
            def callback(x, y):
                global sample_x, sample_y        
                sample_x.append(x)
                sample_y.append(y)
            
            pv = 0.
            pd = (t[-1]-t[p_wave])*0.1
            burnin = 10000
            total = 50000
            min_partitions = 2
            max_partitions = 10
            percentage = 0.01
            xsamples = int(percentage*nn)
            ysamples = int(percentage*nn)
            credible_interval = 0.95
            
            results = rjmcmc.regression_part1d_natural(dataset,
                                                       pv,
                                                       pd,
                                                       burnin,
                                                       total,
                                                       min_partitions,
                                                       max_partitions,
                                                       xsamples,
                                                       ysamples,
                                                       credible_interval,
                                                       callback)
            
            partitions = np.array(results.partitions())[burnin:]
            counts = np.bincount(partitions)
            if False:#np.any(counts > 0.2*total):
                 continue
            else:
                 i += 1
            partitions_mode = np.argmax(counts)
            partitions_mode_locations = np.where(partitions == partitions_mode)[0]
            misfit = np.array(results.misfit())[burnin:]
            best = partitions_mode_locations[np.argmin(misfit[partitions_mode_locations])] + burnin
            
            window2 = window.copy()
            window2[int(window2.shape[0]/2):] = 1.
            solution = np.hstack((np.zeros(p_wave),
                                  np.interp(t[p_wave:], np.array(sample_x[best])+t[p_wave], sample_y[best])*window2))
            
            vel_corr = fil_vel - solution
            
            plt.figure()
            plt.subplot(211)
            plt.plot(t, fil_vel, t, solution)
            plt.fill_between(t, solution + smooth_std, solution - smooth_std, alpha = 0.5, facecolor='r', edgecolor='k', linestyle='--')
            plt.xlabel('Time [s]')
            plt.ylabel('Velocity [m/s]')
            plt.title('Solution RJMCMC Regression Natural\npd = %0.2f, pv = %0.2f, samples = [%0.2f*npoints]' %(pd, pv, percentage))
            plt.grid()
            plt.subplot(212)
            plt.hist(partitions, bins = np.linspace(min_partitions, max_partitions, max_partitions-min_partitions+1, dtype=int),
                align='right', ec='k')
            plt.xlabel('Number of partitions')
            plt.ylabel('Counts')
            plt.xticks(np.arange(min_partitions, max_partitions+1, int((max_partitions-min_partitions)/2)))
            plt.grid()
            
            plt.figure()
            plt.subplot(311)
            plt.plot(t, np.gradient(vel_corr, dt, edge_order=2))
            plt.ylabel(r'Acceleration [m/$s^2$]')
            plt.title('Corrected baseline in velocity\nRJMCMC Regression Natural')
            plt.grid()
            
            plt.subplot(312)
            plt.plot(t, vel_corr)
            plt.ylabel('Velocity [m/s]')
            plt.grid()
            
            plt.subplot(313)
            plt.plot(t, spin.cumtrapz(vel_corr, dx=dt, initial=0.))
            plt.xlabel('Time [s]')
            plt.ylabel('Displacement [m]')
            plt.grid()
    if k > nsim:
        break

plt.ion()
plt.show()
