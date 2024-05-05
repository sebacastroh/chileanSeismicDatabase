# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 11:02:49 2022

@author: sebacastroh
"""

import os
import sys
import scipy
import pyrjmcmc
import numpy as np
import scipy.stats as spsts
import scipy.signal as spsig
import scipy.integrate as spin

from obspy.signal import filter as flt
import obspy.signal.konnoohmachismoothing as kon

spV = scipy.__version__

def sampling_function(x, boundaries, local_values):
    
    locations = []
    for xi in boundaries[:-1]:
        pos = np.argmax(x >= xi)
        locations.append(pos)
    locations.append(len(x))
    
    y = np.empty_like(x)
    
    for i in range(len(locations)-1):
        p1 = locations[i]
        p2 = locations[i+1]
        y[p1:p2] = local_values[i,0]*x[p1:p2] + local_values[i,1]
    
    return y

def correctRecord(acc, dt, status, p_wave, saveInTemp=False, filename=''):
    if len(acc) == 0 or not status:
        acc_corr = np.empty(0)
        acc_fil = np.empty(0)
        freqs = np.array([np.nan, np.nan])
        
        if saveInTemp:
            np.savez(filename, acc_corr = acc_corr, acc_fil = acc_fil, freqs = freqs)
        
        return acc_corr, acc_fil, freqs

    fsamp = 1./dt
    
    n = len(acc)
    t = np.linspace(0., (n-1)*dt, n)
    
    # Zero completation
    if p_wave > 0:
        new_acc = acc - acc[:p_wave].mean()
    else:
        new_acc = acc.copy()
    
    new_vel = spin.cumtrapz(new_acc, dx=dt, initial=0.)
    if p_wave > 0:
        new_vel -= new_vel[:p_wave].mean()
    
    # Base line correction on velocity
    alpha = 0.2
    fsamp = 1./dt # Hz
    vel_mean = np.convolve(new_vel, np.ones(int(fsamp/2.)+1)/(int(fsamp/2.)+1), 'same')
    vel_mean2 = np.convolve(new_vel**2, np.ones(int(fsamp/2.)+1)/(int(fsamp/2.)+1), 'same')
    vel_dif = vel_mean2 - vel_mean**2
    vel_dif[vel_dif < 0.] = 0.
    std = np.sqrt(vel_dif)
    peaks = spsig.find_peaks(std, distance=int(2./dt))[0]
    smooth_std = np.interp(t, t[peaks], std[peaks])*alpha
    
    energy = spin.cumtrapz(new_acc**2, dx=dt, initial=0.)
    energy /= energy[-1]
    
    y = [p_wave]
    y.extend([np.argmin(np.abs(energy-x/10.)) for x in range(1,11)])
    y = np.array(y)  
    min_partitions = 1  
    max_partitions = int(max(10, np.sum(np.int64((t[y[1:]]-t[y[:-1]])/5.))))
    
    burnin = 10000
    total = 100000
    weights = None
    
    # Pre
    non_zero_vel = new_vel[:p_wave+1].max() > new_vel[:p_wave+1].min()
    if p_wave > 0 and non_zero_vel:
        datasets = [np.column_stack((t[:p_wave+1], new_vel[:p_wave+1], smooth_std[:p_wave+1]))]
        samples = 100
        pd = (t[p_wave])*0.1
        
        results1 = pyrjmcmc.multiple_partition_connected_lines(datasets,
                                                                burnin,
                                                                total,
                                                                min_partitions,
                                                                max_partitions,
                                                                samples,
                                                                pd,
                                                                seed=None,
                                                                weights=weights)
        
        if spV == '1.9.1':
            mode1 = spsts.mode(results1.npartitions, keepdims=True)
        else:
            mode1 = spsts.mode(results1.npartitions)
        
        pos1 = np.where(results1.npartitions == mode1.mode)
        sol1 = pos1[0][results1.misfit[pos1].argmin()]
        
        b1 = np.array(results1.boundaries[sol1])
        y1 = sampling_function(b1, results1.boundaries[sol1], results1.local_parameters[sol1])
        lp1 = results1.local_parameters[sol1].copy()
    
    # Post
    datasets = [np.column_stack((t[p_wave:], new_vel[p_wave:], smooth_std[p_wave:]))]
    samples = 0
    pd = (t[-1]-t[p_wave])*0.1
    
    results2 = pyrjmcmc.multiple_partition_connected_lines(datasets,
                                                            burnin,
                                                            total,
                                                            min_partitions,
                                                            max_partitions,
                                                            samples,
                                                            pd,
                                                            seed=None,
                                                            weights=weights)
    
    if spV == '1.9.1':
        mode2 = spsts.mode(results2.npartitions, keepdims=True)
    else:
        mode2 = spsts.mode(results2.npartitions)
        
    pos2 = np.where(results2.npartitions == mode2.mode)
    sol2 = pos2[0][results2.misfit[pos2].argmin()]
    
    b2 = np.array(results2.boundaries[sol2])
    y2 = sampling_function(b2, results2.boundaries[sol2], results2.local_parameters[sol2])
    lp2 = results2.local_parameters[sol2].copy()
    
    # Adjust solutions
    if p_wave > 0 and non_zero_vel:
        ym = 0.5*(y1[-1] + y2[0])
        
        m1 = (ym - y1[-2])/(b1[-1] - b1[-2])
        n1 = ym - m1*b1[-1]
        lp1[-1] = np.array([m1, n1])
        
        m2 = (y2[1] - ym)/(b2[1] - b2[0])
        n2 = ym - m2*b2[0]
        lp2[0] = np.array([m2, n2])
        
        lp = np.vstack((lp1, lp2))
        
        boundaries = np.hstack((b1, b2[1:]))
    elif p_wave > 0 and not non_zero_vel:
        ym = new_vel[:p_wave+1].max()
        m2 = (y2[1] - ym)/(b2[1] - b2[0])
        n2 = ym - m2*b2[0]
        lp2[0] = np.array([m2, n2])

        lp = np.vstack(([0, ym], lp2))

        boundaries = np.hstack(([0], b2))
    else:
        lp = results2.local_parameters[sol2].copy()
        boundaries = results2.boundaries[sol2]
    
    solution = sampling_function(t, boundaries, lp)
    
    vel_corr = new_vel - solution
    acc_corr = np.gradient(vel_corr, dt, edge_order=2)
    
    if p_wave > 0 and non_zero_vel:
        acc_corr = acc_corr - acc_corr[:p_wave].mean()
        signal = acc_corr[p_wave:]
        noise = acc_corr[:p_wave]
        
        S = np.fft.fft(signal)
        freqS = np.fft.fftfreq(S.size, d=dt)
        posS = np.where(freqS >= 0.)
        
        N = np.fft.fft(noise)
        freqN = np.fft.fftfreq(N.size, d=dt)
        posN = np.where(freqN >= 0.)
        
        Ssmooth = kon.konno_ohmachi_smoothing(np.abs(S[posS]), freqS[posS], normalize=True)
        Nsmooth = kon.konno_ohmachi_smoothing(np.abs(N[posN]), freqN[posN], normalize=True)
        
        freq = np.logspace(-2, 0, 1001)
        Sint = np.interp(freq, freqS[posS], Ssmooth)
        Nint = np.interp(freq, freqN[posN], Nsmooth)
        
        SNR = Sint/Nint
        
        f0 = 0.01
        f1 = freqN[1]
        peaks = spsig.find_peaks(SNR)[0]
        if len(peaks) > 0:
            f2 = freq[peaks[0]]
        else:
            f2 = 0.
        
        pos = np.where(SNR < 3.)[0]
        if len(pos) > 0:
            f3 = freq[pos[-1]]
        else:
            f3 = freq[np.argmin(SNR)]
            
        freq_min = np.max([f0, f1, f2, f3])
        
        freq = np.logspace(0, 2, 1001)
        Sint = np.interp(freq, freqS[posS], Ssmooth)
        Nint = np.interp(freq, freqN[posN], Nsmooth)
        
        SNR = Sint/Nint
        pos = np.where(SNR < 3.)[0]
        if len(pos) > 0:
            freq_max = min(max(freq[pos[0]], 30.), fsamp/2.)
        else:
            freq_max = fsamp/2.
    else:
        freq_min = 0.01
        freq_max = fsamp/2.
    
    order = 4
    zerophase = True
    
    if freq_max >= fsamp/2.:
        acc_fil = flt.highpass(acc_corr, freq_min, fsamp, corners=order, zerophase=zerophase)
    else:
        acc_fil = flt.bandpass(acc_corr, freq_min, freq_max, fsamp, corners=order, zerophase=zerophase)
        
    freqs = np.array([freq_min, freq_max])
    
    if saveInTemp:
        np.savez(filename, acc_corr = acc_corr, acc_fil = acc_fil, freqs = freqs)
    
    return acc_corr, acc_fil, freqs


