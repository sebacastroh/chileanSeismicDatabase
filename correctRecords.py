#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 12:28:55 2019

@author: srcastro
"""
import os
import sys
import pickle
import rjmcmc
import numpy as np
import multiprocessing
import scipy.io as spio
import scipy.signal as spsig
import scipy.fftpack as spfft
import scipy.integrate as spin
from obspy.signal import filter as flt
import obspy.signal.konnoohmachismoothing as kon

def correct_seismogram(acc, p_wave, dt, tmp_filename=None):

    classic = False

    n = len(acc)
    t = np.linspace(0., (n-1)*dt, n)

    # Zero completation and band-pass filtering
    nn = len(acc) - p_wave

    new_acc = acc - acc.mean()

    if classic:
        new_vel = spin.cumtrapz(new_acc, dx=dt, initial=0.)
    else:
        new_vel = spfft.diff(new_acc, -1, period=n*dt)
        new_vel -= new_vel[0]

    # Base line correction on velocity
    alpha = 0.2
    fsamp = 1./dt # Hz
    vel_mean = np.convolve(new_vel, np.ones(int(fsamp/2.)+1)/(int(fsamp/2.)+1), 'same')
    vel_mean2 = np.convolve(new_vel**2, np.ones(int(fsamp/2.)+1)/(int(fsamp/2.)+1), 'same')
    std = np.sqrt(vel_mean2 - vel_mean**2)
    peaks = spsig.find_peaks(std, distance=int(2./dt))[0]
    smooth_std = np.interp(t, t[peaks], std[peaks])*alpha

    energy = spin.cumtrapz(new_acc**2, dx=dt, initial=0.)
    energy /= energy[-1]

    dataset_pre = rjmcmc.dataset1d(*np.c_[t[:p_wave], vel_mean[:p_wave], smooth_std[:p_wave]].T.tolist())
    dataset_post = rjmcmc.dataset1d(*np.c_[t[:nn], vel_mean[p_wave:], smooth_std[p_wave:]].T.tolist())

    datasets = [dataset_pre, dataset_post]
    for i in range(2):
        sample_x = []
        sample_y = []

        def callback(x, y): 
            sample_x.append(x)
            sample_y.append(y)

        pv = 0.
        burnin = 10000
        total = 100000
        min_partitions = 2
        y = [p_wave]
        y.extend([np.argmin(np.abs(energy-x/10.)) for x in range(1,11)])
        y = np.array(y)    
        max_partitions = int(max(10, np.sum(np.int64((t[y[1:]]-t[y[:-1]])/5.))))
        pd = (t[-1]-t[p_wave])*0.1
        percentage = 0.01
        xsamples = int(percentage*nn)
        ysamples = int(percentage*nn)
        credible_interval = 0.95

        results = rjmcmc.regression_part1d_natural(datasets[i],
                                                   pv,
                                                   pd,
                                                   burnin,
                                                   total,
                                                   min_partitions,
                                                   max_partitions+1,
                                                   xsamples,
                                                   ysamples,
                                                   credible_interval,
                                                   callback)

        partitions = np.array(results.partitions())[burnin:]
        counts = np.bincount(partitions)

        partitions_mode = np.argmax(counts)
        partitions_mode_locations = np.where(partitions == partitions_mode)[0]
        misfit = np.array(results.misfit())[burnin:]
        best = partitions_mode_locations[np.argmin(misfit[partitions_mode_locations])] + burnin

        if i == 0:
            solution_pre = np.interp(t[:p_wave], np.array(sample_x[best]), sample_y[best])
        else:
            solution_post = np.interp(t[p_wave:], np.array(sample_x[best])+t[p_wave], sample_y[best])

    m1 = (solution_pre[-2] - solution_pre[-1])/dt
    i = 2
    while True:
        yp = m1*(t[p_wave-i] - t[p_wave]) + solution_pre[-1]
        if yp != solution_pre[-i-1]:
            break
        i += 1
        if p_wave - i <= 0:
            i = p_wave
            break

    m2 = (solution_post[0] - solution_post[1])/dt
    j = 2
    while True:
        yp = m2*(t[p_wave+j] - t[p_wave]) + solution_post[0]
        if yp != solution_post[j]:
            break
        j += 1

        if p_wave + j >= n:
            j = n - p_wave
            break

    solution = np.hstack((solution_pre[:-i],
                          np.interp(t[p_wave-i:p_wave+j], [t[p_wave-i], t[p_wave+j]], [solution_pre[-i], solution_post[j]]),
                          solution_post[j:]))

    vel_corr = new_vel - solution

    if classic:
        acc_corr = np.gradient(vel_corr, dt, edge_order=2)
    else:
        acc_corr = spfft.diff(vel_corr, 1, period=n*dt)

    # signal = acc_corr[p_wave:]
    # noise = acc_corr[:p_wave]
    signal = acc[p_wave:]
    noise = acc[:p_wave]

    S = np.fft.fft(signal)
    freqS = np.fft.fftfreq(S.size, d=dt)
    posS = np.where(freqS >= 0.)

    N = np.fft.fft(noise)
    freqN = np.fft.fftfreq(N.size, d=dt)
    posN = np.where(freqN >= 0.)

    Ssmooth = kon.konno_ohmachi_smoothing(np.abs(S[posS]), freqS[posS], normalize=True)
    Nsmooth = kon.konno_ohmachi_smoothing(np.abs(N[posN]), freqN[posN], normalize=True)

    freq = np.logspace(-3, 0, 1001)
    Sint = np.interp(freq, freqS[posS], Ssmooth)
    Nint = np.interp(freq, freqN[posN], Nsmooth)

    SNR = Sint/Nint

    pos = np.where(SNR < 3.)[0]
    if len(pos) > 0:
        freq_min = max(freq[pos[-1]], freqN[1])
    else:
        freq_min = freqN[1]

    freq = np.logspace(0, 2, 1001)
    Sint = np.interp(freq, freqS[posS], Ssmooth)
    Nint = np.interp(freq, freqN[posN], Nsmooth)

    SNR = Sint/Nint
    pos = np.where(SNR < 3.)[0]
    if len(pos) > 0:
        freq_max = min(max(freq[pos[0]], 30.), fsamp/2.)
    else:
        freq_max = 50.

    order = 4
    zerophase = True
    acc_fil = flt.bandpass(acc_corr, freq_min, freq_max, fsamp, corners=order, zerophase=zerophase)
    # window = spsig.tukey(n, alpha = 0.005)
    # acc_fil *= window

    if classic:
        vel_fil = spin.cumtrapz(acc_fil, dx=dt, initial=0.)
        dis_fil = spin.cumtrapz(vel_fil, dx=dt, initial=0.)
        vel = spin.cumtrapz(acc, dx=dt, initial=0.)
        dis = spin.cumtrapz(vel, dx=dt, initial=0.)
    else:
        vel_fil = spfft.diff(acc_fil, -1, period=n*dt)
        vel_fil -= vel_fil[0]

        dis_fil = spfft.diff(vel_fil, -1, period=n*dt)
        dis_fil -= dis_fil[0]

        vel = spfft.diff(acc, -1, period=n*dt)
        vel -= vel[0]

        dis = spfft.diff(vel, -1, period=n*dt)
        dis -= dis[0]
    
    if tmp_filename is not None:
        np.save(tmp_filename, acc_corr)
        return None
    else:
        return acc_corr
    
import matplotlib.pyplot as plt
   
def plots(t, acc, vel, dis, new_vel, solution, freqS, posS, Ssmooth, freqN, posN, Nsmooth, acc_fil, vel_fil, dis_fil):
    plt.close('all')

    plt.figure()
    plt.subplot(311)
    plt.plot(t, acc)
    plt.grid(which='both')
    plt.subplot(312)
    plt.plot(t, vel)
    plt.grid(which='both')
    plt.subplot(313)
    plt.plot(t, dis)
    plt.grid(which='both')

    plt.suptitle('Original')

    plt.figure()
    plt.plot(t, new_vel)
    plt.plot(t, solution)
    plt.grid(which='both')
    plt.title('Solution proposed in velocity')

    plt.figure()
    plt.loglog(freqS[posS], Ssmooth, label='Event')
    plt.loglog(freqN[posN], Nsmooth, label='Noise (pre-event)')
    plt.grid(which='both')
    plt.title('Smoothed Fourier amplitude spectra (Konno Ohmachi method)')
    plt.legend()

    freq = np.logspace(-3, 2, 1001)
    Sint = np.interp(freq, freqS[posS], Ssmooth)
    Nint = np.interp(freq, freqN[posN], Nsmooth)

    # SNR = Sint/Nint
    plt.figure()
    plt.loglog(freq, Sint/Nint)
    plt.grid(which='both')
    plt.title('SNR')

    plt.figure()
    plt.subplot(311)
    plt.plot(t, acc_fil)
    plt.grid(which='both')
    plt.subplot(312)
    plt.plot(t, vel_fil)
    plt.grid(which='both')
    plt.subplot(313)
    plt.plot(t, dis_fil)
    plt.grid(which='both')

    plt.suptitle('Corrected')

def correctStation(station, p_wave):
    acc_1 = station.acc_1
    acc_2 = station.acc_2
    acc_3 = station.acc_3
    dt = station.dt
    
    if np.allclose(acc_1, np.zeros_like(acc_1)):
        new_acc_1 = acc_1.copy()
    else:
        new_acc_1 = correct_seismogram(acc_1, p_wave, dt)
        
    if np.allclose(acc_2, np.zeros_like(acc_2)):
        new_acc_2 = acc_2.copy()
    else:
        new_acc_2 = correct_seismogram(acc_2, p_wave, dt)
        
    if np.allclose(acc_3, np.zeros_like(acc_3)):
        new_acc_3 = acc_3.copy()
    else:
        new_acc_3 = correct_seismogram(acc_3, p_wave, dt)
    
    new_station = {}
    for item in station.__dict__.items():
        if item[0] == 'acc_1':
            new_station['acc_1'] = new_acc_1
        elif item[0] == 'acc_2':
            new_station['acc_2'] = new_acc_2
        elif item[0] == 'acc_3':
            new_station['acc_3'] = new_acc_3
        else:
            new_station[item[0]] = item[1]
            
    return new_station

def correctStationSingleInput(inp):
    return correctStation(*inp)
    
def correctEvent(filename, p_waves, old_path, new_path):
    
    stations = spio.loadmat(os.path.join(old_path, filename),
                            squeeze_me=True,
                            struct_as_record=False)
    
    n = len(p_waves)
    valid_keys = []
    inputs = []
    
    for i in range(n):
        key = 'st%0.2i' %i
        if p_waves[key][2] == 'Valid':
            valid_keys.append(key)
            inputs.append([stations[key], p_waves[key][1]])
    
    if len(valid_keys) == 0:
        return 'No valid records'
    else:
        try:
            if len(inputs) == 1:
                new_stations = [correctStationSingleInput(inp) for inp in inputs]
            else:
                pool = multiprocessing.Pool(processes=len(inputs))
                new_stations = pool.map(correctStationSingleInput, inputs)
                pool.close()
            # new_stations = [correctStationSingleInput(inp) for inp in inputs]
        except:
            return 'Error'
        
        new_file = {}
        for i, new_station in enumerate(new_stations):
            new_key = 'st%0.2i' %i
            new_file[new_key] = new_station
        spio.savemat(os.path.join(new_path, filename), new_file)
        if n > 1:
            return '%i of %i records corrected' %(len(valid_keys), n)
        else:
            return '%i of %i record corrected' %(len(valid_keys), n)
        
def correctEventSingleInput(inp):
    return correctEvent(*inp)

def correctEventShell(filename, save_path):
    if os.path.exists(os.path.join(save_path, filename)):
        return None
    else:
        name_pkl = 'p_waves.pkl'
        with open(name_pkl, 'rb') as f:
            results = pickle.load(f)
        a = {result[0] : {} for result in results}
        {a[result[0]].update({result[1] : result[2:]}) for result in results}
        p_waves = a[filename]
        old_path = os.path.join(os.getcwd(), 'events_mat_uncorrected')
        
        return correctEvent(filename, p_waves, old_path, save_path)

if len(sys.argv) > 1:
    filename = sys.argv[1]
    save_path = sys.argv[2]
    result = correctEventShell(filename, save_path)
    if len(sys.argv) > 3:
        save_output = sys.argv[3]
        if save_output == 'True':
            with open(filename + '.txt', 'w') as f:
                f.write(filename + ' - ' + result)
