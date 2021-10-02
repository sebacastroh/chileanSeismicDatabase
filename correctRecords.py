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
import scipy.integrate as spin
from obspy.signal import filter as flt

def correct_seismogram(acc, p_wave, dt, tmp_filename=None):
    
    n = len(acc)
    t = np.linspace(0., (n-1)*dt, n)

    # Zero completation and band-pass filtering
    # Tmax = 100.
    nn = len(acc) - p_wave
    # mm = int(Tmax/dt)
    # window = spsig.tukey(nn, alpha = 0.005)

    # new_acc = np.hstack((np.zeros(mm), acc[p_wave:]*window))
    new_acc = acc - acc.mean()
    new_acc = new_acc - new_acc[:p_wave].mean()

    freq_min = 0.01 # Hz
    freq_max = 20. # Hz
    fsamp = 1./dt # Hz
    order = 4
    zerophase = False

    # fil_acc = np.hstack((np.zeros(p_wave),
    #         flt.bandpass(new_acc, freq_min, freq_max, fsamp,
    #                        corners=order, zerophase=zerophase)[mm:]))
    fil_acc = flt.bandpass(new_acc, freq_min, freq_max, fsamp,
                            corners=order, zerophase=zerophase)

    fil_vel = spin.cumtrapz(fil_acc, dx=dt, initial=0.)

    # Base line correction on velocity
    alpha = 0.2
    vel_mean = np.convolve(fil_vel, np.ones(int(fsamp/2.)+1)/(int(fsamp/2.)+1), 'same')
    vel_mean2 = np.convolve(fil_vel**2, np.ones(int(fsamp/2.)+1)/(int(fsamp/2.)+1), 'same')
    std = np.sqrt(vel_mean2 - vel_mean**2)
    peaks = spsig.find_peaks(std, distance=int(2./dt))[0]
    smooth_std = np.interp(t, t[peaks], std[peaks])*alpha
    # smooth_std[:p_wave] = 0.

    # beta = 0.1
    energy = spin.cumtrapz(new_acc**2, dx=dt, initial=0.)
    energy /= energy[-1]
    # p = np.where((energy >= beta/2.) & (energy <= 1.- beta/2.))[0]
    # mask = np.ones(len(energy), dtype=bool)
    # mask[p] = 0
    # smooth_std[p] = np.max(smooth_std[mask])

    dataset_pre = rjmcmc.dataset1d(*np.c_[t[:p_wave], vel_mean[:p_wave], smooth_std[:p_wave]].T.tolist())
    dataset_post = rjmcmc.dataset1d(*np.c_[t[:nn], vel_mean[p_wave:], smooth_std[p_wave:]].T.tolist())
    # dataset = rjmcmc.dataset1d(*np.c_[t, vel_mean, smooth_std].T.tolist())

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
        # xsamples = int(percentage*n)
        # ysamples = int(percentage*n)
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

    solution = np.hstack((solution_pre, solution_post))

    vel_corr = fil_vel - solution
    acc_corr = np.gradient(vel_corr, dt, edge_order=2)
    # dis_corr = spin.cumtrapz(vel_corr, dx=dt, initial=0.)

    # vel = spin.cumtrapz(acc, dx=dt, initial=0.)
    # dis = spin.cumtrapz(vel, dx=dt, initial=0.)

    # plt.figure()
    # plt.subplot(311)
    # plt.plot(t, acc)
    # plt.subplot(312)
    # plt.plot(t, vel)
    # plt.subplot(313)
    # plt.plot(t, dis)
    # plt.suptitle('Original')

    # plt.figure()
    # plt.plot(t, fil_vel)
    # plt.plot(t, solution)
    # plt.title('Solution proposed')

    # plt.figure()
    # plt.subplot(311)
    # plt.plot(t, acc_corr)
    # plt.subplot(312)
    # plt.plot(t, vel_corr)
    # plt.subplot(313)
    # plt.plot(t, dis_corr)
    # plt.suptitle('Result')


    # solution = np.hstack((np.zeros(p_wave),
    #                       np.interp(t[p_wave:], np.array(sample_x[best])+t[p_wave], sample_y[best])))


    # vo = fil_vel[0]
    # vp = solution[p_wave]
    # tp = t[p_wave]

    # X = (t[:p_wave]**2 - tp*t[:p_wave])
    # y = (fil_vel[:p_wave] - vo - (vp-vo)/tp*t[:p_wave])

    # a = np.dot(X,y)/np.dot(X,X)
    # b = (vp - vo - a*tp**2)/tp
    # c = vo

    # solution = np.interp(t, np.array(sample_x[best]), sample_y[best])
    # solution3 = np.hstack((np.interp(t[:p_wave], [t[0], t[p_wave]], [fil_vel[0], solution[p_wave]]),
    #                        np.interp(t[p_wave:], np.array(sample_x[best]+t[p_wave]), sample_y[best])))

    # solution2 = np.hstack((a*t[:p_wave]**2 + b*t[:p_wave] + c,
    #                        np.interp(t[p_wave:], np.array(sample_x[best]+t[p_wave]), sample_y[best])))



    # vel_corr = fil_vel - solution
    # vel_corr1 = fil_vel - solution2
    # vel_corr2 = fil_vel - solution3

    # if p_wave != 0:
    #     vel_corr[p_wave:] *= window
    # acc_corr = np.gradient(vel_corr, dt, edge_order=2)

    # fil_acc = flt.bandpass(acc_corr, freq_min, freq_max, fsamp,
    #                        corners=order, zerophase=zerophase)

    # vel_corr = spin.cumtrapz(acc_corr, dx=dt, initial=0.)
    # dis_corr = spin.cumtrapz(vel_corr, dx=dt, initial=0.)
    
    if tmp_filename is not None:
        np.save(tmp_filename, acc_corr)
        return None
    else:
        return acc_corr

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
