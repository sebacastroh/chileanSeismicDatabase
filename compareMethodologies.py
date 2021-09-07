#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 12:18:56 2020

@author: srcastro
"""
import os
import sys
import pickle
import random
import rjmcmc
import numpy as np
from scipy import io
from scipy import signal
from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from obspy.signal import filter as flt
import matplotlib
import matplotlib.backends.backend_pdf

with open('p_waves.pkl', 'rb') as fopen:
    p_waves = pickle.load(fopen)
    
valids = []
for i,p_wave in enumerate(p_waves):
    waveform, key, station_name, pos, status, mechanism, k = p_wave
    if status == 'Valid' and mechanism == 'Automatic' and k != -1:
        valids.append(i)

path = os.path.join(os.getcwd(), 'events_mat_uncorrected')
plt.ioff()
def createFigure(n):
    p_wave = p_waves[n]
    waveform, key, station_name, pos, status, mechanism, k = p_wave

    fig0 = plt.figure(figsize=(16,7.75))
    ax = fig0.add_axes([0, 0, 1, 1])
    p = patches.Rectangle(
    (0.25, 0.25), 0.5, 0.5,
    fill=False, transform=ax.transAxes, clip_on=False
    )
    ax.add_patch(p)
    ax.text(0.5, 0.5, waveform + '\n' + station_name,
        horizontalalignment='center',
        verticalalignment='center',
        fontsize=20, color='red',
        transform=ax.transAxes)
    ax.set_axis_off()

    if status == 'Valid' and mechanism == 'Automatic' and k != -1:        
        fig = plt.figure(figsize=(16,7.75))
        fig2 = plt.figure(figsize=(16,7.75))        
        filename = os.path.join(path, waveform) 
        
        stations = io.loadmat(filename, struct_as_record=False, squeeze_me=True)
        stations.pop('__header__')
        stations.pop('__version__')
        stations.pop('__globals__')
        for station in stations.values():
            if station.name == station_name:
                X = station.acc_3
                dt = station.dt
                fsamp = int(1./dt)
                break
        
        N = len(X)
        L = 800
        AIlim = 0.5
        seconds = 2.
        ntrials = 5
        
        t = np.linspace(0.,(N-1)/fsamp, N)
        
        if X.mean() != 0.:
            Z = (X/X.mean())
        else:
            Z = X.copy()
            
        arias = integrate.cumtrapz(Z**2, dx=1./fsamp, initial=0.)
        arias /= arias.max()
        
        energy_pos = np.where(arias <= AIlim)[0][-1]
        
        if len(Z) <= L:
            L = int(len(Z)/3)
        
        if energy_pos <= L:
            energy_pos = len(arias)
        
        ycorr = np.ones(len(arias))
        ycorr[L-1:] = 0.
        
        indexes = np.array([np.arange(x, x+L) for x in range(N-L+1)])
        for k in range(ntrials):
            Y = np.random.rand(len(arias))
            ycorr[L-1:] += 1./np.sqrt(np.sum((Z[indexes] - Y[indexes])**2, axis=1)/(L-1))
        
        ycorr[L-1:energy_pos] = ycorr[L-1:energy_pos]/ntrials
        ycorr[L-1:energy_pos] = ycorr[L-1:energy_pos]/ycorr[L-1:energy_pos].max()
        xcorr = np.arange(len(arias))
        
        ysmooth = np.convolve(ycorr, np.ones(int(fsamp/2.)+1)/(int(fsamp/2.)+1), 'same')
        mean = ysmooth[L-1:energy_pos].mean()
        
        if ysmooth[energy_pos-1] > mean:
            aux = np.where(ysmooth[:energy_pos] <= mean)[0][-1]
            p1 = np.where(ysmooth[:aux] >= mean)[0][-1]
        else:
            aux = energy_pos
            p1 = np.where(ysmooth[:energy_pos] >= mean)[0][-1]
            
        p2 = np.argmin(np.abs(xcorr - (p1 - seconds*fsamp)))
        if p2 < 0:
            p2 = 0
        
        p3 = np.argmin(ysmooth[p2:aux]) + p2
        if p3 == p2:
            p3 += 1
            
        peaks, _ = signal.find_peaks(ysmooth[p2:p3], distance = int(fsamp))
        
        peaks2, _ = signal.find_peaks(-(ysmooth[p3] - ysmooth[p2:p3])/(xcorr[p3] - xcorr[p2:p3]), distance = int(fsamp))
        
        peaks3, _ = signal.find_peaks(-(ysmooth[p2+1:p3+1] - ysmooth[p2:p3])/(xcorr[p2+1:p3+1] - xcorr[p2:p3]), distance = int(fsamp))
        
        p5 = np.hstack((peaks, peaks2, peaks3))
        p5 = np.unique(p5)
        
        p5 += p2
        
        p4 = np.argmin(np.divide(arias[p5], arias[p5+int(fsamp)], out=np.ones(len(p5)), where=arias[p5+int(fsamp)]!=0.))
        
        pos = p5[p4]
        
        i = random.randint(1,3)
        
        X = getattr(station, 'acc_%i' %i).copy()
        
        a1 = fig.add_subplot(4,3,1)
        a1.plot(t, X/9.81)
        a1.set_ylabel('a(t) [g]', fontsize=8)
        
        V = integrate.cumtrapz(X, dx=dt, initial=0.)
        a2 = fig.add_subplot(4,3,2)
        a2.plot(t, V)
        a2.set_ylabel('v(t) [m/s]', fontsize=8)
        
        D = integrate.cumtrapz(V, dx=dt, initial=0.)
        a3 = fig.add_subplot(4,3,3)
        a3.plot(t,D)
        a3.set_ylabel('d(t) [m]', fontsize=8)
        
        Tmax = 100.
        nn = len(X) - pos
        mm = int(Tmax/dt)
        
        nn = len(X) - pos
        alpha = 0.005
        window = signal.tukey(nn, alpha=alpha)
        
        new_tr = X.copy()
        new_tr[pos:] = new_tr[int(pos):]*window
        new_tr[:pos] = 0
        
        new_acc = np.hstack((np.zeros(mm), X[pos:]*window))
                
        freq_min = 0.01 # Hz
        freq_max = 20. # Hz
        fsamp = 1./dt # Hz
        order = 4
        zerophase = False
        
        fil_acc = np.hstack((np.zeros(pos),
                flt.bandpass(new_acc, freq_min, freq_max, fsamp,
                               corners=order, zerophase=zerophase)[mm:]))
        
        fil_vel = integrate.cumtrapz(fil_acc, dx=dt, initial=0.)
        fil_dis = integrate.cumtrapz(fil_vel, dx=dt, initial=0.)
        
        a4 = fig.add_subplot(4,3,4)
        a4.plot(t, fil_acc/9.81)
        a4.set_ylabel(r'$\hat{a}(t)$ [g]', fontsize=8)
        
        a5 = fig.add_subplot(4,3,5)
        a5.plot(t, fil_vel)
        a5.set_ylabel(r'$\hat{v}(t)$ [m/s]', fontsize=8)
        
        a6 = fig.add_subplot(4,3,6)
        a6.plot(t, fil_dis)
        a6.set_ylabel(r'$\hat{d}(t)$ [m]', fontsize=8)
        
        # Base line correction on velocity
        alphas = [1., 1., 0.2]
        min_y = np.inf
        max_y = -min_y
        for j in range(3):
            alpha = alphas[j]
            vel_mean = np.convolve(fil_vel, np.ones(int(fsamp/2.)+1)/(int(fsamp/2.)+1), 'same')
            vel_mean2 = np.convolve(fil_vel**2, np.ones(int(fsamp/2.)+1)/(int(fsamp/2.)+1), 'same')
            std = np.sqrt(vel_mean2 - vel_mean**2)
            peaks = signal.find_peaks(std, distance=int(2./dt))[0]
            smooth_std = np.interp(t, t[peaks], std[peaks])*alpha
            smooth_std[:pos] = 0.
            
            beta = 0.1
            energy = integrate.cumtrapz(fil_acc**2, dx=dt, initial=0.)
            energy /= energy[-1]
            if j > 0:
                p = np.where((energy >= beta/2.) & (energy <= 1.- beta/2.))[0]
                mask = np.ones(len(energy), dtype=bool)
                mask[p] = 0
                smooth_std[p] = np.max(smooth_std[mask])
            
            a_sigma = fig.add_subplot(4,3,7+j)
            a_sigma.plot(t, smooth_std)
            a_sigma.set_ylabel(r'$\sigma(t)$ [m/s] - Method %i' %(j+1), fontsize=8)
            
            min_y = min(list(a_sigma.get_ylim()) + [min_y])
            max_y = max(list(a_sigma.get_ylim()) + [max_y])
            
            dataset = rjmcmc.dataset1d(*np.c_[t[:nn], fil_vel[pos:], smooth_std[pos:]].T.tolist())
            
            sample_x = []
            sample_y = []
            
            def callback(x, y): 
                sample_x.append(x)
                sample_y.append(y)
            
            pv = 0.
            burnin = 10000
            total = 100000
            min_partitions = 2
            y = [pos]
            y.extend([np.argmin(np.abs(energy-x/10.)) for x in range(1,11)])
            y = np.array(y)    
            max_partitions = max(10, np.sum(np.int64((t[y[1:]]-t[y[:-1]])/5.)))    
            pd = (t[-1]-t[pos])*0.1
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
            
            solution = np.hstack((np.zeros(pos),
                                  np.interp(t[pos:], np.array(sample_x[best])+t[pos], sample_y[best])))
            
            a_vel = fig.add_subplot(4,3,10+j)
            a_vel.plot(t, fil_vel, t, solution)
            a_vel.set_ylabel(r'$\hat{v}(t)$ [m/s]', fontsize=8)
            a_vel.set_xlabel('Time [s]', fontsize=8)
            
            vel_corr = fil_vel - solution
            vel_corr[pos:] *= window
            acc_corr = np.gradient(vel_corr, dt, edge_order=2)
            dis_corr = integrate.cumtrapz(vel_corr, dx=dt, initial=0.)
            
            a_corr = fig2.add_subplot(3,3,1+j)
            a_corr.plot(t, acc_corr/9.81)
            a_corr.set_ylabel(r'$\tilde{a}(t)$ [g]', fontsize=8)

            v_corr = fig2.add_subplot(3,3,4+j)
            v_corr.plot(t, vel_corr)
            v_corr.set_ylabel(r'$\tilde{v}(t)$ [m/s]', fontsize=8)

            d_corr = fig2.add_subplot(3,3,7+j)
            d_corr.plot(t, dis_corr)
            d_corr.set_ylabel(r'$\tilde{d}(t)$ [m]', fontsize=8)
            d_corr.set_xlabel('Time [s]', fontsize=8)
        
        [axe.ticklabel_format(axis='y', style='sci', scilimits=(0,0)) for axe in fig.axes]
        [axe.set_ylim(min_y, max_y) for axe in fig.axes[6:11:2]]
        fig.subplots_adjust(wspace=0.5, hspace=0.5)

        [axe.ticklabel_format(axis='y', style='sci', scilimits=(0,0)) for axe in fig2.axes]
        fig2.subplots_adjust(wspace=0.5, hspace=0.5)
        
    #return (fig0,fig,fig2)
    pdf = matplotlib.backends.backend_pdf.PdfPages("case%i.pdf" %n)
    pdf.savefig(fig0, fmt='pdf', bbox_inches='tight', pad_inches=0)
    pdf.savefig(fig, fmt='pdf', bbox_inches='tight', pad_inches=0)
    pdf.savefig(fig2, fmt='pdf', bbox_inches='tight', pad_inches=0)
    pdf.close()

n = int(sys.argv[1])
createFigure(n)