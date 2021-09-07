# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 11:04:20 2020

@author: sebac
"""
import os
import pickle
import random
import rjmcmc
import numpy as np
from scipy import io
from scipy import signal
from scipy import integrate
import matplotlib.pyplot as plt
from obspy.signal import filter as flt

use_selected = True

if use_selected:
    #p_waves = [['20150921_4.5Mw_30.66S_71.41W_38KM.mat', 'st00', u'C11O', 10400, 'Valid', 'Automatic', 0.23500000000000001]]
    p_waves = [['20171013_4Mw_20.29S_69.14W_91KM.mat', 'st02', u'T11A', 8835, 'Valid', 'Automatic', 0.10000000000000001]]
else:
    with open('p_waves.pkl', 'rb') as fopen:
        p_waves = pickle.load(fopen)
        random.shuffle(p_waves)

path = os.path.join(os.getcwd(), 'events_mat_uncorrected')

plt.ion()

solutions = []

for p_wave in p_waves:
    waveform, key, station_name, pos, status, mechanism, k = p_wave
    if status == 'Valid' and mechanism == 'Automatic' and k != -1:
        plt.close('all')
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
        
        fig1 = plt.figure(figsize=(10.62, 9.82))
        
        a1 = fig1.add_subplot(3,3,1)
        a1.plot(t, station.acc_1/9.81, label='Component: ' + station.component_1)
        
        a2 = fig1.add_subplot(3,3,2)
        a2.plot(t, station.acc_2/9.81, label='Component: ' + station.component_2)
        
        a3 = fig1.add_subplot(3,3,3)
        a3.plot(t, station.acc_3/9.81, label='Component: ' + station.component_3)
        
        a4 = fig1.add_subplot(3,3,4)
        a4.plot(t, integrate.cumtrapz(station.acc_1, dx=dt, initial=0.), label='Component: ' + station.component_1)
        
        a5 = fig1.add_subplot(3,3,5)
        a5.plot(t, integrate.cumtrapz(station.acc_2, dx=dt, initial=0.), label='Component: ' + station.component_2)
        
        a6 = fig1.add_subplot(3,3,6)
        a6.plot(t, integrate.cumtrapz(station.acc_3, dx=dt, initial=0.), label='Component: ' + station.component_3)
        
        a7 = fig1.add_subplot(3,3,7)
        a7.plot(t, integrate.cumtrapz(integrate.cumtrapz(station.acc_1, dx=dt, initial=0.), dx=dt, initial=0.), label='Component: ' + station.component_1)
        
        a8 = fig1.add_subplot(3,3,8)
        a8.plot(t, integrate.cumtrapz(integrate.cumtrapz(station.acc_2, dx=dt, initial=0.), dx=dt, initial=0.), label='Component: ' + station.component_2)
        
        a9 = fig1.add_subplot(3,3,9)
        a9.plot(t, integrate.cumtrapz(integrate.cumtrapz(station.acc_3, dx=dt, initial=0.), dx=dt, initial=0.), label='Component: ' + station.component_3)
        
        a1.set_ylabel('Acceleration [g]')
        a4.set_ylabel('Velocity [m/s]')
        a7.set_ylabel('Displacement [m]')
        
        a7.set_xlabel('Time [s]')
        a8.set_xlabel('Time [s]')
        a9.set_xlabel('Time [s]')
        
        a1.set_title('Component ' + station.component_1)
        a2.set_title('Component ' + station.component_2)
        a3.set_title('Component ' + station.component_3)
        
        [ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0)) for ax in fig1.axes]
        
        ymin = min(a1.get_ylim() + a2.get_ylim() + a3.get_ylim())
        ymax = max(a1.get_ylim() + a2.get_ylim() + a3.get_ylim())
        
        ymax = max(abs(ymin), ymax)
        ymin = -ymax
        
        a1.set_ylim(ymin, ymax)
        a2.set_ylim(ymin, ymax)
        a3.set_ylim(ymin, ymax)
        
        ymin = min(a4.get_ylim() + a5.get_ylim() + a6.get_ylim())
        ymax = max(a4.get_ylim() + a5.get_ylim() + a6.get_ylim())
        
        ymax = max(abs(ymin), ymax)
        ymin = -ymax
        
        a4.set_ylim(ymin, ymax)
        a5.set_ylim(ymin, ymax)
        a6.set_ylim(ymin, ymax)
        
        ymin = min(a7.get_ylim() + a8.get_ylim() + a9.get_ylim())
        ymax = max(a7.get_ylim() + a8.get_ylim() + a9.get_ylim())
        
        ymax = max(abs(ymin), ymax)
        ymin = -ymax
        
        a7.set_ylim(ymin, ymax)
        a8.set_ylim(ymin, ymax)
        a9.set_ylim(ymin, ymax)
        
        fig1.subplots_adjust(hspace=0.3*1.2)
        fig1.savefig('./figures/original_record.pdf', fmt='pdf', bbox_inches='tight', pad_inches=0)
        
        fig2 = plt.figure(figsize=(10.62, 9.82))
        
        if X.mean() != 0.:
            Z = (X/X.mean())
        else:
            Z = X.copy()
        
        a1 = fig2.add_subplot(3,2,1)
        a1.plot(t, X/9.81)
        a1.set_xlabel('Time [s]')
        a1.set_ylabel('$X(t)$ [g]')
        
        min_x, max_x = a1.get_xlim()
        
        a2 = fig2.add_subplot(3,2,2)
        a2.plot(t, Z)
        a2.set_xlabel('Time [s]')
        a2.set_ylabel('$Z(t)$ [g]')
        a2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        
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
        
        a3 = fig2.add_subplot(3,2,3)
        a3.plot(t, ycorr)
        a3.set_xlabel('Time [s]')
        a3.set_ylabel(r'$\gamma(t)$')
        min_y, max_y = a3.get_ylim()
        a3.plot([t[L-1], t[L-1]], [min_y, ycorr[L-1]], 'k--')
        a3.set_ylim([min_y, max_y])
        
        extraticks = [t[L-1]]
        extraticklabels = [r'$t_1$']
        
        a3.set_xticks(list(a3.get_xticks()) + extraticks)
        a3.set_xticklabels(['%i' %num for num in a3.get_xticks()[:-1]] + extraticklabels)
        a3.set_xlim([min_x, max_x])
        
        ysmooth = np.convolve(ycorr, np.ones(int(fsamp/2.)+1)/(int(fsamp/2.)+1), 'same')
        mean = ysmooth[L-1:energy_pos].mean()
        
        a4 = fig2.add_subplot(3,2,5)
        a4.plot(t, ysmooth, label='')
        a4.plot([min_x, max_x], [mean, mean], '--', label=r'$\bar{\gamma}$')
        a4.set_xlim([min_x, max_x])
        a4.set_xlabel('Time [s]')
        a4.set_ylabel(r'$\hat{\gamma}(t)$')
        
        a5 = fig2.add_subplot(3,2,4)
        a5.plot(t, arias, label='')
        a5.plot([min_x, max_x], [arias[energy_pos], arias[energy_pos]], '--', label=r'$\bar{I}_A$')
        a5.plot([t[energy_pos], t[energy_pos]], [arias[energy_pos], min_y], 'k--')
        a5.set_xlim([min_x, max_x])
        a5.set_ylim([min_y, max_y])
        a5.set_xlabel('Time [s]')
        a5.set_ylabel(r'$I_A(t)$')
        a5.legend()
        
        extraticks = [t[energy_pos]]
        extraticklabels = [r'$t_2$']
        
        a5.set_xticks(list(a5.get_xticks()) + extraticks)
        a5.set_xticklabels(['%i' %num for num in a5.get_xticks()[:-1]] + extraticklabels)
        a5.set_xlim([min_x, max_x])
        
        if ysmooth[energy_pos-1] > mean:
            aux = np.where(ysmooth[:energy_pos] <= mean)[0][-1]
            p1 = np.where(ysmooth[:aux] >= mean)[0][-1]
        else:
            aux = energy_pos
            p1 = np.where(ysmooth[:energy_pos] >= mean)[0][-1]
        
        a4.plot([t[p1], t[p1]], [min_y, mean], 'k--', label='')
        a4.plot([t[aux], t[aux]], [min_y, mean], 'k--', label='')
        a4.set_ylim([min_y, max_y])
        a4.legend()
        
        extraticks = [t[p1], t[aux]]
        extraticklabels = [r'$t_a$', r'$t_b$']
        
        a4.set_xticks(list(a4.get_xticks()) + extraticks)
        a4.set_xticklabels(['%i' %num for num in a4.get_xticks()[:-2]] + extraticklabels)
        a4.set_xlim([min_x, max_x])
        
        p2 = np.argmin(np.abs(xcorr - (p1 - seconds*fsamp)))
        if p2 < 0:
            p2 = 0
        
        p3 = np.argmin(ysmooth[p2:aux]) + p2
        if p3 == p2:
            p3 += 1
        
        #a4.plot([t[p1], t[p1]], [0., 1.], '--')
        
        a6 = fig2.add_subplot(3,2,6)
        a6.plot(t, ysmooth)
        a6.plot([t[p2], t[p2]], [min_y, max_y], 'k--')
        a6.plot([t[p3], t[p3]], [min_y, max_y], 'k--')
        
        
        extraticks = [t[p2], t[p3]]
        extraticklabels = [r'$\hat{t}_a$', r'$\hat{t}_b$']
        
        a6.set_xticks(list(a6.get_xticks()) + extraticks)
        a6.set_xticklabels(['%i' %num for num in a6.get_xticks()[:-2]] + extraticklabels)
        a6.set_xlim([min_x, max_x])
        a6.set_ylim([min_y, max_y])
        a6.set_ylabel(r'$\hat{\gamma}(t)$')
        
        a1.set_title('(a)')
        a2.set_title('(b)')
        a3.set_title('(c)')
        a4.set_title('(e)')
        a5.set_title('(d)')
        a6.set_title('(f)')
        
        a4.set_xlabel('Time [s]')
        a6.set_xlabel('Time [s]')
        
        fig2.subplots_adjust(hspace=0.3*1.5)
        fig2.savefig('./figures/p_wave.pdf', fmt='pdf', bbox_inches='tight', pad_inches=0)
        
        top = fig2.subplotpars.top
        bottom = fig2.subplotpars.bottom
        left = fig2.subplotpars.left
        right = fig2.subplotpars.right
        hspace = fig2.subplotpars.hspace
        wspace = fig2.subplotpars.wspace
        
        fig2.show()
        
        fig3 = plt.figure(figsize=(10.62, 9.82))
        #    fig.suptitle('Event name: ' + station.event_name + '\nStation: ' + station.name + ' - Vertical component', fontsize=16)
        
        peaks, _ = signal.find_peaks(ysmooth[p2:p3], distance = int(fsamp))
        
        peaks2, _ = signal.find_peaks(-(ysmooth[p3] - ysmooth[p2:p3])/(xcorr[p3] - xcorr[p2:p3]), distance = int(fsamp))
        
        peaks3, _ = signal.find_peaks(-(ysmooth[p2+1:p3+1] - ysmooth[p2:p3])/(xcorr[p2+1:p3+1] - xcorr[p2:p3]), distance = int(fsamp))
        #    peaks3 = np.hstack([np.arange(x-int(fsamp), x+-int(fsamp)) for x in peaks2])
        
        ta = int(t[p2])+1
        tb = int(t[p3])
        
        a1 = fig3.add_subplot(3,2,1)
        a1.plot(t[p2:p3], ysmooth[p2:p3])
        a1.plot(t[p2:p3][peaks], ysmooth[p2:p3][peaks], 'o')
        
        a2 = fig3.add_subplot(3,2,3)
        a2.plot(t[p2:p3], -(ysmooth[p3] - ysmooth[p2:p3])/(xcorr[p3] - xcorr[p2:p3]))
        a2.plot(t[p2:p3][peaks2], (-(ysmooth[p3] - ysmooth[p2:p3])/(xcorr[p3] - xcorr[p2:p3]))[peaks2], 'o')
        a2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        
        a3 = fig3.add_subplot(3,2,2)
        a3.plot(t[p2:p3], -(ysmooth[p2+1:p3+1] - ysmooth[p2:p3])/(xcorr[p2+1:p3+1] - xcorr[p2:p3]))
        a3.plot(t[p2:p3][peaks3], (-(ysmooth[p2+1:p3+1] - ysmooth[p2:p3])/(xcorr[p2+1:p3+1] - xcorr[p2:p3]))[peaks3], 'o')
        a3.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        
        p5 = np.hstack((peaks, peaks2, peaks3))
        p5 = np.unique(p5)
        
        p5 += p2
        
        p4 = np.argmin(np.divide(arias[p5], arias[p5+int(fsamp)], out=np.ones(len(p5)), where=arias[p5+int(fsamp)]!=0.))
        
        a4 = fig3.add_subplot(3,2,4)
        a4.plot(t[p2:p3], np.divide(arias[p2:p3], arias[p2+int(fsamp):p3+int(fsamp)], out=np.ones(p3-p2), where=arias[p2+int(fsamp):p3+int(fsamp)]!=0.))
        a4.plot(t[p2:p3][p5-p2], np.divide(arias[p5], arias[p5+int(fsamp)], out=np.arange(len(p5))[::-1]+1., where=arias[p5+int(fsamp)]!=0.), 'o')
           
        pos = p5[p4]
        values = np.divide(arias[p5], arias[p5+int(fsamp)], out=np.ones(len(p5)), where=arias[p5+int(fsamp)]!=0.)
        a4.plot(t[pos], values[p4], 'o')
        
        a4.annotate("", xy=((t[p3]-t[pos])/3.+t[pos], values[p4]), xytext=(t[pos], values[p4]), arrowprops=dict(arrowstyle="->"))
        a4.text(x=(t[p3]-t[pos])/3.+t[pos], y=values[p4], s='P wave arrival', va='center')
        
        a5 = fig3.add_subplot(3,1,3)
        a5.plot(t, X/9.81)
        a5.plot(t[pos], X[pos]/9.81, 'o')
        
        extraticks = [t[pos]]
        extraticklabels = [r'$t_P$']
        
        min_x, max_x = a5.get_xlim()
        a5.set_xticks(list(a5.get_xticks()) + extraticks)
        a5.set_xticklabels(['%i' %num for num in a5.get_xticks()[:-1]] + extraticklabels)
        a5.set_xlim([min_x, max_x])
        
        min_y_w, max_y_w = a5.get_ylim()
        a5.plot([t[pos], t[pos]], [X[pos]/9.81, min_y_w], 'k--')
        a5.set_ylim(min_y_w, max_y_w)
        
        extraticks = [t[p2], t[p3]]
        extraticklabels = [r'$\hat{t}_a$', r'$\hat{t}_b$']
        
        ticks = np.linspace(ta, tb, 6)
        if (tb-ta)%5 == 0:
            ticks = np.int64(ticks)
        
        ticks = ticks.tolist()
        
        min_x, max_x = a1.get_xlim()
        a1.set_xticks(ticks[1:-1] + extraticks)
        a2.set_xticks(ticks[1:-1] + extraticks)
        a3.set_xticks(ticks[1:-1] + extraticks)
        a4.set_xticks(ticks[1:-1] + extraticks)
        if (tb-ta)%5 == 0:
            a1.set_xticklabels(['%i' %num for num in ticks[1:-1]] + extraticklabels)
            a2.set_xticklabels(['%i' %num for num in ticks[1:-1]] + extraticklabels)
            a3.set_xticklabels(['%i' %num for num in ticks[1:-1]] + extraticklabels)
            a4.set_xticklabels(['%i' %num for num in ticks[1:-1]] + extraticklabels)
        else:
            a1.set_xticklabels(['%0.1f' %num for num in ticks[1:-1]] + extraticklabels)
            a2.set_xticklabels(['%0.1f' %num for num in ticks[1:-1]] + extraticklabels)
            a3.set_xticklabels(['%0.1f' %num for num in ticks[1:-1]] + extraticklabels)
            a4.set_xticklabels(['%0.1f' %num for num in ticks[1:-1]] + extraticklabels)
        
        a1.set_xlim([min_x, max_x])
        a2.set_xlim([min_x, max_x])
        a3.set_xlim([min_x, max_x])
        a4.set_xlim([min_x, max_x])
        a4.set_ylim([min_y, max_y])

        a1.set_title('(a)')
        a2.set_title('(c)')
        a3.set_title('(b)')
        a4.set_title('(d)')
        a5.set_title('(e)')
        
        a1.set_xlabel('Time [s]')
        a2.set_xlabel('Time [s]')
        a3.set_xlabel('Time [s]')
        a4.set_xlabel('Time [s]')
        a5.set_xlabel('Time [s]')
        
        a1.xaxis.set_label_coords(0.5, -0.19814062281804212)
        a2.xaxis.set_label_coords(0.5, -0.19814062281804212)
        a3.xaxis.set_label_coords(0.5, -0.19814062281804212)
        a4.xaxis.set_label_coords(0.5, -0.19814062281804212)
        a5.xaxis.set_label_coords(0.5, -0.19814062281804212)
        
        a1.set_ylabel(r'$\hat{\gamma}(t)$')
        a3.set_ylabel(r'd$\hat{\gamma}(t)$/d$t$')
        a2.set_ylabel(r'Secant of $\hat{\gamma}(t)$')
        a4.set_ylabel(r'$AIR(t)$')
        a5.set_ylabel(r'$X(t)$ [g]')
        
        fig3.subplots_adjust(left, bottom, right, top, wspace, hspace)
        fig3.savefig('./figures/peaks.pdf', fmt='pdf', bbox_inches='tight', pad_inches=0)
        
        fig3.show()
        
        for i in range(3):
            X = getattr(station, 'acc_%i' %(i+1)).copy()
            Tmax = 100.
            nn = len(X) - pos
            mm = int(Tmax/dt)
            
            nn = len(X) - pos
            alpha = 0.005
            window = signal.tukey(nn, alpha=alpha)
            
            new_tr = X.copy()
            new_tr[pos:] = new_tr[int(pos):]*window
            new_tr[:pos] = 0
            
            if i == 0:
                fig4 = plt.figure(figsize=(10.62, 9.82))
                
                a1 = fig4.add_subplot(3,2,1)
                a1.plot(t[pos:], window)
                a1.set_xlabel('Time [s]')
                a1.set_ylabel(r'$w(t-t_p)$')
                a1.set_title('(a)')
                
                a2 = fig4.add_subplot(3,2,2)
                a2.plot(t, new_tr/9.81)
                a2.set_title('(b)')
                a2.set_xlabel('Time [s]')
                a2.set_ylabel(r'$X_1(t)$ [g]')
                a2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                
            new_acc = np.hstack((np.zeros(mm), X[pos:]*window))
                
            freq_min = 0.01 # Hz
            freq_max = 20. # Hz
            fsamp = 1./dt # Hz
            order = 4
            zerophase = False
            
            fil_acc = np.hstack((np.zeros(pos),
                    flt.bandpass(new_acc, freq_min, freq_max, fsamp,
                                   corners=order, zerophase=zerophase)[mm:]))
            
            if i == 0:
                a4 = fig4.add_subplot(3,2,4)
                a4.plot(t, fil_acc/9.81)
                a4.set_title('(d)')
                a4.set_xlabel('Time [s]')
                a4.set_ylabel(r'$X_2(t)$ [g]')
                a4.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    
            z,p,k = signal.butter(4, [freq_min/fsamp, freq_max/fsamp], btype='band', output='zpk')
            w, h = signal.freqz_zpk(z,p,k, worN=np.hstack((0., np.logspace(-3., 2, 1000), 200))*np.pi*2./fsamp)
            
            if i == 0:
                a3 = fig4.add_subplot(3,2,3)        
                a3.semilogx(fsamp/2.* w/np.pi, np.abs(h))
                a3.set_xlabel('Frequency [Hz]')
                a3.set_ylabel(r'4th order Butterworth filter')
                a3.set_xlim([5e-4, 150.])
                a3.set_title('(c)')
     
                fig4.subplots_adjust(left, bottom, right, top, wspace, hspace)
                fig4.savefig('./figures/filtering_acc.pdf', fmt='pdf', bbox_inches='tight', pad_inches=0)
                
                fig4.show()
            
            fil_vel = integrate.cumtrapz(fil_acc, dx=dt, initial=0.)
            
            if i == 0:
                fig5 = plt.figure(figsize=(10.62, 9.82))
                [fig5.add_subplot(3,3, j+1) for j in range(6)]
            
            ax = fig5.axes[i]
            ax.plot(t, fil_vel)
            
            ax.set_xlabel('Time [s]')
            
            if i == 0:
                ax.set_ylabel(r'$V(t)$ [m/s]')
            
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            
#            fig5.subplots_adjust(left, bottom, right, top, wspace, hspace)
#            fig5.savefig('./figures/velocity_uncorr.pdf', fmt='pdf', bbox_inches='tight', pad_inches=0)
            
#            fig4.show()
            
            if not use_selected:
                inp = input('Use this record? 1: Yes 2: No\n')
            else:
                inp = 1
    
            if inp == 2:
                fig1.clf()
                fig2.clf()
                fig3.clf()
                fig4.clf()
                continue
    
            # Base line correction on velocity
            alpha = 1.#0.2
            vel_mean = np.convolve(fil_vel, np.ones(int(fsamp/2.)+1)/(int(fsamp/2.)+1), 'same')
            vel_mean2 = np.convolve(fil_vel**2, np.ones(int(fsamp/2.)+1)/(int(fsamp/2.)+1), 'same')
            std = np.sqrt(vel_mean2 - vel_mean**2)
            peaks = signal.find_peaks(std, distance=int(2./dt))[0]
            smooth_std = np.interp(t, t[peaks], std[peaks])*alpha
            smooth_std[:pos] = 0.
            
#            beta = 0.1
            energy = integrate.cumtrapz(fil_acc**2, dx=dt, initial=0.)
            energy /= energy[-1]
#            p = np.where((energy >= beta/2.) & (energy <= 1.- beta/2.))[0]
#            mask = np.ones(len(energy), dtype=bool)
#            mask[p] = 0
#            smooth_std[p] = np.max(smooth_std[mask])
            
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
            max_partitions = 150#max(10, np.sum(np.int64((t[y[1:]]-t[y[:-1]])/5.)))    
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
            
            if i == 0:
                plt.figure()
                plt.hist(partitions, range=[min_partitions, max_partitions])
            
            solution = np.hstack((np.zeros(pos),
                                  np.interp(t[pos:], np.array(sample_x[best])+t[pos], sample_y[best])))
            
            solutions.append(solution)
            
            vel_corr = fil_vel - solution
            vel_corr[pos:] *= window
            acc_corr = np.gradient(vel_corr, dt, edge_order=2)
            
#            fig5 = plt.figure(figsize=(10.62, 9.82))
#            a1 = fig5.add_subplot(311)
#            a1.plot(t, fil_vel)
            ax.plot(t, solution)
            
            ay = fig5.axes[i+3]
            ay.plot(t, vel_corr)
            
#            a1.set_xlabel('Time [s]')
#            a1.set_ylabel(r'$V(t)$ [m/s]')
            
            ay.set_xlabel('Time [s]')
            
            if i == 0:
                ay.set_ylabel(r'$\hat{V}(t)$ [m/s]')
            
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            ay.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            
            ax.set_title('Component ' + getattr(station, 'component_%i' %(i+1)))
            
            if i == 2:
                
                ylims = []
                [ylims.extend(axes.get_ylim()) for axes in fig5.axes[:3]]
                ylim = max(-min(ylims), max(ylims))
                [fig5.axes[j].set_ylim(-ylim, ylim) for j in range(3)]
                
                ylims = []
                [ylims.extend(axes.get_ylim()) for axes in fig5.axes[3:6]]
                ylim = max(-min(ylims), max(ylims))
                [fig5.axes[j].set_ylim(-ylim, ylim) for j in range(3,6)]
                
                fig5.subplots_adjust(left, bottom, right, top, wspace, hspace)
                fig5.savefig('./figures/velocity_corr.pdf', fmt='pdf', bbox_inches='tight', pad_inches=0)
            
            if i == 0:
                fig6 = plt.figure(figsize=(10.62, 9.82))
                [fig6.add_subplot(3,3,j+1) for j in range(9)]
            
            a1 = fig6.axes[i]
            a1.plot(t, acc_corr/9.81)
            
            a2 = fig6.axes[i+3]
            a2.plot(t, vel_corr)
            
            a3 = fig6.axes[i+6]
            a3.plot(t, integrate.cumtrapz(vel_corr, dx=dt, initial=0.))
            
            a1.set_title('Component ' + getattr(station, 'component_%i' %(i+1)))
            
            if i == 0:
                a1.set_ylabel(r'$\hat{X}(t)$ [g]')
                a2.set_ylabel(r'$\hat{V}(t)$ [m/s]')
                a3.set_ylabel(r'$\hat{D}(t)$ [m]')
            
            a3.set_xlabel('Time [s]')
            
            
            a1.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            a2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            a3.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            
            if i == 2:
                
                ylims = []
                [ylims.extend(axes.get_ylim()) for axes in fig6.axes[:3]]
                ylim = max(-min(ylims), max(ylims))
                [fig6.axes[j].set_ylim(-ylim, ylim) for j in range(3)]
                
                ylims = []
                [ylims.extend(axes.get_ylim()) for axes in fig6.axes[3:6]]
                ylim = max(-min(ylims), max(ylims))
                [fig6.axes[j].set_ylim(-ylim, ylim) for j in range(3,6)]
                
                ylims = []
                [ylims.extend(axes.get_ylim()) for axes in fig6.axes[6:9]]
                ylim = max(-min(ylims), max(ylims))
                [fig6.axes[j].set_ylim(-ylim, ylim) for j in range(6,9)]
                
                fig6.subplots_adjust(left, bottom, right, top, wspace, hspace)
                fig6.savefig('./figures/seismogram_corrected.pdf', fmt='pdf', bbox_inches='tight', pad_inches=0)
        
        #if not use_selected:
        #    inp = input('Use this waveform? 1: Yes 2: No\n')
        #else:
        #    inp = 1
        
        if inp == 1:
            break
        else:
            fig1.clf()
            fig2.clf()
            fig3.clf()
            fig4.clf()


np.save('solutions.npy', np.array(solutions))
#"""