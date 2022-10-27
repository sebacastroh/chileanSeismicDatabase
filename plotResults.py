# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 23:53:06 2022

@author: sebac
"""
import os
import sys
import numpy as np
import scipy.integrate as spin
import matplotlib.pyplot as plt

if len(sys.argv) >= 1:
    eventId = sys.argv[1]

    def plots(acc, acc_corr, acc_fil, p_wave, dt, station_id, channel_id):
        vel_fil = spin.cumtrapz(acc_fil, dx=dt, initial=0.)
        if p_wave > 0:
            vel_fil -= vel_fil[:p_wave].mean()
        
        dis_fil = spin.cumtrapz(vel_fil, dx=dt, initial=0.)
        if p_wave > 0:
            dis_fil -= dis_fil[:p_wave].mean()
        
        vel = spin.cumtrapz(acc, dx=dt, initial=0.)
        if p_wave > 0:
            vel -= vel[:p_wave].mean()
        
        dis = spin.cumtrapz(vel, dx=dt, initial=0.)
        if p_wave > 0:
            dis -= dis[:p_wave].mean()
        
        vel_corr = spin.cumtrapz(acc_corr, dx=dt, initial=0.)
        if p_wave > 0:
            vel_corr -= vel_corr[:p_wave].mean()
        
        dis_corr = spin.cumtrapz(vel_corr, dx=dt, initial=0.)
        vel = spin.cumtrapz(acc, dx=dt, initial=0.)
        if p_wave > 0:
            dis_corr -= dis_corr[:p_wave].mean()
        
        
        n = len(acc)
        
        t = np.linspace(0., (n-1)*dt, n)
        
        # Figures
        
        lw = 0.25
        
        amin = min(acc.min(), acc_fil.min(), acc_corr.min())
        amax = max(acc.max(), acc_fil.max(), acc_corr.max())
        acen = max(abs(amin), abs(amax))
        
        vmin = min(vel.min(), vel_fil.min(), vel_corr.min())
        vmax = max(vel.max(), vel_fil.max(), vel_corr.max())
        vcen = max(abs(vmin), abs(vmax))
        
        dmin = min(dis.min(), dis_fil.min(), dis_corr.min())
        dmax = max(dis.max(), dis_fil.max(), dis_corr.max())
        dcen = max(abs(dmin), abs(dmax))
        
        plt.figure(figsize=[19.2 ,  9.77])
        
        plt.subplot(3,3,1)
        plt.plot(t, acc, lw=lw, c='k')
        plt.ylim([-acen, acen])
        plt.title('Uncorrected')
        
        plt.subplot(3,3,4)
        plt.plot(t, vel, lw=lw, c='k')
        plt.ylim([-vcen, vcen])
        
        plt.subplot(3,3,7)
        plt.plot(t, dis, lw=lw, c='k')
        plt.ylim([-dcen, dcen])
        
        
        plt.subplot(3,3,2)
        plt.plot(t, acc_corr, lw=lw, c='k')
        plt.ylim([-acen, acen])
        plt.title('Corrected')
        
        plt.subplot(3,3,5)
        plt.plot(t, vel_corr, lw=lw, c='k')
        plt.ylim([-vcen, vcen])
        
        plt.subplot(3,3,8)
        plt.plot(t, dis_corr, lw=lw, c='k')
        plt.ylim([-dcen, dcen])
        
        
        plt.subplot(3,3,3)
        plt.plot(t, acc_fil, lw=lw, c='k')
        plt.ylim([-acen, acen])
        plt.title('Filtered')
        
        plt.subplot(3,3,6)
        plt.plot(t, vel_fil, lw=lw, c='k')
        plt.ylim([-vcen, vcen])
        
        plt.subplot(3,3,9)
        plt.plot(t, dis_fil, lw=lw, c='k')
        plt.ylim([-dcen, dcen])
        
        plt.suptitle('Station %s - Channel %s' %(station_id, channel_id))
    
    def plotCombination(event, station_id, channel_id):
        station = event.get(station_id)
        acc = station['acc_uncorrected_' + channel_id]
        acc_corr = station['acc_corrected_' + channel_id]
        acc_fil = station['acc_filtered_' + channel_id]
        p_wave = station['p_wave']['pos']
        dt = station['dt']
        plots(acc, acc_corr, acc_fil, p_wave, dt, station_id, channel_id)
        
    plot = True
    
    with np.load(os.path.join('seismicDatabase', 'npz', eventId + '.npz'), allow_pickle=True) as f:
        event = {}
        for key, value in f.items():
            event[key] = value.item() 
            if not key.startswith('st'):
                continue
        
            station = f.get(key).item()
            
            for i in range(3):
                channel = '%i' %(i+1)
                acc = station['acc_uncorrected_' + channel]
                acc_corr = station['acc_corrected_' + channel]
                acc_fil = station['acc_filtered_' + channel]
                p_wave = station['p_wave']['pos']
                dt = station['dt']
                
                if plot and len(acc_corr):
                    plots(acc, acc_corr, acc_fil, p_wave, dt, key, channel)
                    plt.pause(0.5)
                    response = input('Presione enter para continuar, escriba C para finalizar: ')
                    plt.close('all')
                    
                    if response == 'C':
                        plot = False
else:
    sys.exit('Expected 1 argument, got 0')
