# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 23:53:06 2022

@author: sebac
"""
import os
import numpy as np
import scipy.integrate as spin
import matplotlib.pyplot as plt

eventId = '20150916_8.4M_31.55S_71.86W_11.0KM'

plt.close('all')

st = 'st00'
channel = 'acc_3'

with np.load(os.path.join('databaseUncorrected', 'npz', eventId + '.npz'), allow_pickle=True) as f:
    station = f.get(st).item()
    acc = station[channel]
    p_wave = station['p_wave']['pos']
    dt = station['dt']

with np.load(os.path.join('databaseCorrected', 'npz', eventId + '.npz'), allow_pickle=True) as f:
    station = f.get(st).item()
    acc_corr = station[channel]
    
with np.load(os.path.join('databaseFiltered', 'npz', eventId + '.npz'), allow_pickle=True) as f:
    station = f.get(st).item()
    acc_fil = station[channel]


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

plt.figure()

plt.subplot(3,1,1)
plt.plot(t, acc, lw=lw, c='k')
plt.ylim([-acen, acen])

plt.subplot(3,1,2)
plt.plot(t, vel, lw=lw, c='k')
plt.ylim([-vcen, vcen])

plt.subplot(3,1,3)
plt.plot(t, dis, lw=lw, c='k')
plt.ylim([-dcen, dcen])

plt.suptitle('Uncorrected')

plt.savefig('uncorrected.pdf', bbox_inches='tight')

plt.figure()

plt.subplot(3,1,1)
plt.plot(t, acc_corr, lw=lw, c='k')
plt.ylim([-acen, acen])

plt.subplot(3,1,2)
plt.plot(t, vel_corr, lw=lw, c='k')
plt.ylim([-vcen, vcen])

plt.subplot(3,1,3)
plt.plot(t, dis_corr, lw=lw, c='k')
plt.ylim([-dcen, dcen])

plt.suptitle('Corrected')

plt.savefig('corrected.pdf', bbox_inches='tight')

plt.figure()

plt.subplot(3,1,1)
plt.plot(t, acc_fil, lw=lw, c='k')
plt.ylim([-acen, acen])

plt.subplot(3,1,2)
plt.plot(t, vel_fil, lw=lw, c='k')
plt.ylim([-vcen, vcen])

plt.subplot(3,1,3)
plt.plot(t, dis_fil, lw=lw, c='k')
plt.ylim([-dcen, dcen])

plt.suptitle('Filtered')

plt.savefig('filtered.pdf', bbox_inches='tight')
