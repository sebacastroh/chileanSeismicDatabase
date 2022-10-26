# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 17:13:41 2022

@author: sebac
"""
import pyrjmcmc
import numpy as np

import matplotlib.pyplot as plt
import scipy.io as spio
import scipy.stats as spsts
import scipy.integrate as spin
import scipy.signal as spsig
import obspy.signal.konnoohmachismoothing as kon
from obspy.signal import filter as flt


plt.close('all')

def sampling_function(x, boundaries, local_values):
    
    y = np.zeros_like(x)
    for i in range(len(boundaries)-1):
        pos = np.where((x >= boundaries[i]) & (x < boundaries[i+1]))
        y[pos] = np.polyval(local_values[i], x[pos])
    y[-1] = np.polyval(local_values[-1], x[-1])
    
    return y

data = spio.loadmat('20190120_6.7M_30.28S_71.36W_50KM.mat', struct_as_record=False, squeeze_me=True)

station = data['st01']
acc = station.acc_3
p_wave = 4863
dt = station.dt


classic = True
fsamp = 1./dt

n = len(acc)
t = np.linspace(0., (n-1)*dt, n)

# Zero completation and band-pass filtering
nn = len(acc) - p_wave

new_acc = acc - acc[:p_wave].mean()


new_vel = spin.cumtrapz(new_acc, dx=dt, initial=0.)
new_vel -= new_vel[:p_wave].mean()


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

y = [p_wave]
y.extend([np.argmin(np.abs(energy-x/10.)) for x in range(1,11)])
y = np.array(y)  
min_partitions = 1  
max_partitions = 10#int(max(10, np.sum(np.int64((t[y[1:]]-t[y[:-1]])/5.))))

burnin = 1000
total = 100000

# Pre

datasets = [np.column_stack((t[:p_wave], vel_mean[:p_wave], smooth_std[:p_wave]))]
samples = 100
pd = (t[p_wave])*0.1

results = multiple_partition_connected_lines(datasets, burnin, total, min_partitions, max_partitions, samples, pd, seed=None)

mode = spsts.mode(results.npartitions, keepdims=True)
pos = np.where(results.npartitions == mode.mode)
sol = pos[0][results.misfit[pos].argmin()]

y_pre = sampling_function(t[:p_wave], results.boundaries[sol], results.local_parameters[sol])

# fig = plt.figure(figsize=(15, 8))

# a1 = fig.add_subplot(2,2,1)

# for di,dataset in enumerate(datasets):
#     a1.plot(dataset[:,0], dataset[:,1], label='Dataset %i' %(di+1))
#     a1.plot(results.x, results.y[sol], label='Selected model %i' %(di+1))
#     a1.plot(results.x, results.y[sol], 'o', label='')

# a1.legend()
# a1.set_xlabel('x-axis')
# a1.set_ylabel('y-axis')

# a2 = fig.add_subplot(2,2,3)
# a2.bar(results.boundaries_histogram[0], results.boundaries_histogram[1], ec='k', width=t[p_wave]/samples)
# a2.set_xlabel('x-axis')
# a2.set_ylabel('Boundary frequency')

# a3 = fig.add_subplot(2,2,2)
# a3.bar(results.partitions_histogram[0], results.partitions_histogram[1], ec='k')
# a3.set_xlabel('Number of partitions')
# a3.set_ylabel('Partition frequency')

# a4 = fig.add_subplot(2,2,4)
# a4.semilogy(results.misfit)
# a4.set_xlabel('Iteration')
# a4.set_ylabel('Loglikelihood')

# fig.show()

# Post
datasets = [np.column_stack((t[p_wave:], vel_mean[p_wave:], smooth_std[p_wave:]))]
samples = 100
pd = (t[-1]-t[p_wave])*0.1
results = pyrjmcmc.multiple_partition_connected_lines(datasets, burnin, total, min_partitions, max_partitions, samples, pd, seed=None)

mode = spsts.mode(results.npartitions, keepdims=True)
pos = np.where(results.npartitions == mode.mode)
sol = pos[0][results.misfit[pos].argmin()]

y_post = sampling_function(t[p_wave:], results.boundaries[sol], results.local_parameters[sol])

# fig = plt.figure(figsize=(15, 8))

# a1 = fig.add_subplot(2,2,1)

# for di,dataset in enumerate(datasets):
#     a1.plot(dataset[:,0], dataset[:,1], label='Dataset %i' %(di+1))
#     a1.plot(results.x, results.y[sol], label='Selected model %i' %(di+1))
#     a1.plot(results.x, results.y[sol], 'o', label='')

# a1.legend()
# a1.set_xlabel('x-axis')
# a1.set_ylabel('y-axis')

# a2 = fig.add_subplot(2,2,3)
# a2.bar(results.boundaries_histogram[0], results.boundaries_histogram[1], ec='k', width=(t[-1]-t[p_wave])/samples)
# a2.set_xlabel('x-axis')
# a2.set_ylabel('Boundary frequency')

# a3 = fig.add_subplot(2,2,2)
# a3.bar(results.partitions_histogram[0], results.partitions_histogram[1], ec='k')
# a3.set_xlabel('Number of partitions')
# a3.set_ylabel('Partition frequency')

# a4 = fig.add_subplot(2,2,4)
# a4.semilogy(results.misfit)
# a4.set_xlabel('Iteration')
# a4.set_ylabel('Loglikelihood')

# fig.show()

m1 = (y_pre[-2] - y_pre[-1])/dt
i = 2
while True:
    yp = m1*(t[p_wave-i] - t[p_wave]) + y_pre[-1]
    if yp != y_pre[-i-1]:
        break
    i += 1
    if p_wave - i <= 0:
        i = p_wave
        break

m2 = (y_post[0] - y_post[1])/dt
j = 2
while True:
    yp = m2*(t[p_wave+j] - t[p_wave]) + y_post[0]
    if yp != y_post[j]:
        break
    j += 1

    if p_wave + j >= n:
        j = n - p_wave
        break

solution = np.hstack((y_pre[:-i],
                      np.interp(t[p_wave-i:p_wave+j], [t[p_wave-i], t[p_wave+j]], [y_pre[-i], y_post[j]]),
                      y_post[j:]))

vel_corr = new_vel - solution
acc_corr = np.gradient(vel_corr, dt, edge_order=2)

acc_corr = acc_corr - acc_corr[:p_wave].mean()
signal = acc_corr[p_wave:]
noise = acc_corr[:p_wave]
# signal = acc[p_wave:]
# noise = acc[:p_wave]

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
# peaks = spsig.find_peaks(SNR)[0]
# freq_min = max(freq[np.argmin(SNR[peaks[0]:peaks[-1]+1])+peaks[0]], freqN[1])

freq = np.logspace(0, 2, 1001)
Sint = np.interp(freq, freqS[posS], Ssmooth)
Nint = np.interp(freq, freqN[posN], Nsmooth)

SNR = Sint/Nint
pos = np.where(SNR < 3.)[0]
if len(pos) > 0:
    freq_max = min(max(freq[pos[0]], 30.), fsamp/2.)
else:
    freq_max = fsamp/2.

order = 4
zerophase = True

if freq_max >= fsamp/2.:
    acc_fil = flt.highpass(acc_corr, freq_min, fsamp, corners=order, zerophase=zerophase)
else:
    acc_fil = flt.bandpass(acc_corr, freq_min, freq_max, fsamp, corners=order, zerophase=zerophase)
   
vel_fil = spin.cumtrapz(acc_fil, dx=dt, initial=0.)
vel_fil -= vel_fil[:p_wave].mean()
dis_fil = spin.cumtrapz(vel_fil, dx=dt, initial=0.)
dis_fil -= dis_fil[:p_wave].mean()
vel = spin.cumtrapz(acc, dx=dt, initial=0.)
vel -= vel[:p_wave].mean()
dis = spin.cumtrapz(vel, dx=dt, initial=0.)
dis -= dis[:p_wave].mean()

lw = 0.25

plt.figure()
plt.subplot(3,1,1)
plt.plot(t, acc, lw=lw, c='k')
plt.subplot(3,1,2)
plt.plot(t, vel, lw=lw, c='k')
plt.subplot(3,1,3)
plt.plot(t, dis, lw=lw, c='k')
plt.ylim([dis.min(), dis.max()])

plt.figure()
plt.subplot(3,1,1)
plt.plot(t, acc_fil, lw=lw, c='k')
plt.subplot(3,1,2)
plt.plot(t, vel_fil, lw=lw, c='k')
plt.subplot(3,1,3)
plt.plot(t, dis_fil, lw=lw, c='k')
plt.ylim([dis.min(), dis.max()])

