#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:35:19 2019

@author: srcastro
"""
import numpy as np
from scipy import integrate
from scipy import signal
from scipy import stats
from obspy.signal import filter as flt
import rjmcmc

def PWaveDetection(X, fsamp, WN=800, AIlim=0.5, seconds=2.,
                   ntrials=5, return_similarity=False):
    """
    Detects the arrival of the P wave of a given waveform analyzing a
    similarity coefficient between the waveform and a white noise.
    
    Parameters
    ----------
    X: numpy_array
        Waveform
    
    fsamp: int
        Sampling frequency of the waveform
    
    WN: (optional) int
        Size of the window. By default is 800
    
    AIlim: (optional) float
        Arias Intensity threshold. By default is 0.5
    
    seconds: (optional) float
        Seconds before the mean limit to start the search of the P wave.
        By default is 2.
    
    ntrials: (optional) int
        Number of iterations. By default is 5
    
    return_similarity: (optional) boolean
        If true, returns similarity coefficients. By default is False
        
    Returns
    -------
    p_wave_pos: int
        Index of the arrival of the P Wave in X
        
    similarity: numpy_array
        The similarity coefficients of the waveform. Only returned if
        return_similarity is set True
    """
    
    if X.mean() != 0.:
        Z = (X/X.mean())
    else:
        Z = X.copy()
    
    arias = integrate.cumtrapz(Z**2, dx=1./fsamp, initial=0.)
    arias /= arias.max()
    
    energy_pos = np.where(arias <= AIlim)[0][-1]
    
    if len(Z) <= WN:
        WN = int(len(Z)/3)
    
    if energy_pos <= WN:
        energy_pos = len(arias)
    
    ycorr = np.ones(len(arias))
    ycorr[WN:energy_pos] = 0.
    
    indexes = np.array([np.arange(x, x+WN) for x in range(energy_pos-WN)])
    for k in range(ntrials):
        Y = np.random.rand(len(arias))
        ycorr[WN:energy_pos] += 1./np.sqrt(np.sum((Z[indexes] - Y[indexes])**2, axis=1)/(WN-1))
    
    ycorr[WN:energy_pos] = ycorr[WN:energy_pos]/ntrials
    ycorr[WN:energy_pos] = ycorr[WN:energy_pos]/ycorr[WN:energy_pos].max()
    xcorr = np.arange(len(arias))
    
    ysmooth = np.convolve(ycorr, np.ones(int(fsamp/2.)+1)/(int(fsamp/2.)+1), 'same')
    
    mean = ysmooth[WN:energy_pos].mean()
    
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
#    peaks3 = np.hstack([np.arange(x-int(fsamp), x+-int(fsamp)) for x in peaks2])
    
    p5 = np.hstack((peaks, peaks2, peaks3))
    p5 = np.unique(p5)
    
    if len(p5) == 0:
        p5 = np.array([p2])
    else:
        p5 += p2
        
    p4 = np.argmin(np.divide(arias[p5], arias[p5+int(fsamp)], out=np.arange(len(p5))[::-1]+1., where=arias[p5+int(fsamp)]!=0.))
    
    pos = p5[p4]
    
    if return_similarity:
        similarity = np.ones_like(X)
        similarity[WN:energy_pos] = ysmooth[WN:energy_pos]
        similarity = np.vstack((np.arange(len(X)), similarity)).T        
        return pos, similarity
    else:
        return pos

sample_x = None
sample_curves = []
sample_i = 0
sample_rate = 1

def sampler_cb(x, y):
    global sample_x, sample_curves, sample_i, sample_rate

    if sample_i == 0:
        sample_x = x

    if sample_i % sample_rate == 0:
        sample_curves.append(y)

    sample_i = sample_i + 1

def RecordCorrection(record, dt, pos, alpha=0.005):
    global sample_x, sample_curves, sample_i, sample_rate
    nn = len(record) - pos
    window = signal.tukey(nn, alpha=alpha)
    
    new_tr = record.copy()
    new_tr[pos:] = new_tr[int(pos):]*window
    new_tr[:pos] = 0
    
    # Step 3: Band-pass filtering on velocity
    vel = integrate.cumtrapz(new_tr, dx=dt)
    n = len(vel)
    
    t = np.linspace(0, (n-1)*dt, n)
    min_std = 0.1*np.max(np.abs(vel))
    
    # Computing std from envelope
    #Filter needs to be adjusted but this is a good start...
    fil_vel = flt.bandpass(vel, .01, 99., 1./dt, corners=3, zerophase=True)
    # Step 4: Baseline correction
    # Envelope of filtered data
    tr_env= flt.envelope(fil_vel)
    
    std_n = tr_env.max()*np.ones(len(fil_vel))
    std_n[np.where(std_n < min_std)] = min_std
    
    data = rjmcmc.dataset1d(*map(list,[t, fil_vel, std_n]))
    
    # This is our callback function which samples the curves generated 
    # during the analysis
    #
    sample_x = None
    sample_curves = []
    sample_i = 0
    sample_rate = 1
    
    pd = 0.1 * t.max()
    burnin = 10000
    total = 50000
    max_partitions = 20
    max_order = 1
 
    results = rjmcmc.regression_part1d_sampled(data, sampler_cb, 
                                       pd, 
                                       burnin, 
                                       total, 
                                       max_partitions,
                                       max_order)
    

    #################################  SEGUNDA PARTE #####################################
    mode = stats.mode(results.partitions())[0][0]
    
    pos1 = np.where(results.partitions()==mode)[0]
    
    m_misfit = []
    for i in range(len(pos1)):
        m_misfit.append(results.misfit()[pos1[i]])
    
    best = np.min(m_misfit)
    pos_b = np.where(m_misfit==best)[0][0]
    
    m_curves = []
    for i in range(len(pos1)):
        m_curves.append(sample_curves[pos1[i]])
    
    best_curv = m_curves[pos_b]
    
    correc = np.interp(t, sample_x, best_curv)
    vel1 = fil_vel - correc
    
    # Step 5: Derivation of velocity
    window2 = signal.tukey(n, alpha=alpha)
    velf = vel1*window2
    acc = np.gradient(velf, dt, edge_order=2)
    
    return acc