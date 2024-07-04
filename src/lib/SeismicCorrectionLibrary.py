#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:35:19 2019

@author: srcastro
"""
import numpy as np
from scipy import integrate
from scipy import signal

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
    
    arias = integrate.cumulative_trapezoid(Z**2, dx=1./fsamp, initial=0.)
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
