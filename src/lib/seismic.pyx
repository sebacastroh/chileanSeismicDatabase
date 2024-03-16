# -*- coding: utf-8 -*-
"""
Created on Wed May 16 08:26:01 2018

@author: srcastro
"""
import numpy as np
cimport numpy as np
import cython
from cython.parallel import prange
from libc.stdlib cimport malloc, free

cdef extern from "math.h":
    double pow(double, double)
    double fabs(double)  nogil
    double M_PI
    double pi "M_PI"
    double sin(double) nogil
    double cos(double) nogil

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef Spectrum(double[::1] ax, double dt, double[::1] T, double xi):
    '''
    Calculates the acceleration spectrum using the Newmark's method
    '''
    cdef int i, j, npts, nT, first
    cdef double beta, gamma, wn, max_disp
    cdef double a1, a2, a3, a4, b1, b2, c1, c2
    cdef double u1, up1, upp1, u2, up2, upp2
    cdef np.ndarray Sa
    
    npts = len(ax)
    nT   = len(T)
    
    Sa = np.zeros(nT)
    # Average constant acceleration scheme for Newmark's method
    gamma = 0.5
    beta  = 0.25

    if T[0] == 0:
        for j in range(npts):
            if fabs(ax[j]) > Sa[0]:
                Sa[0] = fabs(ax[j])
        first = 0
    else:
        first = -1

    b1 = (1. - gamma)*dt
    b2 = gamma*dt

    c1 = (0.5 - beta)*dt**2
    c2 = beta*dt**2

    for i in range(first, nT-1):
        wn = 2.*pi/T[i+1]

        a1 = 1. + 2.*xi*wn*b2 + wn**2*c2
        a2 = 2.*xi*wn*b1 + wn**2*c1
        a3 = 2.*xi*wn + wn**2*dt
        a4 = wn**2

        u1   = 0.
        up1  = 0.
        upp1 = -ax[0]
        
        max_disp = 0.

        for j in range(npts-1):
            upp2 = (-ax[j+1] - upp1*a2 - up1*a3 - u1*a4)/a1
            up2  = up1 + b1*upp1 + b2*upp2
            u2   = u1 + up1*dt + c1*upp1 + c2*upp2

            u1   = u2
            up1  = up2
            upp1 = upp2

            u2_p = fabs(u2)

            if u2_p > max_disp:
                max_disp = u2_p
                Sa[i+1]  = u2_p*a4

    return Sa

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)    
cdef int _Spectrum(double *ax, double dt, double[::1] T, double xi, int npts, int nT, double *Sa, int pos) nogil:
    '''
    Calculates the acceleration spectrum using the Newmark's method
    '''
    cdef int i, j, first
    cdef double beta, gamma, wn, max_disp
    cdef double a1, a2, a3, a4, b1, b2, c1, c2
    cdef double u1, up1, upp1, u2, up2, upp2

    # Average constant acceleration scheme for Newmark's method
    gamma = 0.5
    beta  = 0.25
    
    if T[0] == 0:
        Sa[pos] = 0.
        for j in xrange(npts):
            if fabs(ax[j]) > Sa[pos]:
                Sa[pos] = fabs(ax[j])
        first = 0
    else:
        first = -1

    b1 = (1. - gamma)*dt
    b2 = gamma*dt

    c1 = (0.5 - beta)*dt**2
    c2 = beta*dt**2

    for i in xrange(first, nT-1):
        wn = 2.*pi/T[i+1]

        a1 = 1. + 2.*xi*wn*b2 + wn**2*c2
        a2 = 2.*xi*wn*b1 + wn**2*c1
        a3 = 2.*xi*wn + wn**2*dt
        a4 = wn**2

        u1   = 0.
        up1  = 0.
        upp1 = -ax[0]
        
        max_disp = 0.

        for j in xrange(npts-1):
            upp2 = (-ax[j+1] - upp1*a2 - up1*a3 - u1*a4)/a1
            up2  = up1 + b1*upp1 + b2*upp2
            u2   = u1 + up1*dt + c1*upp1 + c2*upp2

            u1   = u2
            up1  = up2
            upp1 = upp2

            u2_p = fabs(u2)

            if u2_p > max_disp:
                max_disp    = u2_p
                Sa[pos+i+1] = u2_p*a4

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)  
cpdef SpectraRot(double[::1] ax, double[::1] ay, double dt, double[::1] T, double xi, int nTheta):
    
    cdef double theta, s, c
    cdef int i, j, n, nT
    cdef double *acc
    cdef double *thisSa
    cdef np.ndarray Sa
    
    n  = min(len(ax), len(ay))
    nT = len(T)
    
    thisSa = <double *>malloc(nT * nTheta * sizeof(double))
    acc    = <double *>malloc(n * sizeof(double))

    for i in prange(nTheta, nogil=True):
        theta = pi*i/(nTheta - 1)
        s     = sin(theta)
        c     = cos(theta)
        for j in range(n):
            acc[j] = ax[j]*c + ay[j]*s
        
        _Spectrum(acc, dt, T, xi, n, nT, thisSa, i*nT)
        
    free(acc)
        
    Sa = np.empty((nTheta, nT))
    for i in range(nTheta):
        for j in range(nT):
            Sa[i,j] = thisSa[i*nT+j]
    free(thisSa)
    
    return Sa