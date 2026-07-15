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

    for i in prange(nTheta, nogil=True):
        theta = pi*i/(nTheta - 1)
        acc   = <double *>malloc(n * sizeof(double))
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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)    
cdef int _Spectrum_Combined(double *ax, double *ay, double s, double c, double dt, double[::1] T, double[::1] xi, int npts, int nT, int nXi, double *Sa, int pos) nogil:

    cdef int i, j, t, first, idx, idx_data
    cdef double beta, gamma, wn
    cdef double b1, b2, c1, c2
    cdef double u2, up2, upp2
    cdef double acc, u2_p
    cdef double *u        = <double *>malloc(nT * nXi * sizeof(double))
    cdef double *up       = <double *>malloc(nT * nXi * sizeof(double))
    cdef double *upp      = <double *>malloc(nT * nXi * sizeof(double))
    cdef double *max_disp = <double *>malloc(nT * nXi * sizeof(double))
    cdef double *a1 = <double *>malloc(nT * nXi * sizeof(double))
    cdef double *a2 = <double *>malloc(nT * nXi * sizeof(double))
    cdef double *a3 = <double *>malloc(nT * nXi * sizeof(double))
    cdef double *a4 = <double *>malloc(nT * nXi * sizeof(double))

    gamma = 0.5
    beta  = 0.25

    b1 = (1. - gamma) * dt
    b2 = gamma * dt

    c1 = (0.5 - beta) * (dt**2)
    c2 = beta * (dt**2)

    if T[0] == 0:
        Sa[pos] = 0.
        for t in range(npts):
            acc = ax[t]*c + ay[t]*s
            if fabs(acc) > Sa[pos]:
                Sa[pos] = fabs(acc)
        for j in range(nXi):
            idx = pos + j
            Sa[idx] = Sa[pos]
        
        first = 1
    else:
        first = 0

    acc = ax[0]*c + ay[0]*s
    for i in range(first, nT):
        wn = 2. * pi / T[i]
        for j in range(nXi):
            idx = i * nXi + j

            u[idx]        = 0.0
            up[idx]       = 0.0
            upp[idx]      = -acc
            max_disp[idx] = 0.0

            a1[idx] = 1. + 2.*xi[j]*wn*b2 + (wn**2)*c2
            a2[idx] = 2.*xi[j]*wn*b1 + (wn**2)*c1
            a3[idx] = 2.*xi[j]*wn + (wn**2)*dt
            a4[idx] = wn**2

    for t in range(npts-1):
        acc = ax[t+1]*c + ay[t+1]*s
        for i in range(first, nT):
            for j in range(nXi):
                idx = i * nXi + j

                upp2 = (-acc - upp[idx]*a2[idx] - up[idx]*a3[idx] - u[idx]*a4[idx])/a1[idx]
                up2  = up[idx] + b1*upp[idx] + b2*upp2
                u2   = u[idx] + up[idx]*dt + c1*upp[idx] + c2*upp2

                u[idx]   = u2
                up[idx]  = up2
                upp[idx] = upp2

                u2_p = fabs(u2)
                if u2_p > max_disp[idx]:
                    max_disp[idx] = u2_p

    for i in range(first, nT):
        for j in range(nXi):
            idx = i * nXi + j
            Sa[pos+idx] = max_disp[idx]*a4[idx]

    free(u)
    free(up)
    free(upp)
    free(max_disp)
    free(a1)
    free(a2)
    free(a3)
    free(a4)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)  
cpdef SpectraRotFull(double[::1] ax, double[::1] ay, double dt, double[::1] T, double[::1] xi, int nTheta):

    cdef double theta, s, c
    cdef int j, n, nT, nXi, offset, x, t, p
    cdef double *acc
    cdef double *thisSa
    cdef np.ndarray Sa

    n   = min(len(ax), len(ay))
    nT  = len(T)
    nXi = len(xi)

    thisSa = <double *>malloc(nXi * nT * nTheta * sizeof(double))

    for j in prange(nTheta, nogil=True):
        theta = pi*j/(nTheta - 1)
        s     = sin(theta)
        c     = cos(theta)

        offset = j*nXi*nT
        _Spectrum_Combined(&ax[0], &ay[0], s, c, dt, T, xi, n, nT, nXi, thisSa, offset)

    Sa = np.empty((nXi, nTheta, nT))
    for x in range(nXi):
        for t in range(nTheta):
            for p in range(nT):
                Sa[x, t, p] = thisSa[(t * nXi * nT) + (p * nXi) + x]
    free(thisSa)

    return Sa

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)  
cpdef SpectraMultiXi(double[::1] ax, double dt, double[::1] T, double[::1] xi):

    cdef int x, p, n, nT, nXi
    cdef double *thisSa
    cdef np.ndarray Sa

    n   = len(ax)
    nT  = len(T)
    nXi = len(xi)

    thisSa = <double *>malloc(nXi * nT * sizeof(double))
    
    _Spectrum_Combined(&ax[0], &ax[0], 1., 0., dt, T, xi, n, nT, nXi, thisSa, 0)

    Sa = np.empty((nXi, nT))
    for x in range(nXi):
        for p in range(nT):
            Sa[x, p] = thisSa[(p * nXi) + x]
    free(thisSa)

    return Sa
