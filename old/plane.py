#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 12:42:02 2019

@author: srcastro
"""
import math
import numpy as np

#%% Functions
def LatLonDepth2XYZNum(lat, lon, depth):
    a = 6378.137 # km
    f = 1./298.257223563
    b = a*(1. - f)
    d = math.sqrt(a**2*math.cos(lat*math.pi/180.)**2 + b**2*math.sin(lat*math.pi/180.)**2)
    
    xo = (a**2/d + depth)*math.cos(lat*math.pi/180.)*math.cos(lon*math.pi/180.)
    yo = (a**2/d + depth)*math.cos(lat*math.pi/180.)*math.sin(lon*math.pi/180.)
    zo = (b**2/d + depth)*math.sin(lat*math.pi/180.)
    
    xyz = np.array([xo, yo, zo])
    
    return xyz

def LatLonDepth2XYZNumPy(lat, lon, depth):
    a = 6378.137 # km
    f = 1./298.257223563
    b = a*(1. - f)
    d = np.sqrt(a**2*np.cos(lat*np.pi/180.)**2 + b**2*np.sin(lat*np.pi/180.)**2)
    
    xo = (a**2/d + depth)*np.cos(lat*np.pi/180.)*np.cos(lon*np.pi/180.)
    yo = (a**2/d + depth)*np.cos(lat*np.pi/180.)*np.sin(lon*np.pi/180.)
    zo = (b**2/d + depth)*np.sin(lat*np.pi/180.)
    
    xyz = np.c_[xo, yo, zo]
    
    return xyz

Mw = 7.1
L = 10**(-2.9 + 0.63*Mw)
W = 10**(-0.86 + 0.35*Mw)
strike = 20.*np.pi/180.
dip = 16.*np.pi/180.
rake = 112.*np.pi/180.

n = np.array([-np.sin(dip)*np.sin(strike), np.sin(dip)*np.cos(strike), -np.cos(dip)])

depth = 33.8
lat = -35.2
lon = -72.22

xyz = LatLonDepth2XYZNum(lat, lon, depth)