#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 17:54:03 2019

@author: srcastro
"""
import os
import numpy as np
import unidecode
import chardet
import scipy.io
import geopy.distance
from zipfile import ZipFile
import pandas as pd
import pickle

from stations_renadic import stations_renadic
from stations_csn import stations_csn

import distanceSurfaceFault

licensing = 'This SIBER-RISK Strong Motion Database is made available '
licensing += 'under the Creative Commons Attribution-NonCommercial-'
licensing += 'ShareAlike 4.0 International Public License: '
licensing += 'https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode. '
licensing += 'Any rights in individual contents of the database are '
licensing += 'licensed under the Creative Commons Attribution-'
licensing += 'NonCommercial-ShareAlike 4.0 International Public License: '
licensing += 'https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode'

slab = None#np.load('sam_slab2.npz')

def write_dataframe(stations):
    table = []
    for station in stations.values():
        table.append([station.event_name, station.start_time, station.magnitude,
                  station.hypocenter_lat, station.hypocenter_lon,
                  station.event_type, station.depth, station.name,
                  station.lat, station.lon, station.Rrup, station.Vs30,
                  station.azimut, np.max(np.abs(np.hstack((station.acc_1, station.acc_2, station.acc_3))))])
    
    df = pd.DataFrame(np.array(table), columns=['Earthquake Name',
                      'Start time record', 'Magnitude [Mw]',
                      'Hypocenter latitude', 'Hypocenter longitude',
                      'Event type', 'Depth [km]', 'Station name',
                      'Station latitude', 'Station longitude',
                      'Distance of rupture [km]', 'Vs30 [m/s]',
                      'Azimut [o]', 'PGA (from 3 components) [g]'])
    
    return df

def write_dataframe2(stations):
    table = []
    for station in stations.values():
        table.append([station['event_name'], station['start_time'], station['magnitude'],
                  station['hypocenter_lat'], station['hypocenter_lon'],
                  station['event_type'], station['depth'], station['name'],
                  station['lat'], station['lon'], station['Rrup'], station['Vs30'],
                  station['azimut'], np.max(np.abs(np.c_[station['acc_1'],
                                                   station['acc_2'],
                                                   station['acc_3']]))])
    
    df = pd.DataFrame(np.array(table), columns=['Earthquake Name',
                      'Start time record', 'Magnitude [Mw]',
                      'Hypocenter latitude', 'Hypocenter longitude',
                      'Event type', 'Depth [km]', 'Station name',
                      'Station latitude', 'Station longitude',
                      'Distance of rupture [km]', 'Vs30 [m/s]',
                      'Azimut [o]', 'PGA (from 3 components) [g]'])
    
    return df

def transformRenadicv1(GMPaths, event_path, save_path, replace_old=False):
    """
    Stores information of several ground motion files with the "Renadic" format 
    in a pickle file.
    
    A new .pickle file is stored in the same path where the original file is
    located.

    Parameters
    ----------
    GMPaths : list of strings
        Paths of the files containing each ground motions to transform.
    fileType : string
        Type of output file. Options: pickle or mat.
    Variables stored
    ----------
    - dt : time step
    - chx : acceleration of the first horizontal component (cm/s/s)
    - chy : acceleration of the second horizontal component (cm/s/s)
    - chz : acceleration of the vertical component (cm/s/s)
    """
    global lat_slab, lon_slab, dep_slab
    
    if isinstance(GMPaths, str):
        GMPaths = [GMPaths]
    
    # Get information of the event
    event = GMPaths[0].split('/')[-2]
    info = event.split('_')

    matfile = os.path.join(save_path, event + '.mat')

    if os.path.exists(matfile) and not replace_old:
        stations = scipy.io.loadmat(matfile, squeeze_me=True, struct_as_record=False)
        stations.pop('__version__')
        stations.pop('__globals__')
        stations.pop('__header__')
        return write_dataframe(stations), False
    
    magnitude = float(info[1][:-2])
    
    try:
        hypocenter_lat = -float(info[2][:-1])
        hypocenter_lon = -float(info[3][:-1])
        depth = float(info[4][:-2])
        details = True
    except:
        details = False
    
    n_registros = len(GMPaths)
    stations = {}
    
    for i in np.arange(n_registros):
        station = {}
        # Read the text file
        GM_txt = open(GMPaths[i],'r')
        GM_lines = GM_txt.readlines()
        GM_txt.close()
        
        # Get station name
        encoding = chardet.detect(GM_lines[5].strip())['encoding']
        station_name = unidecode.unidecode(GM_lines[5].strip().decode(encoding))
        pos = station_name.find('S/N')
        station_name = station_name[:pos-1]
        
        if station_name.endswith(' '):
            station_name = station_name[:-1]
            
        if station_name.endswith('(SMA-1)'):
            station_name = u'ANTOFAGASTA UCN'
        elif station_name == u'IQUIQUE E. CHIPANA' or station_name == u'IQUQUE ESCUELA CHIPANA':
            station_name = u'IQUIQUE ESCUELA CHIPANA'
        elif station_name == u'STGO CENTRO':
            station_name = u'SANTIAGO CENTRO'
        elif station_name == u'VINA DEL MAR - EL SALTO' or station_name == u'VINA EL SALTO':
            station_name = u'VINA DEL MAR EL SALTO'
        elif station_name == u'VINA EL SALTO CERRO':
            station_name = u'VINA DEL MAR EL SALTO CERRO'
        elif station_name == u'ETNA':
            station_name = u'PICA'
        elif station_name == 'QDR':
            if event[:4] == '2009':
                station_name = u'ALTO HOSPICIO'
            else:
                station_name = u'VALLENAR'
        
        num_points = int(GM_lines[10].split()[4])
        dt = float(GM_lines[10].split()[-2].replace('=',''))/num_points
        
#         Extract accelerations
        lines_to_read = (num_points*2-1)/10 + 1
        chx = np.array([])
        chy = np.array([])
        chz = np.array([])
        baseCh1 = 27
        baseCh2 = 27*2 + lines_to_read + 1
        baseCh3 = 27*3 + lines_to_read*2 + 2
        if GM_lines[6].split()[2] == 'V' or GM_lines[6].split()[2] == 'Z':
            component1 = GM_lines[baseCh2-27+6].split()[2]
            component2 = GM_lines[baseCh3-27+6].split()[2]
            component3 = GM_lines[6].split()[2]
            for j in np.arange(lines_to_read):
                try:
                    chx = np.hstack([chx,np.array(GM_lines[baseCh2+j].split()[1::2]).astype(np.float)])
                    chy = np.hstack([chy,np.array(GM_lines[baseCh3+j].split()[1::2]).astype(np.float)])
                    chz = np.hstack([chz,np.array(GM_lines[baseCh1+j].split()[1::2]).astype(np.float)])
                except:
                    thischx = []
                    thischy = []
                    thischz = []
                    stringsX = GM_lines[baseCh2+j].split()
                    stringsY = GM_lines[baseCh3+j].split()
                    stringsZ = GM_lines[baseCh1+j].split()
                    for k in np.arange(1, len(stringsX)):
                        if len(stringsX[k])<8:
                            thischx.append(float(stringsX[k]))
                            thischy.append(float(stringsY[k]))
                            thischz.append(float(stringsZ[k]))
                        else:
                            thischx.append(float(stringsX[k][0:-7]))
                            thischy.append(float(stringsY[k][0:-7]))
                            thischz.append(float(stringsZ[k][0:-7]))
                    chx = np.hstack([chx, thischx])
                    chy = np.hstack([chy, thischy])
                    chz = np.hstack([chz, thischz])
        elif GM_lines[baseCh2-27+6].split()[2] == 'V' or GM_lines[baseCh2-27+6].split()[2] == 'Z':
            component1 = GM_lines[6].split()[2]
            component2 = GM_lines[baseCh3-27+6].split()[2]
            component3 = GM_lines[baseCh2-27+6].split()[2]
            for j in np.arange(lines_to_read):
                try:
                    chx = np.hstack([chx,np.array(GM_lines[baseCh1+j].split()[1::2]).astype(np.float)])
                    chy = np.hstack([chy,np.array(GM_lines[baseCh3+j].split()[1::2]).astype(np.float)])
                    chz = np.hstack([chz,np.array(GM_lines[baseCh2+j].split()[1::2]).astype(np.float)])
                except:
                    if len(GM_lines[baseCh1+j].split())<10:
                        thischx = []
                        thischy = []
                        thischz = []
                        stringsX = GM_lines[baseCh1+j].split()
                        stringsY = GM_lines[baseCh3+j].split()
                        stringsZ = GM_lines[baseCh2+j].split()
                        for k in np.arange(1, len(stringsX)):
                            if len(stringsX[k])<8:
                                thischx.append(float(stringsX[k]))
                                thischy.append(float(stringsY[k]))
                                thischz.append(float(stringsZ[k]))
                            else:
                                thischx.append(float(stringsX[k][0:-7]))
                                thischy.append(float(stringsY[k][0:-7]))
                                thischz.append(float(stringsZ[k][0:-7]))
                        chx = np.hstack([chx, thischx])
                        chy = np.hstack([chy, thischy])
                        chz = np.hstack([chz, thischz])

        elif GM_lines[baseCh3-27+6].split()[2] == 'V' or GM_lines[baseCh3-27+6].split()[2] == 'Z':
            component1 = GM_lines[6].split()[2]
            component2 = GM_lines[baseCh2-27+6].split()[2]
            component3 = GM_lines[baseCh3-27+6].split()[2]
            for j in np.arange(lines_to_read):
                try:
                    chx = np.hstack([chx,np.array(GM_lines[baseCh1+j].split()[1::2]).astype(np.float)])
                    chy = np.hstack([chy,np.array(GM_lines[baseCh2+j].split()[1::2]).astype(np.float)])
                    chz = np.hstack([chz,np.array(GM_lines[baseCh3+j].split()[1::2]).astype(np.float)])
                except:
                    thischx = []
                    thischy = []
                    thischz = []
                    stringsX = GM_lines[baseCh1+j].split()
                    stringsY = GM_lines[baseCh2+j].split()
                    stringsZ = GM_lines[baseCh3+j].split()
                    for k in np.arange(1, len(stringsX)):
                        if len(stringsX[k])<8:
                            thischx.append(float(stringsX[k]))
                            thischy.append(float(stringsY[k]))
                            thischz.append(float(stringsZ[k]))
                        else:
                            thischx.append(float(stringsX[k][0:-7]))
                            thischy.append(float(stringsY[k][0:-7]))
                            thischz.append(float(stringsZ[k][0:-7]))
                    chx = np.hstack([chx, thischx])
                    chy = np.hstack([chy, thischy])
                    chz = np.hstack([chz, thischz])

        else:
            raise Exception('The vertical component of groun motion '+str(i+1)+' is not defined by "V" or "Z"')
        
        acc_component1 = chx*9.81/10. # m/s/s
        acc_component2 = chy*9.81/10. # m/s/s
        acc_component3 = chz*9.81/10. # m/s/s
        
        station['name'] = station_name
        station['lat'] = stations_renadic[station_name][0]
        station['lon'] = stations_renadic[station_name][1]
        station['azimut'] = stations_renadic[station_name][4]
        
        if details:
            if magnitude < 6.:
                distance = np.sqrt(geopy.distance.geodesic((hypocenter_lat,
                                                            hypocenter_lon),
                                                            (stations_renadic[station_name][0],
                                                             stations_renadic[station_name][1])).km**2 + depth**2)
                station['Rrup'] = distance
            else:
                distance = 'To determine'
                
                station['Rrup'] = distance
        else:
            station['Rrup'] = 'N/A'
            
        station['Vs30'] = stations_renadic[station_name][2]
        
        station['acc_1'] = acc_component1.copy()
        station['component_1'] = component1
        
        station['acc_2'] = acc_component2.copy()
        station['component_2'] = component2
        
        station['acc_3'] = acc_component3.copy()
        station['component_3'] = component3
        
        start_time = GM_lines[3].strip()[14:].replace('/', '-')
        start_time = start_time.split(' ')[1]
        pos = start_time.find(':')
        if pos == 1:
            start_time = '0' + start_time
        start_time = event[:4] + '-' + event[4:6] + '-' + event[6:8] + ' ' + start_time
        
        if start_time.find('.') == -1:
            start_time += '.000'
        
        station['start_time'] = start_time
        station['dt'] = dt
        
        station['event_name'] = event
        station['magnitude'] = magnitude
        if details:
            station['hypocenter_lat'] = hypocenter_lat
            station['hypocenter_lon'] = hypocenter_lon
            station['depth'] = depth
            
            dep_pos = np.nanargmin(np.sqrt((slab['lon'] - hypocenter_lon)**2 + (slab['lat'] - hypocenter_lat)**2))
            # La profundidad del slab es negativa
            difference = -slab['dep'][dep_pos] - depth
            
            if difference > 20:
                station['event_type'] = 'crustal'
            elif difference < -20:
                station['event_type'] = 'interface'
            else:
                station['event_type'] = 'intraslab'
        else:
            station['hypocenter_lat'] = 'N/A'
            station['hypocenter_lon'] = 'N/A'
            station['depth'] = 'N/A'
            station['event_type'] = 'N/A'
        station['licensing'] = licensing
        station['units'] = 'Acceleration: m/s/s; Magnitude: Mw; dt: s; Depth: km; Vs30: m/s, Rrup: km'
        
        stations['st%0.2i' %i] = station
    
    scipy.io.savemat(matfile, stations)
    
    df = write_dataframe2(stations)
    
    return df, True

def transformCSN(GMPaths, event_path, save_path, replace_old=False):
    if isinstance(GMPaths, str):
        GMPaths = [GMPaths]
    
    # Get information of the event
    event = GMPaths[0].split('/')[-2]
    info = event.split('_')
    
    matfile = os.path.join(save_path, event + '.mat')
    
    if os.path.exists(matfile) and not replace_old:
        stations = scipy.io.loadmat(matfile, squeeze_me=True, struct_as_record=False)
        stations.pop('__version__')
        stations.pop('__globals__')
        stations.pop('__header__')
        return write_dataframe(stations), False
    
    magnitude = float(info[1][:-2])
    
    hypocenter_lat = -float(info[2][:-1])
    hypocenter_lon = -float(info[3][:-1])
    depth = float(info[4][:-2])
    
    stations_names = list(set([filename.split('/')[-1][:-4] for filename in GMPaths]))
    stations = {}
    for i,station_name in enumerate(stations_names):
        ch = []
        times = []
        components = []
        station = {}
        filename = event_path + '/' + station_name + '.zip'
        with ZipFile(filename, 'r') as zipObj:
            zipObj.extractall(event_path)
        filenames = os.listdir(event_path)
        for filename in filenames:
            if filename.endswith('.txt'):
                with open(event_path + '/' + filename, 'r') as fopen:
                    for line in fopen:
                        if line.startswith('#'):
                            if line.startswith('# Tiempo'):
                                this_start_time = line.strip().split()[-1]
                                pos = this_start_time.find('T') + 1
                                time = this_start_time[pos:-1].split(':')
                                time = int(time[0])*3600 + int(time[1])*60 + float(time[2])
                            elif line.startswith('# Tasa'):
                                dt = 1./float(line.strip().split()[-2])
                            elif line.startswith('# Esta'):
                                this_station_name = line.strip().split()[2]
                                component = line.strip().split()[-1]
                            elif line.startswith('# Lat'):
                                lat = float(line.strip().split()[2])
                                lon = float(line.strip().split()[4])
                        else:
                            break
                ch.append(np.loadtxt(event_path + '/' + filename))
                times.append(time)
                components.append(component)
                os.remove(event_path + '/' + filename)
                
        min_time = min(times)
        max_time = max(times)
        
        if magnitude < 6.:
            distance = np.sqrt(geopy.distance.geodesic((hypocenter_lat, hypocenter_lon), (lat, lon)).km**2 + depth**2)
            station['Rrup'] = distance
        else:
            fault_plane_properties = '/home/srcastro/projects/correctSeismicDatabase/newMethodology/fault_plane_properties.pkl'
            with open(fault_plane_properties, 'rb') as f:
                data = pickle.load(f)
            if event in data.keys():
                strike, dip, rake = data[event]
                distance = distanceSurfaceFault.MinDistToFault(hypocenter_lat, hypocenter_lon, depth, strike, dip, magnitude, lat, lon)
                station['Rrup'] = distance
            else:
                station['Rrup'] = 'To determine'
        
        station['name'] = this_station_name
        station['dt'] = dt
        
        if abs(max_time - min_time) >= dt:
            sizes = np.int64((np.array(times) - min_time)/dt)
            ch[0] = np.hstack((np.zeros(sizes[0]), ch[0]))
            ch[1] = np.hstack((np.zeros(sizes[1]), ch[1]))
            ch[2] = np.hstack((np.zeros(sizes[2]), ch[2]))
        
        sizes = np.array([channel.shape[0] for channel in ch])
        sizes = sizes.max()-sizes
        
        ch[0] = np.hstack((ch[0], np.zeros(sizes[0])))
        ch[1] = np.hstack((ch[1], np.zeros(sizes[1])))
        ch[2] = np.hstack((ch[2], np.zeros(sizes[2])))
        
        for j,component in enumerate(components):
            if component.endswith('E') or component.endswith('X'):
                station['acc_1'] = ch[j].copy()
                station['component_1'] = component
            elif component.endswith('N') or component.endswith('Y'):
                station['acc_2'] = ch[j].copy()
                station['component_2'] = component
            else:
                station['acc_3'] = ch[j].copy()
                station['component_3'] = component
                
        start_time = event[:4] + '-' + event[4:6] + '-' + event[6:8] + ' '
        hours = int(min_time/3600.)
        minutes = int((min_time/3600. - hours)*60)
        seconds = (min_time- (hours*3600 + minutes*60))
        start_time += '%0.2i:%0.2i:%09.6f' %(hours, minutes, seconds)        
        station['start_time'] = start_time
        
        try:
            station['Vs30'] = stations_csn[this_station_name][2]
        except:
            station['Vs30'] = 'N/A'
                
        try:
            station['lat'] = stations_csn[this_station_name][0]
            station['lon'] = stations_csn[this_station_name][1]
        except:
            station['lat'] = lat
            station['lon'] = lon
            
        try:
            station['azimut'] = stations_csn[this_station_name][4]
        except:
            station['azimut'] = 'N/A'
        
        station['event_name'] = event
        station['magnitude'] = magnitude
        station['hypocenter_lat'] = hypocenter_lat
        station['hypocenter_lon'] = hypocenter_lon
        station['depth'] = depth
        
        dep_pos = np.nanargmin(np.sqrt((slab['lon'] - hypocenter_lon)**2 + (slab['lat'] - hypocenter_lat)**2))
        # La profundidad del slab es negativa
        difference = -slab['dep'][dep_pos] - depth
        
        if difference > 20:
            station['event_type'] = 'crustal'
        elif difference < -20:
            station['event_type'] = 'interface'
        else:
            station['event_type'] = 'intraslab'
        
        station['licensing'] = licensing
        station['units'] = 'Acceleration: m/s/s; Magnitude: Mw; dt: s; Depth: km; Vs30: m/s, Rrup: km'
        
        stations['st%0.2i' %i] = station
    
    scipy.io.savemat(matfile, stations)
    
    df = write_dataframe2(stations)
    
    return df, True

def writeMat(window, widget, event, events_path, save_path, replace_old = False):
    
    stations = sorted(os.listdir(os.path.join(events_path, event)))
    stations = [os.path.join(events_path, event, station) for station in stations]
    
    event_path = os.path.join(events_path, event)
    
    if int(event[:4]) == 1985:
        filenames = [filename for filename in stations if filename.endswith('.zip')]
        df, msg = transformCSN(filenames, event_path, save_path, replace_old)
        
    elif int(event[:4]) <= 2010: # Renadic v1
        df, msg = transformRenadicv1(stations, event_path, save_path, replace_old)
        
    else: # CSN
        filenames = [filename for filename in stations if filename.endswith('.zip')]
        df, msg = transformCSN(filenames, event_path, save_path, replace_old)
        
    if msg:
        widget.insert('end', 'Event ' + event + ' transformed\n')
        widget.see('end')
        window.update_idletasks()
    return df
