# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 09:29:28 2022

@author: sebac
"""
import os
import json
import numpy as np
import pandas as pd

def updateFlatFile(window, widget, basePath):

    if not os.path.exists(os.path.join(basePath, 'data', 'seismicDatabase', 'npz')):
        widget.insert('end', 'No existen registros almacenados para registrar en la base de datos.\n')
        window.update_idletasks()
        return

    filenames = sorted(os.listdir(os.path.join(basePath, 'data', 'seismicDatabase', 'npz')))
    if len(filenames) == 0:
        widget.insert('end', 'No existen registros almacenados para registrar en la base de datos.\n')
        window.update_idletasks()
        return

    table = []
    
    with open(os.path.join(basePath, 'data', 'p_waves.json')) as f:
        p_waves = json.load(f)

    for filename in filenames:
        p_wave_info = p_waves.get(filename[:-4])

        with np.load(os.path.join(basePath, 'data', 'seismicDatabase', 'npz', filename), allow_pickle=True) as stations:
            for key in sorted(stations.keys()):
                if not key.startswith('st'):
                    continue

                station = stations.get(key).item()

                Rrup = station.get('Rrup')
                Rjb  = station.get('Rjb')
                if isinstance(Rrup, str):
                    Rrup = -1
                    Rjb = -1    
                    
                Vs30 = station.get('Vs30')
                if isinstance(Vs30, str):
                    Vs30 = -1
                
                azimut = station.get('azimut')
                if isinstance(azimut, str):
                    azimut = -1

                if p_wave_info is not None:
                    station_info = p_wave_info.get(station.get('station_code'))
                    if station_info is None:
                        corrected = False
                    elif station_info.get('corrected'):
                        corrected = True
                    else:
                        corrected = False
                else:
                    corrected = False

                table.append([
                    filename[:-4],                                          # Earthquake name
                    '-'.join([filename[:4], filename[4:6], filename[6:8]]), # Earthquake date
                    station.get('starttime'),
                    station.get('magnitude'),
                    station.get('hypocenter_lat'),
                    station.get('hypocenter_lon'),
                    station.get('depth'),
                    station.get('event_type'),
                    station.get('station_name'),
                    station.get('station_code'),
                    station.get('station_lat'),
                    station.get('station_lon'),
                    station.get('dt'),
                    station.get('Rhypo'),
                    station.get('Repi'),
                    station.get('Rrup'),
                    station.get('Rjb'),
                    station.get('vs30'),
                    station.get('azimuth'),
                    station.get('hvsr'),
                    corrected,
                    station.get('last_update')
                ])
    
    df = pd.DataFrame(np.array(table), columns=['Earthquake Name', 'Earthquake date',
                      'Start time record', 'Magnitude [Mw]',
                      'Hypocenter latitude', 'Hypocenter longitude',
                      'Depth [km]', 'Event type', 'Station name',
                      'Station code',
                      'Station latitude', 'Station longitude', 'Station dt [s]',
                      'Hypocentral distance [km]',
                      'Epicentral distance [km]',
                      'Rupture distance [km]',
                      'Joyner-Boore distance [km]', 'Vs30 [m/s]',
                      'Azimut [o]', 'HVSR', 'Corrected records', 'Last update'])
    
    for col in df.columns[[3,4,5,6,10,11,12,13,14,15,16,17,18]]:
        df[col] = df[col].astype(float) 

    df['Earthquake date'] = pd.to_datetime(df['Earthquake date'], format='%Y-%m-%d')
    df['Start time record'] = pd.to_datetime(df['Start time record'], format='%Y-%m-%dT%H:%M:%S.%fZ')
    df['Last update'] = pd.to_datetime(df['Last update'], format='%Y-%m-%dT%H:%M:%S.%f')

    df.to_excel(os.path.join(basePath, 'data', 'flatFile.xlsx'), index=False)
    df.to_csv(os.path.join(basePath, 'data', 'flatFile.csv'), index=False)
    
    widget.insert('end', 'Flat file actualizado.\n')
    window.update_idletasks()