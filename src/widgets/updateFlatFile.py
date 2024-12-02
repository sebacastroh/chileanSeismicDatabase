# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 09:29:28 2022

@author: sebac
"""
import os
import json
import datetime
import numpy as np
import pandas as pd

def updateFlatFile(window, widget, basePath, dataPath, draftPath):

    if not os.path.exists(os.path.join(dataPath, 'seismicDatabase', 'npz')):
        widget.insert('end', 'No existen registros almacenados para registrar en la base de datos.\n')
        window.update_idletasks()
        return

    if not os.path.exists(os.path.join(basePath, 'data', 'p_waves.json')):
        widget.insert('end', 'No se ha encontrado el archivo p_waves.json necesario.\n')
        window.update_idletasks()
        return

    columns=['Earthquake Name', 'Earthquake date',
        'Start time record', 'Magnitude [Mw]',
        'Hypocenter latitude', 'Hypocenter longitude',
        'Depth [km]', 'Event type', 'Station name',
        'Station code',
        'Station latitude', 'Station longitude', 'Station dt [s]',
        'Hypocentral distance [km]',
        'Epicentral distance [km]',
        'Rupture distance [km]',
        'Joyner-Boore distance [km]', 'Vs30 [m/s]',
        'Azimuth [o]', 'HVSR', 'Corrected records', 'Last update']

    if os.path.exists(os.path.join(draftPath, 'flatFile.csv')):
        df = pd.read_csv(os.path.join(draftPath, 'flatFile.csv'))
    elif os.path.exists(os.path.join(dataPath, 'flatFile.csv')):
        df = pd.read_csv(os.path.join(dataPath, 'flatFile.csv'))
    else:
        df = pd.DataFrame([], columns=columns)

    with open(os.path.join(basePath, 'data', 'p_waves.json')) as f:
        p_waves = json.load(f)

    table = []
    for event_id, p_wave in p_waves.items():
        station_codes = []
        for scode, sinfo in p_wave.items():
            subdf = df[(df['Earthquake Name'] == event_id) & (df['Station code'] == scode)]

            if len(subdf) == 0:
                station_codes.append(scode)
            else:
                lastUpdateCSV  = datetime.datetime.strptime(subdf['Last update'].iloc[0], '%Y-%m-%d %H:%M:%S.%f')
                lastUpdateJSON = datetime.datetime.strptime(sinfo['updated'], '%Y-%m-%dT%H:%M:%S.%f')
                if lastUpdateCSV != lastUpdateJSON:
                    station_codes.append(scode)
                    df.drop(labels=subdf.index[0], inplace=True)

        if len(station_codes) == 0:
            continue

        if os.path.exists(os.path.join(draftPath, 'seismicDatabase', 'npz', event_id + '.npz')):
            filename = os.path.join(draftPath, 'seismicDatabase', 'npz', event_id + '.npz')
        else:
            filename = os.path.join(dataPath, 'seismicDatabase', 'npz', event_id + '.npz')

        with np.load(filename, allow_pickle=True) as stations:
            for key in sorted(stations.keys()):
                if not key.startswith('st'):
                    continue

                station = stations.get(key).item()

                if station.get('station_code') not in station_codes:
                    continue

                Rrup = station.get('Rrup')
                Rjb  = station.get('Rjb')
                if isinstance(Rrup, str):
                    Rrup = -1
                    Rjb = -1    
                    
                Vs30 = station.get('vs30')
                if isinstance(Vs30, str):
                    Vs30 = -1
                
                azimuth = station.get('azimuth')
                if isinstance(azimuth, str):
                    azimuth = -1

                corrected = p_wave.get(station.get('station_code')).get('corrected')

                table.append([
                    event_id,
                    '-'.join([event_id[:4], event_id[4:6], event_id[6:8]]), # Earthquake date
                    station.get('starttime').replace('T', ' ').replace('Z', ''),
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
                    Rrup,
                    Rjb,
                    Vs30,
                    azimuth,
                    station.get('hvsr'),
                    corrected,
                    station.get('last_update').replace('T', ' ').replace('Z', '')
                ])
    
    if len(table) > 0:
        new_rows = pd.DataFrame(table, columns=columns)
        df = pd.concat([df, new_rows], ignore_index=True)
        df.sort_values(by=['Earthquake Name', 'Station code'], inplace=True)
        
        for col in df.columns[[3,4,5,6,10,11,12,13,14,15,16,17,18]]:
            df[col] = df[col].astype(float)

        df['Earthquake date']   = pd.to_datetime(df['Earthquake date'], format='%Y-%m-%d')
        df['Start time record'] = pd.to_datetime(df['Start time record'], format='%Y-%m-%d %H:%M:%S.%f')
        df['Last update']       = pd.to_datetime(df['Last update'], format='%Y-%m-%d %H:%M:%S.%f')

        df.to_excel(os.path.join(draftPath, 'flatFile.xlsx'), index=False)
        df.to_csv(os.path.join(draftPath, 'flatFile.csv'), index=False)
        df.to_csv(os.path.join(basePath, 'data', 'flatFile - backup.csv'), index=False)
    
    widget.insert('end', 'Flat file actualizado.\n')
    window.update_idletasks()
