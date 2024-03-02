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

    if os.path.exists(os.path.join(basePath, 'data', 'flatFile.csv')):
        df = pd.read_csv(os.path.join(basePath, 'data', 'flatFile.csv'))
    else:
        df = pd.DataFrame([], columns=columns)

    table = []
    
    with open(os.path.join(basePath, 'data', 'p_waves.json')) as f:
        p_waves = json.load(f)

    for filename in filenames:
        subdf = df[df['Earthquake Name'] == filename[:-4]]
        p_wave_info = p_waves.get(filename[:-4])

        updateRows = False

        if len(subdf) == 0:
            updateRows = True
        else:
            if p_wave_info is None:
                updateRows = True
            else:
                stations_codes_json = list(p_wave_info.keys())
                station_codes_csv   = subdf['Station code'].tolist()

                for station_code_json in stations_codes_json:
                    if station_code_json not in station_codes_csv:
                        updateRows = True
                        break
                    else:
                        row = subdf[subdf['Station code'] == station_code_json]
                        lastUpdateCSV  = datetime.datetime.strptime(row['Last update'].iloc[0], '%Y-%m-%d %H:%M:%S.%f')
                        lastUpdateJSON = datetime.datetime.strptime(p_wave_info[station_code_json]['updated'], '%Y-%m-%dT%H:%M:%S.%f')

                        if lastUpdateCSV != lastUpdateJSON:
                            updateRows = True
                            break

        if not updateRows:
            continue

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
                
                azimuth = station.get('azimuth')
                if isinstance(azimuth, str):
                    azimuth = -1

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
    
    new_rows = pd.DataFrame(table, columns=columns)
    df = pd.concat([df, new_rows], ignore_index=True)
    
    for col in df.columns[[3,4,5,6,10,11,12,13,14,15,16,17,18]]:
        df[col] = df[col].astype(float) 

    df['Earthquake date'] = pd.to_datetime(df['Earthquake date'], format='%Y-%m-%d')
    df['Start time record'] = pd.to_datetime(df['Start time record'], format='%Y-%m-%d %H:%M:%S.%f')
    df['Last update'] = pd.to_datetime(df['Last update'], format='%Y-%m-%d %H:%M:%S.%f')

    df.to_excel(os.path.join(basePath, 'data', 'flatFile.xlsx'), index=False)
    df.to_csv(os.path.join(basePath, 'data', 'flatFile.csv'), index=False)
    
    widget.insert('end', 'Flat file actualizado.\n')
    window.update_idletasks()
