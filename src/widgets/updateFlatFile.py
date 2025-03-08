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
    elif os.path.exists(os.path.join(basePath, 'data', 'flatFile - backup.csv')):
        df = pd.read_csv(os.path.join(basePath, 'data', 'flatFile - backup.csv'))
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

        # Load original flatfile
        if os.path.exists(os.path.join(draftPath, 'flatFile.csv')):
            df_old = pd.read_csv(os.path.join(draftPath, 'flatFile.csv'), parse_dates=['Earthquake date', 'Start time record', 'Last update'])
            df_old_path = os.path.join(draftPath, 'flatFile.csv')
        elif os.path.exists(os.path.join(dataPath, 'flatFile.csv')):
            df_old = pd.read_csv(os.path.join(dataPath, 'flatFile.csv'), parse_dates=['Earthquake date', 'Start time record', 'Last update'])
            df_old_path = os.path.join(dataPath, 'flatFile.csv')
        else:
            df_old = pd.DataFrame([], columns=columns)
            df_old_path = None

        # Save new/modified lines
        df_merge = pd.merge(df_old, df, how='right', on=['Earthquake Name', 'Station code'])
        new_rows = df[df_merge.apply(lambda row: row['Last update_x'] != row['Last update_y'], axis=1)]

        new_rows.to_csv(os.path.join(basePath, 'tmp', 'new_rows.csv'), index=False)

        with open(os.path.join(basePath, 'tmp', 'new_rows.csv')) as f:
            new_lines = f.read().split('\n')[1:]
            
        os.remove(os.path.join(basePath, 'tmp', 'new_rows.csv'))

        # Insert new lines in original flatfile
        if df_old_path is not None:
            with open(df_old_path) as f:
                new_csv = f.readline()
                for line in f:
                    event_id     = line.split(',')[0]
                    station_code = line.split(',')[9]
                    skip = False

                    if len(new_lines) > 0:
                        new_event_id     = new_lines[0].split(',')[0]
                        new_station_code = new_lines[0].split(',')[9]

                    while len(new_lines) > 0 and event_id >= new_event_id and station_code >= new_station_code:
                        new_csv += new_lines[0].strip() + '\n'

                        if event_id == new_event_id and station_code == new_station_code:
                            skip = True

                        new_lines.pop(0)

                        if len(new_lines) > 0:
                            new_event_id     = new_lines[0].split(',')[0]
                            new_station_code = new_lines[0].split(',')[9]

                    if skip:
                        continue

                    new_csv += line

            for new_line in new_lines:
                new_csv += new_line.strip() + '\n'

        # Save flatfile
        df.to_excel(os.path.join(draftPath, 'flatFile.xlsx'), index=False)

        if df_old_path is None:
            df.to_csv(os.path.join(draftPath, 'flatFile.csv'), index=False)
            df.to_csv(os.path.join(basePath, 'data', 'flatFile - backup.csv'), index=False)
        else:
            with open(os.path.join(draftPath, 'flatFile.csv'), 'w') as f:
                f.write(new_csv)

            with open(os.path.join(basePath, 'data', 'flatFile - backup.csv'), 'w') as f:
                f.write(new_csv)

    widget.insert('end', 'Flat file actualizado.\n')
    window.update_idletasks()
