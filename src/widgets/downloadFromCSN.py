# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 17:37:11 2024

@author: sebac
"""
import io
import os
import time
import json
import obspy
import urllib
import zipfile
import datetime
import numpy as np
import pandas as pd

DEFAULT_INDENT = 2
SORT_KEYS      = True

def is_json(val):
    try:
        json.dumps(val)
        return True
    except (TypeError, OverflowError):
        return False

def downloadNewEvents(window, widget, basePath, dataPath):

    with open(os.path.join(basePath, 'data', 'eventLists', 'registry.json')) as f:
        registry = json.load(f)

    records = []
    for event_id, stations in registry.items():
        for station_id, status in stations.items():
            if not isinstance(status, bool):
                records.append([event_id, station_id])

    failed_files = []
    for r, record in enumerate(records):        
        event_id, station_id = record

        if not os.path.exists(os.path.join(basePath, 'data', 'rawEvents', 'mseed', event_id + '.zip')):
            url = 'https://evtdb.csn.uchile.cl/raw/' + event_id
            urllib.request.urlretrieve(url, filename=os.path.join(basePath, 'data', 'rawEvents', 'mseed', event_id + '.zip'))
            time.sleep(0.5)

        if not os.path.exists(os.path.join(basePath, 'data', 'rawEvents', 'txt', event_id + '_' + station_id + '.zip')):
            try:
                url = 'https://evtdb.csn.uchile.cl/write/{event_id}/{station_id}'.format(
                    event_id   = event_id,
                    station_id = station_id
                )
                urllib.request.urlretrieve(url, filename=os.path.join(basePath, 'data', 'rawEvents', 'txt', event_id + '_' + station_id + '.zip'))
                time.sleep(0.5)
            except:
                failed_files.append([event_id, station_id])

    failed_files = pd.DataFrame(failed_files, columns=['Identificador', 'Estación'])
    # failed_files.to_csv(os.path.join('data', 'failed_files.csv'), index=False)

    current_event = None
    data = {}
    save = False
    for event_id, stations in registry.items():
        if event_id != current_event:            
            if save:
                np.savez_compressed(os.path.join(basePath, 'data', 'rawEvents', current_event + '.npz'), **data)
                widget.insert('end', 'Archivo creado o modificado {element_id}.npz\n'.format(element_id=current_event))
                widget.see('end')
                window.update_idletasks()
                save = False
            
            current_event = event_id
            
            filename = os.path.join(basePath, 'data', 'rawEvents', event_id + '.npz')
            data = {}
            
            if os.path.exists(filename):
                with np.load(filename, allow_pickle=True) as f:
                    for key, value in f.items():
                        data[key] = value.item()
        
        for station_id, status in stations.items():
            if status:
                continue
            
            if not os.path.exists(os.path.join(basePath, 'data', 'rawEvents', 'txt', event_id + '_' + station_id + '.zip')):
                continue
            
            save = True
            with zipfile.ZipFile(os.path.join(basePath, 'data', 'rawEvents', 'mseed', event_id + '.zip')) as zf:
                filenames = zf.namelist()

                evt_files   = []
                mseed_files = []

                for filename in filenames:
                    if filename.endswith('.evt'):
                        evt_files.append(filename)
                    else:
                        mseed_files.append(filename)

                metadata = False
                for evt_file in evt_files:
                    if station_id in evt_file:
                        metadata  = True
                        seed      = io.BytesIO(zf.read(evt_file))
                        stream    = obspy.read(seed)
                        seed.close()
                        break

                if not metadata:
                    for mseed_file in mseed_files:
                        seed   = io.BytesIO(zf.read(mseed_file))
                        stream = obspy.read(seed)
                        seed.close()
                        for trace in stream.traces:
                            if station_id == trace.meta['station'].strip():
                                metadata = True
                                break
                        if metadata:
                            break

            record = {}

            with zipfile.ZipFile(os.path.join(basePath, 'data', 'rawEvents', 'txt', event_id + '_' + station_id + '.zip')) as zf:
                filenames = zf.namelist()
                for filename in filenames:
                    txt = io.BytesIO(zf.read(filename))
                    acc = np.loadtxt(txt, skiprows=6)
                    txt.seek(0)
                    lines = []
                    while True:
                        line = txt.readline().decode().strip()
                        if not line.startswith('#'):
                            break
                        lines.append(line)
                    txt.close()
                    
                    npts = len(acc)

                    starttime = pd.Timestamp(lines[0].strip().split()[-1])
                    dt = 1./float(lines[1].strip().split()[4])
                    endtime = starttime + datetime.timedelta(seconds=(npts-1)*dt)
                    
                    t = pd.to_datetime(np.linspace(starttime.value, endtime.value, npts)).strftime("%Y-%m-%d %H:%M:%S.%f")
                    
                    channel = filename.split('-')[-1][:-4]
                    
                    latitude  = None
                    longitude = None
                    
                    if metadata:
                        for trace in stream.traces:
                            this_station = trace.meta['station'].strip()
                            if station_id != this_station:
                                continue
                            
                            this_channel = trace.meta['channel']
                            if this_channel != channel:
                                if this_channel.isdigit():
                                    if trace.meta.get('kinemetrics_evt') is not None:
                                        this_channel = trace.meta.kinemetrics_evt['chan_id']
                                        if this_channel != channel:
                                            continue
                                    else:
                                        continue
                                else:
                                    continue
                            
                            kinemetrics = trace.meta.get('kinemetrics_evt')
                            if kinemetrics is not None:
                                latitude  = trace.meta['kinemetrics_evt']['latitude']
                                longitude = trace.meta['kinemetrics_evt']['longitude']
                            else:
                                latitude  = float(lines[4].split()[2])
                                longitude = float(lines[4].split()[4])
                            break
                    
                    if metadata:
                        meta = dict(trace.meta)
                    else:
                        meta = {}
                    
                    record[channel] = {'x': np.array(t, dtype=np.datetime64),
                                        'y': acc.copy(),
                                        'metadata': meta,
                                        'location': {'lon': longitude,
                                                    'lat': latitude},
                                        'labels': {'x': 'Fecha (UTC)',
                                                'y': 'Aceleración m/s/s'}
                                        }
                    
                    try:
                        record[channel]['metadata']['starttime'] = str(record[channel]['m']['starttime'])
                        record[channel]['metadata']['endtime']   = str(record[channel]['m']['endtime'])
                    except:
                        record[channel]['metadata']['starttime'] = starttime.strftime('%Y-%m-%dT%H:%M:%S.%fZ')
                        record[channel]['metadata']['endtime']   = endtime.strftime('%Y-%m-%dT%H:%M:%S.%fZ')
            
                    if metadata and 'mseed' in record[channel]['metadata'].keys():
                        record[channel]['metadata']['mseed'] = dict(record[channel]['metadata']['mseed'])
                    elif metadata:
                        record[channel]['metadata']['kinemetrics_evt'] = dict(record[channel]['metadata']['kinemetrics_evt'])
                        for key in record[channel]['metadata']['kinemetrics_evt'].keys():
                            if not is_json(record[channel]['metadata']['kinemetrics_evt'][key]):
                                record[channel]['metadata']['kinemetrics_evt'][key] = str(record[channel]['metadata']['kinemetrics_evt'][key])
                
                data[station_id] = record.copy()

    widget.insert('end', 'Proceso completado.\n')
    window.update_idletasks()
