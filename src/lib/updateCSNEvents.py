# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 11:50:36 2024

@author: sebac
"""
import os
import json
import time
import urllib
import pandas as pd

DEFAULT_INDENT = 2
SORT_KEYS      = True

def updateCSNEvents(basePath, dataPath, filename):

    new_events = []
    all_done = False
    i = 1

    with open(os.path.join(basePath, 'data', 'eventLists', 'registry.json')) as f:
        registry = json.load(f)

    while not all_done:
        time.sleep(1)
        url = 'http://evtdb.csn.uchile.cl/?page=%i' %i
        try:
            urllib.request.urlretrieve(url, filename=os.path.join(basePath, 'tmp', 'download.wget'))
        except:
            all_done = True
            break
        next_line = ''
        saved = False
        with open(os.path.join(basePath, 'tmp', 'download.wget'), 'r') as fopen:
            for line in fopen:
                if line.startswith('<title>500 Internal Server Error</title>'):
                    all_done = True
                    break
                
                if line.find('"/event/') >= 0:
                    event = line
                    pos = event.find('/event/')
                    uid = event[pos+7:-3]
                    next_line = 'date'
                    continue
                
                # if not saved:
                if next_line == 'date':
                    date = line.strip()
                    next_line = ''

                elif line.find('<td class="latitude">') >= 0:
                    next_line = 'latitude'

                elif next_line == 'latitude':
                    lat = line.strip()
                    next_line = ''

                elif line.find('<td class="longitude">') >= 0:
                    next_line = 'longitude'

                elif next_line == 'longitude':
                    lon = line.strip()
                    next_line = ''

                elif line.find('<td class="depth">') >= 0:
                    next_line = 'depth'

                elif next_line == 'depth':
                    depth = line.strip()
                    next_line = ''

                elif line.find('<td class="magnitude">') >= 0:
                    next_line = 'magnitude'

                elif next_line == 'magnitude':
                    mag = line.strip()
                    row = [date, float(lat), float(lon), float(depth), float(mag), '', uid]
                    next_line = ''
                    new_events.append(row)

        if not all_done:
            i += 1
            os.remove(os.path.join('tmp', 'download.wget'))
            
    new_events = pd.DataFrame(new_events, columns=['Fecha', 'Latitud', 'Longitud', 'Profundidad', 'Magnitud', 'Estaciones', 'Identificador'])

    for r, row in new_events.iterrows():
        time.sleep(1)
        event_id = row['ID']
        url = 'http://evtdb.csn.uchile.cl/event/' + event_id
        urllib.request.urlretrieve(url, filename=os.path.join('tmp', 'download.wget'))
        
        if registry.get(event_id) is None:
            registry[event_id] = {}

        stations = []
        with open(os.path.join('tmp', 'download.wget'), 'r') as fopen:
            for line in fopen:
                if line.find('/write/') >= 0:
                    pos = line.find("/write/")
                    sl = line[pos:].split("/")
                    evt = sl[2]
                    sta = sl[3][:-2].strip()
                    stations.append(sta)

                    if not isinstance(registry[event_id].get(sta), bool):
                        registry[event_id][sta] = None
                    
        stations = '; '.join(sorted(set(stations)))
        new_events.iloc[r,5] = stations
        
        os.remove(os.path.join('tmp', 'download.wget'))

    new_events.to_csv(os.path.join(basePath, 'data', 'eventLists', filename + '.csv'), index=False)

    with open(os.path.join(basePath, 'data', 'eventLists', 'registry.json'), 'w') as f:
        json.dump(registry, f, indent=DEFAULT_INDENT, sort_keys=SORT_KEYS)
