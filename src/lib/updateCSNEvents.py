# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 11:50:36 2024

@author: sebac
"""
import os
import json
import time
import urllib
import datetime
import requests
import pandas as pd

DEFAULT_INDENT = 2
SORT_KEYS      = True

def updateCSNEvents(window, widget, basePath, dataPath, draftPath, filename, tmp_file=None, start_event=None):

    new_events = []

    with open(os.path.join(basePath, 'data', 'eventLists', 'registry.json')) as f:
        registry = json.load(f)

    if tmp_file is not None:
        widget.insert('end', 'Cargando archivo temporal %s.\n' %tmp_file)
        widget.see('end')
        window.update_idletasks()
        new_events = pd.read_csv(os.path.join(basePath, 'tmp', tmp_file))
    else:
        widget.insert('end', 'Revisión de eventos disponibles en sitio web https://evtdb.csn.uchile.cl/.\n')
        widget.see('end')
        window.update_idletasks()

        url = 'https://evtdb.csn.uchile.cl/events'
        payload = {
            'min_date' : '2012-01-01',
            'max_date' : datetime.datetime.now().strftime('%Y-%m-%d'),
            'min_lat'  : '-90',
            'max_lat'  : '90',
            'min_lon'  : '-180',
            'max_lon'  : '180',
            'min_depth': '0',
            'max_depth': '9999',
            'min_mag'  : '0',
            'max_mag'  : '9999',
            'filter'   : 'Buscar'
        }

        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        }

        next_line = ''
        try:
            response = requests.post(url, data=payload, headers=headers)
            response.raise_for_status()
            website = response.text.split('\n')
            for line in website:
                if line.find('"/event/') >= 0:
                    event = line
                    pos = event.find('/event/')
                    uid = event[pos+7:-2]
                    next_line = 'date'
                    continue
                
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
        except:
            widget.insert('end', '\n¡Ha ocurrido un error al descargar la lista actualizada de eventos!\n')
            widget.see('end')
            window.update_idletasks()

            return False

    new_events = pd.DataFrame(new_events, columns=['Fecha (UTC)', 'Latitud', 'Longitud', 'Profundidad [km]', 'Magnitud [*]', 'Estaciones', 'Identificador'])
    new_events = new_events.sort_values(by=['Fecha (UTC)', 'Identificador']).reset_index(drop=True)

    widget.insert('end', '\nRevisión de estaciones dentro de eventos\n')
    widget.see('end')
    window.update_idletasks()

    if start_event is not None:
        start_event_pos = new_events[new_events['Identificador'] == start_event].iloc[0].name
    else:
        start_event_pos = None

    for r, row in new_events.iterrows():
        event_id = row['Identificador']
        url = 'http://evtdb.csn.uchile.cl/event/' + event_id

        if start_event_pos is not None and r < start_event_pos:
            if registry.get(event_id) is None:
                registry[event_id] = {}

            stations = row['Estaciones'].split('; ')
            for station in stations:
                if not isinstance(registry[event_id].get(station), bool):
                    registry[event_id][station] = None

            continue
        time.sleep(1)
        widget.insert('end', 'Obteniendo estaciones del evento %s (%i/%i)\n' %(event_id, r+1, len(new_events)))
        widget.see('end')
        window.update_idletasks()

        max_retries = 3
        success = False
        for attempt in range(max_retries):
            try:
                urllib.request.urlretrieve(url, filename=os.path.join('tmp', 'download.wget'))
                success = True
                break
            except:
                time.sleep(1.5)
                continue

        if not success:
            timestamp = datetime.datetime.utcnow().strftime('%Y%m%d%H%M%S')
            new_events.to_csv(os.path.join(basePath, 'tmp', filename + '_' + timestamp + '.csv'), index=False)

            widget.insert('end', '\n¡Ha ocurrido un error al descargar la página del evento %s!.\n' %event_id)
            widget.insert('end', 'Se han guardado los resultados parciales en la ruta %s.\n' %(os.path.join(basePath, 'tmp', filename + '_' + timestamp + '.csv')))
            widget.see('end')
            window.update_idletasks()

            return False
        
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

    old_events = pd.read_csv(os.path.join(basePath, 'data', 'eventLists', filename + '.csv'))
    removed_rows = []
    for r, row in old_events.iterrows():
        identifier     = row['Identificador']
        new_event_rows = new_events[new_events['Identificador'] == identifier]

        if len(new_event_rows) == 0:
            removed_rows.append(row)
            continue
        
        new_stations = new_event_rows.iloc[0]['Estaciones'].split('; ')

        changed = False
        for station in row['Estaciones'].split('; '):
            if station not in new_stations:
                new_stations.append(station)
                changed = True

        if changed:
            new_events.loc[new_event_rows.iloc[0].name, 'Estaciones'] = '; '.join(sorted(set(new_stations)))

    if len(removed_rows) > 0:
        new_events = pd.concat([new_events, pd.DataFrame(removed_rows)], ignore_index=True)

    new_events.sort_values(by=['Fecha (UTC)', 'Identificador'], inplace=True)
    new_events.to_csv(os.path.join(basePath, 'data', 'eventLists', filename + '.csv'), index=False)

    with open(os.path.join(basePath, 'data', 'eventLists', 'registry.json'), 'w') as f:
        json.dump(registry, f, indent=DEFAULT_INDENT, sort_keys=SORT_KEYS)

    return True
