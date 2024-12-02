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
import pandas as pd

DEFAULT_INDENT = 2
SORT_KEYS      = True

def updateCSNEvents(window, widget, basePath, dataPath, draftPath, filename, tmp_file=None, start_event=''):

    new_events = []
    all_done = False
    i = 1

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

    while not all_done:
        if tmp_file is not None:
            break

        time.sleep(1)
        url = 'http://evtdb.csn.uchile.cl/?page=%i' %i
        try:
            urllib.request.urlretrieve(url, filename=os.path.join(basePath, 'tmp', 'download.wget'))
        except:
            all_done = True
            break

        widget.insert('end', 'Página número %0.3i completada.\n' %i)
        widget.see('end')
        window.update_idletasks()

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
            
    new_events = pd.DataFrame(new_events, columns=['Fecha (UTC)', 'Latitud', 'Longitud', 'Profundidad [km]', 'Magnitud [*]', 'Estaciones', 'Identificador'])

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

        if r < start_event_pos:
            if registry.get(event_id) is None:
                registry[event_id] = {}

            stations = row['Estaciones'].split('; ')
            for station in stations:
                if not isinstance(registry[event_id].get(station), bool):
                    registry[event_id][station] = None

            continue
        time.sleep(1)
        widget.insert('end', 'Obteniendo estaciones del evento %s\n' %event_id)
        widget.see('end')
        window.update_idletasks()

        try:
            urllib.request.urlretrieve(url, filename=os.path.join('tmp', 'download.wget'))
        except:
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
