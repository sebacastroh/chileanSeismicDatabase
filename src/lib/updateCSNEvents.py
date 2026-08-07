# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 11:50:36 2024

@author: sebac
"""
import os
import re
import json
import time
import datetime
import requests
import pandas as pd

DEFAULT_INDENT = 2
SORT_KEYS      = True

def updateCSNEvents(window, widget, basePath, dataPath, draftPath, filename, tmp_file=None, start_event=None):

    new_events = []

    with open(os.path.join(basePath, 'data', 'eventLists', 'registry.json')) as f:
        registry = json.load(f)

    session = requests.Session()
    session.headers.update({
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
        'Accept-Encoding': 'gzip, deflate'
    })

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

        try:
            response = session.post(url, data=payload, timeout=30)
            response.raise_for_status()

            pattern = re.compile(
                r'href="/event/([a-f0-9]+)".*?>\s*(.*?)\s*</a>.*?'
                r'<td class="latitude">\s*(.*?)\s*</td>.*?'
                r'<td class="longitude">\s*(.*?)\s*</td>.*?'
                r'<td class="depth">\s*(.*?)\s*</td>.*?'
                r'<td class="magnitude">\s*(.*?)\s*</td>',
                re.DOTALL
            )
            
            matches = pattern.findall(response.text)
            for uid, date, lat, lon, depth, mag in matches:
                new_events.append([date, float(lat), float(lon), float(depth), float(mag), '', uid])

        except Exception as e:
            widget.insert('end', f'\n¡Ha ocurrido un error al descargar la lista actualizada de eventos!: {e}\n')
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

    station_regex = re.compile(r'/write/[^/]+/([^"/]+)')

    total_events = len(new_events)
    stations_list_column = []

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
            
            stations_list_column.append(row['Estaciones'])
            continue

        widget.insert('end', f'Obteniendo estaciones del evento {event_id} ({r+1}/{total_events})\n')
        widget.see('end')
        window.update_idletasks()

        html_content = None
        for attempt in range(3):
            try:
                res = session.get(url, timeout=10)
                if res.status_code == 200:
                    html_content = res.text
                    break
            except requests.RequestException:
                time.sleep(2)

        if html_content is None:
            timestamp = datetime.datetime.utcnow().strftime('%Y%m%d%H%M%S')
            new_events.to_csv(os.path.join(basePath, 'tmp', filename + '_' + timestamp + '.csv'), index=False)

            widget.insert('end', '\n¡Ha ocurrido un error al descargar la página del evento %s!.\n' %event_id)
            widget.insert('end', 'Se han guardado los resultados parciales en la ruta %s.\n' %(os.path.join(basePath, 'tmp', filename + '_' + timestamp + '.csv')))
            widget.see('end')
            window.update_idletasks()

            return False

        if event_id not in registry:
            registry[event_id] = {}

        found_stations = sorted(list(set(station_regex.findall(html_content))))
        
        for sta in found_stations:
            if not isinstance(registry[event_id].get(sta), bool):
                registry[event_id][sta] = None

        stations_str = '; '.join(found_stations)
        stations_list_column.append(stations_str)

    new_events['Estaciones'] = stations_list_column

    old_events = pd.read_csv(os.path.join(basePath, 'data', 'eventLists', filename + '.csv'))
    
    old_events_map = dict(zip(old_events['Identificador'], old_events['Estaciones']))
    new_events_map = dict(zip(new_events['Identificador'], new_events['Estaciones']))
    
    for ident, old_sta_str in old_events_map.items():
        if ident in new_events_map:
            old_set = set(str(old_sta_str).split('; ')) if pd.notna(old_sta_str) else set()
            new_set = set(str(new_events_map[ident]).split('; ')) if pd.notna(new_events_map[ident]) else set()
            
            combined = sorted(list(old_set.union(new_set)))
            new_events.loc[new_events['Identificador'] == ident, 'Estaciones'] = '; '.join(combined)
        else:
            row_to_add = old_events[old_events['Identificador'] == ident]
            new_events = pd.concat([new_events, row_to_add], ignore_index=True)

    new_events.sort_values(by=['Fecha (UTC)', 'Identificador'], inplace=True)
    new_events.to_csv(os.path.join(basePath, 'data', 'eventLists', filename + '.csv'), index=False)

    with open(os.path.join(basePath, 'data', 'eventLists', 'registry.json'), 'w') as f:
        json.dump(registry, f, indent=DEFAULT_INDENT, sort_keys=SORT_KEYS)

    return True
