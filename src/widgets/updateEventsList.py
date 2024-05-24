# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 20:56:28 2022

@author: sebac
"""
import os
import json
import numpy as np
import pandas as pd
import lib.pytrend as pytrend

def updateEventsList(window, widget, basePath, dataPath):

    if not os.path.exists(os.path.join(basePath, 'tmp')):
        os.mkdir(os.path.join(basePath, 'tmp'))

    with open(os.path.join(basePath, 'secrets', 'pdd-secrets.json')) as f:
        credentials = json.load(f)

    session = pytrend.itrend_developer_tools()

    session.set_credentials(
        access_key_id     = credentials['access_key_id'],
        secret_access_key = credentials['secret_access_key'],
    )

    session.download_file('7c473f40d7862735', 'csv', filename = os.path.join(basePath, 'data', 'eventLists', '7c473f40d7862735.csv'))
    session.download_file('k2dimntt23gm0b7t', 'csv', filename = os.path.join(basePath, 'data', 'eventLists', 'k2dimntt23gm0b7t.csv'))
    session.download_file('ZW5TFERBT8B0GMQ' , 'csv', filename = os.path.join(basePath, 'data', 'eventLists', 'ZW5TFERBT8B0GMQ.csv'))

    df1 = pd.read_csv(os.path.join(basePath, 'data', 'eventLists', '7c473f40d7862735.csv'), dtype={'Identificador': str}) # 1985
    df2 = pd.read_csv(os.path.join(basePath, 'data', 'eventLists', 'k2dimntt23gm0b7t.csv'), dtype={'Identificador': str}) # 1994 a 2010
    df3 = pd.read_csv(os.path.join(basePath, 'data', 'eventLists', 'ZW5TFERBT8B0GMQ.csv') , dtype={'Identificador': str}) # 2012 a la fecha

    df1['Fuente'] = '7c473f40d7862735'
    df2['Fuente'] = 'k2dimntt23gm0b7t'
    df3['Fuente'] = 'ZW5TFERBT8B0GMQ'

    df = pd.concat([df1, df2, df3])
    events = []

    for r, row in df.iterrows():
        event_id = row['Fecha (UTC)'][:10].replace('-', '')
        if np.isnan(row['Magnitud [*]']):
            event_id += '_magM'
        else:
            event_id += '_' + '%0.1fM' %row['Magnitud [*]']
        if np.isnan(row['Latitud']):
            event_id += '_latS'
        else:
            event_id += '_' + '%0.2fS' %np.abs(row['Latitud'])
        if np.isnan(row['Longitud']):
            event_id += '_lonW'
        else:
            event_id += '_' + '%0.2fW' %np.abs(row['Longitud'])
        if np.isnan(row['Profundidad [km]']):
            event_id += '_depthKM'
        else:
            event_id += '_' + '%0.1fKM' %np.abs(row['Profundidad [km]'])
        events.append(event_id)

    df['ID'] = events
    df.to_csv(os.path.join(basePath, 'data', 'events.csv'), index=False)

    widget.insert('end', 'Base de datos de registros sísmicos disponibles actualizada!\n')
    window.update_idletasks()
