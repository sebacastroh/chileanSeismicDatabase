# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 20:56:28 2022

@author: sebac
"""
import os
import numpy as np
import pandas as pd
from lib.updateCSNEvents import updateCSNEvents

def updateEventsList(window, widget, basePath, dataPath, draftPath):

    if not os.path.exists(os.path.join(basePath, 'tmp')):
        os.mkdir(os.path.join(basePath, 'tmp'))

    filename_COSMOS  = '7c473f40d7862735'
    filename_RENADIC = 'k2dimntt23gm0b7t'
    filename_CSN     = 'ZW5TFERBT8B0GMQ'

    status = updateCSNEvents(window, widget, basePath, dataPath, draftPath, filename_CSN)

    if not status:
        return

    df1 = pd.read_csv(os.path.join(basePath, 'data', 'eventLists', filename_COSMOS  + '.csv'), dtype={'Identificador': str}) # 1985
    df2 = pd.read_csv(os.path.join(basePath, 'data', 'eventLists', filename_RENADIC + '.csv'), dtype={'Identificador': str}) # 1994 a 2010
    df3 = pd.read_csv(os.path.join(basePath, 'data', 'eventLists', filename_CSN     + '.csv'), dtype={'Identificador': str}) # 2012 a la fecha

    df1['Fuente'] = filename_COSMOS
    df2['Fuente'] = filename_RENADIC
    df3['Fuente'] = filename_CSN

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
