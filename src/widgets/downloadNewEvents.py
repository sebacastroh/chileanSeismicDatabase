# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 16:00:26 2023

@author: sebacastroh
"""
import os
import json
import pandas as pd
import lib.pytrend as pytrend

def downloadNewEvents(window, widget, basePath):
    
    if not os.path.exists(os.path.join(basePath, 'data', 'events.csv')):
        widget.insert('end', 'No existe el archivo "events.csv", primero ejecute la función para descargar la base de datos actualizada\n')
        window.update_idletasks()
        return

    if not os.path.exists(os.path.join(basePath, 'data', 'rawEvents')):
        os.mkdir(os.path.join(basePath, 'data', 'rawEvents'))

    df = pd.read_csv(os.path.join(basePath, 'data', 'events.csv'), dtype={'Identificador': str})
    
    with open(os.path.join(basePath, 'secrets', 'pdd-secrets.json')) as f:
        credentials = json.load(f)
        
    session = pytrend.itrend_developer_tools()
    
    session.set_credentials(
        access_key_id     = credentials['access_key_id'],
        secret_access_key = credentials['secret_access_key'],
    )
    
    for r, row in df.iterrows():
        rawName = row['Identificador']
        sdbName = row['ID']
        if os.path.exists(os.path.join(basePath, 'data', 'rawEvents', rawName + '.npz')) \
        or os.path.exists(os.path.join(basePath, 'data', 'seismicDatabase', 'npz', sdbName + '.npz')):
            continue
        
        dataset_id = row['Fuente']
        element_id = row['Identificador']
        
        widget.insert('end', 'Descargando archivo {element_id}... '.format(element_id=element_id))
        widget.see('end')
        window.update_idletasks()

        session.download_file(dataset_id, 'npz', element_id, filename = os.path.join(basePath, 'data', 'rawEvents', element_id + '.npz'))
    
    widget.insert('end', 'Done!\n')
    window.update_idletasks()
