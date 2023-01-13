# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 16:00:26 2023

@author: sebacastroh
"""
import os
import json
import pytrend
import pandas as pd

def downloadNewEvents(window, widget):
    
    df = pd.read_csv('siberrisk.csv', dtype={'Identificador': str})
    
    with open('pdd-secrets.json') as f:
        credentials = json.load(f)
        
    session = pytrend.itrend_developer_tools()
    
    session.set_credentials(
        access_key_id = credentials['access_key_id'],
        secret_access_key = credentials['secret_access_key'],
    )
    
    for r, row in df.iterrows():
        name = row['ID']
        if os.path.exists(os.path.join('seismicDatabase', 'npz', name + '.npz')):
            continue
        
        dataset_id = row['Fuente']
        element_id = row['Identificador']
        
        session.download_file(dataset_id, 'npz', element_id, os.path.join('rawData', element_id + '.npz'))
    
    widget.insert('end', 'Done!\n')
    window.update_idletasks()
