# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 10:14:02 2022

@author: sebacastroh
"""
import os
import sys
import datetime
import subprocess
import numpy as np
import scipy.io as spio

if len(sys.argv) >= 1:
    event_id = sys.argv[1]
    
    nStations = 0
    with np.load(os.path.join('seismicDatabase', 'npz', event_id + '.npz'), allow_pickle=True) as f:
        for key in f.keys():
            if key.startswith('st'):
                nStations += 1
    
    cwd = os.getcwd()
    filepath = os.path.join(cwd, '02_correctStation.py')
    cmd = 'python ' + filepath + ' ' + event_id
    for i in range(nStations):
        subprocess.run(cmd + ' ' + 'st%0.2i' %i + ' True', shell=True)
    
    with np.load(os.path.join('seismicDatabase', 'npz', event_id + '.npz'), allow_pickle=True) as f:
        event = {}
        for key, value in f.items():
            event[key] = value.item()
    
    for i in range(nStations):
        station_id = 'st%0.2i' %i
        for j in range(3):
            channel_id = 'acc_uncorrected_%i.npz' %(j+1)
            filename = '-'.join([event_id, station_id, channel_id])
            if os.path.exists(os.path.join('tmp', filename)):
                with np.load(os.path.join('tmp', filename)) as f:
                    acc_corr = f['acc_corr']
                    acc_fil = f['acc_fil']
                    freqs = f['freqs']
                
                event[station_id]['acc_corrected_%i' %(j+1)] = acc_corr.copy()
                event[station_id]['acc_filtered_%i' %(j+1)] = acc_fil.copy()
                event[station_id]['corner_freqs_%i' %(j+1)] = freqs.copy()
                event[station_id]['last_update'] = datetime.datetime.now().isoformat()
                
                os.remove(os.path.join('tmp', filename))
                    
    np.savez_compressed(os.path.join('seismicDatabase', 'npz', event_id), **event)
    spio.savemat(os.path.join('seismicDatabase', 'mat', event_id + '.mat'), event, do_compression=True)
else:
    sys.exit('Expected 1 argument, got 0')
