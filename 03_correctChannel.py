# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 17:13:41 2022

@author: sebac
"""
import os
import sys
import numpy as np
import automaticCorrection

if len(sys.argv) >= 4:
    event_id = sys.argv[1]
    station_id = sys.argv[2]
    channel_id = sys.argv[3]
    saveInTemp = True
    
    if len(sys.argv) >= 5:
        if sys.argv[4] == 'False':
            saveInTemp = False    
    
    with np.load(os.path.join('seismicDatabase', 'npz', event_id + '.npz'), allow_pickle=True) as f:
        acc_uncorrected = f.get(station_id).item().get(channel_id)
        dt = f.get(station_id).item().get('dt')
        status = f.get(station_id).item().get('p_wave').get('status')
        p_wave = f.get(station_id).item().get('p_wave').get('pos')
    
    filename = '-'.join([event_id, station_id, channel_id])
    acc_corrected, acc_filtered, corner_freqs = automaticCorrection.correctRecord(acc_uncorrected, dt, status, p_wave, saveInTemp, filename)
else:
    sys.exit('Expected at least 3 arguments, got %i' %(len(sys.argv) - 1))
