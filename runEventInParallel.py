# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 12:11:37 2022

@author: sebacastroh
"""
import os
import sys
import subprocess
import numpy as np
import multiprocessing

def callCorrectStation(key):
    event_id, station_id = key
    cwd = os.getcwd()
    filepath = os.path.join(cwd, '02_correctStation.py')
    cmd = 'python ' + filepath + ' ' + event_id
    subprocess.run(cmd + ' ' + station_id + ' True', shell=True)

if len(sys.argv) >= 1:
    event_id = sys.argv[1]
    nStations = 0
    keys = []
    with np.load(os.path.join('seismicDatabase', 'npz', event_id + '.npz'), allow_pickle=True) as f:
        for key in f.keys():
            if key.startswith('st'):
                keys.append([event_id, key])
                nStations += 1
    
    processes = min(nStations, 50)
    
    if __name__ == '__main__':
        pool = multiprocessing.Pool(processes=processes)
        pool.map(callCorrectStation, keys)
        pool.close()
else:
    sys.exit('Expected 1 argument, got 0')
