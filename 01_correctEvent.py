# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 10:14:02 2022

@author: sebacastroh
"""
import os
import sys
import subprocess
import numpy as np

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
        break
else:
    sys.exit('Expected 1 argument, got 0')
