# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 10:36:19 2022

@author: sebacastroh
"""

import os
import sys
import subprocess

if len(sys.argv) >= 3:
    event_id = sys.argv[1]
    station_id = sys.argv[2]
    saveInTemp = True
    
    if len(sys.argv) >= 4:
        if sys.argv[3] == 'False':
            saveInTemp = False    
    
    cwd = os.getcwd()
    filepath = os.path.join(cwd, '03_correctRecord.py')
    cmd = 'python ' + filepath + ' ' + event_id + ' ' + station_id
    for i in range(3):
        subprocess.run(cmd + ' ' + 'acc_uncorrected_%i' %(i+1) + ' ' + str(saveInTemp))
else:
    sys.exit('Expected at least 2 arguments, got %i' %(len(sys.argv) - 1))
