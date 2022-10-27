# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 12:22:18 2022

@author: sebacastroh
"""
import os
import sys
import subprocess
import multiprocessing

def callCorrectChannel(key):
    event_id, station_id, channel_id = key
    cwd = os.getcwd()
    filepath = os.path.join(cwd, '03_correctChannel.py')
    cmd = 'python ' + filepath + ' ' + event_id + ' ' + station_id + ' ' + channel_id + ' True'
    subprocess.run(cmd, shell=True)

if len(sys.argv) >= 3:
    event_id = sys.argv[1]
    station_id = sys.argv[2]
    
    keys = [[event_id, station_id, 'acc_uncorrected_1'],
            [event_id, station_id, 'acc_uncorrected_2'],
            [event_id, station_id, 'acc_uncorrected_3']]
    
    processes = 3
    
    if __name__ == '__main__':
        pool = multiprocessing.Pool(processes=processes)
        pool.map(callCorrectChannel, keys)
        pool.close()
else:
    sys.exit('Expected 2 argument, got %i' %(len(sys.argv)-1))
