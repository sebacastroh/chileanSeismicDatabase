#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 21:21:58 2019

@author: srcastro
"""

import os
import subprocess
import multiprocessing

cwd = os.getcwd()
old_path = 'events_mat_uncorrected'
new_path = os.path.join(cwd, 'events_mat_corrected_v2')

cmd = 'python ' + os.path.join(cwd, 'correctRecords.py ')

def run(filename):
    global cmd
    subprocess.call(cmd + filename + ' ' + new_path, shell=True)

filenames = os.listdir(os.path.join(cwd, old_path))

these_filenames = []
for filename in sorted(filenames):
    if not os.path.exists(os.path.join(new_path, filename)):
        these_filenames.append(filename)

if __name__ == '__main__':
    pool = multiprocessing.Pool(processes=56)
    pool.map(run, these_filenames)
    pool.close()