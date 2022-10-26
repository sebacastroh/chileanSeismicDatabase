#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 08:48:03 2020

@author: srcastro
"""

import os
import zipfile as zp

database_path = '/home/srcastro/projects/correctSeismicDatabase/newMethodology/events_mat_corrected_v2/'
zipfile = 'SIBER-RISK_database_matlab.zip'

# database_path = '/home/srcastro/projects/correctSeismicDatabase/newMethodology/events_pkl_corrected_v2/'
# zipfile = 'SIBER-RISK_database_python.zip'

all_files = os.listdir(database_path)
files_to_add = []
for filename in all_files:
    if int(filename[:8]) >= 20210701:
    # if filename.startswith('20210527_4.3'):
        files_to_add.append(filename)

with zp.ZipFile(zipfile, 'a', zp.ZIP_DEFLATED) as Zp:
    for filename in sorted(files_to_add):
        Zp.write(os.path.join(database_path, filename), filename)
