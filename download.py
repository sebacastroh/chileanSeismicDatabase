#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 21:40:28 2019

@author: srcastro
"""

import os
import wget
import shutil
import urllib.request

def download(window, widget):
    
    if os.path.exists('download.wget'):
        os.remove('download.wget')
    
    all_events = []
    all_done = False
    i = 1
    
    while not all_done:
        url = 'http://evtdb.csn.uchile.cl/?page=%i' %i
        widget.insert('end', 'Checking web ' + url + '\n')
        window.update_idletasks()
        wget.download(url)
        next_line = ''
        with open('download.wget', 'r') as fopen:
            for line in fopen:
                if line.startswith('<title>500 Internal Server Error</title>'):
                    all_done = True
                    break
                
                if line.find('"/event/') >= 0:
                    event = line
                    next_line = 'date'
                    
                elif next_line == 'date':
                    date = line.strip().split()[0].replace('-','')
                    next_line = ''
    
                elif line.find('<td class="latitude">') >= 0:
                    next_line = 'latitude'
    
                elif next_line == 'latitude':
                    lat = line.strip().split()[0].replace('-','') + 'S'
                    next_line = ''
    
                elif line.find('<td class="longitude">') >= 0:
                    next_line = 'longitude'
    
                elif next_line == 'longitude':
                    lon = line.strip().split()[0].replace('-','') + 'W'
                    next_line = ''
    
                elif line.find('<td class="depth">') >= 0:
                    next_line = 'depth'
    
                elif next_line == 'depth':
                    depth = line.strip().split()[0] + 'KM'
                    next_line = ''
    
                elif line.find('<td class="magnitude">') >= 0:
                    next_line = 'magnitude'
    
                elif next_line == 'magnitude':
                    mag = line.strip().split()[0] + 'Mw'
                    directory = date + '_' + mag + '_' + lat + '_' + lon + '_' + depth
                    next_line = ''
    
                    if os.path.exists('./events/' + directory):
                        all_done = True
                        break
                    else:
                        all_events.append(event)
    
        if not all_done:
            i += 1
            os.remove('download.wget')
    
    os.remove('download.wget')
    
    if len(all_events) == 0:
        widget.insert('end', '\nThere are 0 events to add\n')
        window.update_idletasks()
    else:
        widget.insert('end', '\nThere are %i events to check\n' %len(all_events))
        window.update_idletasks()
        
        basedir = os.getcwd()
        
        directory = './phase2'
        if not os.path.exists(directory):
            os.makedirs(directory)
        else:    
            shutil.rmtree(directory)
            os.makedirs(directory)
        
        os.chdir('./phase2')
        for event in all_events:
            pos = event.find('/event/')
            if pos >= 0:
                evt = event[pos:-3]
                url = 'http://evtdb.csn.uchile.cl' + evt
                wget.download(url)
                
        os.chdir(basedir)
    
        files = os.listdir('./phase2/')
        
        if not os.path.exists('./events'):
            os.makedirs('./events')
        
        os.chdir('./events')
        basedir2 = os.getcwd()
    
        def download_zip(f):
            directory = ''
            with open(os.path.join(basedir, 'phase2', f), 'r') as fopen:
                for line in fopen:
                    if line.find('Evento del') >= 0:
                        pos  = line.find('Evento del')
                        date = line[pos:].split()[2].replace('-','')
                        
                    elif line.find('<strong>Latitud') >= 0:
                        pos = line.find('<strong>')
                        lat = line[pos:].split()[1]
                        pos = lat.find('<')
                        lat = lat[1:pos] + 'S'
                        
                    elif line.find('<strong>Longitud') >= 0:
                        pos = line.find('<strong>')
                        lon = line[pos:].split()[1]
                        pos = lon.find('<')
                        lon = lon[1:pos] + 'W'
    
                    elif line.find('<strong>Profundidad') >= 0:
                        pos = line.find('<strong>')
                        pro = line[pos:].split()[1]
                        pos = pro.find('<')
                        pro = pro[:pos] + 'KM'
                        
                    elif line.find('<strong>Magnitud') >= 0:
                        pos = line.find('<strong>')
                        mag = line[pos:].split()[1]
                        pos = mag.find('<')
                        mag = mag[:pos] + 'Mw'
                        directory = date + '_' + mag + '_' + lat + '_' + lon + '_' + pro
                        
                    elif line.find('/write/') >= 0:
                        pos = line.find("/write/")
                        if  pos >= 0:
                            sl = line[pos:].split("/")
                            evt = sl[2]
                            sta = sl[3][:-2]
                            url = "http://evtdb.csn.uchile.cl/write/{}/{}".format(evt, sta)
        
                            if not os.path.exists(directory):
                                os.makedirs(directory)
        
                            os.chdir(directory)
                            
                            if not os.path.exists(sta + '.zip'):
                                widget.insert('end', 'Downloading event ' + directory + ', station ' + sta + '\n')
                                window.update_idletasks()
                                try:
                                    wget.download(url, sta + '.zip')
                                except:
                                    widget.insert('end', 'Event ' + directory + ', station ' + sta + ' cannot be downloaded\n')
                                # urllib.request.urlretrieve(url, filename=sta + '.zip')
                            else:
                                widget.insert('end', 'Event ' + directory + ', station ' + sta + ' already downloaded\n')
                                window.update_idletasks()
                            os.chdir(basedir2)
        
        for filename in files:
            download_zip(filename)
        os.chdir(basedir)
        shutil.rmtree(directory)
        widget.insert('end', 'Done!\n')
        window.update_idletasks()
