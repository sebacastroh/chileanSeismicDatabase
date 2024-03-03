# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 15:57:27 2022

@author: sebac
"""
import os
import json
import numpy as np
import pandas as pd
import lib.SeismicCorrectionLibrary as SeismicCorrectionLibrary

import tkinter
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure

first    = True
action   = None
stations = []
p_waves  = {}

DEFAULT_INDENT = 2
SORT_KEYS      = True

def automatic_p_wave(window, widget, basePath=None):
    global stations, p_waves

    # Registered P-Waves
    if os.path.exists(os.path.join(basePath, 'data', 'p_waves.json')):
        with open(os.path.join(basePath, 'data', 'p_waves.json')) as f:
            p_waves = json.load(f)
    else:
        p_waves = {}

    stations = []
    for event_id, stations_info in p_waves.items():
        for station_code, station_info in stations_info.items():
            if station_info.get('status') is None:
                stations.append([event_id, station_code])

    nStations = len(stations)
    disable   = False

    if nStations > 0:
        widget.insert('end', 'Eventos con estaciones pendientes\n\n' + '\n'.join([event_id + ' - ' + station_code for event_id, station_code in stations]) + '\n\n')

    if nStations == 0:
        widget.insert('end', 'No hay eventos pendientes por procesar.')
        disable = True

    widget.see('end')
    window.update_idletasks()

    return disable

def detect_p_wave(masterWindow, basePath):
    global stations, p_waves

    action = 'continue'

    # Detect P-Wave on events with pending stations
    old_event_id = None
    for event_station in stations:
        event_id, station_code = event_station

        if not os.path.exists(os.path.join(basePath, 'data', 'seismicDatabase', 'npz', event_id + '.npz')):
            continue

        if event_id != old_event_id:
            old_event_id = event_id
            with np.load(os.path.join(basePath, 'data', 'seismicDatabase', 'npz', event_id + '.npz'), allow_pickle=True) as f:
                event = {}
                for key, value in f.items():
                    event[key] = value.item()

        event_id = event.get('event_id')
        for st, station in event.items():
            if not st.startswith('st'):
                continue

            if station.get('station_code') == station_code:
                action = plot_p_wave(masterWindow, event_id, station, basePath)
                break

        if action != 'continue':
            break

    # Stop in case user pressed quit
    if action != 'continue':
        return

    # Show end message
    if action == 'continue':
        pop_up(basePath)

def plot_p_wave(masterWindow, event_id, station, basePath):
    global first, p_waves, action

    first = True
    action = None

    window = tkinter.Toplevel()
    window.attributes('-fullscreen', True)
    window.wm_title("P Wave Detection")

    station_code = station.get('station_code')

    record_x = station.get('acc_uncorrected_1')
    record_y = station.get('acc_uncorrected_2')
    record_z = station.get('acc_uncorrected_3')
    dt       = station.get('dt')
    
    n_x        = len(record_x)
    n_y        = len(record_y)
    n_z        = len(record_z)

    n = np.max([n_x, n_y, n_z])
    t = np.linspace(0., (n-1)*dt, n)

    record_x = np.hstack((record_x, np.zeros(n-n_x)))
    record_y = np.hstack((record_y, np.zeros(n-n_y)))
    record_z = np.hstack((record_z, np.zeros(n-n_z)))

    component_1 = station.get('component_1')
    component_2 = station.get('component_2')
    component_3 = station.get('component_3')

    try:
        pos, sim = SeismicCorrectionLibrary.PWaveDetection(record_z, 1./dt,
                              ntrials=5,
                              return_similarity=True)
    except:
        pos = 0
        sim = np.vstack((np.arange(len(t)), np.ones(len(t)))).T

    this_pos = [pos]

    lw = 1
    fig = Figure(figsize=(5, 4), dpi=100)

    a1 = fig.add_subplot(411)
    a1.plot(t, record_x, label=component_1, lw=lw)
    a1.plot(t[pos], record_x[pos], 'ro', label='')
    a1.set_title(event_id + ' - ' + station_code)
    a1.legend()
    a1.grid()
    a1.set_xlim(-t[-1]*0.1, t[-1]*1.1)
    a1.set_ylim(-np.max(np.abs(record_x))*1.1, np.max(np.abs(record_x))*1.1)

    a2 = fig.add_subplot(412)
    a2.plot(t, record_y, label=component_2, lw=lw)
    a2.plot(t[pos], record_y[pos], 'ro', label='')
    a2.legend()
    a2.grid()
    a2.set_xlim(-t[-1]*0.1, t[-1]*1.1)
    a2.set_ylim(-np.max(np.abs(record_y))*1.1, np.max(np.abs(record_y))*1.1)

    a3 = fig.add_subplot(413)
    a3.plot(t, record_z, label=component_3, lw=lw)
    a3.plot(t[pos], record_z[pos], 'ro', label='')
    a3.legend()
    a3.grid()
    a3.set_xlim(-t[-1]*0.1, t[-1]*1.1)
    a3.set_ylim(-np.max(np.abs(record_z))*1.1, np.max(np.abs(record_z))*1.1)

    a4 = fig.add_subplot(414)
    a4.plot(t, sim[:,1], lw=lw)
    a4.plot(t[pos], sim[pos,1], 'ro')
    a4.set_xlim(-t[-1]*0.1, t[-1]*1.1)
    a4.set_ylim(-0.1, 1.1)
    a4.grid()

    a1.set_picker(True)
    a2.set_picker(True)
    a3.set_picker(True)
    a4.set_picker(True)

    def onclick(event):
        global first
        if event.dblclick:
            first = False
            if event.inaxes == a1:
                xdata = a1.lines[0].get_xdata()
                ydata = a1.lines[0].get_ydata()
            elif event.inaxes == a2:
                xdata = a2.lines[0].get_xdata()
                ydata = a2.lines[0].get_ydata()
            elif event.inaxes == a3:
                xdata = a3.lines[0].get_xdata()
                ydata = a3.lines[0].get_ydata()
            elif event.inaxes == a4:
                xdata = a4.lines[0].get_xdata()
                ydata = a4.lines[0].get_ydata()
            
            pos = np.argmin(np.sqrt((xdata-event.xdata)**2 + (ydata-event.ydata)**2))
            this_pos[0] = pos
            
            xdata1 = a1.lines[0].get_xdata()
            ydata1 = a1.lines[0].get_ydata()
            a1.lines[1].set_data((xdata1[pos], ydata1[pos]))
            
            xdata2 = a2.lines[0].get_xdata()
            ydata2 = a2.lines[0].get_ydata()
            a2.lines[1].set_data((xdata2[pos], ydata2[pos]))
            
            xdata3 = a3.lines[0].get_xdata()
            ydata3 = a3.lines[0].get_ydata()
            a3.lines[1].set_data((xdata3[pos], ydata3[pos]))
            
            xdata4 = a4.lines[0].get_xdata()
            ydata4 = a4.lines[0].get_ydata()
            a4.lines[1].set_data((xdata4[pos], ydata4[pos]))
            fig.canvas.draw()

    canvas = FigureCanvasTkAgg(fig, master=window)  # A tk.DrawingArea.
    canvas.draw()
    canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
    
    toolbar = NavigationToolbar2Tk(canvas, window)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
    
    canvas.mpl_connect('button_press_event', onclick)

    def _valid():
        global first, p_waves, action
        
        if p_waves.get(event_id) is None:
            p_waves[event_id] = {}

        if first:
            method = 'Automatic'
        else:
            method = 'Manual'

        p_waves[event_id][station_code]['pos']    = int(this_pos[0])
        p_waves[event_id][station_code]['status'] = True
        p_waves[event_id][station_code]['method'] = method

        action = 'continue'
        window.destroy()

    def _discarded():
        global first, p_waves, action
        
        if p_waves.get(event_id) is None:
            p_waves[event_id] = {}

        if first:
            method = 'Automatic'
        else:
            method = 'Manual'

        p_waves[event_id][station_code]['pos']    = int(this_pos[0])
        p_waves[event_id][station_code]['status'] = False
        p_waves[event_id][station_code]['method'] = method

        action = 'continue'
        window.destroy()

    def _omit():
        global action
        action = 'continue'
        window.destroy()

    def _save():
        global action
        with open(os.path.join(basePath, 'data', 'p_waves.json'), 'w') as f:
            json.dump(p_waves, f, indent=DEFAULT_INDENT, sort_keys=SORT_KEYS)
        action = 'continue'
        
    def _save_and_quit():
        global action
        with open(os.path.join(basePath, 'data', 'p_waves.json'), 'w') as f:
            json.dump(p_waves, f, indent=DEFAULT_INDENT, sort_keys=SORT_KEYS)
        action = 'quit'
        window.destroy()

    def _quit():
        global action
        action = 'quit'
        window.destroy()

    button = tkinter.Button(master=window, text="Registrar onda P",
                            command=_valid)
    button.pack(side=tkinter.LEFT)
    
    button2 = tkinter.Button(master=window, text="Descartar registro",
                             command=_discarded)
    button2.pack(side=tkinter.LEFT)

    button3 = tkinter.Button(master=window, text="Dejar para después",
                             command=_omit)
    button3.pack(side=tkinter.LEFT)
    
    button4 = tkinter.Button(master=window, text="Guardar",
                             command=_save)
    button4.pack(side=tkinter.LEFT)
    
    button5 = tkinter.Button(master=window, text="Guardar y salir",
                             command=_save_and_quit)
    button5.pack(side=tkinter.LEFT)
    
    button6 = tkinter.Button(master=window, text="Salir sin guardar",
                             command=_quit)
    button6.pack(side=tkinter.LEFT)

    masterWindow.wait_window(window)

    return action

def pop_up(basePath):
    global p_waves

    window = tkinter.Toplevel()
    window.wm_title("Atención")

    def _save():
        global p_waves
        with open(os.path.join('data', 'p_waves.json'), 'w') as f:
            json.dump(p_waves, f, indent=DEFAULT_INDENT, sort_keys=SORT_KEYS)
        window.destroy()

    l = tkinter.Label(window, text='Ya se han revisado todos los registros')
    l.grid(row=0, column=0)

    b1 = tkinter.Button(window, text='Cerrar sin guardar', command=window.destroy)
    b1.grid(row=1, column=1)

    b2 = tkinter.Button(window, text='Guardar', command=_save)
    b2.grid(row=1, column=0)
