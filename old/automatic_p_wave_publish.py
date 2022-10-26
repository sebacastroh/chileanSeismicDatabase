# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 15:57:27 2022

@author: sebac
"""

import os
import json
import numpy as np
import scipy.io as sio
import SeismicCorrectionLibrary

import tkinter
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure

first = True
events = {}
event_name = ''
key = ''
currentStation = 0
nStations = 0
event = None

#%%
def find_p_wave(filename):
    p_wave_path = '/home/srcastro/projects/correctSeismicDatabase/newMethodology/p_wave/'
    mat_files_wd = '/home/srcastro/projects/correctSeismicDatabase/newMethodology/events_mat_uncorrected/'
    event = sio.loadmat(mat_files_wd + filename,
                        struct_as_record=False, squeeze_me=True)
    results = []
    
    for key,value in event.items():
        output = 'It has to be checked'
        method = 'To be determined'
        difference = -1.
        if key.startswith('st'):
            fsamp = int(1./value.dt)
            X  = value.acc_3.copy()
            
            if True:
                pos = SeismicCorrectionLibrary.PWaveDetection(X, fsamp)
                
                txtfile = os.path.join(p_wave_path, filename[:-4], value.name + '.txt')
                
                if os.path.exists(txtfile):
                    manual_pos = int(np.loadtxt(txtfile))
                    difference = np.abs(manual_pos - pos)*value.dt
                    if difference < 1.:
                        output = 'Valid'
                        method = 'Automatic'
            # except:
                # print(filename, key)
                # pos = 0
            
            results.append([filename, key, value.name, pos, output, method, difference])
    
    return results
   
def pop_up():
    win = tkinter.Toplevel()
    win.wm_title("Warning")
    
    l = tkinter.Label(win, text='All records have been checked')
    l.grid(row=0, column=0)
    
    b = tkinter.Button(win, text="Close", command=win.destroy)
    b.grid(row=1, column=0)    

def main_window(window, mat_files_wd, p_waves, filename):
    global first, events, event_name, key, nStations, currentStation, event
    matfiles = sorted(os.listdir(mat_files_wd))
    
    key = 'st00'
    nStations = 0
    
    to_check = False
    for matfile in matfiles:
        if matfile[:-4] not in p_waves.keys():
            event_name = matfile[:-4]
            event = sio.loadmat(os.path.join(mat_files_wd, matfile),
                                struct_as_record=False, squeeze_me=True)
            
            nStations = len(event) - 3
            
            p = 1
            while event.get(key) is None:
                key = 'st%0.2i' %p
                p += 1
            
            record_x = event[key].acc_1/9.81
            record_y = event[key].acc_2/9.81
            record_z = event[key].acc_3/9.81
            
            component_1 = event[key].component_1
            component_2 = event[key].component_2
            component_3 = event[key].component_3
            
            dt = event['st00'].dt
            
            t = np.linspace(0., (len(record_z)-1)*dt, len(record_z))
            currentStation = 1
            to_check = True
            break
        # i[0] += 1
    
    if not to_check:
        window.wm_title("P Wave Detection")
        
        l = tkinter.Label(window, text='All records have been checked')
        l.grid(row=0, column=0)
    
        b = tkinter.Button(window, text="Close", command=window.destroy)
        b.grid(row=1, column=0)
    else:
        first = True
        window.wm_title("P Wave Detection")
        
        try:
            pos, sim = SeismicCorrectionLibrary.PWaveDetection(record_z, 1./dt,
                                  ntrials=1,
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
        a1.set_title(event['st00'].event_id + ' - %i/%i' %(currentStation, nStations))
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
    
        def update_plot(t, sim, pos, event_name, record_x, record_y, record_z, component_1, component_2, component_3):
            global first            
            first = True
            a1.lines[0].set_data((t, record_x))
            a1.lines[1].set_data((t[pos], record_x[pos]))
            
            a2.lines[0].set_data((t, record_y))
            a2.lines[1].set_data((t[pos], record_y[pos]))
            
            a3.lines[0].set_data((t, record_z))
            a3.lines[1].set_data((t[pos], record_z[pos]))
            
            a4.lines[0].set_data((t, sim[:,1]))
            a4.lines[1].set_data((t[pos], sim[pos,1]))
            
            a1.set_xlim(-t[-1]*0.1, t[-1]*1.1)
            a1.set_ylim(-np.max(np.abs(record_x))*1.1, np.max(np.abs(record_x))*1.1)
            
            a2.set_xlim(-t[-1]*0.1, t[-1]*1.1)
            a2.set_ylim(-np.max(np.abs(record_y))*1.1, np.max(np.abs(record_y))*1.1)
            
            a3.set_xlim(-t[-1]*0.1, t[-1]*1.1)
            a3.set_ylim(-np.max(np.abs(record_z))*1.1, np.max(np.abs(record_z))*1.1)
            
            a4.set_xlim(-t[-1]*0.1, t[-1]*1.1)
            a4.set_ylim(-0.1, 1.1)
            
            a1.set_label(component_1)
            a2.set_label(component_2)
            a3.set_label(component_3)
            
            a1.set_title(event_name)
            
            a1.legend([component_1])
            a2.legend([component_2])
            a3.legend([component_3])
            fig.canvas.draw()
    
        def _valid():
            global first, events, event_name, key, nStations, currentStation, event
            if events.get(event_name) is None:
                events[event_name] = {}
            
            if first:
                method = 'Automatic'
            else:
                method = 'Manual'
            
            events[event_name][key] = {
              'pos': int(this_pos[0]),
              'status': True,
              'method': method
            } 
            
            currentStation += 1
            
            if currentStation <= nStations:
                to_check = True
                p = int(key[2:])
                key = 'st%0.2i' %(p+1)
            else:
                p_waves[event_name] = events[event_name]
                to_check = False
                for matfile in matfiles:
                    if p_waves.get(matfile[:-4]) is None:
                        event_name = matfile[:-4]
                        event = sio.loadmat(os.path.join(mat_files_wd, matfile),
                                            struct_as_record=False, squeeze_me=True)
                        
                        nStations = len(event) - 3
                        currentStation = 1
                        to_check = True
                        key = 'st00'
                        break
            
                        
            if not to_check:
                pop_up()
            else:
                p = int(key[2:])
                while event.get(key) is None:
                    p += 1
                    key = 'st%0.2i' %p
                    
                record_x = event[key].acc_1/9.81
                record_y = event[key].acc_2/9.81
                record_z = event[key].acc_3/9.81
                
                component_1 = event[key].component_1
                component_2 = event[key].component_2
                component_3 = event[key].component_3
                
                dt = event[key].dt
                t = np.linspace(0., (len(record_z)-1)*dt, len(record_z))
                    
                try:
                    pos, sim = SeismicCorrectionLibrary.PWaveDetection(record_z, 1./dt,
                                          return_similarity=True)
                except:
                    pos = 0
                    sim = np.vstack((np.arange(len(t)), np.ones(len(t)))).T
                
                this_pos[0] = pos
                update_plot(t, sim, pos, event[key].event_id + ' - %i/%i' %(currentStation, nStations), record_x, record_y, record_z, component_1, component_2, component_3)
            # i[0] += k
        
        def _discarded(): 
            global first, events, event_name, key, nStations, currentStation, event
            if events.get(event_name) is None:
                events[event_name] = {}
            
            if first:
                method = 'Automatic'
            else:
                method = 'Manual'
            
            events[event_name][key] = {
              'pos': int(this_pos[0]),
              'status': False,
              'method': method
            } 
            
            currentStation += 1
            
            if currentStation <= nStations:
                to_check = True
                p = int(key[2:])
                key = 'st%0.2i' %(p+1)
            else:
                p_waves[event_name] = events[event_name]
                to_check = False
                for matfile in matfiles:
                    if p_waves.get(matfile[:-4]) is None:
                        event_name = matfile[:-4]
                        event = sio.loadmat(os.path.join(mat_files_wd, matfile),
                                            struct_as_record=False, squeeze_me=True)
                        
                        nStations = len(event) - 3
                        currentStation = 1
                        to_check = True
                        key = 'st00'
                        break
            
                        
            if not to_check:
                pop_up()
            else:
                p = int(key[2:])
                while event.get(key) is None:
                    p += 1
                    key = 'st%0.2i' %p
                    
                record_x = event[key].acc_1/9.81
                record_y = event[key].acc_2/9.81
                record_z = event[key].acc_3/9.81
                
                component_1 = event[key].component_1
                component_2 = event[key].component_2
                component_3 = event[key].component_3
                
                dt = event[key].dt
                t = np.linspace(0., (len(record_z)-1)*dt, len(record_z))
                    
                try:
                    pos, sim = SeismicCorrectionLibrary.PWaveDetection(record_z, 1./dt,
                                          return_similarity=True)
                except:
                    pos = 0
                    sim = np.vstack((np.arange(len(t)), np.ones(len(t)))).T
                
                this_pos[0] = pos
                update_plot(t, sim, pos, event[key].event_id + ' - %i/%i' %(currentStation, nStations), record_x, record_y, record_z, component_1, component_2, component_3)
        
        def _save():
            with open(filename, 'w') as f:
                json.dump(p_waves, f)
            
        def _save_and_quit():
            with open(filename, 'w') as f:
                json.dump(p_waves, f)
                    
            window.quit()     # stops mainloop
            window.destroy()
            
        def _quit():            
            window.quit()     # stops mainloop
            window.destroy()  # this is necessary on Windows to prevent
                              # Fatal Python Error: PyEval_RestoreThread: NULL tstate
    
    
        
        canvas = FigureCanvasTkAgg(fig, master=window)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
        
        toolbar = NavigationToolbar2Tk(canvas, window)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
        
        canvas.mpl_connect('button_press_event', onclick)
        
        button = tkinter.Button(master=window, text="Accept",
                                command=_valid)
        button.pack(side=tkinter.LEFT)
        
        button2 = tkinter.Button(master=window, text="Discard",
                                 command=_discarded)
        button2.pack(side=tkinter.LEFT)
        
        button4 = tkinter.Button(master=window, text="Save",
                                 command=_save)
        button4.pack(side=tkinter.LEFT)
        
        button5 = tkinter.Button(master=window, text="Save and Quit",
                                 command=_save_and_quit)
        button5.pack(side=tkinter.LEFT)
        
        button6 = tkinter.Button(master=window, text="Quit without saving",
                                 command=_quit)
        button6.pack(side=tkinter.LEFT)
        
        tkinter.mainloop()