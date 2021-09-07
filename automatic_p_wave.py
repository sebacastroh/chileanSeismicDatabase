#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 10:17:12 2019

@author: srcastro
"""
import os
import pickle
import numpy as np
import scipy.io as sio
import SeismicCorrectionLibrary

import tkinter
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure

first = True

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

def main_window(window, results, mat_files_wd, name_pkl):
    global first
    filenames = sorted(os.listdir(mat_files_wd))
    
    unique_events = np.unique(np.array(results)[:,0])
#    last_event_considered = unique_events[-1]
#    position = filenames.index(last_event_considered) + 1
    
    to_consider = []
#    for filename in filenames[position:]:
#        to_consider.append(filename)
        
    for filename in filenames:
        if filename not in unique_events:
            to_consider.append(filename)
    
    for filename in to_consider:
        this_results = find_p_wave(filename)
        results.extend(this_results)
    
    i = [-1]
    to_check = False
    for result in results:
        if result[4] == 'It has to be checked':
            event = sio.loadmat(os.path.join(mat_files_wd, result[0]),
                                struct_as_record=False, squeeze_me=True)
            
            record_x = event[result[1]].acc_1/9.81
            record_y = event[result[1]].acc_2/9.81
            record_z = event[result[1]].acc_3/9.81
            
#            station_name = event[result[1]].name
            
            component_1 = event[result[1]].component_1
            component_2 = event[result[1]].component_2
            component_3 = event[result[1]].component_3
            
            dt = event[result[1]].dt
            pos = result[3]
            
            t = np.linspace(0., (len(record_z)-1)*dt, len(record_z))
            i[0] += 1
            to_check = True
            break
        i[0] += 1
    
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
        a1.set_title(event[result[1]].event_name)
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
            first = False
            
            if event.dblclick:
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
            global first
            results[i[0]][3] = this_pos[0]
            results[i[0]][4] = 'Valid'
            
            if first:
                results[i[0]][5] = 'Automatic'
            else:
                results[i[0]][5] = 'Manual'
            
            i[0] += 1
            k = -1
            
            if len(results[i[0]:]) == 0:
                pop_up()
            
            for result in results[i[0]:]:
                k += 1
                if result[4] == 'It has to be checked':
                    event = sio.loadmat(os.path.join(mat_files_wd, result[0]),
                                        struct_as_record=False, squeeze_me=True)
                    record_x = event[result[1]].acc_1/9.81
                    record_y = event[result[1]].acc_2/9.81
                    record_z = event[result[1]].acc_3/9.81
                    
#                    station_name = event[result[1]].name
                    
                    component_1 = event[result[1]].component_1
                    component_2 = event[result[1]].component_2
                    component_3 = event[result[1]].component_3
                    
                    dt = event[result[1]].dt
                    t = np.linspace(0., (len(record_z)-1)*dt, len(record_z))
                    
                    try:
                        pos, sim = SeismicCorrectionLibrary.PWaveDetection(record_z, 1./dt,
                                              return_similarity=True)
                    except:
                        pos = 0
                        sim = np.vstack((np.arange(len(t)), np.ones(len(t)))).T
                    
                    this_pos[0] = pos
                    update_plot(t, sim, pos, event[result[1]].event_name, record_x, record_y, record_z, component_1, component_2, component_3)
                    break
            i[0] += k
        
        def _discarded():            
            results[i[0]][3] = -1
            results[i[0]][4] = 'Discarded'
            results[i[0]][5] = 'Not considered'
            
            i[0] += 1
            k = -1
            
            if len(results[i[0]:]) == 0:
                pop_up()
            
            for result in results[i[0]:]:
                k += 1
                if result[4] == 'It has to be checked':
                    event = sio.loadmat(os.path.join(mat_files_wd, result[0]),
                                        struct_as_record=False, squeeze_me=True)
                    
                    record_x = event[result[1]].acc_1/9.81
                    record_y = event[result[1]].acc_2/9.81
                    record_z = event[result[1]].acc_3/9.81
                    
#                    station_name = event[result[1]].name
                    
                    component_1 = event[result[1]].component_1
                    component_2 = event[result[1]].component_2
                    component_3 = event[result[1]].component_3
                    
                    dt = event[result[1]].dt
                    t = np.linspace(0., (len(record_z)-1)*dt, len(record_z))
                    
                    try:
                        pos, sim = SeismicCorrectionLibrary.PWaveDetection(record_z, 1./dt,
                                              return_similarity=True)
                    except:
                        pos = 0
                        sim = np.vstack((np.arange(len(t)), np.ones(len(t)))).T
                    
                    this_pos[0] = pos
                    update_plot(t, sim, pos, event[result[1]].event_name, record_x, record_y, record_z, component_1, component_2, component_3)
                    break
            i[0] += k
        
        def _check_later():
            i[0] += 1
            k = -1
            
            if len(results[i[0]:]) == 0:
                pop_up()
            
            for result in results[i[0]:]:
                k += 1
                if result[4] == 'It has to be checked':
                    event = sio.loadmat(os.path.join(mat_files_wd, result[0]),
                                        struct_as_record=False, squeeze_me=True)
                    
                    record_x = event[result[1]].acc_1/9.81
                    record_y = event[result[1]].acc_2/9.81
                    record_z = event[result[1]].acc_3/9.81
                    
#                    station_name = event[result[1]].name
                    
                    component_1 = event[result[1]].component_1
                    component_2 = event[result[1]].component_2
                    component_3 = event[result[1]].component_3
                    
                    dt = event[result[1]].dt
                    t = np.linspace(0., (len(record_z)-1)*dt, len(record_z))
                    
                    try:
                        pos, sim = SeismicCorrectionLibrary.PWaveDetection(record_z, 1./dt,
                                              return_similarity=True)
                    except:
                        pos = 0
                        sim = np.vstack((np.arange(len(t)), np.ones(len(t)))).T
                    
                    this_pos[0] = pos
                    update_plot(t, sim, pos, event[result[1]].event_name, record_x, record_y, record_z, component_1, component_2, component_3)
                    break
            i[0] += k
        
        def _save():            
            with open(name_pkl, 'wb') as f:
                pickle.dump(results, f)  
            
        def _save_and_quit():
            with open(name_pkl, 'wb') as f:
                pickle.dump(results, f)
                    
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
        
        button3 = tkinter.Button(master=window, text="Check later",
                                 command=_check_later)
        button3.pack(side=tkinter.LEFT)
        
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