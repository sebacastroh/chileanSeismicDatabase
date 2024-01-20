#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 21:02:53 2019

@author: srcastro
"""
import os
import sys
import json
import pickle
import numpy as np
import pandas as pd
import scipy.io as spio

import pyperclip
import webbrowser
import subprocess
import multiprocessing

import tkinter
from PIL import ImageTk, Image

sys.path.append(os.path.abspath(os.path.join('.', 'lib', 'pyrjmcmc')))
sys.path.append(os.path.abspath(os.path.join('.', 'lib')))

# Widgets
from widgets.updateEventsList  import updateEventsList
from widgets.downloadNewEvents import downloadNewEvents
from widgets.transformRecords  import transformRecords
from widgets.automatic_p_wave  import automatic_p_wave, detect_p_wave
from widgets.correctRecords    import correctRecords, applyCorrection
from widgets.updateFlatFile   import updateFlatFile

import copyreg
from types import MethodType

version = 2

def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)
    
def run_correction(inp):    
    cwd = os.getcwd()
    cmd = 'python ' + os.path.join(cwd, 'correctRecords.py ')
    subprocess.call(cmd + inp[0] + ' ' + inp[3] + ' True', shell=True)
    with open(inp[0] + '.txt', 'r') as f:
        result = f.readline()
    os.remove(inp[0] + '.txt')    
    output = 'Status: ' + result
        
    return output

class TkThread:
    def __init__(self, basePath = None):
        self.root = tkinter.Tk()
        self.basePath = basePath
    
    def run_tk(self):
        self.root.wm_title("SIBER-RISK Strong Motion Database")
        self.root.geometry("400x600")
        self.root.configure(background='white')
        
        logo = ImageTk.PhotoImage(file="assets/logo.png")
        self.root.iconphoto(True, logo)
        
        img = ImageTk.PhotoImage(Image.open("assets/logo_fondecyt.png"))
        
        imglabel = tkinter.Label(self.root, image=img, borderwidth=0)
        
        imglabel.pack(side = tkinter.TOP)
        
        button_update = tkinter.Button(master=self.root, text="Actualizar lista de eventos",
                                    command=self._updateEvents)
        button_update.pack(side=tkinter.TOP, pady=10)
        
        button_download = tkinter.Button(master=self.root, text="Descargar nuevos registros",
                                    command=self._download)
        button_download.pack(side=tkinter.TOP, pady=10)
        
        # button_properties = tkinter.Button(master=self.root, text="Plane fault properties",
        #                             command=self._planeFaultProperties)
        # button_properties.pack(side=tkinter.TOP)
        
        button_matfiles = tkinter.Button(master=self.root, text="Transformar registros a formato paper",
                                    command=self._transformRecords)
        button_matfiles.pack(side=tkinter.TOP, pady=10)
        
        button_pwave = tkinter.Button(master=self.root, text="Detección de onda P",
                                 command=self._pwave_detection)
        button_pwave.pack(side=tkinter.TOP, pady=10)
        
        button_correction = tkinter.Button(master=self.root, text="Corrección de línea base",
                                 command=self._correction)
        button_correction.pack(side=tkinter.TOP, pady=10)
        
        button_flatfile = tkinter.Button(master=self.root, text="Actualizar Flat file",
                                 command=self._updateFlatFile)
        button_flatfile.pack(side=tkinter.TOP, pady=10)
        
        button_quit = tkinter.Button(master=self.root, text="Salir",
                                 command=self._quit)
        button_quit.pack(side=tkinter.TOP, pady=10)
        
        self.root.mainloop()

    def _updateEvents(self):
        window = tkinter.Toplevel(self.root)
        toolbar = tkinter.Frame(window)
        toolbar.pack(side="top", fill="x")
        text = tkinter.Text(toolbar, wrap="word")
        text.pack(side="top", fill="both", expand=True)
        text.tag_configure("stderr", foreground="#b22222")
        
        def _close():
            window.destroy()
        
        button_quit = tkinter.Button(master=window, text="Close",
                                  command=_close)
        button_quit.pack(side=tkinter.BOTTOM)
        window.update_idletasks()
        updateEventsList(window, text, self.basePath)

    def _download(self):
        window = tkinter.Toplevel(self.root)
        toolbar = tkinter.Frame(window)
        toolbar.pack(side="top", fill="x")
        text = tkinter.Text(toolbar, wrap="word")
        text.pack(side="top", fill="both", expand=True)
        text.tag_configure("stderr", foreground="#b22222")
        
        def _close():
            window.destroy()
        
        button_quit = tkinter.Button(master=window, text="Close",
                                  command=_close)
        button_quit.pack(side=tkinter.BOTTOM)
        window.update_idletasks()
        downloadNewEvents(window, text, self.basePath)
        
    def _planeFaultProperties(self):
        
        fault_plane_properties = 'fault_plane_properties_2.pkl'
        with open(fault_plane_properties, 'rb') as f:
            data = pickle.load(f)
        
        filenames = sorted(os.listdir(os.path.join(os.getcwd(), 'events')))
        these_filenames = []
        for filename in filenames:
            p1 = filename.index('_') + 1
            p2 = filename.index('Mw')
            magnitude = float(filename[p1:p2])
            if magnitude >= 7.1 and filename not in data.keys():
                these_filenames.append(filename)
        
        if len(these_filenames) > 0:
            
            window = tkinter.Toplevel(self.root)
            title = tkinter.Label(window, text="Event:")
            event = tkinter.Label(window, text= these_filenames[0])
            title.grid(row=0, column=0)
            event.grid(row=0, column=1)
            tkinter.Label(window, text="Strike").grid(row=1)
            tkinter.Label(window, text="Dip").grid(row=2)
            tkinter.Label(window, text="Rake").grid(row=3)
            
            e1 = tkinter.Entry(window)
            e2 = tkinter.Entry(window)
            e3 = tkinter.Entry(window)
            
            e1.grid(row=1, column=1)
            e2.grid(row=2, column=1)
            e3.grid(row=3, column=1)
            
            e1.delete(0,tkinter.END)
            e1.insert(0,'Insert value')
            
            e2.delete(0,tkinter.END)
            e2.insert(0,'Insert value')
            
            e3.delete(0,tkinter.END)
            e3.insert(0,'Insert value')
            i = [0]
            
            def _get_suggested_link():
                base_link = "https://www.globalcmt.org/cgi-bin/globalcmt-cgi-bin/CMT5/form?itype=ymd&yr=YYYY&mo=MM&day=DD&oyr=1976&omo=1&oday=1&jyr=1976&jday=1&ojyr=1976&ojday=1&otype=nd&nday=1&lmw=MWmin&umw=MWmax&lms=0&ums=10&lmb=0&umb=10&llat=-90&ulat=90&llon=-180&ulon=180&lhd=0&uhd=1000&lts=-9999&uts=9999&lpe1=0&upe1=90&lpe2=0&upe2=90&list=0"
                filename = str(event.cget('text'))
                year = filename[:4]
                month = filename[4:6]
                day = filename[6:8]
                p1 = filename.index('_') + 1
                p2 = filename.index('Mw')
                magnitude = float(filename[p1:p2])
                MWmin = '%0.1f' %(magnitude - 0.5)
                MWmax = '%0.1f' %(magnitude + 0.5)
            
                suggested_link = base_link.replace('YYYY', year).replace('MM', month).replace('DD', day).replace('MWmin', MWmin).replace('MWmax', MWmax)
                
                webbrowser.open_new_tab(suggested_link)
                
            def _copy_suggested_link():
                base_link = "https://www.globalcmt.org/cgi-bin/globalcmt-cgi-bin/CMT5/form?itype=ymd&yr=YYYY&mo=MM&day=DD&oyr=1976&omo=1&oday=1&jyr=1976&jday=1&ojyr=1976&ojday=1&otype=nd&nday=1&lmw=MWmin&umw=MWmax&lms=0&ums=10&lmb=0&umb=10&llat=-90&ulat=90&llon=-180&ulon=180&lhd=0&uhd=1000&lts=-9999&uts=9999&lpe1=0&upe1=90&lpe2=0&upe2=90&list=0"
                filename = str(event.cget('text'))
                year = filename[:4]
                month = filename[4:6]
                day = filename[6:8]
                p1 = filename.index('_') + 1
                p2 = filename.index('Mw')
                magnitude = float(filename[p1:p2])
                MWmin = '%0.1f' %(magnitude - 0.5)
                MWmax = '%0.1f' %(magnitude + 0.5)
            
                suggested_link = base_link.replace('YYYY', year).replace('MM', month).replace('DD', day).replace('MWmin', MWmin).replace('MWmax', MWmax)
                
                pyperclip.copy(suggested_link)
                
            
            def _nextPlaneFault():
                strike = e1.get()
                dip = e2.get()
                rake = e3.get()
                if strike != 'Insert value' and dip != 'Insert value' and rake != 'Insert value':
                    strike = float(strike)
                    dip = float(dip)
                    rake = float(rake)
                    
                    e1.delete(0,tkinter.END)
                    e1.insert(0,'Insert value')
                    
                    e2.delete(0,tkinter.END)
                    e2.insert(0,'Insert value')
                    
                    e3.delete(0,tkinter.END)
                    e3.insert(0,'Insert value')
                    
                    data[these_filenames[i[0]]] = [strike, dip, rake]
                
                if i[0] == len(these_filenames)-1:
                    pass
                else:
                    i[0] += 1
                    event['text'] = these_filenames[i[0]]
            
            def _close():
                with open('fault_plane_properties_2.pkl', 'wb') as f:
                    pickle.dump(data, f)
                window.destroy()
            
            button_link = tkinter.Button(master=window, text='Open link',
                                            command=_get_suggested_link)
            button_link.grid(row=4, column=0, sticky=tkinter.W, pady=4)
            
            button_copylink = tkinter.Button(master=window, text='Copy link',
                                            command=_copy_suggested_link)
            button_copylink.grid(row=4, column=1, sticky=tkinter.W, pady=4)
            
            button_newProp = tkinter.Button(master=window, text='Save and continue',
                                            command=_nextPlaneFault)
            button_newProp.grid(row=4, column=2, sticky=tkinter.W, pady=4)
            
            button_quit = tkinter.Button(master=window, text="Close",
                                     command=_close)
            button_quit.grid(row=4, column=3, sticky=tkinter.W, pady=4)
            window.update_idletasks()
    
    def _transformRecords(self):        
        window = tkinter.Toplevel(self.root)
        toolbar = tkinter.Frame(window)
        toolbar.pack(side="top", fill="x")
        text = tkinter.Text(toolbar, wrap="word")
        text.pack(side="top", fill="both", expand=True)
        text.tag_configure("stderr", foreground="#b22222")
        
        def _close():
            window.destroy()
                           
        button_quit = tkinter.Button(master=window, text="Close",
                                  command=_close)
        button_quit.pack(side=tkinter.BOTTOM)
        window.update_idletasks()
        transformRecords(window, text, self.basePath)
    
    def _pwave_detection(self):
        window = tkinter.Toplevel(self.root)
        toolbar = tkinter.Frame(window)
        toolbar.pack(side="top", fill="x")
        text = tkinter.Text(toolbar, wrap="word")
        text.pack(side="top", fill="both", expand=True)
        text.tag_configure("stderr", foreground="#b22222")

        def _close():
            window.destroy()
                           
        button_quit = tkinter.Button(master=window, text="Cerrar",
                                  command=_close)
        button_quit.pack(side=tkinter.RIGHT)

        def _process():
            button_process['state'] = 'disabled'
            detect_p_wave(window, self.basePath)

        button_process = tkinter.Button(master=window, text="Iniciar detección onda P",
                                  command=_process)
        button_process.pack(side=tkinter.RIGHT)

        window.update_idletasks()
        disable = automatic_p_wave(window, text, self.basePath)
        if disable:
            button_process['state'] = 'disabled'
    
    def _correction(self):
        window = tkinter.Toplevel(self.root)
        
        toolbar = tkinter.Frame(window)
        toolbar.pack(side="top", fill="x")
        text = tkinter.Text(toolbar, wrap="word")
        text.pack(side="top", fill="both", expand=True)
        text.tag_configure("stderr", foreground="#b22222")
        
        def _close():
            window.destroy()
        
        button_quit = tkinter.Button(master=window, text="Salir",
                                 command=_close)
        button_quit.pack(side=tkinter.RIGHT)
        window.update_idletasks()

        def _serial_process():
            button_serial_process['state']   = 'disabled'
            button_parallel_process['state'] = 'disabled'
            applyCorrection(window, text, self.basePath, False)

        def _parallel_process():
            button_serial_process['state']   = 'disabled'
            button_parallel_process['state'] = 'disabled'
            applyCorrection(window, text, self.basePath, True)

        button_serial_process = tkinter.Button(master=window, text="Ejecución en serie",
                                  command=_serial_process)
        button_serial_process.pack(side=tkinter.RIGHT)

        button_parallel_process = tkinter.Button(master=window, text="Ejecución en paralelo",
                                  command=_parallel_process)
        button_parallel_process.pack(side=tkinter.RIGHT)

        window.update_idletasks()
        enable = correctRecords(window, text, self.basePath)
        
        if not enable:
            button_serial_process['state']   = 'disabled'
            button_parallel_process['state'] = 'disabled'
        
    def _updateFlatFile(self):
        
        window = tkinter.Toplevel(self.root)
        toolbar = tkinter.Frame(window)
        toolbar.pack(side="top", fill="x")
        text = tkinter.Text(toolbar, wrap="word")
        text.pack(side="top", fill="both", expand=True)
        text.tag_configure("stderr", foreground="#b22222")
        
        def _close():
            window.destroy()
        
        button_quit = tkinter.Button(master=window, text="Close",
                                 command=_close)
        button_quit.pack(side=tkinter.BOTTOM)
        
        updateFlatFile(window, text, self.basePath)
        
    
    def _quit(self):
        self.root.quit()     # stops mainloop
        self.root.destroy()  # this is necessary on Windows to prevent
                        # Fatal Python Error: PyEval_RestoreThread: NULL tstate

if __name__ == '__main__':

    if not os.path.exists('data'):
        os.mkdir('data')

    basePath = os.path.abspath('.')

    copyreg.pickle(MethodType, _pickle_method, _unpickle_method)
    tk_thread = TkThread(basePath)
    tk_thread.run_tk()
