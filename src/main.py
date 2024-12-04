#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 21:02:53 2019

@author: srcastro
"""
import os
import sys
import subprocess

import tkinter
from PIL import ImageTk, Image

sys.path.append(os.path.abspath(os.path.join('.', 'lib', 'pyrjmcmc')))
sys.path.append(os.path.abspath(os.path.join('.', 'lib')))

# Widgets
from widgets.updateEventsList     import updateEventsList
from widgets.downloadFromCSN      import downloadNewEvents
from widgets.transformRecords     import transformRecords
from widgets.automatic_p_wave     import automatic_p_wave, detect_p_wave
from widgets.correctRecords       import correctRecords, applyCorrection
from widgets.updateFlatFile       import updateFlatFile
from widgets.updateSpectralValues import updateSpectralValues
from widgets.updateDistances      import updateDistances
from widgets.draftToMain          import draftToMain

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
    def __init__(self, basePath = None, dataPath = None, draftPath = None):
        self.root = tkinter.Tk()
        self.basePath = basePath
        self.dataPath = dataPath
        self.draftPath = draftPath
    
    def run_tk(self):
        self.root.wm_title("SIBER-RISK Strong Motion Database")
        self.root.geometry("400x700")
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
        
        button_matfiles = tkinter.Button(master=self.root, text="Transformar registros a formato paper",
                                    command=self._transformRecords)
        button_matfiles.pack(side=tkinter.TOP, pady=10)
        
        button_pwave = tkinter.Button(master=self.root, text="Detección de onda P",
                                 command=self._pwave_detection)
        button_pwave.pack(side=tkinter.TOP, pady=10)

        button_pwave_postponed = tkinter.Button(master=self.root, text="Detección de onda P registros pospuestos",
                                 command=self._pwave_detection_postponed)
        button_pwave_postponed.pack(side=tkinter.TOP, pady=10)
        
        button_correction = tkinter.Button(master=self.root, text="Corrección de línea base",
                                 command=self._correction)
        button_correction.pack(side=tkinter.TOP, pady=10)
        
        button_flatfile = tkinter.Button(master=self.root, text="Actualizar Flat file",
                                 command=self._updateFlatFile)
        button_flatfile.pack(side=tkinter.TOP, pady=10)

        button_spectral = tkinter.Button(master=self.root, text="Actualizar espectros",
                                 command=self._updateSpectralValues)
        button_spectral.pack(side=tkinter.TOP, pady=10)
        
        button_distances = tkinter.Button(master=self.root, text="Actualizar distancias",
                                 command=self._updateDistances)
        button_distances.pack(side=tkinter.TOP, pady=10)

        button_draft2main = tkinter.Button(master=self.root, text="Traspasar borrador a producción",
                                 command=self._draftToMain)
        button_draft2main.pack(side=tkinter.TOP, pady=10)
        
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
        updateEventsList(window, text, self.basePath, self.dataPath, self.draftPath)

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
        downloadNewEvents(window, text, self.basePath, self.dataPath, self.draftPath)
    
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
        transformRecords(window, text, self.basePath, self.dataPath, self.draftPath)
    
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
            detect_p_wave(window, self.basePath, self.dataPath, self.draftPath)

        button_process = tkinter.Button(master=window, text="Iniciar detección onda P",
                                  command=_process)
        button_process.pack(side=tkinter.RIGHT)

        window.update_idletasks()
        disable = automatic_p_wave(window, text, self.basePath, self.dataPath, self.draftPath)
        if disable:
            button_process['state'] = 'disabled'
    
    def _pwave_detection_postponed(self):
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
            detect_p_wave(window, self.basePath, self.dataPath, self.draftPath)

        button_process = tkinter.Button(master=window, text="Iniciar detección onda P registros pospuestos",
                                  command=_process)
        button_process.pack(side=tkinter.RIGHT)

        window.update_idletasks()
        disable = automatic_p_wave(window, text, self.basePath, self.dataPath, self.draftPath, str)
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
            button_quit['state'] = 'disabled'
            applyCorrection(window, text, self.basePath, self.dataPath, self.draftPath, False)
            button_quit['state'] = 'normal'

        def _parallel_process():
            button_serial_process['state']   = 'disabled'
            button_parallel_process['state'] = 'disabled'
            button_quit['state'] = 'disabled'
            applyCorrection(window, text, self.basePath, self.dataPath, self.draftPath, True)
            button_quit['state'] = 'normal'

        button_serial_process = tkinter.Button(master=window, text="Ejecución en serie",
                                  command=_serial_process)
        button_serial_process.pack(side=tkinter.RIGHT)

        button_parallel_process = tkinter.Button(master=window, text="Ejecución en paralelo",
                                  command=_parallel_process)
        button_parallel_process.pack(side=tkinter.RIGHT)

        window.update_idletasks()
        disable = correctRecords(window, text, self.basePath, self.dataPath, self.draftPath)
        
        if disable:
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
        
        updateFlatFile(window, text, self.basePath, self.dataPath, self.draftPath)
        
    def _updateSpectralValues(self):
        
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
        
        updateSpectralValues(window, text, self.basePath, self.dataPath, self.draftPath)
        
    def _updateDistances(self):
        
        window = tkinter.Toplevel(self.root)
        window.geometry("800x350")
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
        
        updateDistances(window, text, self.basePath, self.dataPath, self.draftPath)
    
    def _draftToMain(self):
        
        window = tkinter.Toplevel(self.root)
        window.geometry("800x350")
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
        
        draftToMain(window, text, self.basePath, self.dataPath, self.draftPath)

    def _quit(self):
        self.root.quit()     # stops mainloop
        self.root.destroy()  # this is necessary on Windows to prevent
                        # Fatal Python Error: PyEval_RestoreThread: NULL tstate

if __name__ == '__main__':

    basePath  = os.path.abspath('.')
    dataPath  = os.path.abspath(os.path.join('.', '..', 'data'))
    draftPath = os.path.abspath(os.path.join('.', '..', 'draft'))

    copyreg.pickle(MethodType, _pickle_method, _unpickle_method)
    tk_thread = TkThread(basePath, dataPath, draftPath)
    tk_thread.run_tk()
