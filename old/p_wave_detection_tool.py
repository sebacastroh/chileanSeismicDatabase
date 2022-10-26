# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 14:12:03 2022

@author: sebac
"""
import os
import json
import copyreg
from types import MethodType

import automatic_p_wave_publish

import tkinter
from PIL import ImageTk, Image

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


class TkThread:
    def __init__(self):
        self.root = tkinter.Tk()
    
    def run_tk(self):
        self.root.wm_title("SIBER-RISK Strong Motion Database")
        self.root.geometry("400x300")
        self.root.configure(background='white')
        
        logo = ImageTk.PhotoImage(file="logo.png")
        self.root.iconphoto(True, logo)
        
        img = ImageTk.PhotoImage(Image.open("logo_fondecyt.png"))
        
        imglabel = tkinter.Label(self.root, image=img, borderwidth=0)
        
        imglabel.pack(side = tkinter.TOP)
        
        button_pwave = tkinter.Button(master=self.root, text="P-Wave detection",
                                 command=self._pwave_detection)
        button_pwave.pack(side=tkinter.TOP)
        
        button_quit = tkinter.Button(master=self.root, text="Quit",
                                 command=self._quit)
        button_quit.pack(side=tkinter.TOP)
        
        self.root.mainloop()
        
    def _pwave_detection(self):
        window = tkinter.Toplevel(self.root)
        
        filename = 'published_p_waves.json'
        if os.path.exists(filename):
            with open(filename) as f:
                p_waves = json.load(f)
        else:
            p_waves = {}
        
        mat_files_wd = os.path.join(os.getcwd(), 'events', 'matUncorrected')
        
        automatic_p_wave_publish.main_window(window, mat_files_wd, p_waves, filename)
        
    def _quit(self):
        self.root.quit()     # stops mainloop
        self.root.destroy()  # this is necessary on Windows to prevent
                        # Fatal Python Error: PyEval_RestoreThread: NULL tstate

if __name__ == '__main__':
    copyreg.pickle(MethodType, _pickle_method, _unpickle_method)
    tk_thread = TkThread()
    tk_thread.run_tk()