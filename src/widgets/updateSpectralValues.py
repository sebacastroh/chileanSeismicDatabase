# -*- coding: utf-8 -*-
"""
Created on Tus Mar 26 11:45:34 2024

@author: sebac
"""
import os
import seismic
import numpy as np
import pandas as pd

def updateSpectralValues(window, widget, basePath, dataPath, drafPath):

    if not os.path.exists(os.path.join(dataPath, 'seismicDatabase', 'npz')):
        widget.insert('end', 'No existen registros almacenados para registrar en la base de datos.\n')
        window.update_idletasks()
        return

    if not os.path.exists(os.path.join(basePath, 'data', 'p_waves.json')):
        widget.insert('end', 'No se ha encontrado el archivo p_waves.json necesario.\n')
        window.update_idletasks()
        return

    if not os.path.exists(os.path.join(dataPath, 'flatFile.csv')):
        widget.insert('end', 'No se ha encontrado el archivo flatFile.csv necesario.\n')
        window.update_idletasks()
        return

    Tn = np.array([
        0.   , 0.01 , 0.02 , 0.022, 0.025,  0.029, 0.03 , 0.032, 0.035, 0.036,\
        0.04 , 0.042, 0.044, 0.045, 0.046,  0.048, 0.05 , 0.055, 0.06 , 0.065,\
        0.067, 0.07 , 0.075, 0.08 , 0.085,  0.09 , 0.095, 0.1  , 0.11 , 0.12 ,\
        0.13 , 0.133, 0.14 , 0.15 , 0.16 ,  0.17 , 0.18 , 0.19 , 0.2  , 0.22 ,\
        0.24 , 0.25 , 0.26 , 0.28 , 0.29 ,  0.3  , 0.32 , 0.34 , 0.35 , 0.36 ,\
        0.38 , 0.4  , 0.42 , 0.44 , 0.45 ,  0.46 , 0.48 , 0.5  , 0.55 , 0.6  ,\
        0.65 , 0.667, 0.7  , 0.75 , 0.8  ,  0.85 , 0.9  , 0.95 , 1.   , 1.1  ,\
        1.2  , 1.3  , 1.4  , 1.5  , 1.6  ,  1.7  , 1.8  , 1.9  , 2.   , 2.2  ,\
        2.4  , 2.5  , 2.6  , 2.8  , 3.   ,  3.2  , 3.4  , 3.5  , 3.6  , 3.8  ,\
        4.   , 4.2  , 4.4  , 4.6  , 4.8  ,  5.   , 5.5  , 6.   , 6.5  , 7.   ,\
        7.5  , 8.   , 8.5  , 9.   , 9.5  , 10.
    ])
    n  = len(Tn)

    xis = [0.02, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2, 0.3, 0.5]
    g   = 9.81

    spectrum_names = ['component_1', 'component_2', 'component_3',
        'geometric_mean', 'rotd0', 'rotd50', 'rotd100']

    columns = ['Earthquake Name', 'Station code'] + Tn.tolist()

    if os.path.exists(os.path.join(draftPath, 'flatFile.csv')):
        flatfile = pd.read_csv(os.path.join(draftPath, 'flatFile.csv'), parse_dates=['Earthquake date', 'Start time record', 'Last update'])
    else:
        flatfile = pd.read_csv(os.path.join(dataPath, 'flatFile.csv'), parse_dates=['Earthquake date', 'Start time record', 'Last update'])
    flatfile = flatfile[flatfile['Corrected records']]

    if not os.path.exists(os.path.join(draftPath, 'spectralValues')):
        os.mkdir(os.path.join(draftPath, 'spectralValues'))

    if not os.path.exists(os.path.join(dataPath, 'spectralValues')):
        os.mkdir(os.path.join(dataPath, 'spectralValues'))

    if os.path.exists(os.path.join(draftPath, 'spectralValues', 'computed.xlsx')):
        computed = pd.read_excel(os.path.join(draftPath, 'spectralValues', 'computed.xlsx'))
    elif os.path.exists(os.path.join(dataPath, 'spectralValues', 'computed.xlsx')):
        computed = pd.read_excel(os.path.join(dataPath, 'spectralValues', 'computed.xlsx'))
    else:
        computed = pd.DataFrame([], columns=['Earthquake Name', 'Station code', 'Component 1', 'Component 2', 'Component 3', 'Last update'])

    pending = computed.merge(flatfile, how='right', on=['Earthquake Name', 'Station code'])
    pending = pending[pending['Last update_x'] != pending['Last update_y']][['Earthquake Name', 'Station code']]
    pending.drop_duplicates(keep=False, inplace=True)

    if len(pending) == 0:
        widget.insert('end', 'No hay espectros nuevos que calcular.\n')
        window.update_idletasks()
        return

    new_rows = []

    indices = pending.merge(computed.reset_index(), how='inner', on=['Earthquake Name', 'Station code'])['index'].tolist()

    for k, xi in enumerate(xis):
        widget.insert('end', 'Calculando espectro para xi = %0.2f.\n' %xi)
        widget.see('end')
        window.update_idletasks()

        foldername = os.path.join(dataPath, 'spectralValues', 'xi_%0.2f' %xi)
        if not os.path.exists(foldername):
            os.mkdir(foldername)

        foldernameDraft = os.path.join(draftPath, 'spectralValues', 'xi_%0.2f' %xi)
        if not os.path.exists(foldernameDraft):
            os.mkdir(foldernameDraft)

        spectral_values = [[] for i in range(7)]

        old_event_id = None
        for r, row in pending.iterrows():
            event_id     = row['Earthquake Name']
            station_code = row['Station code']

            if event_id != old_event_id:
                if os.path.exists(os.path.join(draftPath, 'seismicDatabase', 'npz', event_id + '.npz')):
                    filename = os.path.join(draftPath, 'seismicDatabase', 'npz', event_id + '.npz')
                else:
                    filename = os.path.join(dataPath, 'seismicDatabase', 'npz', event_id + '.npz')

                with np.load(filename, allow_pickle=True) as f:
                    data = {}
                    for key, value in f.items():
                        data[key] = value.item()

                old_event_id = event_id

            for st, station in data.items():
                if not st.startswith('st'):
                    continue

                if station['station_code'] != station_code:
                    continue

                if k == 0:
                    new_rows.append([event_id, station_code, station['component_1'], station['component_2'], station['component_3'], station['last_update']])

                component_1 = station['acc_filtered_1']/g
                component_2 = station['acc_filtered_2']/g
                component_3 = station['acc_filtered_3']/g

                compute_1   = False
                compute_2   = False
                compute_3   = False
                compute_rot = False

                if len(component_1) > 0:
                    compute_1 = True

                if len(component_2) > 0:
                    compute_2 = True

                if len(component_3) > 0:
                    compute_3 = True

                if compute_1 and compute_2:
                    compute_rot = True

                dt = station['dt']

                if compute_rot:
                    spectra       = seismic.SpectraRot(component_1, component_2, dt, Tn, xi, 181)
                    spectrum_1    = spectra[0]
                    spectrum_2    = spectra[90]
                    spectrum_g    = np.sqrt(spectra[0]*spectra[90])
                    spectrum_r0   = np.min(spectra, axis=0)
                    spectrum_r50  = np.median(spectra, axis=0)
                    spectrum_r100 = np.max(spectra, axis=0)
                else:
                    spectrum_1    = np.array([np.nan for i in range(n)])
                    spectrum_2    = np.array([np.nan for i in range(n)])
                    spectrum_g    = np.array([np.nan for i in range(n)])
                    spectrum_r0   = np.array([np.nan for i in range(n)])
                    spectrum_r50  = np.array([np.nan for i in range(n)])
                    spectrum_r100 = np.array([np.nan for i in range(n)])
                    if compute_1:
                        spectrum_1 = seismic.Spectrum(component_1, dt, Tn, xi)
                    if compute_2:
                        spectrum_2 = seismic.Spectrum(component_2, dt, Tn, xi)

                if compute_3:
                    spectrum_3 = seismic.Spectrum(component_3, dt, Tn, xi)
                else:
                    spectrum_3 = np.array([np.nan for i in range(n)])

                spectral_values[0].append(spectrum_1)
                spectral_values[1].append(spectrum_2)
                spectral_values[2].append(spectrum_3)
                spectral_values[3].append(spectrum_g)
                spectral_values[4].append(spectrum_r0)
                spectral_values[5].append(spectrum_r50)
                spectral_values[6].append(spectrum_r100)

        for i, spectrum_name in enumerate(spectrum_names):
            widget.insert('end', 'Guardando espectro %s.\n' %spectrum_name)
            widget.see('end')
            window.update_idletasks()

            new_spectrum_values = pd.concat([pending.reset_index(drop=True), pd.DataFrame(spectral_values[i], columns=columns[2:])], axis=1)

            if os.path.exists(os.path.join(foldernameDraft, spectrum_name + '.xlsx')):
                spectrum_values = pd.read_excel(os.path.join(foldernameDraft, spectrum_name + '.xlsx'))
            elif os.path.exists(os.path.join(foldername, spectrum_name + '.xlsx')):
                spectrum_values = pd.read_excel(os.path.join(foldername, spectrum_name + '.xlsx'))
            else:
                spectrum_values = pd.DataFrame([], columns=columns)

            filename = os.path.join(foldernameDraft, spectrum_name + '.xlsx')

            spectrum_values.drop(indices, inplace=True)
            spectrum_values = pd.concat([spectrum_values, new_spectrum_values], ignore_index=True)
            spectrum_values.sort_values(by=['Earthquake Name', 'Station code'], inplace=True)

            spectrum_values.to_excel(filename, index=False)

    widget.insert('end', 'Guardando tabla índice.\n\n')
    widget.see('end')
    window.update_idletasks()

    computed.drop(indices, inplace=True)
    computed = pd.concat([computed, pd.DataFrame(new_rows, columns=['Earthquake Name', 'Station code', 'Component 1', 'Component 2', 'Component 3', 'Last update'])], ignore_index=True)
    computed.sort_values(by=['Earthquake Name', 'Station code'], inplace=True)

    computed.to_excel(os.path.join(draftPath, 'spectralValues', 'computed.xlsx'), index=False)

    widget.insert('end', 'Espectros actualizados.\n')
    widget.see('end')
    window.update_idletasks()
