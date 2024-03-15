import os
import json
import datetime
import numpy as np
import multiprocessing
import scipy.io as spio
import lib.automaticCorrection as automaticCorrection

p_waves    = {}
to_correct = []
base_path  = ''
data_path  = ''

DEFAULT_INDENT = 2
SORT_KEYS      = True

def correctRecords(window, widget, basePath, dataPath):
    """
    Función que determina las estaciones que deben ser corregidasde acuerdo al archivo ```siberrisk.csv```.

    :param window: Ventana de ```tkinter``` donde se presentan los resultados.
    :type window: tkinter.Toplevel

    :param widget: Texto de ```tkinter``` para presentar los avances del método.
    :type widget: tkinter.Text

    :param basePath: Ruta base de ejecución del código.
    :type basePath: str

    :return: Booleano que indica si hay estaciones a corregir (True) o no (False)
    :rtype: boolean
    """
    global p_waves, to_correct

    # Actualización log
    widget.insert('end', 'Inicio proceso de corrección de línea base\n\n')
    widget.see('end')
    window.update_idletasks()

    # Registered P-Waves
    if os.path.exists(os.path.join(basePath, 'data', 'p_waves.json')):
        with open(os.path.join(basePath, 'data', 'p_waves.json')) as f:
            p_waves = json.load(f)
    else:
        p_waves = {}

    # Select events and stations to correct
    to_correct = []
    widget.insert('end', 'Eventos y estaciones a corregir\n\n')
    widget.see('end')
    window.update_idletasks()
    to_correct = []
    for event_id, p_wave in p_waves.items():
        for scode, sinfo in p_wave.items():
             if sinfo['status'] and not sinfo['corrected']:
                to_correct.append([event_id, scode])
                widget.insert('end', event_id + ' - ' + scode + '\n')
                widget.see('end')
                window.update_idletasks()

    if len(to_correct) == 0:
        widget.insert('end', 'No hay registros que corregir.\n')
        widget.see('end')
        window.update_idletasks()
        return False

    widget.insert('end', '\nNúmero total de registros a corregir: %i\n' %len(to_correct))
    widget.see('end')
    window.update_idletasks()

    return True

def applyCorrection(window, widget, basePath, dataPath, parallel):
    """
    Función que aplica la corrección de línea base a los registros identificados.

    :param window: Ventana de ```tkinter``` donde se presentan los resultados.
    :type window: tkinter.Toplevel

    :param widget: Texto de ```tkinter``` para presentar los avances del método.
    :type widget: tkinter.Text

    :param basePath: Ruta base de ejecución del código.
    :type basePath: str

    :param parallel: Booleano que indica si la corrección de las estaciones se ejecuta en paralelo.
    Este proceso puede requerir un alto espacio de memoria en el disco. Por defecto es ```False```.
    :type parallel: bool
    """
    global p_waves, to_correct, base_path, data_path
    base_path = basePath
    data_path = dataPath

    widget.insert('end', '\nIniciando proceso de corrección.\n')
    widget.see('end')
    window.update_idletasks()

    if parallel:
        widget.insert('end', '\nEjecución en paralelo, el avance no se mostrará en pantalla.\n')
        widget.see('end')
        window.update_idletasks()

        if not os.path.exists(os.path.join(basePath, 'tmp')):
            os.mkdir(os.path.join(basePath, 'tmp'))

        pool = multiprocessing.Pool(processes=8)
        pool.map(parallel_run, to_correct)
        pool.close()

        widget.insert('end', '\nOperación terminada, se procede a guardar los resultados.\n')
        widget.see('end')
        window.update_idletasks()

    save = False
    current_event_id = None
    for event_id, station_code in to_correct:
        if current_event_id != event_id:
            if save:
                np.savez_compressed(os.path.join(dataPath, 'seismicDatabase', 'npz', current_event_id), **data)
                spio.savemat(os.path.join(dataPath, 'seismicDatabase', 'mat', current_event_id + '.mat'), data, do_compression=True)
                
                with open(os.path.join(basePath, 'data', 'p_waves.json'), 'w') as f:
                    json.dump(p_waves, f, indent=DEFAULT_INDENT, sort_keys=SORT_KEYS)

            current_event_id = event_id
            with np.load(os.path.join(dataPath, 'seismicDatabase', 'npz', event_id + '.npz'), allow_pickle=True) as f:
                data = {}
                for key, value in f.items():
                    data[key] = value.item()

            save = False

        for key, station in data.items():
            if not key.startswith('st'):
                continue

            if station.get('station_code') == station_code:
                dt     = station.get('dt')
                status = p_waves[event_id][station_code]['status']
                p_wave = p_waves[event_id][station_code]['pos']

                if not parallel:
                    widget.insert('end', '\nCorrigiendo '  + event_id + ' - ' + station_code + '.\n')
                    widget.see('end')
                    window.update_idletasks()

                for i in range(3):
                    if parallel:
                        temp = np.load(os.path.join(basePath, 'tmp', event_id + '_' + station_code + '_%i.npz' %(i+1)))
                        acc_corr = temp['acc_corr']
                        acc_fil  = temp['acc_fil']
                        freqs    = temp['freqs']
                        os.remove(os.path.join(basePath, 'tmp', event_id + '_' + station_code + '_%i.npz' %(i+1)))
                    else:
                        acc = station.get('acc_uncorrected_%i' %(i+1))
                        acc_corr, acc_fil, freqs = automaticCorrection.correctRecord(acc, dt, status, p_wave, saveInTemp=False, filename='')

                    data[key]['acc_corrected_%i' %(i+1)] = acc_corr.copy()
                    data[key]['acc_filtered_%i' %(i+1)]  = acc_fil.copy()
                    data[key]['corner_freqs_%i' %(i+1)]  = freqs.copy()

                    if not parallel:
                        widget.insert('end', 'Componente %i completada.\n' %(i+1))
                        widget.see('end')
                        window.update_idletasks()

                updated = datetime.datetime.now().isoformat()
                if status:
                    data[key]['p_wave']      = p_wave
                    data[key]['last_update'] = updated

                p_waves[event_id][station_code]['corrected'] = True
                p_waves[event_id][station_code]['updated']   = updated
                save = True
                break

    if save:
        np.savez_compressed(os.path.join(dataPath, 'seismicDatabase', 'npz', current_event_id), **data)
        spio.savemat(os.path.join(dataPath, 'seismicDatabase', 'mat', current_event_id + '.mat'), data, do_compression=True)
        
        with open(os.path.join(basePath, 'data', 'p_waves.json'), 'w') as f:
            json.dump(p_waves, f, indent=DEFAULT_INDENT, sort_keys=SORT_KEYS)

    widget.insert('end', '\nProceso finalizado.')
    widget.see('end')
    window.update_idletasks()


def parallel_run(combination):
    global base_path, data_path
    basePath = base_path
    dataPath = data_path

    event_id, station_code = combination
    with np.load(os.path.join(dataPath, 'data', 'seismicDatabase', 'npz', event_id + '.npz'), allow_pickle=True) as f:
        data = {}
        for key, value in f.items():
            data[key] = value.item()

    for key, station in data.items():
        if not key.startswith('st'):
            continue
        if station.get('station_code') == station_code:
            dt     = station.get('dt')
            status = p_waves[event_id][station_code]['status']
            p_wave = p_waves[event_id][station_code]['pos']

            for i in range(3):
                acc = station.get('acc_uncorrected_%i' %(i+1))
                acc_corr, acc_fil, freqs = automaticCorrection.correctRecord(acc, dt, status, p_wave, saveInTemp=True,
                    filename=os.path.join(basePath, 'tmp', event_id + '_' + station_code + '_%i.npz' %(i+1)))
