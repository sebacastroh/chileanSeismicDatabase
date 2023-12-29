import os
import json
import numpy as np
import pandas as pd

def correctRecords(window, widget, basePath):
    """
    Función que ejecuta la corrección de línea base de los registros listados en el archivo ```siberrisk.csv```.

    :param window: Ventana de ```tkinter``` donde se presentan los resultados.
    :type window: tkinter.Toplevel

    :param widget: Texto de ```tkinter``` para presentar los avances del método.
    :type widget: tkinter.Text

    :param basePath: Ruta base de ejecución del código.
    :type basePath: str
    """
    # Actualización log
    widget.insert('end', 'Inicio proceso de corrección de línea base\n\n')
    widget.see('end')
    window.update_idletasks()

    # Lista de eventos
    df = pd.read_csv(os.path.join(basePath, 'data', 'events.csv'), dtype={'Identificador': str})

    # Registered P-Waves
    if os.path.exists(os.path.join(basePath, 'data', 'p_waves.json')):
        with open(os.path.join(basePath, 'data', 'p_waves.json')) as f:
            p_waves = json.load(f)
    else:
        p_waves = {}

    events_with_p_waves = list(p_waves.keys())

    # List events with registered P-Waves
    subdf = df[df['ID'].isin(events_with_p_waves)]

    # Select events and stations to correct
    to_correct = []
    widget.insert('end', 'Eventos y estaciones a corregir\n\n')
    widget.see('end')
    window.update_idletasks()
    for r, row in subdf.iterrows():
        event_id       = row['ID']
        stations_codes = row['Estaciones'].split('; ')

        for station_code in stations_codes:
            station_info = p_waves.get(event_id).get(station_code)
            if station_info is None:
                continue
            elif station_info.get('status') and station_info.get('corrected'):
                continue
            else:
                to_correct.append([event_id, station_code])

                widget.insert('end', event_id + ' - ' + station_code + '\n')
                widget.see('end')
                window.update_idletasks()

    if len(to_correct) == 0:
        widget.insert('end', 'No hay registros que corregir.\n')
        widget.see('end')
        window.update_idletasks()
        return

    widget.insert('end', '\nIniciando proceso de corrección.\n\n')
    widget.see('end')
    window.update_idletasks()
