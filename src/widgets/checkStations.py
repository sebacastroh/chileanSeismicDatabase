import os
import json
import pandas as pd

def checkStations(window, widget, basePath, dataPath, draftPath):
    # Registered stations
    with open(os.path.join(basePath, 'data', 'stationsInfo.json')) as f:
        sinfo = json.load(f)
    registeredStations = list(sinfo.keys())
    registeredStations = [registeredStation.replace(' ', '').upper() for registeredStation in registeredStations]

    # I
    events = pd.read_csv(os.path.join(basePath, 'data', 'events.csv'), parse_dates=True)
    foundStations = list(events['Estaciones'].str.split('; ').explode().unique())

    # Missing data
    missingStations = []
    for station in foundStations:
        if station.upper() not in registeredStations:
            missingStations.append(station)

    if len(missingStations) == 0:
        return

    # Find events with those missing stations
    stations  = events['Estaciones'].str.split('; ').explode()
    events_id = events.loc[stations.index, 'Identificador']

    df1 = pd.DataFrame({'Estación': stations.to_list(), 'Identificador': events_id.to_list()})
    df2 = pd.DataFrame({'Estación': missingStations})

    merged = df1.merge(df2, on='Estación')
    merged.sort_values(['Estación','Identificador'], inplace=True)
