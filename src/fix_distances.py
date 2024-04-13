import os
import sys
import pickle
import pyperclip
import webbrowser
import subprocess

import json
import pyproj
import datetime
import numpy as np
import pandas as pd
import scipy.io as spio
import lib.computeDistances as computeDistances

basePath = os.path.abspath('.')
dataPath = os.path.abspath(os.path.join('.', '..', 'data'))

DEFAULT_INDENT = 2
SORT_KEYS      = True

flatfile = pd.read_csv(os.path.join(dataPath, 'flatFile.csv'), parse_dates=['Earthquake date', 'Start time record', 'Last update'])

## Registros con distancias no calculadas pese a tener la información
filtered = flatfile[(~flatfile['Magnitude [Mw]'].isna()) &
                    (~flatfile['Hypocenter latitude'].isna()) &
                    (~flatfile['Hypocenter longitude'].isna()) &
                    (~flatfile['Depth [km]'].isna()) &
                    (~flatfile['Station latitude'].isna()) &
                    (~flatfile['Station longitude'].isna()) &
                    (flatfile['Hypocentral distance [km]'].isna() |
                     flatfile['Epicentral distance [km]'].isna() |
                     flatfile['Rupture distance [km]'].isna() |
                     flatfile['Joyner-Boore distance [km]'].isna() |
                     (flatfile['Hypocentral distance [km]'] <= 0) |
                     (flatfile['Epicentral distance [km]'] <= 0) |
                     (flatfile['Rupture distance [km]'] <= 0) |
                     (flatfile['Joyner-Boore distance [km]'] <= 0))]

slab = np.load(os.path.join(basePath, 'data', 'sam_slab2.npz'))
geod = pyproj.Geod(ellps='WGS84')

distances = []

ffm_files = os.listdir(os.path.join(basePath, 'data', 'ffm'))

save = False
current_event_id = None
for r, row in filtered.iterrows():
    event_id     = row['Earthquake Name']
    station_code = row['Station code']

    if current_event_id != event_id:
        # if save:
        #     np.savez_compressed(os.path.join(dataPath, 'seismicDatabase', 'npz', current_event_id), **data)
        #     spio.savemat(os.path.join(dataPath, 'seismicDatabase', 'mat', current_event_id + '.mat'), data, do_compression=True)

        #     with open(os.path.join(basePath, 'data', 'p_waves.json'), 'w') as f:
        #         json.dump(p_waves, f, indent=DEFAULT_INDENT, sort_keys=SORT_KEYS)

        current_event_id = event_id
        # with np.load(os.path.join(dataPath, 'seismicDatabase', 'npz', event_id + '.npz'), allow_pickle=True) as f:
        #     data = {}
        #     for key, value in f.items():
        #         data[key] = value.item()

        # event_lat = data['st00']['hypocenter_lat']
        # event_lon = data['st00']['hypocenter_lon']
        # event_dep = data['st00']['depth']
        # event_mag = data['st00']['magnitude']

        event_lat = row['Hypocenter latitude']
        event_lon = row['Hypocenter longitude']
        event_dep = row['Depth [km]']
        event_mag = row['Magnitude [Mw]']

        hypocenter = computeDistances.LatLonDepth2XYZNum(event_lat, event_lon, event_dep)
        save = False

    # for st, station in data.items():
    #     if not st.startswith('st'):
    #         continue

    #     if station.get('station_code') != station_code:
    #         continue

    #     
    #     stationLat = station.get('station_lat')
    #     stationLon = station.get('station_lon')

    stationLat = row['Station latitude']
    stationLon = row['Station longitude']

    Rhypo = None
    Repi  = None
    Rrup  = None
    Rjb   = None
    station_xyz = computeDistances.LatLonDepth2XYZNum(stationLat, stationLon, 0.)

    if event_mag < 7.1:
        Rhypo = np.sqrt(np.sum((hypocenter - station_xyz)**2))
        azimuth1, azimuth2, Repi = geod.inv(stationLon, stationLat, event_lon, event_lat)
        Repi /= 1000.
        Rrup = Rhypo
        Rjb = Repi
    else:
        if event_id + '.geojson' in ffm_files:
            with open(os.path.join(basePath, 'data', 'ffm', event_id + '.geojson')) as f:
                ffm = json.load(f)

            max_slip = -np.inf
            min_slip =  np.inf
            for feature in ffm['features']:
                slip = feature['properties']['slip']
                max_slip = max(max_slip, slip)
                min_slip = min(min_slip, slip)
            
            dslip = 0.15*(max_slip - min_slip) + min_slip

            Rhypo = np.sqrt(np.sum((hypocenter - station_xyz)**2))
            azimuth1, azimuth2, Repi = geod.inv(stationLon, stationLat, event_lon, event_lat)
            Repi /= 1000.

            Rrup = np.inf
            Rjb = np.inf
            for feature in ffm['features']:
                slip = feature['properties']['slip']
                if slip < dslip:
                    continue

                points = np.array(feature['geometry']['coordinates'][0])[:-1]
                sub_fault = computeDistances.LatLonDepth2XYZNumPy(points[:,1], points[:,0], points[:,2]/1000.)
                dist = np.sqrt(np.sum((sub_fault - station_xyz)**2, axis=1)).min()

                Rrup = min(Rrup, dist)

                if computeDistances.inPolygon(station_xyz[:2], sub_fault[:,:2]):
                    Rjb = 0.
                elif Rjb != 0.:
                    azimuth1, azimuth2, dist = geod.inv(stationLon*np.ones(len(points)),\
                                                        stationLat*np.ones(len(points)),\
                                                        points[:,0], points[:,1])
                    Rjb = min(Rjb, dist.min()/1000.)
            else:
                pass

    distances.append([event_id, station_code, Rhypo, Repi, Rrup, Rjb])
