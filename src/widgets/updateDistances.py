import os
import json
import pyproj
import datetime
import numpy as np
import pandas as pd
import scipy.io as spio
import lib.computeDistances as computeDistances

DEFAULT_INDENT = 2
SORT_KEYS      = True

def updateDistances(window, widget, basePath, dataPath):
    """
    Función que intenta recalcular las cuatro distancias definidas para cada estación en un evento.

    Para esto, filtra el archivo flatFile.csv y si alguna de las cuatro distancias no es un valor positivo, intenta calcular las distancias.

    En caso de poder calcularlas, se actualiza la base de datos, p_waves.json y computed.xlsx (este último porque no se requiere actualizar el espectro).

    :param window: Ventana de ```tkinter``` donde se presentan los resultados.
    :type window: tkinter.Toplevel

    :param widget: Texto de ```tkinter``` para presentar los avances del método.
    :type widget: tkinter.Text

    :param basePath: Ruta base de ejecución del código.
    :type basePath: str

    :param dataPath: Ruta de acceso a la base de datos.
    :type dataPath: str
    """

    flatfile = pd.read_csv(os.path.join(basePath, 'data', 'flatFile - backup.csv'), parse_dates=['Earthquake date', 'Start time record', 'Last update'])
    if os.path.exists(os.path.join(dataPath, 'spectralValues', 'computed.xlsx')):
        computed = pd.read_excel(os.path.join(dataPath, 'spectralValues', 'computed.xlsx'))
    else:
        computed = None

    candidates = flatfile[~((flatfile['Hypocentral distance [km]']  >= 0)
                          & (flatfile['Epicentral distance [km]']   >= 0)
                          & (flatfile['Rupture distance [km]']      >= 0)
                          & (flatfile['Joyner-Boore distance [km]'] >= 0))]

    candidates = candidates[~((candidates['Magnitude [Mw]'].isna()      )
                            | (candidates['Hypocenter latitude'].isna() )
                            | (candidates['Hypocenter longitude'].isna())
                            | (candidates['Depth [km]'].isna()          )
                            | (candidates['Station latitude'].isna()    )
                            | (candidates['Station longitude'].isna()   ))]

    geod = pyproj.Geod(ellps='WGS84')

    with open(os.path.join(basePath, 'data', 'p_waves.json')) as f:
        p_waves = json.load(f)

    with open(os.path.join(basePath, 'data', 'fault_plane_properties.json')) as f:
        fault_properties = json.load(f)

    data     = {}
    modified = False

    for c, candidate in candidates.iterrows():
        event_id     = candidate['Earthquake Name']
        station_code = candidate['Station code']

        magnitude = candidate['Magnitude [Mw]']
        event_lat = candidate['Hypocenter latitude']
        event_lon = candidate['Hypocenter longitude']
        event_dep = candidate['Depth [km]']

        stationLat = candidate['Station latitude']
        stationLon = candidate['Station longitude']

        hypocenter  = computeDistances.LatLonDepth2XYZNum(event_lat, event_lon, event_dep)
        station_xyz = computeDistances.LatLonDepth2XYZNum(stationLat, stationLon, 0.)

        widget.insert('end', 'Analizando estación %s de evento %s... ' %(station_code, event_id))
        widget.see('end')
        window.update_idletasks()

        if magnitude < 7.1:
            Rhypo = np.sqrt(np.sum((hypocenter - station_xyz)**2))
            azimuth1, azimuth2, Repi = geod.inv(stationLon, stationLat, event_lon, event_lat)
            Repi /= 1000.
            Rrup  = Rhypo
            Rjb   = Repi
        else:
            Rhypo = np.sqrt(np.sum((hypocenter - station_xyz)**2))
            azimuth1, azimuth2, Repi = geod.inv(stationLon, stationLat, event_lon, event_lat)
            Repi /= 1000.

            if os.path.exists(os.path.join(basePath, 'data', 'ffm', event_id + '.geojson')):
                with open(os.path.join(basePath, 'data', 'ffm', event_id + '.geojson')) as f:
                    ffm = json.load(f)

                max_slip = -np.inf
                min_slip =  np.inf
                for feature in ffm['features']:
                    slip = feature['properties']['slip']
                    max_slip = max(max_slip, slip)
                    min_slip = min(min_slip, slip)

                dslip = 0.15*(max_slip - min_slip) + min_slip

                Rrup = np.inf
                Rjb  = np.inf
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

            elif fault_properties.get(event_id) is not None:
                strike, dip, rake = fault_properties.get(event_id)

                L = 10**(-2.9  + 0.63*magnitude) # Allen et al 2017, Table 2
                W = 10**(-0.86 + 0.35*magnitude) # Allen et al 2017, Table 2
                subfaults = computeDistances.generateSubFaults(event_lat, event_lon, event_dep, strike, dip, L, W)

                faultXYZ = computeDistances.LatLonDepth2XYZNumPy(subfaults[:,1],\
                                                                 subfaults[:,0],\
                                                                 subfaults[:,2])

                Rrup = np.sqrt(np.sum((faultXYZ - station_xyz)**2, axis=1)).min()

                nx = int(L)
                ny = int(W)
                polygon = subfaults[[0,ny-1,nx*ny-1,nx*ny-ny]]

                if computeDistances.inPolygon([stationLon, stationLat], polygon):
                    Rjb = 0.
                else:
                    azimuth1, azimuth2, dist = geod.inv(stationLon*np.ones(nx*ny),\
                                                        stationLat*np.ones(nx*ny),\
                                                        subfaults[:,0], subfaults[:,1])

                    Rjb = dist.min()/1000.
            else:
                Rrup = -1
                Rjb  = -1

        if data.get('event_id') != event_id:
            if modified:
                np.savez_compressed(os.path.join(dataPath, 'seismicDatabase', 'npz', data.get('event_id') + '.npz'), **data)
                spio.savemat(os.path.join(dataPath, 'seismicDatabase', 'mat', data.get('event_id') + '.mat'), data, do_compression=True)

            with np.load(os.path.join(dataPath, 'seismicDatabase', 'npz', event_id + '.npz'), allow_pickle=True) as f:
                data = {}
                for key, value in f.items():
                    data[key] = value.item()
            modified = False

        nStations = len(data) - 4
        for i in range(nStations):
            st = 'st%0.2i' %i
            station = data[st]

            if station['station_code'] != station_code:
                continue

            if distant_distances(station['Rhypo'], Rhypo) or distant_distances(station['Repi'], Repi) \
                or distant_distances(station['Rrup'], Rrup) or distant_distances(station['Rjb'], Rjb):

                data[st]['magnitude']      = magnitude
                data[st]['hypocenter_lat'] = event_lat
                data[st]['hypocenter_lon'] = event_lon
                data[st]['depth']          = event_dep
                data[st]['station_lat']    = stationLat
                data[st]['station_lon']    = stationLon

                data[st]['Rhypo'] = Rhypo
                data[st]['Repi']  = Repi
                data[st]['Rrup']  = Rrup if Rrup >= 0 else 'Finite fault model required'
                data[st]['Rjb']   = Rjb  if Rjb  >= 0 else 'Finite fault model required'

                update = datetime.datetime.now().isoformat()
                data[st]['last_update']                    = update
                p_waves[event_id][station_code]['updated'] = update

                if computed is not None:
                    index = computed[(computed['Earthquake Name'] == event_id) & (computed['Station code'] == station_code)].index
                    computed.loc[index, 'Last update'] = update

                flatfile.loc[c, 'Last update'] = update
                flatfile.loc[c, 'Hypocentral distance [km]']  = Rhypo
                flatfile.loc[c, 'Epicentral distance [km]']   = Repi
                flatfile.loc[c, 'Rupture distance [km]']      = Rrup
                flatfile.loc[c, 'Joyner-Boore distance [km]'] = Rjb

                modified = True

                widget.insert('end', 'Distancias actualizadas.\n')
                widget.see('end')
                window.update_idletasks()
            else:
                widget.insert('end', 'No se registraron cambios.\n')
                widget.see('end')
                window.update_idletasks()

    if modified:
        np.savez_compressed(os.path.join(dataPath, 'seismicDatabase', 'npz', event_id + '.npz'), **data)
        spio.savemat(os.path.join(dataPath, 'seismicDatabase', 'mat', event_id + '.mat'), data, do_compression=True)

    with open(os.path.join(basePath, 'data', 'p_waves.json'), 'w') as f:
        json.dump(p_waves, f, indent=DEFAULT_INDENT, sort_keys=SORT_KEYS)

    flatfile.to_excel(os.path.join(dataPath, 'flatFile.xlsx'), index=False)
    flatfile.to_csv(os.path.join(dataPath, 'flatFile.csv'), index=False)
    flatfile.to_csv(os.path.join(basePath, 'data', 'flatFile - backup.csv'), index=False)

    if computed is not None:
        computed.to_excel(os.path.join(dataPath, 'spectralValues', 'computed.xlsx'), index=False)

    widget.insert('end', '\nProceso finalizado.')
    widget.see('end')
    window.update_idletasks()

def distant_distances(d1, d2):

    if isinstance(d1, str):
        d1 = -1

    if isinstance(d2, str):
        d2 = -1

    if abs(d1 - d2) < 1e-8:
        return False
    
    return True

