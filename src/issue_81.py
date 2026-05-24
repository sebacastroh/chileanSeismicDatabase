import io
import os
import json
import zipfile
import datetime
import requests
import numpy as np
import pandas as pd
import scipy.io as spio

basePath  = os.path.abspath('.')
draftPath = os.path.abspath(os.path.join('.', '..', 'draft'))
dataPath  = os.path.abspath(os.path.join('.', '..', 'data'))

DEFAULT_INDENT = 2
SORT_KEYS      = True

df = pd.read_excel(os.path.join(dataPath, 'spectralValues', 'xi_0.05', 'component_1.xlsx'))
to_fix = df[df[0] > 1][['Earthquake Name', 'Station code']]

events = pd.read_csv(os.path.join(basePath, 'data', 'events.csv'))
events.rename(columns={'ID': 'Earthquake Name'}, inplace=True)

to_fix = to_fix.merge(events[['Earthquake Name', 'Identificador']], on='Earthquake Name', how='left')
failed_files = []

with open(os.path.join(basePath, 'data', 'p_waves.json')) as f:
    p_waves = json.load(f)

for r, row in to_fix.iterrows():
    event_id        = row['Identificador']
    station_id      = row['Station code']
    earthquake_name = row['Earthquake Name']

    zipfilename = os.path.join(basePath, 'tmp', event_id + '_' + station_id + '.zip')

    if not os.path.exists(zipfilename):
        url = f'https://evtdb.csn.uchile.cl/write/{event_id}/{station_id}'
        response = requests.get(url)
        try:
            response.raise_for_status()
            with open(zipfilename, 'wb') as f:
                f.write(response.content)
        except:
            failed_files.append([event_id, station_id])
            continue

    if os.path.exists(os.path.join(draftPath, 'seismicDatabase', 'npz', f'{earthquake_name}.npz')):
        npz_filename = os.path.join(draftPath, 'seismicDatabase', 'npz', f'{earthquake_name}.npz')
    else:
        npz_filename = os.path.join(dataPath, 'seismicDatabase', 'npz', f'{earthquake_name}.npz')

    with np.load(npz_filename, allow_pickle=True) as f:
        data = {}
        for key, value in f.items():
            data[key] = value.item()

    for st, record in data.items():
        if not st.startswith('st'):
            continue
        if record['station_code'] == station_id:
            break

    acc_1 = np.empty(0)
    acc_2 = np.empty(0)
    acc_3 = np.empty(0)

    x1 = np.inf
    x2 = np.inf
    x3 = np.inf

    with zipfile.ZipFile(zipfilename) as zf:
        filenames = zf.namelist()
        for filename in filenames:
            txt = io.BytesIO(zf.read(filename))
            acc = np.loadtxt(txt, skiprows=6)
            txt.seek(0)
            lines = []
            while True:
                line = txt.readline().decode().strip()
                if not line.startswith('#'):
                    break
                lines.append(line)
            txt.close()
            npts = len(acc)
            channel = filename.split('-')[-1][:-4]

            if data[st]['component_1'] == channel:
                x1 = datetime.datetime.fromisoformat(lines[0][20:])
                acc_1 = acc.copy()
            elif data[st]['component_2'] == channel:
                x2 = datetime.datetime.fromisoformat(lines[0][20:])
                acc_2 = acc.copy()
            elif data[st]['component_3'] == channel:
                x3 = datetime.datetime.fromisoformat(lines[0][20:])
                acc_3 = acc.copy()

    xini = np.min([x1, x2, x3])
    if xini < x1:
        delta = x1 - xini
        n     = int(delta.total_seconds()/data[st]['dt'])
        if n > 0:
            acc_1 = np.hstack((np.zeros(n), acc_1))

    if xini < x2:
        delta = x2 - xini
        n     = int(delta.total_seconds()/data[st]['dt'])
        if n > 0:
            acc_2 = np.hstack((np.zeros(n), acc_2))

    if xini < x3:
        delta = x3 - xini
        n     = int(delta.total_seconds()/data[st]['dt'])
        if n > 0:
            acc_3 = np.hstack((np.zeros(n), acc_3))

    data[st]['acc_uncorrected_1'] = acc_1
    data[st]['acc_corrected_1']   = np.empty(0)
    data[st]['acc_filtered_1']    = np.empty(0)

    data[st]['acc_uncorrected_2'] = acc_2
    data[st]['acc_corrected_2']   = np.empty(0)
    data[st]['acc_filtered_2']    = np.empty(0)

    data[st]['acc_uncorrected_3'] = acc_3
    data[st]['acc_corrected_3']   = np.empty(0)
    data[st]['acc_filtered_3']    = np.empty(0)

    updated = datetime.datetime.now().isoformat()

    data[st]['last_update'] = updated

    p_waves[earthquake_name][station_id] = {
        "status": None,
        "corrected": False,
        'updated': updated
    }

    np.savez_compressed(os.path.join(draftPath, 'seismicDatabase', 'npz', earthquake_name), **data)
    spio.savemat(os.path.join(draftPath, 'seismicDatabase', 'mat', earthquake_name + '.mat'), data, do_compression=True)
    
    with open(os.path.join(basePath, 'data', 'p_waves.json'), 'w') as f:
        json.dump(p_waves, f, indent=DEFAULT_INDENT, sort_keys=SORT_KEYS)
