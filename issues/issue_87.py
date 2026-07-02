import os
import json
import numpy as np
import pandas as pd
import scipy.io as spio

basePath  = os.path.abspath(os.path.join('.', '..', 'src'))
dataPath  = os.path.abspath(os.path.join('.', '..', 'data'))
draftPath = os.path.abspath(os.path.join('.', '..', 'draft'))

DEFAULT_INDENT = 2
SORT_KEYS      = True

event_id = '20100227_8.8M_36.17S_73.14W_30.1KM'

with np.load(os.path.join(dataPath, 'seismicDatabase', 'npz', f'{event_id}.npz'), allow_pickle=True) as f:
    data = {}
    for key, value in f.items():
        data[key] = value.item()

with open(os.path.join(basePath, 'data', 'stationsInfo.json')) as f:
    sinfo = json.load(f)

fix_station_names = pd.read_csv(os.path.join('extras', 'issue_87_corrected_station_codes.csv'), encoding='utf8')

for st, station in data.items():
    if not st.startswith('st'):
        continue

    old_station_code = station['station_code']
    properties = fix_station_names[fix_station_names['Old Station Code'] == old_station_code].iloc[0]

    data[st]['station_code'] = properties['Station Code']
    data[st]['station_name'] = properties['Station Name']
    data[st]['vs30']         = sinfo[properties['Station Code']][2]

    sinfo[properties['Station Code']][5] = properties['Station Name']

np.savez_compressed(os.path.join(draftPath, 'seismicDatabase', 'npz', event_id), **data)
spio.savemat(os.path.join(draftPath, 'seismicDatabase', 'mat', event_id + '.mat'), data, do_compression=True)

with open(os.path.join(basePath, 'data', 'stationsInfo.json'), 'w') as f:
    json.dump(sinfo, f, indent=DEFAULT_INDENT, ensure_ascii=False)

flatfile = ''
with open(os.path.join(dataPath, 'flatFile.csv')) as f:
    lines = f.readlines()

for line in lines:
    if not line.startswith(event_id):
        flatfile += line
        continue

    elements = line.split(',')

    old_station_code = elements[9]
    properties = fix_station_names[fix_station_names['Old Station Code'] == old_station_code].iloc[0]

    elements[8]  = properties['Station Name']
    elements[9]  = properties['Station Code']
    elements[17] = '%i' %sinfo[properties['Station Code']][2]

    new_line = ','.join(elements)
    flatfile += new_line

with open(os.path.join(basePath, 'data', 'flatFile - backup.csv'), 'w', encoding='utf8') as f:
    f.write(flatfile)

with open(os.path.join(draftPath, 'flatFile.csv'), 'w', encoding='utf8') as f:
    f.write(flatfile)

df = pd.read_csv(os.path.join(basePath, 'data', 'flatFile - backup.csv'))
df.to_excel(os.path.join(draftPath, 'flatFile.xlsx'), index=False)

with open(os.path.join(basePath, 'data', 'p_waves.json')) as f:
    p_waves = json.load(f)

old_station_codes = list(p_waves[event_id].keys())

for old_station_code in old_station_codes:
    properties = fix_station_names[fix_station_names['Old Station Code'] == old_station_code].iloc[0]

    p_waves[event_id][properties['Station Code']] = p_waves[event_id][old_station_code].copy()
    p_waves[event_id].pop(old_station_code, None)

with open(os.path.join(basePath, 'data', 'p_waves.json'), 'w') as f:
    json.dump(p_waves, f, indent=DEFAULT_INDENT, sort_keys=SORT_KEYS)
