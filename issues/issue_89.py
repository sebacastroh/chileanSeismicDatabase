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

flatfile = pd.read_csv(os.path.join(basePath, 'data', 'flatFile - backup.csv'))
bastiasFlatfile = pd.read_csv(os.path.join('extras', 'issue_89_Flatfile_chileanSiteStation_v7.csv'))

with open(os.path.join(basePath, 'data', 'stationsInfo.json')) as f:
    sinfo = json.load(f)

events = flatfile[(flatfile['Vs30 [m/s]'].isna()) | (flatfile['Vs30 [m/s]'] < 0)]

stationCodes = events['Station code'].unique().tolist()

valid = []
for stationCode in stationCodes:
    station = bastiasFlatfile[bastiasFlatfile['CodeSta'] == stationCode]
    
    if len(station) == 0:
        continue
    
    vs30 = station.iloc[0]['PreferedVs30']
    
    valid.append((stationCode, vs30))
    sinfo[stationCode][2] = vs30

with open(os.path.join(basePath, 'data', 'stationsInfo.json'), 'w') as f:
    json.dump(sinfo, f, indent=DEFAULT_INDENT, ensure_ascii=False)
    
valid  = pd.DataFrame(valid, columns=['Station code', 'New Vs30'])
to_fix = events.reset_index().merge(valid, on='Station code', how='inner')

currentEvent = None
save = False
data = {}
for r, row in to_fix.iterrows():
    event_id = row['Earthquake Name']
    
    if currentEvent != event_id:
        if save:
            np.savez_compressed(os.path.join(draftPath, 'seismicDatabase', 'npz', currentEvent), **data)
            spio.savemat(os.path.join(draftPath, 'seismicDatabase', 'mat', currentEvent + '.mat'), data, do_compression=True)
            
            currentEvent = event_id
        
        with np.load(os.path.join(dataPath, 'seismicDatabase', 'npz', f'{event_id}.npz'), allow_pickle=True) as f:
            data = {}
            for key, value in f.items():
                data[key] = value.item()
            save = True
        
    for st, station in data.items():
        if not st.startswith('st'):
            continue
        
        if station['station_code'] == row['Station code']:
            data[st]['vs30'] = row['New Vs30']
            break

np.savez_compressed(os.path.join(draftPath, 'seismicDatabase', 'npz', event_id), **data)
spio.savemat(os.path.join(draftPath, 'seismicDatabase', 'mat', event_id + '.mat'), data, do_compression=True)

flatfile = ''
with open(os.path.join(basePath, 'data', 'flatFile - backup.csv')) as f:
    lines = f.readlines()

indices = to_fix['index'].tolist()

for i, line in enumerate(lines):
    if i-1 not in indices:
        flatfile += line
        continue

    elements = line.split(',')
    elements[-5] = '%i' %to_fix[to_fix['index'] == i-1].iloc[0]['New Vs30']

    new_line = ','.join(elements)
    flatfile += new_line

with open(os.path.join(basePath, 'data', 'flatFile - backup.csv'), 'w', encoding='utf8') as f:
    f.write(flatfile)

# with open(os.path.join(draftPath, 'flatFile.csv'), 'w', encoding='utf8') as f:
    # f.write(flatfile)

# df = pd.read_csv(os.path.join(basePath, 'data', 'flatFile - backup.csv'))
# df.to_excel(os.path.join(draftPath, 'flatFile.xlsx'), index=False)