import os
import numpy as np
import pandas as pd
import scipy.io as spio

basePath  = os.path.abspath(os.path.join('.', '..', 'src'))
dataPath  = os.path.abspath(os.path.join('.', '..', 'data'))
draftPath = os.path.abspath(os.path.join('.', '..', 'draft'))

flatfile = pd.read_csv(os.path.join(basePath, 'data', 'flatFile - backup.csv'))

to_change = [
 ['20190315_6.3M_17.85S_65.91W_353.0KM', 'Undetermined'],
 ['20150916_8.4M_31.55S_71.86W_11.0KM', 'interface'],
 ['20140401_8.2M_19.57S_70.91W_39.0KM', 'interface']
]

events = flatfile[
    (flatfile['Event type'] == 'interface') &
    (flatfile['Depth [km]'] > 80) &
    (flatfile['Earthquake Name'] != to_change[0][0]) &
    (flatfile['Earthquake Name'] != to_change[1][0]) &
    (flatfile['Earthquake Name'] != to_change[2][0])].drop_duplicates(subset='Earthquake Name')['Earthquake Name'].tolist()

for event in events:
    to_change.append([event, 'intraslab'])

for event_id, event_type in to_change:
    filename = os.path.join(draftPath, 'seismicDatabase', 'npz', f'{event_id}.npz')
    if not os.path.exists(filename):
        filename = os.path.join(dataPath, 'seismicDatabase', 'npz', f'{event_id}.npz')

    with np.load(filename, allow_pickle=True) as f:
        data = {}
        for key, value in f.items():
            data[key] = value.item()

    for st, station in data.items():
        if not st.startswith('st'):
            continue

        data[st]['event_type'] = event_type

    np.savez_compressed(os.path.join(draftPath, 'seismicDatabase', 'npz', event_id), **data)
    spio.savemat(os.path.join(draftPath, 'seismicDatabase', 'mat', event_id + '.mat'), data, do_compression=True)

# Update flatfile
to_change = pd.DataFrame(to_change, columns=['Earthquake Name', 'New event type'])
to_fix    = flatfile.reset_index().merge(to_change, on='Earthquake Name', how='inner')
indices   = to_fix['index'].tolist()

flatfile = ''
with open(os.path.join(basePath, 'data', 'flatFile - backup.csv')) as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if i-1 not in indices:
        flatfile += line
        continue

    elements = line.split(',')
    elements[7] = to_fix[to_fix['index'] == i-1].iloc[0]['New event type']
    
    new_line = ','.join(elements)
    flatfile += new_line

with open(os.path.join(basePath, 'data', 'flatFile - backup.csv'), 'w', encoding='utf8') as f:
    f.write(flatfile)

with open(os.path.join(draftPath, 'flatFile.csv'), 'w', encoding='utf8') as f:
    f.write(flatfile)

df = pd.read_csv(os.path.join(basePath, 'data', 'flatFile - backup.csv'))
df.to_excel(os.path.join(draftPath, 'flatFile.xlsx'), index=False)
