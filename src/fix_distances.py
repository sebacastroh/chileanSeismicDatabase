import os
import sys
import pickle
import pyperclip
import webbrowser
import subprocess
import pandas as pd

basePath = os.path.abspath('.')
dataPath = os.path.abspath(os.path.join('.', '..', 'data'))

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

filtered[filtered['Magnitude [Mw]'] >= 7.1]['Earthquake Name'].unique()
