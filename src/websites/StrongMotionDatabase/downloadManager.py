import os
import sys
import json
import random
import string
import datetime
import numpy as np
import pandas as pd

###################
## Bokeh modules ##
###################
from bokeh.events import ButtonClick, DocumentReady

from bokeh.io import curdoc

from bokeh.layouts import grid

from bokeh.models import CustomJS, DatePicker
from bokeh.models.widgets import Button, NumericInput, Select

from jinja2 import Environment, FileSystemLoader

# Paths
currentDir = os.path.dirname(__file__)
dataPath   = os.path.abspath(os.path.join(currentDir, '..', '..', '..', 'data'))
srcPath    = os.path.abspath(os.path.join(currentDir, '..', '..'))

################
##  Settings  ##
################
pwdith  = 1
pheight = 250
formats = ['Matlab (MAT-binary, *.mat)', 'Python (Numpy, *.npz)']
lettersandnumbers  = string.ascii_lowercase + string.digits

#########################
##  Filter components  ##
#########################
filter_since = DatePicker(title='Show events since', sizing_mode='stretch_width', width=pwdith)

filter_until = DatePicker(title='Show events until', sizing_mode='stretch_width', width=pwdith)

filter_eType = Select(title='Event type',sizing_mode='stretch_width', width=pwdith)

filter_minMw = NumericInput(title='Show events larger or equal than', mode='float', sizing_mode='stretch_width', width=pwdith)

filter_maxMw = NumericInput(title='Show events smaller or equal than', mode='float', sizing_mode='stretch_width', width=pwdith)

filter_sCode = Select(title='Recorded by station', sizing_mode='stretch_width', width=pwdith)

###############
##  Buttons  ##
###############
button_filter   = Button(label='Apply filters', button_type='primary',
    sizing_mode='stretch_width', width=pwdith, align='end')

button_download = Button(label='Download', button_type='success',
    sizing_mode='stretch_width', width=pwdith, align='end')

####################################
##  Custom Javascript Functions   ##
####################################

# Para chequear si se seleccionaron todas las filas
# table.rows({selected: true})[0].length == n

# # Para obtener los ids de las seleccionadas
# for i in range(m):
#     $('#earthquakes').DataTable().rows('.selected').data()[i][7]

Download = CustomJS(args=dict(n=0, extension='mat'),code="""
var table = $('#earthquakes').DataTable();
var events = [];
var nSelected = table.rows({selected: true})[0].length;
if (nSelected == 0) {
    alert("You must select at least one event");
    return;
} else if (nSelected == n) {
    for (let i=0; i < n; i++) {
        events.push(table.rows(i).data()[0][7]);
    }
    table.rows().data()
} else {
    var selected = table.rows('.selected').data();
    for (let i=0; i < nSelected; i++) {
        events.push(selected[i][7]);
    }
}

var url = '';
var path = 'data/seismicDatabase/';
for (let i=0; i < nSelected; i++) {
    url = '%s' + path + extension + '/' + events[i] + '.' + extension;
    console.log(url);
}
"""%sys.argv[1])

CreateTable = CustomJS(args=dict(data=[]), code="""
var table = document.querySelector('#earthquakes');
new DataTable(table, {
    searching: false,
    pagingType: 'simple_numbers',
    columnDefs: [
        {
            orderable: false,
            render: DataTable.render.select(),
            targets: 0
        },
        {
            targets: [ -1 ],
            "visible": false
        }
    ],
    select: {
        style: 'multi',
        selector: 'td:first-child'
    },
    order: [[1, 'des']],
    columns: [
        { title: '' },
        { title: 'Date' },
        { title: 'Magnitude' },
        { title: 'Latitude' },
        { title: 'Longitude' },
        { title: 'Depth' },
        { title: 'Event type' },
        { title: '' }
    ],
    data: data,
});
""")

UpdateTable = CustomJS(args=dict(data=[], n=0), code="""
if (n == 0) {
    alert('Zero events meet your criteria! Please try again');
} else {
    if (n == 1) {
        var message = 'A total of 1 event meets your criteria';
    } else {
        var message = 'A total of '.concat(n.toString());
        message = message.concat(' events meet your criteria');
    }
    alert(message);
}

var table = $('#earthquakes').DataTable();
table.clear();
table.rows.add( data ).draw();
""")

########################
##  Filter functions  ##
########################
def filter_events():
    global flatfile

    since = filter_since.value
    since = datetime.date(int(since[:4]), int(since[5:7]), int(since[8:]))
    
    until = filter_until.value
    until = datetime.date(int(until[:4]), int(until[5:7]), int(until[8:]))
    
    minMag = filter_minMw.value
    maxMag = filter_maxMw.value
    
    eType = filter_eType.value
    sCode = filter_sCode.value
    
    conditions = (flatfile['Magnitude [Mw]']  >= minMag) & \
                 (flatfile['Magnitude [Mw]']  <= maxMag) & \
                 (flatfile['Earthquake date'] >= since) & \
                 (flatfile['Earthquake date'] <= until)

    if eType != 'Any':
        conditions &= (flatfile['Event type'].str.lower() == eType.lower())

    if sCode != 'Any':
        conditions &= (flatfile['Station code'] == sCode)

    EarthquakeNames = flatfile[conditions]['Earthquake Name'].str.cat(flatfile[conditions]['Event type'], sep='_')
    EarthquakeNames = EarthquakeNames.drop_duplicates().values.tolist()

    rows = []
    for earthquakeName in EarthquakeNames:
        date, mag, lat, lon, depth, etype = earthquakeName.split('_')
        event_id = '_'.join([date, mag, lat, lon, depth])
        date  = '-'.join([date[:4], date[4:6], date[6:]])
        mag   = mag.replace('mag','')[:-1]
        lat   = '-' + lat.replace('lat','')[:-1]
        lon   = '-' + lon.replace('lon','')[:-1]
        depth = depth.replace('depth','')[:-2]
        etype = etype.capitalize()
        rows.append([None, date, mag, lat, lon, depth, etype, event_id])
    
    UpdateTable.args = dict(data=rows, n=len(EarthquakeNames))
    Download.args    = dict(n=len(EarthquakeNames), extension='mat')
    
    button_filter.name = ''.join(random.choice(lettersandnumbers) for i in range(10))
    
########################
## Javascript events  ##
########################
button_filter.js_on_change('name', UpdateTable)
button_filter.on_click(filter_events)
button_download.js_on_event(ButtonClick, Download)
curdoc().js_on_event(DocumentReady, CreateTable)

################
##  Database  ##
################
flatfile = pd.read_csv(os.path.join(dataPath, 'flatFile.csv'))

# Dates
flatfile['Earthquake date'] = pd.to_datetime(flatfile['Earthquake date']).dt.date
min_date = flatfile['Earthquake date'].min()
max_date = flatfile['Earthquake date'].max()

# Earthquake types
eTypes = ['Any'] + list(map(lambda x: x.capitalize(), np.sort(flatfile['Event type'].unique()).tolist()))

# Magnitudes
min_mag = flatfile['Magnitude [Mw]'].min()
max_mag = flatfile['Magnitude [Mw]'].max()

# Station codes
sCodes = ['Any'] + np.sort(flatfile['Station code'].unique()).tolist()

# Seismic events
seismic_events = list(reversed(flatfile['Earthquake Name'].unique().tolist()))
station_codes  = np.sort(flatfile[flatfile['Earthquake Name'] == seismic_events[0]]['Station code']).tolist()

#############################
##  Initial configuration  ##
#############################
filter_since.update(min_date=min_date, max_date=max_date, value=min_date)
filter_until.update(min_date=min_date, max_date=max_date, value=max_date)
filter_eType.update(value=eTypes[0], options=eTypes)
filter_minMw.update(value=min_mag)
filter_maxMw.update(value=max_mag)
filter_sCode.update(value=sCodes[0], options=sCodes)

EarthquakeNames = flatfile['Earthquake Name'].str.cat(flatfile['Event type'], sep='_')
EarthquakeNames = EarthquakeNames.drop_duplicates().values.tolist()

rows = []
for earthquakeName in EarthquakeNames:
    date, mag, lat, lon, depth, etype = earthquakeName.split('_')
    event_id = '_'.join([date, mag, lat, lon, depth])
    date  = '-'.join([date[:4], date[4:6], date[6:]])
    mag   = mag.replace('mag','')[:-1]
    lat   = '-' + lat.replace('lat','')[:-1]
    lon   = '-' + lon.replace('lon','')[:-1]
    depth = depth.replace('depth','')[:-2]
    etype = etype.capitalize()
    rows.append([None, date, mag, lat, lon, depth, etype, event_id])

CreateTable.args = dict(data=rows)
Download.args    = dict(n=len(EarthquakeNames), extension='mat')

######################
##  Export website  ##
######################
_env = Environment(loader=FileSystemLoader('StrongMotionDatabase'))
FILE = _env.get_template('siberrisk_downloadmanager.html')
curdoc().template = FILE

distribution = grid([[filter_since, filter_minMw, filter_eType, button_filter],
                     [filter_until, filter_maxMw, filter_sCode, button_download]], sizing_mode='stretch_width')

curdoc().add_root(distribution)
curdoc().title = 'Strong Motion Database Download Manager' + u'\u2013' + ' SIBER-RISK'
