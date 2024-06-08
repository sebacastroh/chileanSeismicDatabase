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
from bokeh.models.widgets import Button, Div, NumericInput, RadioGroup, Select

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

############
##  Divs  ##
############
div_format   = Div(text='<h2>1. Choose the format of the files</h2>', sizing_mode='stretch_width', width=pwdith)

div_filter   = Div(text='<h2>2. Filter the events table (optional)</h2>', sizing_mode='stretch_width', width=pwdith)

div_select   = Div(text='<h2>3. Select the events to download</h2>', sizing_mode='stretch_width', width=pwdith)

div_download = Div(text='<h2>4. Download your events</h2>', sizing_mode='stretch_width', width=pwdith)

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
button_format   = RadioGroup(labels=formats)

button_filter   = Button(label='Apply filters', button_type='primary',
    sizing_mode='stretch_width', width=pwdith, align='end')

button_download = Button(label='Download', button_type='success',
    width=pwdith, align='end', name='download', styles={'width': '10%', 'margin-left': '45%'})

####################################
##  Custom Javascript Functions   ##
####################################
Download = CustomJS(args=dict(n=0, button_format=button_format, sizes_npz={}, sizes_mat={}),code="""
async function downloadFiles() {

    if (button_format.active == null) {
        alert("You must choose a format to download");
        return;
    } else if (button_format.active == 0) {
        var extension = 'mat';
    } else {
        var extension = 'npz';
    }

    var table = $('#earthquakes').DataTable();
    var events = [];
    var nSelected = table.rows({selected: true})[0].length;
    if (nSelected == 0) {
        alert("You must select at least one event");
        return;
    } else if (nSelected == n) {
        for (let i=0; i < n; i++) {
            events.push(table.rows(i).data()[0][8]);
        }
        table.rows().data()
    } else {
        var selected = table.rows('.selected').data();
        for (let i=0; i < nSelected; i++) {
            events.push(selected[i][8]);
        }
    }

    // Get modal and elements
    var button = document.getElementsByClassName("closeDownloadModal")[0];
    button.style.display = "none";

    var progress = document.getElementsByClassName("downloadProgress")[0];
    progress.textContent = "Downloading event 0" + " of " + nSelected.toString();

    var modal = document.getElementById("myModal");
    modal.style.display = "block";

    var url = '';
    var path = 'data/seismicDatabase/';

    var totalSize = 0;
    for (let i=0; i < nSelected; i++) {
        if (button_format.active == 0) {
            totalSize += sizes_mat[events[i]];
        } else {
            totalSize += sizes_npz[events[i]];
        }
    }

    var partial_size = 0;
    for (let i=0; i < nSelected; i++) {
        if (button_format.active == 0) {
            partial_size += sizes_mat[events[i]];
        } else {
            partial_size += sizes_npz[events[i]];
        }
        progress.textContent = "Downloading event " + (i+1).toString() + " of " + nSelected.toString() + " (" + partial_size.toString() + "MB / " + total_size.toString() + "MB)";
        url = '%s' + path + extension + '/' + events[i] + '.' + extension;
        let data = await fetch(url);
        let content = await data.blob();
        let a = document.createElement('a');
        a.href = window.URL.createObjectURL(content);
        a.download = events[i] + '.' + extension;
        a.click();
        a.remove();
    }
    progress.textContent = "Download complete";

    button.onclick = function() {
      modal.style.display = "none";
    }
    button.style.display = "block";
}
downloadFiles();
"""%sys.argv[1])

CreateTable = CustomJS(args=dict(data=[]), code="""
var table = document.querySelector('#earthquakes');
new DataTable(table, {
    searching: false,
    pagingType: 'simple_numbers',
    columnDefs: [
        {
            width: '2%',
            orderable: false,
            render: DataTable.render.select(),
            targets: 0
        },
        {
            className: "dt-type-date",
            targets: 7
        },
        {
            width: '13.5%',
            targets: [1, 2, 3, 4, 5, 6]
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
        { title: 'Last update' },
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
    LastUpdates     = flatfile[conditions].groupby(['Earthquake Name'])['Last update'].max().tolist()

    rows = []
    for earthquakeName, lastUpdate in zip(EarthquakeNames, LastUpdates):
        date, mag, lat, lon, depth, etype = earthquakeName.split('_')
        event_id = '_'.join([date, mag, lat, lon, depth])
        date  = '-'.join([date[:4], date[4:6], date[6:]])
        mag   = mag.replace('mag','')[:-1]
        lat   = ('-' + lat).replace('-lat','')[:-1]
        lon   = ('-' + lon).replace('-lon','')[:-1]
        depth = depth.replace('depth','')[:-2]
        etype = etype.capitalize()
        rows.append([None, date, mag, lat, lon, depth, etype, lastUpdate, event_id])
    
    UpdateTable.args = dict(data=rows, n=len(EarthquakeNames))
    Download.args    = dict(n=len(EarthquakeNames), button_format=button_format, sizes_npz=sizes_npz, sizes_mat=sizes_mat)
    
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
LastUpdates     = flatfile.groupby(['Earthquake Name'])['Last update'].max().tolist()

rows = []
for earthquakeName, lastUpdate in zip(EarthquakeNames, LastUpdates):
    date, mag, lat, lon, depth, etype = earthquakeName.split('_')
    event_id = '_'.join([date, mag, lat, lon, depth])
    date  = '-'.join([date[:4], date[4:6], date[6:]])
    mag   = mag.replace('mag','')[:-1]
    lat   = ('-' + lat).replace('-lat','')[:-1]
    lon   = ('-' + lon).replace('-lon','')[:-1]
    depth = depth.replace('depth','')[:-2]
    etype = etype.capitalize()
    rows.append([None, date, mag, lat, lon, depth, etype, lastUpdate, event_id])

sizes_npz = {earthquakeName: os.path.getsize(os.path.join(dataPath, 'seismicDatabase', 'npz', earthquakeName + '.npz'))/1024**2 for earthquakeName in EarthquakeNames}
sizes_mat = {earthquakeName: os.path.getsize(os.path.join(dataPath, 'seismicDatabase', 'mat', earthquakeName + '.mat'))/1024**2 for earthquakeName in EarthquakeNames}

CreateTable.args = dict(data=rows)
Download.args    = dict(n=len(EarthquakeNames), button_format=button_format, sizes_npz=sizes_npz, sizes_mat=sizes_mat)

######################
##  Export website  ##
######################
_env = Environment(loader=FileSystemLoader('StrongMotionDatabase'))
FILE = _env.get_template('siberrisk_downloadmanager.html')
curdoc().template = FILE

distribution = grid([[filter_since, filter_minMw, filter_eType, button_filter],
                     [filter_until, filter_maxMw, filter_sCode, None]], sizing_mode='stretch_width')

curdoc().add_root(div_format)
curdoc().add_root(button_format)
curdoc().add_root(div_filter)
curdoc().add_root(distribution)
curdoc().add_root(div_select)
curdoc().add_root(div_download)
curdoc().add_root(button_download)
curdoc().title = 'Strong Motion Database Download Manager' + u'\u2013' + ' SIBER-RISK'
