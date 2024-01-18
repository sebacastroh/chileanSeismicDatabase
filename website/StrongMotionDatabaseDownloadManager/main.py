import os
import random
import string
import datetime
import subprocess
import numpy as np
import pandas as pd
import zipfile as zp

import json
from jinja2 import Environment, PackageLoader, Markup

from bokeh.io import curdoc
from bokeh.models import CheckboxGroup, ColumnDataSource, CustomJS, DatePicker
from bokeh.models.widgets import TextInput, Button, PreText, Select
from bokeh.layouts import grid

path = os.getcwd()
cwd = os.path.join(path, 'StrongMotionDatabaseDownloadManager')
database_path = '/home/srcastro/projects/correctSeismicDatabase/newMethodology/events_fType_corrected_v2'

with open('/home/siberrisk/private.txt', 'r') as fopen:
    password = fopen.readline().strip()

echo = subprocess.Popen(['echo', password], stdout=subprocess.PIPE)
command = 'sudo -S python ' + os.path.join(cwd, 'remove_file.py') + ' '

lettersandnumbers = string.ascii_lowercase + string.digits

pwdith = 1

#filenames = os.listdir(database_path.replace('#folder#', 'mat')
users_registered = '/home/srcastro/projects/strongMotionDatabaseWeb/seismic_db_users.csv'

df = pd.read_csv(users_registered, encoding='latin1')
users = list(df['email'])
passwords = list(df['password'])

flatfile = '/home/srcastro/projects/correctSeismicDatabase/newMethodology/flatFile_corrected_v2.csv'
tdata = pd.read_csv(flatfile, index_col=0)

events, indices = np.unique(tdata['Earthquake Name'], return_index=True)
dates = pd.to_datetime(tdata.iloc[indices,1].str.split().map(lambda lst : lst[0][:10])).dt.date
magnitudes = tdata.iloc[indices,2]
eventTypes = tdata.iloc[indices,5]

N = len(events)

eventsSince = DatePicker(min_date=dates.iloc[0],
                         max_date=dates.iloc[-1],
                         value=dates.iloc[0],
                         title="Since", sizing_mode='stretch_width', width=pwdith)

eventsUntil = DatePicker(min_date=dates.iloc[0],
                         max_date=dates.iloc[-1],
                         value=dates.iloc[-1],
                         title="Until", sizing_mode='stretch_width', width=pwdith)

minMagnitude = TextInput(title="Larger or equal than", value="4.0", sizing_mode='stretch_width', width=pwdith)

maxMagnitude = TextInput(title="Smaller or equal than", value="9.0", sizing_mode='stretch_width', width=pwdith)

eventType    = Select(title='Event type', value='All', options=['All', 'Crustal', 'Interface', 'Intraslab'], sizing_mode='stretch_width', width=pwdith)

fileType = Select(title='File format', value='Matlab (MAT-binary, *.mat)', options=['Matlab (MAT-binary, *.mat)', 'Python (pickle, *.pkl)'], sizing_mode='stretch_width', width=pwdith)

checkbox_all = CheckboxGroup(labels=['Select all'])

filter_button = Button(label="Apply filters", button_type="success", sizing_mode='stretch_width', width=pwdith, align='end')

download_button = Button(label="Download selected", button_type="success", sizing_mode='stretch_width', width=pwdith, align='end')

status = PreText(text='', background='white', sizing_mode='stretch_width', width=pwdith)

checkbox_group = CheckboxGroup(labels=events.tolist(), height=int(18.4*len(events))+1, height_policy='min')

def _selectAll(attrname, old, new):
    if checkbox_all.active:
        checkbox_group.active = list(range(len(checkbox_group.labels)))
    else:
        checkbox_group.active = []

def _updateList(event):
    since = eventsSince.value
    since = datetime.date(int(since[:4]), int(since[5:7]), int(since[8:]))
    until = eventsUntil.value
    until = datetime.date(int(until[:4]), int(until[5:7]), int(until[8:]))
    minMag = float(minMagnitude.value)
    maxMag = float(maxMagnitude.value)
    
    eType = eventType.value
    if eType == 'All':
        pos = np.where((dates >= since) & (dates <= until) \
                       & (magnitudes >= minMag) & (magnitudes <= maxMag))
    else:
        eType = eType.lower()
        pos = np.where((dates >= since) & (dates <= until) \
                       & (magnitudes >= minMag) & (magnitudes <= maxMag) \
                       & (eventTypes==eType))

    new_labels = events[pos].tolist()

    checkbox_all.active = []
    checkbox_group.active = []
    checkbox_group.labels = new_labels
    checkbox_group.height = int(18.4*len(new_labels))+1
    
def _downloadList(attr, old, new):
    global status, downloadCode

    if status.text == 'Preparing download...':
        n = len(checkbox_group.active)

        if fileType.value == 'Matlab (MAT-binary, *.mat)':
            fType = 'mat'
        else:
            fType = 'pkl'
        
        if n == 0:
            downloadCode.args = dict(n=0, filename='', text='You must select at least one file', download=0, path='')
            status.text = 'You must select at least one file'
        elif n == 1:
            filename  = checkbox_group.labels[checkbox_group.active[0]] + '.' + fType
            dpath = 'events_' + fType + '_corrected_v2/'
            downloadCode.args = dict(n=1, filename=filename, text='Download in progress', download=1, path=dpath)
            status.text = 'Download in progress'
        elif n > 1 and n < N:
            files_to_download = [checkbox_group.labels[i] for i in checkbox_group.active]
            dpath = database_path.replace('fType', fType)
            random_name = ''.join(random.choice(lettersandnumbers) for i in range(10))
            zipname = 'SIBER-RISK_database_' + random_name + '.zip'
            with zp.ZipFile(os.path.join(path, 'tmp', zipname), 'w', zp.ZIP_DEFLATED) as Zp:
                for filename in files_to_download:
                    Zp.write(os.path.join(dpath, filename) + '.' + fType, filename + '.' + fType)
            downloadCode.args = dict(n=n, filename=zipname, text='Download in progress', download=1, path='')
            status.text = 'Download in progress'
            this_command = command + os.path.join(path, 'tmp', zipname)
            subprocess.Popen(this_command, shell=True, stdin=echo.stdout, stdout=subprocess.PIPE)
        elif n == N:
            if fType == 'mat':
                zipname = 'SIBER-RISK_database_matlab.zip'
            else:
                zipname = 'SIBER-RISK_database_python.zip'
            downloadCode.args = dict(n=N, filename=zipname, text='Download in progress', download=1, path='')
            status.text = 'Download in progress'

        downloadCode.args = dict(n=0, filename='', text='', download=0, path='')


validated = ColumnDataSource(data=dict(val=[0]))

downloadCode = CustomJS(args=dict(n=0, filename='', text='', download=0, path=''), code="""
if (download == 1) {
    var link = document.createElement('a');
    link.setAttribute('download', filename);
    if (n == 1) {
        link.href = path + filename;
    } else {
        link.href = 'tmp/'.concat(filename);
    }
    document.body.appendChild(link);
    link.click();
    link.remove();
}""")

status.js_on_change('text', downloadCode)
status.on_change('text', _downloadList)

login = CustomJS(args=dict(validated=validated, status=status, users=users, passwords=passwords),\
        code="""
const data = validated.data;

if (status.text == 'Preparing download...') {
    alert('A download is being prepared, please wait until is finished');
} else {
    if (data['val'][0] == 1) {
        status.text = 'Preparing download...';
        status.change.emit();
    } else {
        var email = prompt('E-mail:', '');
        var pass  = prompt('Code:', '');
        
        for (var i = 0; i < users.length; i++) {
            if (email == users[i] && pass == passwords[i]) {
                data['val'][0] = 1;
                status.text = 'Preparing download...';
                status.change.emit();
                validated.change.emit();
                break;
            }
        }
        if (data['val'][0] == 0) {
            alert('Username and/or password wrong!');
            status.text = 'Invalid username and/or password';
            status.change.emit();
        }
    }
}""")

download_button.js_on_click(login)


checkbox_all.on_change('active', _selectAll)
filter_button.on_click(_updateList)

distribution = grid([[eventsSince, minMagnitude, eventType, filter_button],
                     [eventsUntil, maxMagnitude, fileType, download_button],
                     [checkbox_all, None, None, status],
                     [checkbox_group]], sizing_mode='stretch_width')


# Export website
_env = Environment(loader=PackageLoader('bokeh.core', '_templates'))
_env.filters['json'] = lambda obj: Markup(json.dumps(obj))

FILE = _env.get_template("siberrisk_downloadmanager.html")

curdoc().template = FILE

curdoc().add_root(distribution)
curdoc().title = 'Strong Motion Database Download Manager ' + u'\u2013' + ' SIBER-RISK'
