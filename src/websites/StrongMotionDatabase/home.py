import os
import sys
import json
from jinja2 import Environment, FileSystemLoader
import jinja2
Markup = jinja2.utils.markupsafe.Markup

# Bokeh libraries
import bokeh
from bokeh.plotting import figure
from bokeh.layouts import layout, column
from bokeh.models import ColumnDataSource, CustomJS, DatePicker, Tabs, TabPanel
from bokeh.models.widgets import TextInput, Select, DataTable, TableColumn, Button, Div, CheckboxGroup, Toggle
from bokeh.tile_providers import Vendors
from bokeh.io import curdoc
from bokeh.palettes import Spectral6
from bokeh.events import ButtonClick
from datetime import datetime, timedelta, date

# Paths
currentDir = os.path.dirname(__file__)
libPath    = os.path.abspath(os.path.join(currentDir, '..', '..', 'lib'))
dataPath   = os.path.abspath(os.path.join(currentDir, '..', '..', 'data'))

if not libPath in sys.path:
    sys.path.append(libPath)

# Seismic libraries
import seismic

# System libraries
import subprocess

# Numerical libraries
import scipy.fftpack as spf
import scipy.integrate as spi
import scipy.io as spio
import numpy as np
import pyproj

# User and password libraries
import pandas as pd
import random
import string
import unicodedata
import datetime

lettersandnumbers = string.ascii_lowercase + string.digits

# Mail libraries
#import smtplib
#from email.mime.multipart import MIMEMultipart
#from email.mime.text import MIMEText

# Paths
matfiles = os.path.join(dataPath, 'seismicDatabase', 'mat')
flatfile = os.path.join(dataPath, 'flatFile.csv') #'/home/srcastro/projects/correctSeismicDatabase/newMethodology/flatFile_corrected_v2.csv'
users_registered = os.path.join(dataPath, 'seismic_db_users.csv')

#%% Read database
# df = pd.read_csv(users_registered, encoding='latin1')
df = pd.DataFrame([], columns=['email', 'password'])
users = list(df['email'])
passwords = list(df['password'])

#%% Website
g = 9.81 # m/s**2
npoints = 250
allowed = False

pwdith = 1

nTheta = 91
nx = 0
ny = 45

#%% List of events and stations
tdata = pd.read_csv(flatfile)
Events    = sorted(os.listdir(matfiles), reverse=True)
Events    = [event[:-4] for event in Events] 

events    = Select(title='Seismic events', value=Events[0], options=Events, sizing_mode='stretch_width', width=pwdith*2)
event     = events.value

Stations  = spio.loadmat(os.path.join(matfiles, event) + '.mat', struct_as_record=False, squeeze_me=True)
Stations.pop('__version__')
Stations.pop('__header__')
Stations.pop('__globals__')
keyValues = [item for item in Stations.items()]

keys   = []
snames = []
for key, value in keyValues:
    if not key.startswith('st'):
        continue
    keys.append(key)
    snames.append(value.station_code)

order  = np.argsort(snames)
snames = np.array(snames)[order].tolist()
keys   = np.array(keys)[order].tolist()

stations  = Select(title='Stations\n', value=snames[0], options=snames, sizing_mode='stretch_width', width=pwdith)
Station   = Stations[keys[0]]

#%% Filter options

year = int(event[:4])
month = int(event[4:6])
day = int(event[6:8])

one_day = timedelta(days=1)

if bokeh.__version__[0] == '1':
    eventsSince = DatePicker(min_date=datetime(1985,  3,  3) + one_day,
            max_date=datetime(year, month, day) + one_day,
            value=datetime(1985, 3, 3) + one_day,
                     title="Show events since", sizing_mode='stretch_width', width=pwdith)

    eventsUntil = DatePicker(min_date=datetime(1985,  3,  3) + one_day,
            max_date=datetime(year, month, day) + one_day,
            value=datetime(year, month, day) + one_day,
                     title="Show events until", sizing_mode='stretch_width', width=pwdith)
else:
    eventsSince = DatePicker(min_date=date(1985,  3,  3),
            max_date=date(year, month, day),
            value=date(1985, 3, 3),
                     title="Show events since", sizing_mode='stretch_width', width=pwdith)
    eventsUntil = DatePicker(min_date=date(1985,  3,  3),
            max_date=date(year, month, day),
            value=date(year, month, day),
                     title="Show events until", sizing_mode='stretch_width', width=pwdith)

minMagnitude = TextInput(title="Show events larger or equal than", value="4.0", sizing_mode='stretch_width', width=pwdith)

maxMagnitude = TextInput(title="Show events smaller or equal than", value="9.0", sizing_mode='stretch_width', width=pwdith)

eventType    = Select(title='Event type', value='All', options=['All', 'Crustal', 'Interface', 'Intraslab'], sizing_mode='stretch_width', width=pwdith*2)

stations_codes = ['Any'] + list(sorted(tdata['Station name'].unique()))
stationSelect = Select(title='Recorded by station', value='Any', options=stations_codes, sizing_mode='stretch_width', width=pwdith*2)

fileType = Select(title='File format', value='Matlab (MAT-binary, *.mat)', options=['Matlab (MAT-binary, *.mat)', 'Python (pickle, *.pkl)'], sizing_mode='stretch_width', width=pwdith)

dates = [this_date[:10] for this_date in tdata['Start time record']]
dates = pd.to_datetime(dates).date

def update_inputs(attr, old, new):
    global users, events, passwords
    download_button_callback.args = dict(validated=validated, source=events, users=users, passwords=passwords, extension=fileType)

fileType.on_change('value', update_inputs)


def update_events():
    since = eventsSince.value
    since = date(int(since[:4]), int(since[5:7]), int(since[8:]))
    
    until = eventsUntil.value
    until = date(int(until[:4]), int(until[5:7]), int(until[8:]))
    
    minMag = float(minMagnitude.value)
    maxMag = float(maxMagnitude.value)
    
    eType = eventType.value
    station_name = stationSelect.value
    
    if eType == 'All' and station_name == 'Any':
        EarthquakeNames = tdata[(tdata['Magnitude [Mw]'] >= minMag) & \
                                (tdata['Magnitude [Mw]'] <= maxMag) & \
                                (dates >= since) & \
                                (dates <= until)]['Earthquake Name']
    elif eType != 'All' and station_name == 'Any':
        eType = eType.lower()
        EarthquakeNames = tdata[(tdata['Magnitude [Mw]'] >= minMag) & \
                                (tdata['Magnitude [Mw]'] <= maxMag) & \
                                (dates >= since) & \
                                (dates <= until) & \
                                (tdata['Event type'] == eType)]['Earthquake Name']
    elif eType != 'All' and station_name != 'Any':
        eType = eType.lower()
        EarthquakeNames = tdata[(tdata['Magnitude [Mw]'] >= minMag) & \
                                (tdata['Magnitude [Mw]'] <= maxMag) & \
                                (dates >= since) & \
                                (dates <= until) & \
                                (tdata['Event type'] == eType) & \
                                (tdata['Station name'] == station_name)]['Earthquake Name']
    else:
        EarthquakeNames = tdata[(tdata['Magnitude [Mw]'] >= minMag) & \
                                (tdata['Magnitude [Mw]'] <= maxMag) & \
                                (dates >= since) & \
                                (dates <= until) & \
                                (tdata['Station name'] == station_name)]['Earthquake Name']

        
    EarthquakeNames = sorted(EarthquakeNames.drop_duplicates().values.tolist(), reverse=True)
    
    Alert.args = dict(n=len(EarthquakeNames))
    
    filterButton.name = ''.join(random.choice(lettersandnumbers) for i in range(10))
    if len(EarthquakeNames) > 0:
        events.options = EarthquakeNames
        events.value = EarthquakeNames[0]

filterButton = Button(label="Apply filters", button_type="success", sizing_mode='stretch_width', width=pwdith, align='end')
filterButton.on_click(update_events)

Alert = CustomJS(args=dict(n=0), code="""
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
}""")
filterButton.js_on_change('name', Alert)

#%% Plot options
tooltips_a  = [("(t, a)", "($x, $y)"),]
tooltips_v  = [("(t, v)", "($x, $y)"),]
tooltips_d  = [("(t, d)", "($x, $y)"),]
tooltips_sa = [("(Tn, Sa)", "($x, $y)"),]
tooltips_sv = [("(Tn, Sv)", "($x, $y)"),]
tooltips_f  = [("(f, A)", "($x, $y)"),]
tooltips_ai = [("(t, IA)", "($x, $y)"),]
tooltips_cav = [("(t, CAV)", "($x, $y)"),]

plot_height = 250
#%% Acceleration
# ax
tax = np.linspace(0., (len(Station.acc_uncorrected_1)-1)*Station.dt, len(Station.acc_uncorrected_1))
source_ax = ColumnDataSource(data=dict(x=tax, y=Station.acc_uncorrected_1/g))

ax = figure(title='Acceleration', x_axis_label='Time [s]', y_axis_label='Acceleration [g]',
            tooltips=tooltips_a, height=plot_height, sizing_mode='stretch_width', width=pwdith*4)
ax.line('x', 'y', legend_label=Station.station_name, source=source_ax)
ax.xgrid[0].minor_grid_line_color = '#e5e5e5'
ax.ygrid[0].minor_grid_line_color = '#e5e5e5'

# ay
tay = np.linspace(0., (len(Station.acc_uncorrected_2)-1)*Station.dt, len(Station.acc_uncorrected_2))
source_ay = ColumnDataSource(data=dict(x=tay, y=Station.acc_uncorrected_2/g))

ay = figure(title='Acceleration', x_axis_label='Time [s]', y_axis_label='Acceleration [g]',
            tooltips=tooltips_a, height=plot_height, sizing_mode='stretch_width', width=pwdith*4)
ay.line('x', 'y', legend_label=Station.station_name, source=source_ay)
ay.xgrid[0].minor_grid_line_color = '#e5e5e5'
ay.ygrid[0].minor_grid_line_color = '#e5e5e5'

# az
taz = np.linspace(0., (len(Station.acc_uncorrected_3)-1)*Station.dt, len(Station.acc_uncorrected_3))
source_az = ColumnDataSource(data=dict(x=taz, y=Station.acc_uncorrected_3/g))

az = figure(title='Acceleration', x_axis_label='Time [s]', y_axis_label='Acceleration [g]',
            tooltips=tooltips_a, height=plot_height, sizing_mode='stretch_width', width=pwdith*4)
az.line('x', 'y', legend_label=Station.station_name, source=source_az)
az.xgrid[0].minor_grid_line_color = '#e5e5e5'
az.ygrid[0].minor_grid_line_color = '#e5e5e5'
#%% Velocity
# vx
velocity_x = spi.cumtrapz(Station.acc_uncorrected_1, x=tax, initial=0.)
source_vx = ColumnDataSource(data=dict(x=tax, y=velocity_x))

vx = figure(title='Velocity', x_axis_label='Time [s]', y_axis_label='Velocity [m/s]',
            tooltips=tooltips_v, height=plot_height, sizing_mode='stretch_width', width=pwdith*4)
vx.line('x', 'y', legend_label=Station.station_name, source=source_vx)
vx.xgrid[0].minor_grid_line_color = '#e5e5e5'
vx.ygrid[0].minor_grid_line_color = '#e5e5e5'

# vy
velocity_y = spi.cumtrapz(Station.acc_uncorrected_2, x=tay, initial=0.)
source_vy = ColumnDataSource(data=dict(x=tay, y=velocity_y))

vy = figure(title='Velocity', x_axis_label='Time [s]', y_axis_label='Velocity [m/s]',
            tooltips=tooltips_v, height=plot_height, sizing_mode='stretch_width', width=pwdith*4)
vy.line('x', 'y', legend_label=Station.station_name, source=source_vy)
vy.xgrid[0].minor_grid_line_color = '#e5e5e5'
vy.ygrid[0].minor_grid_line_color = '#e5e5e5'

# vz
velocity_z = spi.cumtrapz(Station.acc_uncorrected_3, x=taz, initial=0.)
source_vz = ColumnDataSource(data=dict(x=taz, y=velocity_z))

vz = figure(title='Velocity', x_axis_label='Time [s]', y_axis_label='Velocity [m/s]',
            tooltips=tooltips_v, height=plot_height, sizing_mode='stretch_width', width=pwdith*4)
vz.line('x', 'y', legend_label=Station.station_name, source=source_vz)
vz.xgrid[0].minor_grid_line_color = '#e5e5e5'
vz.ygrid[0].minor_grid_line_color = '#e5e5e5'
#%% Displacement
# dx
displacement_x = spi.cumtrapz(velocity_x, x=tax, initial=0.)
source_dx = ColumnDataSource(data=dict(x=tax, y=displacement_x))

dx = figure(title='Displacement', x_axis_label='Time [s]', y_axis_label='Displacement [m]',
            tooltips=tooltips_d, height=plot_height, sizing_mode='stretch_width', width=pwdith*4)
dx.line('x', 'y', legend_label=Station.station_name, source=source_dx)
dx.xgrid[0].minor_grid_line_color = '#e5e5e5'
dx.ygrid[0].minor_grid_line_color = '#e5e5e5'

# dy
displacement_y = spi.cumtrapz(velocity_y, x=tay, initial=0.)
source_dy = ColumnDataSource(data=dict(x=tay, y=displacement_y))

dy = figure(title='Displacement', x_axis_label='Time [s]', y_axis_label='Displacement [m]',
            tooltips=tooltips_d, height=plot_height, sizing_mode='stretch_width', width=pwdith*4)
dy.line('x', 'y', legend_label=Station.station_name, source=source_dy)
dy.xgrid[0].minor_grid_line_color = '#e5e5e5'
dy.ygrid[0].minor_grid_line_color = '#e5e5e5'

# dz
displacement_z = spi.cumtrapz(velocity_z, x=taz, initial=0.)
source_dz = ColumnDataSource(data=dict(x=taz, y=displacement_z))

dz = figure(title='Displacement', x_axis_label='Time [s]', y_axis_label='Displacement [m]',
            tooltips=tooltips_d, height=plot_height, sizing_mode='stretch_width', width=pwdith*4)
dz.line('x', 'y', legend_label=Station.station_name, source=source_dz)
dz.xgrid[0].minor_grid_line_color = '#e5e5e5'
dz.ygrid[0].minor_grid_line_color = '#e5e5e5'

panels_stations = [    TabPanel(child=layout([[ax], [vx], [dx]], sizing_mode='stretch_width', width=pwdith*4), title=Station.component_1),
                    TabPanel(child=layout([[ay], [vy], [dy]], sizing_mode='stretch_width', width=pwdith*4), title=Station.component_2),
                    TabPanel(child=layout([[az], [vz], [dz]], sizing_mode='stretch_width', width=pwdith*4), title=Station.component_3)]

tabs_acceleration = Tabs(tabs=panels_stations, sizing_mode='stretch_width', width=pwdith*4)
#%% Inputs
xi = TextInput(title='Damping ratio', value='0.05', sizing_mode='stretch_width', width=pwdith)
ta = TextInput(title='Initial period [s]', value='0.0', sizing_mode='stretch_width', width=pwdith)
tb = TextInput(title='Ending period [s]', value='1.5', sizing_mode='stretch_width', width=pwdith)

checkbox_axis = CheckboxGroup(labels=["Logarithmic X axis",
                          "Logarithmic Y axis"], sizing_mode='stretch_width', width=pwdith)
checkbox_grid = CheckboxGroup(labels=["X axis minor grid",
                          "Y axis minor grid"], sizing_mode='stretch_width', width=pwdith)

Tn_lin = np.linspace(float(ta.value), float(tb.value), npoints)
Tn_log = np.hstack((0., np.logspace(np.log10(0.001), np.log10(float(tb.value)), npoints)))

#%% Spectrum
Spectrum_lin = seismic.SpectraRot(Station.acc_uncorrected_1, Station.acc_uncorrected_2, Station.dt, Tn_lin, float(xi.value), nTheta)
Spectrum_linx = Spectrum_lin[nx]
Spectrum_liny = Spectrum_lin[ny]

source_Sa1_lin = ColumnDataSource(data=dict(x=Tn_lin, y=Spectrum_linx/g))
source_Sa2_lin = ColumnDataSource(data=dict(x=Tn_lin, y=Spectrum_liny/g))
source_SaRotD50_lin = ColumnDataSource(data=dict(x=Tn_lin, y=np.median(Spectrum_lin, axis=0)/g))
source_SaRotD0_lin = ColumnDataSource(data=dict(x=Tn_lin, y=np.min(Spectrum_lin, axis=0)/g))
source_SaRotD100_lin = ColumnDataSource(data=dict(x=Tn_lin, y=np.max(Spectrum_lin, axis=0)/g))
source_Saxy_lin = ColumnDataSource(data=dict(x=Tn_lin, y=np.sqrt(Spectrum_linx*Spectrum_liny)/g))

Spectrum_log = seismic.SpectraRot(Station.acc_uncorrected_1, Station.acc_uncorrected_2, Station.dt, Tn_log, float(xi.value), nTheta)
Spectrum_logx = Spectrum_log[nx]
Spectrum_logy = Spectrum_log[ny]

source_Sa1_log = ColumnDataSource(data=dict(x=Tn_log, y=Spectrum_logx/g))
source_Sa2_log = ColumnDataSource(data=dict(x=Tn_log, y=Spectrum_logy/g))
source_SaRotD50_log = ColumnDataSource(data=dict(x=Tn_log, y=np.median(Spectrum_log, axis=0)/g))
source_SaRotD0_log = ColumnDataSource(data=dict(x=Tn_log, y=np.min(Spectrum_log, axis=0)/g))
source_SaRotD100_log = ColumnDataSource(data=dict(x=Tn_log, y=np.max(Spectrum_log, axis=0)/g))
source_Saxy_log = ColumnDataSource(data=dict(x=Tn_log, y=np.sqrt(Spectrum_logx*Spectrum_logy)/g))

Sa_linx_liny = figure(title="Pseudo-acceleration Response Spectrum",
            x_axis_label='Period [s]', y_axis_label='Spectral Acceleration [g]',
            tooltips=tooltips_sa, sizing_mode='stretch_width', width=pwdith*3)

Sa_linx_liny.line('x', 'y', source=source_Sa1_lin, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
Sa_linx_liny.line('x', 'y', source=source_Sa2_lin, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
Sa_linx_liny.line('x', 'y', source=source_SaRotD50_lin, legend_label='RotD50', line_width=2, color=Spectral6[2])
Sa_linx_liny.line('x', 'y', source=source_SaRotD0_lin, legend_label='RotD0', line_width=2, color=Spectral6[3])
Sa_linx_liny.line('x', 'y', source=source_SaRotD100_lin, legend_label='RotD100', line_width=2, color=Spectral6[4])
Sa_linx_liny.line('x', 'y', source=source_Saxy_lin, legend_label='Geometric Mean', line_width=2, color=Spectral6[5])

Sa_linx_liny.legend.click_policy = "hide"

panel_Sa_linx_liny = TabPanel(child=Sa_linx_liny, title='Sa Spectra')

Sa_linx_logy = figure(title="Pseudo-acceleration Response Spectrum",
            x_axis_label='Period [s]', y_axis_label='Spectral Acceleration [g]',
            tooltips=tooltips_sa, sizing_mode='stretch_width', y_axis_type='log', width=pwdith*3)

Sa_linx_logy.line('x', 'y', source=source_Sa1_lin, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
Sa_linx_logy.line('x', 'y', source=source_Sa2_lin, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
Sa_linx_logy.line('x', 'y', source=source_SaRotD50_lin, legend_label='RotD50', line_width=2, color=Spectral6[2])
Sa_linx_logy.line('x', 'y', source=source_SaRotD0_lin, legend_label='RotD0', line_width=2, color=Spectral6[3])
Sa_linx_logy.line('x', 'y', source=source_SaRotD100_lin, legend_label='RotD100', line_width=2, color=Spectral6[4])
Sa_linx_logy.line('x', 'y', source=source_Saxy_lin, legend_label='Geometric Mean', line_width=2, color=Spectral6[5])

Sa_linx_logy.legend.click_policy = "hide"

panel_Sa_linx_logy = TabPanel(child=Sa_linx_logy, title='Sa Spectrum')

Sa_logx_liny = figure(title="Pseudo-acceleration Response Spectrum",
            x_axis_label='Period [s]', y_axis_label='Spectral Acceleration [g]',
            tooltips=tooltips_sa, sizing_mode='stretch_width', x_axis_type='log', width=pwdith*3)

Sa_logx_liny.line('x', 'y', source=source_Sa1_log, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
Sa_logx_liny.line('x', 'y', source=source_Sa2_log, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
Sa_logx_liny.line('x', 'y', source=source_SaRotD50_log, legend_label='RotD50', line_width=2, color=Spectral6[2])
Sa_logx_liny.line('x', 'y', source=source_SaRotD0_log, legend_label='RotD0', line_width=2, color=Spectral6[3])
Sa_logx_liny.line('x', 'y', source=source_SaRotD100_log, legend_label='RotD100', line_width=2, color=Spectral6[4])
Sa_logx_liny.line('x', 'y', source=source_Saxy_log, legend_label='Geometric Mean', line_width=2, color=Spectral6[5])

Sa_logx_liny.legend.click_policy = "hide"

panel_Sa_logx_liny = TabPanel(child=Sa_logx_liny, title='Sa Spectrum')

Sa_logx_logy = figure(title="Pseudo-acceleration Response Spectrum",
            x_axis_label='Period [s]', y_axis_label='Spectral Acceleration [g]',
            tooltips=tooltips_sa, x_axis_type='log', y_axis_type='log', sizing_mode='stretch_width', width=pwdith*3)

Sa_logx_logy.line('x', 'y', source=source_Sa1_log, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
Sa_logx_logy.line('x', 'y', source=source_Sa2_log, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
Sa_logx_logy.line('x', 'y', source=source_SaRotD50_log, legend_label='RotD50', line_width=2, color=Spectral6[2])
Sa_logx_logy.line('x', 'y', source=source_SaRotD0_log, legend_label='RotD0', line_width=2, color=Spectral6[3])
Sa_logx_logy.line('x', 'y', source=source_SaRotD100_log, legend_label='RotD100', line_width=2, color=Spectral6[4])
Sa_logx_logy.line('x', 'y', source=source_Saxy_log, legend_label='Geometric Mean', line_width=2, color=Spectral6[5])

Sa_logx_logy.legend.click_policy = "hide"

panel_Sa_logx_logy = TabPanel(child=Sa_logx_logy, title='Sa Spectrum')

#%% Combined D-V-A Spectrum
source_Sv1_lin = ColumnDataSource(data=dict(x=Tn_lin, y=Spectrum_linx*Tn_lin/(2.*np.pi)))
source_Sv2_lin = ColumnDataSource(data=dict(x=Tn_lin, y=Spectrum_liny*Tn_lin/(2.*np.pi)))
source_SvRotD50_lin = ColumnDataSource(data=dict(x=Tn_lin, y=np.median(Spectrum_lin, axis=0)*Tn_lin/(2.*np.pi)))
source_SvRotD0_lin = ColumnDataSource(data=dict(x=Tn_lin, y=np.min(Spectrum_lin, axis=0)*Tn_lin/(2.*np.pi)))
source_SvRotD100_lin = ColumnDataSource(data=dict(x=Tn_lin, y=np.max(Spectrum_lin, axis=0)*Tn_lin/(2.*np.pi)))
source_Svxy_lin = ColumnDataSource(data=dict(x=Tn_lin, y=np.sqrt(Spectrum_linx*Spectrum_liny)*Tn_lin/(2.*np.pi)))

source_Sv1_log = ColumnDataSource(data=dict(x=Tn_log, y=Spectrum_logx*Tn_log/(2.*np.pi)))
source_Sv2_log = ColumnDataSource(data=dict(x=Tn_log, y=Spectrum_logy*Tn_log/(2.*np.pi)))
source_SvRotD50_log = ColumnDataSource(data=dict(x=Tn_log, y=np.median(Spectrum_log, axis=0)*Tn_log/(2.*np.pi)))
source_SvRotD0_log = ColumnDataSource(data=dict(x=Tn_log, y=np.min(Spectrum_log, axis=0)*Tn_log/(2.*np.pi)))
source_SvRotD100_log = ColumnDataSource(data=dict(x=Tn_log, y=np.max(Spectrum_log, axis=0)*Tn_log/(2.*np.pi)))
source_Svxy_log = ColumnDataSource(data=dict(x=Tn_log, y=np.sqrt(Spectrum_logx*Spectrum_logy)*Tn_log/(2.*np.pi)))

Sv_linx_liny = figure(title="Combined D-V-A Spectrum",
            x_axis_label='Period [s]', y_axis_label='Spectral Velocity [m/s]',
            tooltips=tooltips_sv, sizing_mode='stretch_width', width=pwdith*3)

#Sv_linx_liny.line('x', 'y', source=source_Sv_lin)
Sv_linx_liny.line('x', 'y', source=source_Sv1_lin, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
Sv_linx_liny.line('x', 'y', source=source_Sv2_lin, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
Sv_linx_liny.line('x', 'y', source=source_SvRotD50_lin, legend_label='RotD50', line_width=2, color=Spectral6[2])
Sv_linx_liny.line('x', 'y', source=source_SvRotD0_lin, legend_label='RotD0', line_width=2, color=Spectral6[3])
Sv_linx_liny.line('x', 'y', source=source_SvRotD100_lin, legend_label='RotD100', line_width=2, color=Spectral6[4])
Sv_linx_liny.line('x', 'y', source=source_Svxy_lin, legend_label='Geometric Mean', line_width=2, color=Spectral6[5])

Sv_linx_liny.legend.click_policy = "hide"

panel_Sv_linx_liny = TabPanel(child=Sv_linx_liny, title='D-V-A Spectrum')

Sv_linx_logy = figure(title="Combined D-V-A Spectrum",
            x_axis_label='Period [s]', y_axis_label='Spectral Velocity [m/s]',
            y_axis_type='log', tooltips=tooltips_sv, sizing_mode='stretch_width', width=pwdith*3)

Sv_linx_logy.line('x', 'y', source=source_Sv1_lin, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
Sv_linx_logy.line('x', 'y', source=source_Sv2_lin, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
Sv_linx_logy.line('x', 'y', source=source_SvRotD50_lin, legend_label='RotD50', line_width=2, color=Spectral6[2])
Sv_linx_logy.line('x', 'y', source=source_SvRotD0_lin, legend_label='RotD0', line_width=2, color=Spectral6[3])
Sv_linx_logy.line('x', 'y', source=source_SvRotD100_lin, legend_label='RotD100', line_width=2, color=Spectral6[4])
Sv_linx_logy.line('x', 'y', source=source_Svxy_lin, legend_label='Geometric Mean', line_width=2, color=Spectral6[5])

Sv_linx_liny.legend.click_policy = "hide"

#Sv_linx_logy.line('x', 'y', source=source_Sv_lin)
panel_Sv_linx_logy = TabPanel(child=Sv_linx_logy, title='D-V-A Spectrum')

Sv_logx_liny = figure(title="Combined D-V-A Spectrum",
            x_axis_label='Period [s]', y_axis_label='Spectral Velocity [m/s]',
            x_axis_type="log", tooltips=tooltips_sv, sizing_mode='stretch_width', width=pwdith*3)

#Sv_logx_liny.line('x', 'y', source=source_Sv_log)
Sv_logx_liny.line('x', 'y', source=source_Sv1_log, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
Sv_logx_liny.line('x', 'y', source=source_Sv2_log, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
Sv_logx_liny.line('x', 'y', source=source_SvRotD50_log, legend_label='RotD50', line_width=2, color=Spectral6[2])
Sv_logx_liny.line('x', 'y', source=source_SvRotD0_log, legend_label='RotD0', line_width=2, color=Spectral6[3])
Sv_logx_liny.line('x', 'y', source=source_SvRotD100_log, legend_label='RotD100', line_width=2, color=Spectral6[4])
Sv_logx_liny.line('x', 'y', source=source_Svxy_log, legend_label='Geometric Mean', line_width=2, color=Spectral6[5])

Sv_logx_liny.legend.click_policy = "hide"

panel_Sv_logx_liny = TabPanel(child=Sv_logx_liny, title='D-V-A Spectrum')

Sv_logx_logy = figure(title="Combined D-V-A Spectrum",
        x_axis_label='Period [s]', y_axis_label='Spectral Velocity [m/s]',
              x_axis_type="log", y_axis_type="log", tooltips=tooltips_sv, sizing_mode='stretch_width', width=pwdith*3)

#Sv_logx_logy.line('x', 'y', source=source_Sv_log)
Sv_logx_logy.line('x', 'y', source=source_Sv1_log, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
Sv_logx_logy.line('x', 'y', source=source_Sv2_log, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
Sv_logx_logy.line('x', 'y', source=source_SvRotD50_log, legend_label='RotD50', line_width=2, color=Spectral6[2])
Sv_logx_logy.line('x', 'y', source=source_SvRotD0_log, legend_label='RotD0', line_width=2, color=Spectral6[3])
Sv_logx_logy.line('x', 'y', source=source_SvRotD100_log, legend_label='RotD100', line_width=2, color=Spectral6[4])
Sv_logx_logy.line('x', 'y', source=source_Svxy_log, legend_label='Geometric Mean', line_width=2, color=Spectral6[5])

Sv_logx_logy.legend.click_policy = "hide"

panel_Sv_logx_logy = TabPanel(child=Sv_logx_logy, title='D-V-A Spectrum')

#%% Fourier spectrum
L = Station.acc_uncorrected_1.shape[0]
Fax = 2.*np.abs(spf.fft(Station.acc_uncorrected_1)*Station.dt)[:L//2]
Fay = 2.*np.abs(spf.fft(Station.acc_uncorrected_2)*Station.dt)[:L//2]
Faz = 2.*np.abs(spf.fft(Station.acc_uncorrected_3)*Station.dt)[:L//2]

Fs = 1./Station.dt
freq = Fs*np.linspace(0., L//2, L//2)/L

source_fx = ColumnDataSource(data=dict(x=freq, y=Fax))
source_fy = ColumnDataSource(data=dict(x=freq, y=Fay))
source_fz = ColumnDataSource(data=dict(x=freq, y=Faz))

Fourier_linx_liny = figure(title='Fourier Spectrum', x_axis_label='Frequency [Hz]',
        y_axis_label='Amplitude [m/s]', tooltips=tooltips_f,
                   sizing_mode='stretch_width', width=pwdith*3)
Fourier_linx_liny.line('x', 'y', source=source_fx, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
Fourier_linx_liny.line('x', 'y', source=source_fy, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
Fourier_linx_liny.line('x', 'y', source=source_fz, legend_label=Station.component_3, line_width=2, color=Spectral6[2])
Fourier_linx_liny.legend.click_policy = "hide"

panel_Fourier_linx_liny = TabPanel(child=Fourier_linx_liny, title='Fourier Spectrum')

Fourier_linx_logy = figure(title='Fourier Spectrum', x_axis_label='Frequency [Hz]',
        y_axis_label='Amplitude [m/s]', tooltips=tooltips_f,
                   y_axis_type='log', sizing_mode='stretch_width', width=pwdith*3)
Fourier_linx_logy.line('x', 'y', source=source_fx, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
Fourier_linx_logy.line('x', 'y', source=source_fy, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
Fourier_linx_logy.line('x', 'y', source=source_fz, legend_label=Station.component_3, line_width=2, color=Spectral6[2])
Fourier_linx_logy.legend.click_policy = "hide"

panel_Fourier_linx_logy = TabPanel(child=Fourier_linx_logy, title='Fourier Spectrum')

Fourier_logx_liny = figure(title='Fourier Spectrum', x_axis_label='Frequency [Hz]',
        y_axis_label='Amplitude [m/s]', tooltips=tooltips_f,
                   x_axis_type='log', sizing_mode='stretch_width', width=pwdith*3)
Fourier_logx_liny.line('x', 'y', source=source_fx, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
Fourier_logx_liny.line('x', 'y', source=source_fy, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
Fourier_logx_liny.line('x', 'y', source=source_fz, legend_label=Station.component_3, line_width=2, color=Spectral6[2])
Fourier_logx_liny.legend.click_policy = "hide"

panel_Fourier_logx_liny = TabPanel(child=Fourier_logx_liny, title='Fourier Spectrum')

Fourier_logx_logy = figure(title='Fourier Spectrum', x_axis_label='Frequency [Hz]',
        y_axis_label='Amplitude [m/s]', tooltips=tooltips_f,
                   x_axis_type='log', y_axis_type='log', sizing_mode='stretch_width', width=pwdith*3)
Fourier_logx_logy.line('x', 'y', source=source_fx, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
Fourier_logx_logy.line('x', 'y', source=source_fy, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
Fourier_logx_logy.line('x', 'y', source=source_fz, legend_label=Station.component_3, line_width=2, color=Spectral6[2])
Fourier_logx_logy.legend.click_policy = "hide"

panel_Fourier_logx_logy = TabPanel(child=Fourier_logx_logy, title='Fourier Spectrum')

#%% Husid plot
t = np.linspace(0., (len(Station.acc_uncorrected_1)-1)*Station.dt, len(Station.acc_uncorrected_1))
iax = np.pi/(2.*g)*spi.cumtrapz(Station.acc_uncorrected_1**2, t, initial=0.)
source_aix = ColumnDataSource(data=dict(x=t, y=iax))

t = np.linspace(0., (len(Station.acc_uncorrected_2)-1)*Station.dt, len(Station.acc_uncorrected_2))
iay = np.pi/(2.*g)*spi.cumtrapz(Station.acc_uncorrected_2**2, t, initial=0.)
source_aiy = ColumnDataSource(data=dict(x=t, y=iay))

t = np.linspace(0., (len(Station.acc_uncorrected_3)-1)*Station.dt, len(Station.acc_uncorrected_3))
iaz = np.pi/(2.*g)*spi.cumtrapz(Station.acc_uncorrected_3**2, t, initial=0.)
source_aiz = ColumnDataSource(data=dict(x=t, y=iaz))

Husid_linx_liny = figure(title='Husid Plot', x_axis_label='Time [s]',
        y_axis_label='Arias intensity [m/s]', tooltips=tooltips_ai,
                 sizing_mode='stretch_width', width=pwdith*3)
Husid_linx_liny.line('x', 'y', source=source_aix, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
Husid_linx_liny.line('x', 'y', source=source_aiy, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
Husid_linx_liny.line('x', 'y', source=source_aiz, legend_label=Station.component_3, line_width=2, color=Spectral6[2])
Husid_linx_liny.legend.click_policy = "hide"

panel_Husid_linx_liny = TabPanel(child=Husid_linx_liny, title='Husid Plot')

Husid_linx_logy = figure(title='Husid Plot', x_axis_label='Time [s]',
        y_axis_label='Arias intensity [m/s]', tooltips=tooltips_ai,
                 y_axis_type='log', sizing_mode='stretch_width', width=pwdith*3)
Husid_linx_logy.line('x', 'y', source=source_aix, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
Husid_linx_logy.line('x', 'y', source=source_aiy, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
Husid_linx_logy.line('x', 'y', source=source_aiz, legend_label=Station.component_3, line_width=2, color=Spectral6[2])
Husid_linx_logy.legend.click_policy = "hide"

panel_Husid_linx_logy = TabPanel(child=Husid_linx_logy, title='Husid Plot')

Husid_logx_liny = figure(title='Husid Plot', x_axis_label='Time [s]',
        y_axis_label='Arias intensity [m/s]', tooltips=tooltips_ai,
                 x_axis_type='log', sizing_mode='stretch_width', width=pwdith*3)
Husid_logx_liny.line('x', 'y', source=source_aix, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
Husid_logx_liny.line('x', 'y', source=source_aiy, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
Husid_logx_liny.line('x', 'y', source=source_aiz, legend_label=Station.component_3, line_width=2, color=Spectral6[2])
Husid_logx_liny.legend.click_policy = "hide"

panel_Husid_logx_liny = TabPanel(child=Husid_logx_liny, title='Husid Plot')

Husid_logx_logy = figure(title='Husid Plot', x_axis_label='Time [s]',
        y_axis_label='Arias intensity [m/s]', tooltips=tooltips_ai,
                 x_axis_type='log', y_axis_type='log', sizing_mode='stretch_width', width=pwdith*3)
Husid_logx_logy.line('x', 'y', source=source_aix, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
Husid_logx_logy.line('x', 'y', source=source_aiy, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
Husid_logx_logy.line('x', 'y', source=source_aiz, legend_label=Station.component_3, line_width=2, color=Spectral6[2])
Husid_logx_logy.legend.click_policy = "hide"

panel_Husid_logx_logy = TabPanel(child=Husid_logx_logy, title='Husid Plot')

#%% Cumulative absolute velocity plot
#t = np.linspace(0., (len(Station.acc_uncorrected_1)-1)*Station.dt, len(Station.acc_uncorrected_1))

t = np.linspace(0., (len(Station.acc_uncorrected_1)-1)*Station.dt, len(Station.acc_uncorrected_1))
cavx = spi.cumtrapz(np.abs(Station.acc_uncorrected_1), t, initial=0.)
source_cavx = ColumnDataSource(data=dict(x=t, y=cavx))

t = np.linspace(0., (len(Station.acc_uncorrected_2)-1)*Station.dt, len(Station.acc_uncorrected_2))
cavy = spi.cumtrapz(np.abs(Station.acc_uncorrected_2), t, initial=0.)
source_cavy = ColumnDataSource(data=dict(x=t, y=cavy))

t = np.linspace(0., (len(Station.acc_uncorrected_3)-1)*Station.dt, len(Station.acc_uncorrected_3))
cavz = spi.cumtrapz(np.abs(Station.acc_uncorrected_3), t, initial=0.)
source_cavz = ColumnDataSource(data=dict(x=t, y=cavz))

CAV_linx_liny = figure(title='Cumulative Absolute Velocity', x_axis_label='Time [s]',
        y_axis_label='CAV [m/s]', tooltips=tooltips_cav,
                 sizing_mode='stretch_width', width=pwdith*3)
CAV_linx_liny.line('x', 'y', source=source_cavx, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
CAV_linx_liny.line('x', 'y', source=source_cavy, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
CAV_linx_liny.line('x', 'y', source=source_cavz, legend_label=Station.component_3, line_width=2, color=Spectral6[2])
CAV_linx_liny.legend.click_policy = "hide"

panel_CAV_linx_liny = TabPanel(child=CAV_linx_liny, title='Cumulative Absolute Velocity')

CAV_linx_logy = figure(title='Cumulative Absolute Velocity', x_axis_label='Time [s]',
        y_axis_label='CAV [m/s]', tooltips=tooltips_cav,
                 y_axis_type='log', sizing_mode='stretch_width', width=pwdith*3)
CAV_linx_logy.line('x', 'y', source=source_cavx, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
CAV_linx_logy.line('x', 'y', source=source_cavy, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
CAV_linx_logy.line('x', 'y', source=source_cavz, legend_label=Station.component_3, line_width=2, color=Spectral6[2])
CAV_linx_logy.legend.click_policy = "hide"

panel_CAV_linx_logy = TabPanel(child=CAV_linx_logy, title='Cumulative Absolute Velocity')

CAV_logx_liny = figure(title='Cumulative Absolute Velocity', x_axis_label='Time [s]',
        y_axis_label='CAV [m/s]', tooltips=tooltips_cav,
                 x_axis_type='log', sizing_mode='stretch_width', width=pwdith*3)
CAV_logx_liny.line('x', 'y', source=source_cavx, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
CAV_logx_liny.line('x', 'y', source=source_cavy, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
CAV_logx_liny.line('x', 'y', source=source_cavz, legend_label=Station.component_3, line_width=2, color=Spectral6[2])
CAV_logx_liny.legend.click_policy = "hide"

panel_CAV_logx_liny = TabPanel(child=CAV_logx_liny, title='Cumulative Absolute Veocity')

CAV_logx_logy = figure(title='Cumulative Absolute Velocity', x_axis_label='Time [s]',
        y_axis_label='CAV [m/s]', tooltips=tooltips_cav,
                 x_axis_type='log', y_axis_type='log', sizing_mode='stretch_width', width=pwdith*3)
CAV_logx_logy.line('x', 'y', source=source_cavx, legend_label=Station.component_1, line_width=2, color=Spectral6[0])
CAV_logx_logy.line('x', 'y', source=source_cavy, legend_label=Station.component_2, line_width=2, color=Spectral6[1])
CAV_logx_logy.line('x', 'y', source=source_cavz, legend_label=Station.component_3, line_width=2, color=Spectral6[2])
CAV_logx_logy.legend.click_policy = "hide"

panel_CAV_logx_logy = TabPanel(child=CAV_logx_logy, title='Cumulative Absolute Velocity')

#%% TabPanels
panels_spectrum = [panel_Sa_linx_liny,
panel_Sv_linx_liny,
panel_Fourier_linx_liny,
panel_Husid_linx_liny,
panel_CAV_linx_liny]
tabs_spectrum = Tabs(tabs=panels_spectrum, sizing_mode='stretch_width', width=pwdith*3)

# Table
table = dict(fields=['Start time',
        'Magnitude [Mw]',
        'Longitude hypocenter',
        'Latitude hypocenter',
        'Event type',
        'Depth [km]',
        'Station name',
        'Longitude station',
        'Latitude station',
        'Hypocentral distance [km]',
        'Epicentral distance [km]',
        'Rupture distance [km]',
        'Joyner-Boore distance [km]',
        'Vs30 [m/s]',
        'Azimut [o]'],
        values=[Station.starttime,
        Station.magnitude,
        Station.hypocenter_lon,
        Station.hypocenter_lat,
        Station.event_type,
        Station.depth,
        Station.station_name,
        Station.station_lon,
        Station.station_lat,
        Station.Rhypo,
        Station.Repi,
        Station.Rrup,
        Station.Rjb,
        Station.vs30,
        Station.azimuth])

source_table = ColumnDataSource(data=table)
columns = [TableColumn(field="fields", title="Parameter"),
TableColumn(field="values", title="Value")]
data_table = DataTable(source=source_table, columns=columns, sizing_mode='stretch_width', width=pwdith*4//3, reorderable=False)

# Map
inProj = pyproj.Proj("EPSG:{0}".format(4326)) # WGS84
outProj = pyproj.Proj("EPSG:{0}".format(3857)) # Mercator

hypo_x, hypo_y = pyproj.transform(inProj, outProj, Station.hypocenter_lon, Station.hypocenter_lat, always_xy=True)
sta_x, sta_y   = pyproj.transform(inProj, outProj, Station.station_lon, Station.station_lat, always_xy=True)

xmean = 0.5*(hypo_x + sta_x)
ymean = 0.5*(hypo_y + sta_y)

dist = max([abs(hypo_x-xmean), abs(sta_x-xmean), abs(hypo_y-ymean), abs(sta_y-ymean)])

xmin = xmean - dist - 30000
xmax = xmean + dist + 30000

ymin = ymean - dist - 30000
ymax = ymean + dist + 30000

p = figure(x_range=(xmin, xmax), y_range=(ymin, ymax),
        x_axis_type="mercator", y_axis_type="mercator", sizing_mode='stretch_width', width=pwdith*4*2//3)
p.add_tile("CartoDB Positron")

source_hypo = ColumnDataSource(
        data=dict(lat=[ hypo_y],
            lon=[ hypo_x]))

source_sta  = ColumnDataSource(
        data=dict(lat=[ sta_y],
            lon=[ sta_x]))

p.circle(x='lon', y='lat', size=15, fill_color='blue', line_color='black', fill_alpha=0.8, source=source_sta, legend_label='Station')
p.star(x='lon', y='lat', size=25, fill_color='yellow', line_color='black', fill_alpha=0.8, source=source_hypo, legend_label='Hypocenter')
p.legend.location = "top_left"

#%% Login and download
username = TextInput(title='E-mail', value='')
password = TextInput(title='Access code', value='')
#login    = Button(label="Submit", button_type="success", align='end')

# Button download
validated = ColumnDataSource(data=dict(val=[0]))

download_button_callback_code = """
if (validated.data['val'][0] == 1) {
    
    if (extension.value == 'Matlab (MAT-binary, *.mat)') {
        var filename = source.value.concat('.mat');
        var path = 'events_mat_corrected_v2/';
    } else {
        var filename = source.value.concat('.pkl');
        var path = 'events_pkl_corrected_v2/';
    }
        
    var link = document.createElement('a');
    link.setAttribute('download', filename);
    link.href = path + filename;
    document.body.appendChild(link);
    link.click();
    link.remove();

    var d = new Date,
    dformat = [d.getMonth()+1,
               d.getDate(),
               d.getFullYear()].join('/')+' '+
              [d.getHours(),
               d.getMinutes(),
               d.getSeconds()].join(':');
    
    var line = email + ',' + dformat + ',' + filename;
     
} else {
    var email = prompt('E-mail:', '');
    var pass  = prompt('Code:', '');

    for (var i = 0; i < users.length; i++){
        if (email == users[i] && pass == passwords[i]){
            validated.data['val'] = [1]
            validated.change.emit();
            
            if (extension.value == 'Matlab (MAT-binary, *.mat)') {
                var filename = source.value.concat('.mat');
                var path = 'events_mat_corrected_v2/';
            } else {
                var filename = source.value.concat('.pkl');
                var path = 'events_pkl_corrected_v2/';
            }
            
            var link = document.createElement('a');
            link.setAttribute('download', filename);
            link.href = path + filename;
            document.body.appendChild(link);
            link.click();
            link.remove();
        }
    }
    if (validated.data['val'][0] == 0) {    
        alert('Username and/or password wrong!');
    }
 }
"""

download_button_callback = CustomJS(args=dict(validated=validated, source=events, users=users, passwords=passwords, extension=fileType),
        code=download_button_callback_code)

download_button = Button(label="Download", button_type="success", sizing_mode='stretch_width', width=pwdith, align='end')

if bokeh.__version__[0] == '1':
    download_button.callback = download_button_callback
else:
    download_button.js_on_event(ButtonClick, download_button_callback)

# Update functions
def update_stations(attrname, old, new):
    global Stations, keys, users, passwords
    
    if Stations['st00'].event_name == events.value:
        print(events.value, Stations['st00'].event_name)
        if stationSelect.value in stations.options:
            stations.value = stationSelect.value
        return None
    
    event = events.value
    Stations  = spio.loadmat(os.path.join(matfiles, event) + '.mat', struct_as_record=False, squeeze_me=True)
    Stations.pop('__version__')
    Stations.pop('__header__')
    Stations.pop('__globals__')


    keyValues = [item for item in Stations.items()]
    snames    = []

    for keyValue in keyValues:
        Station = keyValue[1]
        snames.append(Station.station_name)

    order     = np.argsort(snames)
    keys      = []
    snames    = []
    for pos in order:
        keys.append(keyValues[pos][0])
        snames.append(keyValues[pos][1].name)

    stations.options = snames
    if not stationSelect.value in snames:
        stations.value   = snames[0]
    else:
        stations.value   = ''
        stations.value   = stationSelect.value
        
    # download_button_callback.args = dict(validated=validated, source=event + '.mat', users=users, passwords=passwords)

def update_acceleration(attrname, old, new):
    global Tn_lin, Tn_log, Spectrum_lin, Spectrum_log
    global source_Sa1_lin, source_Sa2_lin, source_SaRotD50_lin, source_SaRotD0_lin, source_SaRotD100_lin, source_Saxy_lin
    global source_Sa1_log, source_Sa2_log, source_SaRotD50_log, source_SaRotD0_log, source_SaRotD100_log, source_Saxy_log
    global source_fx, source_fy, source_fz, source_aix, source_aiy, source_aiz, source_cavx, source_cavy, source_cavz
    global source_ax, source_ay, source_az
    global source_vx, source_vy, source_vz
    global source_dx, source_dy, source_dz
    global tabs_acceleration

    if stations.value == '':
        return None
    pos     = stations.options.index(stations.value)
    Station = Stations[keys[pos]]

    tax = np.linspace(0., (len(Station.acc_uncorrected_1)-1)*Station.dt, len(Station.acc_uncorrected_1))
    source_ax.data = dict(x=tax, y=Station.acc_uncorrected_1/g)
    velocity_x = spi.cumtrapz(Station.acc_uncorrected_1, x=tax, initial=0.)
    source_vx.data = dict(x=tax, y=velocity_x)
    displacement_x = spi.cumtrapz(velocity_x, x=tax, initial=0.)
    source_dx.data = dict(x=tax, y=displacement_x)
    tabs_acceleration.tabs[0].title = Station.component_1

    tay = np.linspace(0., (len(Station.acc_uncorrected_2)-1)*Station.dt, len(Station.acc_uncorrected_2))
    source_ay.data = dict(x=tay, y=Station.acc_uncorrected_2/g)
    velocity_y = spi.cumtrapz(Station.acc_uncorrected_2, x=tay, initial=0.)
    source_vy.data = dict(x=tay, y=velocity_y)
    displacement_y = spi.cumtrapz(velocity_y, x=tay, initial=0.)
    source_dy.data = dict(x=tay, y=displacement_y)
    tabs_acceleration.tabs[1].title = Station.component_2

    taz = np.linspace(0., (len(Station.acc_uncorrected_3)-1)*Station.dt, len(Station.acc_uncorrected_3))
    source_az.data = dict(x=taz, y=Station.acc_uncorrected_3/g)
    velocity_z = spi.cumtrapz(Station.acc_uncorrected_3, x=taz, initial=0.)
    source_vz.data = dict(x=taz, y=velocity_z)
    displacement_z = spi.cumtrapz(velocity_z, x=taz, initial=0.)
    source_dz.data = dict(x=taz, y=displacement_z)
    tabs_acceleration.tabs[2].title = Station.component_3

    for tab in tabs_acceleration.tabs:
        for row in tab.child.children:
            for fig in row.children:
                for item in fig.legend.items:
                    item.update(label=dict(value = Station.station_name))

    if float(ta.value) <= 0.:
        Tn_lin = np.linspace(0., float(tb.value), npoints)
        Tn_log = np.hstack((0., np.logspace(np.log10(0.001), np.log10(float(tb.value)), npoints)))
    else:
        Tn_lin = np.linspace(float(ta.value), float(tb.value), npoints)
        Tn_log = np.logspace(np.log10(float(ta.value)), np.log10(float(tb.value)), npoints)
        
    Spectrum_lin = seismic.SpectraRot(Station.acc_uncorrected_1, Station.acc_uncorrected_2, Station.dt, Tn_lin, float(xi.value), nTheta)
    Spectrum_linx = Spectrum_lin[nx]
    Spectrum_liny = Spectrum_lin[ny]
    
    Spectrum_log = seismic.SpectraRot(Station.acc_uncorrected_1, Station.acc_uncorrected_2, Station.dt, Tn_log, float(xi.value), nTheta)
    Spectrum_logx = Spectrum_log[nx]
    Spectrum_logy = Spectrum_log[ny]

    source_Sa1_lin.data = dict(x=Tn_lin, y=Spectrum_linx/g)
    source_Sa2_lin.data = dict(x=Tn_lin, y=Spectrum_liny/g)
    source_SaRotD50_lin.data = dict(x=Tn_lin, y=np.median(Spectrum_lin, axis=0)/g)
    source_SaRotD0_lin.data = dict(x=Tn_lin, y=np.min(Spectrum_lin, axis=0)/g)
    source_SaRotD100_lin.data = dict(x=Tn_lin, y=np.max(Spectrum_lin, axis=0)/g)
    source_Saxy_lin.data = dict(x=Tn_lin, y=np.sqrt(Spectrum_linx*Spectrum_liny)/g)

    source_Sa1_log.data = dict(x=Tn_log, y=Spectrum_logx/g)
    source_Sa2_log.data = dict(x=Tn_log, y=Spectrum_logy/g)
    source_SaRotD50_log.data = dict(x=Tn_log, y=np.median(Spectrum_log, axis=0)/g)
    source_SaRotD0_log.data = dict(x=Tn_log, y=np.min(Spectrum_log, axis=0)/g)
    source_SaRotD100_log.data = dict(x=Tn_log, y=np.max(Spectrum_log, axis=0)/g)
    source_Saxy_log.data = dict(x=Tn_log, y=np.sqrt(Spectrum_logx*Spectrum_logy)/g)

    source_Sv1_lin.data = dict(x=Tn_lin, y=Spectrum_linx*Tn_lin/(2.*np.pi))
    source_Sv2_lin.data = dict(x=Tn_lin, y=Spectrum_liny*Tn_lin/(2.*np.pi))
    source_SvRotD50_lin.data = dict(x=Tn_lin, y=np.median(Spectrum_lin, axis=0)*Tn_lin/(2.*np.pi))
    source_SvRotD0_lin.data = dict(x=Tn_lin, y=np.min(Spectrum_lin, axis=0)*Tn_lin/(2.*np.pi))
    source_SvRotD100_lin.data = dict(x=Tn_lin, y=np.max(Spectrum_lin, axis=0)*Tn_lin/(2.*np.pi))
    source_Svxy_lin.data = dict(x=Tn_lin, y=np.sqrt(Spectrum_linx*Spectrum_liny)*Tn_lin/(2.*np.pi))

    source_Sv1_log.data = dict(x=Tn_log, y=Spectrum_logx*Tn_log/(2.*np.pi))
    source_Sv2_log.data = dict(x=Tn_log, y=Spectrum_logy*Tn_log/(2.*np.pi))
    source_SvRotD50_log.data = dict(x=Tn_log, y=np.median(Spectrum_log, axis=0)*Tn_log/(2.*np.pi))
    source_SvRotD0_log.data = dict(x=Tn_log, y=np.min(Spectrum_log, axis=0)*Tn_log/(2.*np.pi))
    source_SvRotD100_log.data = dict(x=Tn_log, y=np.max(Spectrum_log, axis=0)*Tn_log/(2.*np.pi))
    source_Svxy_log.data = dict(x=Tn_log, y=np.sqrt(Spectrum_logx*Spectrum_logy)*Tn_log/(2.*np.pi))

    L = Station.acc_uncorrected_1.shape[0]
    Fax = 2.*np.abs(spf.fft(Station.acc_uncorrected_1*Station.dt))[:L//2]
    Fay = 2.*np.abs(spf.fft(Station.acc_uncorrected_2*Station.dt))[:L//2]
    Faz = 2.*np.abs(spf.fft(Station.acc_uncorrected_3*Station.dt))[:L//2]
    Fs = 1./Station.dt
    freq = Fs*np.linspace(0., L//2, L//2)/L
    source_fx.data = dict(x=freq, y=Fax)
    source_fy.data = dict(x=freq, y=Fay)
    source_fz.data = dict(x=freq, y=Faz)

    t = np.linspace(0., (len(Station.acc_uncorrected_1)-1)*Station.dt, len(Station.acc_uncorrected_1))

    iax = np.pi/(2.*g)*spi.cumtrapz(Station.acc_uncorrected_1**2, t, initial=0.)
    source_aix.data  = dict(x=t, y=iax)

    iay = np.pi/(2.*g)*spi.cumtrapz(Station.acc_uncorrected_2**2, t, initial=0.)
    source_aiy.data  = dict(x=t, y=iay)

    iaz = np.pi/(2.*g)*spi.cumtrapz(Station.acc_uncorrected_3**2, t, initial=0.)
    source_aiz.data  = dict(x=t, y=iaz)

    cavx = spi.cumtrapz(np.abs(Station.acc_uncorrected_1), t, initial=0.)
    source_cavx.data  = dict(x=t, y=cavx)

    cavy = spi.cumtrapz(np.abs(Station.acc_uncorrected_2), t, initial=0.)
    source_cavy.data  = dict(x=t, y=cavy)

    cavz = spi.cumtrapz(np.abs(Station.acc_uncorrected_3), t, initial=0.)
    source_cavz.data  = dict(x=t, y=cavz)

    hypo_x, hypo_y = pyproj.transform(inProj, outProj, Station.hypocenter_lon, Station.hypocenter_lat, always_xy=True)
    sta_x, sta_y   = pyproj.transform(inProj, outProj, Station.lon, Station.lat, always_xy=True)

    xmean = 0.5*(hypo_x + sta_x)
    ymean = 0.5*(hypo_y + sta_y)

    dist = max([abs(hypo_x-xmean), abs(sta_x-xmean), abs(hypo_y-ymean), abs(sta_y-ymean)])

    xmin = xmean - dist - 30000
    xmax = xmean + dist + 30000

    ymin = ymean - dist - 30000
    ymax = ymean + dist + 30000

    p.x_range.start = xmin
    p.x_range.end   = xmax
    p.y_range.start = ymin
    p.y_range.end   = ymax

    source_hypo.data = dict(lat=[ hypo_y], lon=[ hypo_x])
    source_sta.data  = dict(lat=[ sta_y], lon=[ sta_x])

    source_table.data = dict(fields=['Event name',
            'Start time',
            'Magnitude [Mw]',
            'Longitude hypocenter',
            'Latitude hypocenter',
            'Event type',
            'Depth [km]',
            'Station name',
            'Longitude station',
            'Latitude station',
            'Hypocentral distance [km]',
            'Epicentral distance [km]',
            'Rupture distance [km]',
            'Joyner-Boore distance [km]',
            'Vs30 [m/s]',
            'Azimut [o]'],
            values=[Station.event_name,
            Station.start_time,
            Station.magnitude,
            Station.hypocenter_lon,
            Station.hypocenter_lat,
            Station.event_type,
            Station.depth,
            Station.station_name,
            Station.lon,
            Station.lat,
            Station.Rhypo,
            Station.Repi,
            Station.Rrup,
            Station.Rjb,
            Station.Vs30,
            Station.azimut])

def update_spectrum(attrname, old, new):
    global Tn_lin, Tn_log, Sprectrum_lin, Spectrum_log

    pos     = stations.options.index(stations.value)
    Station = Stations[keys[pos]]

    if float(ta.value) <= 0.:
        Tn_lin = np.linspace(0., float(tb.value), npoints)
        Tn_log = np.hstack((0., np.logspace(np.log10(0.001), np.log10(float(tb.value)), npoints)))
    else:
        Tn_lin = np.linspace(float(ta.value), float(tb.value), npoints)
        Tn_log = np.logspace(np.log10(float(ta.value)), np.log10(float(tb.value)), npoints)

    Spectrum_lin = seismic.SpectraRot(Station.acc_uncorrected_1, Station.acc_uncorrected_2, Station.dt, Tn_lin, float(xi.value), nTheta)
    Spectrum_linx = Spectrum_lin[nx]
    Spectrum_liny = Spectrum_lin[ny]

    Spectrum_log = seismic.SpectraRot(Station.acc_uncorrected_1, Station.acc_uncorrected_2, Station.dt, Tn_log, float(xi.value), nTheta)
    Spectrum_logx = Spectrum_log[nx]
    Spectrum_logy = Spectrum_log[ny]
    
    source_Sa1_lin.data = dict(x=Tn_lin, y=Spectrum_linx/g)
    source_Sa2_lin.data = dict(x=Tn_lin, y=Spectrum_liny/g)
    source_SaRotD50_lin.data = dict(x=Tn_lin, y=np.median(Spectrum_lin, axis=0)/g)
    source_SaRotD0_lin.data = dict(x=Tn_lin, y=np.min(Spectrum_lin, axis=0)/g)
    source_SaRotD100_lin.data = dict(x=Tn_lin, y=np.max(Spectrum_lin, axis=0)/g)
    source_Saxy_lin.data = dict(x=Tn_lin, y=np.sqrt(Spectrum_linx*Spectrum_liny)/g)
    
    source_Sa1_log.data = dict(x=Tn_log, y=Spectrum_logx/g)
    source_Sa2_log.data = dict(x=Tn_log, y=Spectrum_logy/g)
    source_SaRotD50_log.data = dict(x=Tn_log, y=np.median(Spectrum_log, axis=0)/g)
    source_SaRotD0_log.data = dict(x=Tn_log, y=np.min(Spectrum_log, axis=0)/g)
    source_SaRotD100_log.data = dict(x=Tn_log, y=np.max(Spectrum_log, axis=0)/g)
    source_Saxy_log.data = dict(x=Tn_log, y=np.sqrt(Spectrum_logx*Spectrum_logy)/g)
    
    source_Sv1_lin.data = dict(x=Tn_lin, y=Spectrum_linx*Tn_lin/(2.*np.pi))
    source_Sv2_lin.data = dict(x=Tn_lin, y=Spectrum_liny*Tn_lin/(2.*np.pi))
    source_SvRotD50_lin.data = dict(x=Tn_lin, y=np.median(Spectrum_lin, axis=0)*Tn_lin/(2.*np.pi))
    source_SvRotD0_lin.data = dict(x=Tn_lin, y=np.min(Spectrum_lin, axis=0)*Tn_lin/(2.*np.pi))
    source_SvRotD100_lin.data = dict(x=Tn_lin, y=np.max(Spectrum_lin, axis=0)*Tn_lin/(2.*np.pi))
    source_Svxy_lin.data = dict(x=Tn_lin, y=np.sqrt(Spectrum_linx*Spectrum_liny)*Tn_lin/(2.*np.pi))

#    Sv_linx_liny.x_range.update(start=float(ta.value), end=float(tb.value), bounds=bounds)
#    Sv_linx_logy.x_range.update(start=float(ta.value), end=float(tb.value), bounds=bounds)

    source_Sv1_log.data = dict(x=Tn_log, y=Spectrum_logx*Tn_log/(2.*np.pi))
    source_Sv2_log.data = dict(x=Tn_log, y=Spectrum_logy*Tn_log/(2.*np.pi))
    source_SvRotD50_log.data = dict(x=Tn_log, y=np.median(Spectrum_log, axis=0)*Tn_log/(2.*np.pi))
    source_SvRotD0_log.data = dict(x=Tn_log, y=np.min(Spectrum_log, axis=0)*Tn_log/(2.*np.pi))
    source_SvRotD100_log.data = dict(x=Tn_log, y=np.max(Spectrum_log, axis=0)*Tn_log/(2.*np.pi))
    source_Svxy_log.data = dict(x=Tn_log, y=np.sqrt(Spectrum_logx*Spectrum_logy)*Tn_log/(2.*np.pi))

def update_axis(attrname, old, new):
    if len(checkbox_axis.active) == 0:
        tabs_spectrum.tabs[0] = panel_Sa_linx_liny
        tabs_spectrum.tabs[1] = panel_Sv_linx_liny
        tabs_spectrum.tabs[2] = panel_Fourier_linx_liny
        tabs_spectrum.tabs[3] = panel_Husid_linx_liny
        tabs_spectrum.tabs[4] = panel_CAV_linx_liny
    elif len(checkbox_axis.active) == 1:
        if 0 in checkbox_axis.active:
            tabs_spectrum.tabs[0] = panel_Sa_logx_liny
            tabs_spectrum.tabs[1] = panel_Sv_logx_liny
            tabs_spectrum.tabs[2] = panel_Fourier_logx_liny
            tabs_spectrum.tabs[3] = panel_Husid_logx_liny
            tabs_spectrum.tabs[4] = panel_CAV_logx_liny
        else:
            tabs_spectrum.tabs[0] = panel_Sa_linx_logy
            tabs_spectrum.tabs[1] = panel_Sv_linx_logy
            tabs_spectrum.tabs[2] = panel_Fourier_linx_logy
            tabs_spectrum.tabs[3] = panel_Husid_linx_logy
            tabs_spectrum.tabs[4] = panel_CAV_linx_logy
    else:
        tabs_spectrum.tabs[0] = panel_Sa_logx_logy
        tabs_spectrum.tabs[1] = panel_Sv_logx_logy
        tabs_spectrum.tabs[2] = panel_Fourier_logx_logy
        tabs_spectrum.tabs[3] = panel_Husid_logx_logy
        tabs_spectrum.tabs[4] = panel_CAV_logx_logy 

def update_grid(attrname, old, new):

    if 0 in checkbox_grid.active:
        colorx = '#e5e5e5'
    else:
        colorx = None

    if 1 in checkbox_grid.active:
        colory = '#e5e5e5'
    else:
        colory = None

    Sa_linx_liny.xgrid[0].minor_grid_line_color = colorx
    Sa_linx_liny.ygrid[0].minor_grid_line_color = colory

    Sa_linx_logy.xgrid[0].minor_grid_line_color = colorx
    Sa_linx_logy.ygrid[0].minor_grid_line_color = colory

    Sa_logx_liny.xgrid[0].minor_grid_line_color = colorx
    Sa_logx_liny.ygrid[0].minor_grid_line_color = colory

    Sa_logx_logy.xgrid[0].minor_grid_line_color = colorx
    Sa_logx_logy.ygrid[0].minor_grid_line_color = colory

    Sv_linx_liny.xgrid[0].minor_grid_line_color = colorx
    Sv_linx_liny.ygrid[0].minor_grid_line_color = colory

    Sv_linx_logy.xgrid[0].minor_grid_line_color = colorx
    Sv_linx_logy.ygrid[0].minor_grid_line_color = colory

    Sv_logx_liny.xgrid[0].minor_grid_line_color = colorx
    Sv_logx_liny.ygrid[0].minor_grid_line_color = colory

    Sv_logx_logy.xgrid[0].minor_grid_line_color = colorx
    Sv_logx_logy.ygrid[0].minor_grid_line_color = colory

    Fourier_linx_liny.xgrid[0].minor_grid_line_color = colorx
    Fourier_linx_liny.ygrid[0].minor_grid_line_color = colory

    Fourier_linx_logy.xgrid[0].minor_grid_line_color = colorx
    Fourier_linx_logy.ygrid[0].minor_grid_line_color = colory

    Fourier_logx_liny.xgrid[0].minor_grid_line_color = colorx
    Fourier_logx_liny.ygrid[0].minor_grid_line_color = colory

    Fourier_logx_logy.xgrid[0].minor_grid_line_color = colorx
    Fourier_logx_logy.ygrid[0].minor_grid_line_color = colory

    Husid_linx_liny.xgrid[0].minor_grid_line_color = colorx
    Husid_linx_liny.ygrid[0].minor_grid_line_color = colory

    Husid_linx_logy.xgrid[0].minor_grid_line_color = colorx
    Husid_linx_logy.ygrid[0].minor_grid_line_color = colory

    Husid_logx_liny.xgrid[0].minor_grid_line_color = colorx
    Husid_logx_liny.ygrid[0].minor_grid_line_color = colory

    Husid_logx_logy.xgrid[0].minor_grid_line_color = colorx
    Husid_logx_logy.ygrid[0].minor_grid_line_color = colory

    CAV_linx_liny.xgrid[0].minor_grid_line_color = colorx
    CAV_linx_liny.ygrid[0].minor_grid_line_color = colory

    CAV_linx_logy.xgrid[0].minor_grid_line_color = colorx
    CAV_linx_logy.ygrid[0].minor_grid_line_color = colory

    CAV_logx_liny.xgrid[0].minor_grid_line_color = colorx
    CAV_logx_liny.ygrid[0].minor_grid_line_color = colory

    CAV_logx_logy.xgrid[0].minor_grid_line_color = colorx
    CAV_logx_logy.ygrid[0].minor_grid_line_color = colory


#%% Website creation
events.on_change('value', update_stations)
stations.on_change('value', update_acceleration)
xi.on_change('value', update_spectrum)
ta.on_change('value', update_spectrum)
tb.on_change('value', update_spectrum)
checkbox_axis.on_change('active', update_axis)
checkbox_grid.on_change('active', update_grid)

# Export website
_env = Environment(loader=FileSystemLoader('StrongMotionDatabase'))
FILE = _env.get_template("siberrisk_seismicdatabase.html")
curdoc().template = FILE

div2 = Div(text='<h2>Filter options</h2>', sizing_mode='stretch_width', width=pwdith)

div3 = Div(text='<h2>Strong Motion Database</h2>\nTo download several files at once please click <a href="https://siberrisk.ing.puc.cl/StrongMotionDatabaseDownloadManager" target="_blank">here</a>. Examples to read the database files are available for <a href="examplePython.py" target="_blank">Python</a> and <a href="exampleMatlab.m" target="_blank">Matlab</a>.<br/>', sizing_mode='stretch_width', width=pwdith)

div4 = Div(text='<b>TIP:</b> You can hide/show the curves by pressing their names at the legend box!', sizing_mode='stretch_width', width=pwdith)

inputs_plots = column([ta, tb, xi, checkbox_axis, checkbox_grid, div4], sizing_mode='stretch_width', width=pwdith)

from bokeh.layouts import grid

distribution = grid([[div2, None, None, None],
                     [eventsSince, minMagnitude, eventType, filterButton],
                     [eventsUntil, maxMagnitude, stationSelect, None],
                     [div3],
                     [events, stations, fileType, download_button],
                     [tabs_acceleration],
                     [inputs_plots, tabs_spectrum],
                     [data_table, p, None]], sizing_mode='stretch_width')

distribution.children[14] = (distribution.children[14][0], distribution.children[14][1], distribution.children[14][2], 1, 3)
distribution.children[15] = (distribution.children[15][0], distribution.children[15][1], distribution.children[14][2]+3, 1, 9)
distribution.children[16] = (distribution.children[16][0], distribution.children[16][1], distribution.children[16][2], 1, 4)
distribution.children[17] = (distribution.children[17][0], distribution.children[17][1], distribution.children[16][2]+4, 1, 8)

curdoc().add_root(distribution)
curdoc().title = 'Strong Motion Database ' + u'\u2013' + ' SIBER-RISK'
