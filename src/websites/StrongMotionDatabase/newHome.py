import os
import sys
import random
import string
import datetime
import numpy as np
import pandas as pd
import scipy.fftpack as spf
import scipy.integrate as spi

###################
## Bokeh modules ##
###################
from bokeh.io import curdoc

from bokeh.layouts import column, grid, layout

from bokeh.models import ColumnDataSource, CustomJS, DatePicker, Tabs, TabPanel
from bokeh.models.widgets import Button, CheckboxGroup, DataTable, Div, NumericInput, Select, TableColumn

from bokeh.palettes import Spectral6

from bokeh.plotting import figure

from jinja2 import Environment, FileSystemLoader

# Paths
currentDir = os.path.dirname(__file__)
libPath    = os.path.abspath(os.path.join(currentDir, '..', '..', 'lib'))
dataPath   = os.path.abspath(os.path.join(currentDir, '..', '..', '..', 'data'))

if not libPath in sys.path:
    sys.path.append(libPath)

import seismic

################
##  Settings  ##
################
g          = 9.81 # m/s**2
npoints    = 250
allowed    = False
pwdith     = 1
pheight    = 250
lwidth     = 2
nTheta     = 91
nx         = 0
ny         = 45
formats    = ['Matlab (MAT-binary, *.mat)', 'Python (Numpy, *.npz)']
axis_types = [('linear', 'linear'), ('linear', 'log'), ('log', 'linear'), ('log', 'log')]

########################
##  Global variables  ##
########################
event = {}

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

############
##  Divs  ##
############
div_filter  = Div(text='<h2>Filter options</h2>', sizing_mode='stretch_width', width=pwdith)

div_records = Div(text='<h2>Strong Motion Database</h2>\nTo download several files at once please ' +
    'click <a href="https://siberrisk.ing.puc.cl/StrongMotionDatabaseDownloadManager" target=' +
    '"_blank">here</a>. Examples to read the database files are available for <a href="examplePython.py"' +
    ' target="_blank">Python</a> and <a href="exampleMatlab.m" target="_blank">Matlab</a>.<br/>',
    sizing_mode='stretch_width', width=pwdith)

div_plots   = Div(text='<b>TIP:</b> You can hide/show the curves by pressing their names at the legend box!',
    sizing_mode='stretch_width', width=pwdith)

#########################
##  Filter components  ##
#########################
filter_since = DatePicker(min_date=min_date, max_date=max_date, value=min_date,
    title='Show events since', sizing_mode='stretch_width', width=pwdith)

filter_until = DatePicker(min_date=min_date, max_date=max_date, value=max_date,
    title='Show events until', sizing_mode='stretch_width', width=pwdith)

filter_eType = Select(title='Event type', value=eTypes[0], options=eTypes,
    sizing_mode='stretch_width', width=pwdith*2)

filter_minMw = NumericInput(title='Show events larger or equal than', value=min_mag, mode='float',
    sizing_mode='stretch_width', width=pwdith)

filter_maxMw = NumericInput(title='Show events smaller or equal than', value=max_mag, mode='float',
    sizing_mode='stretch_width', width=pwdith)

filter_sCode = Select(title='Recorded by station', value=sCodes[0], options=sCodes,
    sizing_mode='stretch_width', width=pwdith*2)

#########################
##  Select components  ##
#########################
select_event   = Select(title='Seismic events', value=seismic_events[0], options=seismic_events,
    sizing_mode='stretch_width', width=pwdith*2)

select_station = Select(title='Stations', value=station_codes[0], options=station_codes,
    sizing_mode='stretch_width', width=pwdith)

select_format = Select(title='File format', value=formats[0], options=formats,
    sizing_mode='stretch_width', width=pwdith)

###############
##  Buttons  ##
###############
button_filter   = Button(label='Apply filters', button_type='success',
    sizing_mode='stretch_width', width=pwdith, align='end')

button_download = Button(label='Download', button_type='success',
    sizing_mode='stretch_width', width=pwdith, align='end')

################
## Tool tips  ##
################
tooltips_a   = [('(t, a)', '($x, $y)'),]
tooltips_v   = [('(t, v)', '($x, $y)'),]
tooltips_d   = [('(t, d)', '($x, $y)'),]
tooltips_sa  = [('(Tn, Sa)', '($x, $y)'),]
tooltips_sv  = [('(Tn, Sv)', '($x, $y)'),]
tooltips_f   = [('(f, A)', '($x, $y)'),]
tooltips_ai  = [('(t, IA)', '($x, $y)'),]
tooltips_cav = [('(t, CAV)', '($x, $y)'),]

##########################
##  Acceleration plots  ##
##########################
plots_acc = [figure(title='Acceleration', x_axis_label='Time [s]', y_axis_label='Acceleration [g]',
    tooltips=tooltips_a, height=pheight, sizing_mode='stretch_width', width=pwdith*4) for i in range(3)]

######################
##  Velocity plots  ##
######################
plots_vel = [figure(title='Velocity', x_axis_label='Time [s]', y_axis_label='Velocity [m/s]',
    tooltips=tooltips_v, height=pheight, sizing_mode='stretch_width', width=pwdith*4) for i in range(3)]

##########################
##  Displacement plots  ##
##########################
plots_dis = [figure(title='Displacement', x_axis_label='Time [s]', y_axis_label='Displacement [m]',
    tooltips=tooltips_d, height=pheight, sizing_mode='stretch_width', width=pwdith*4) for i in range(3)]

#########################
##  Sa Spectrum plots  ##
#########################
plots_sa = [figure(title='Pseudo-acceleration Response Spectrum',
    x_axis_label='Period [s]', y_axis_label='Spectral Acceleration [g]',
    tooltips=tooltips_sa, sizing_mode='stretch_width', width=pwdith*3,
    x_axis_type=axis_type[0], y_axis_type=axis_type[1]) for axis_type in axis_types]

##########################
##  DVA Spectrum plots  ##
##########################
plots_dva = [figure(title='Combined D-V-A Spectrum',
    x_axis_label='Period [s]', y_axis_label='Spectral Velocity [m/s]',
    tooltips=tooltips_sv, sizing_mode='stretch_width', width=pwdith*3,
    x_axis_type=axis_type[0], y_axis_type=axis_type[1]) for axis_type in axis_types]

##############################
##  Fourier Spectrum plots  ##
##############################
plots_fourier = [figure(title='Fourier Spectrum',
    x_axis_label='Frequency [Hz]', y_axis_label='Amplitude [m/s]',
    tooltips=tooltips_f, sizing_mode='stretch_width', width=pwdith*3,
    x_axis_type=axis_type[0], y_axis_type=axis_type[1]) for axis_type in axis_types]

###################
##  Husid plots  ##
###################
plots_husid = [figure(title='Husid Plot',
    x_axis_label='Time [s]', y_axis_label='Arias intensity [m/s]',
    tooltips=tooltips_ai, sizing_mode='stretch_width', width=pwdith*3,
    x_axis_type=axis_type[0], y_axis_type=axis_type[1]) for axis_type in axis_types]

#################
##  CAV plots  ##
#################
plots_cav = [figure(title='Cumulative Absolute Plot',
    x_axis_label='Time [s]', y_axis_label='CAV [m/s]',
    tooltips=tooltips_cav, sizing_mode='stretch_width', width=pwdith*3,
    x_axis_type=axis_type[0], y_axis_type=axis_type[1]) for axis_type in axis_types]

###############
##  Sources  ##
###############
sources_acc     = [ColumnDataSource(data=dict(x=[], y=[])) for i in range(3)]
sources_vel     = [ColumnDataSource(data=dict(x=[], y=[])) for i in range(3)]
sources_dis     = [ColumnDataSource(data=dict(x=[], y=[])) for i in range(3)]
sources_sa      = [ColumnDataSource(data=dict(x=[], y=[])) for i in range(6)]
sources_dva     = [ColumnDataSource(data=dict(x=[], y=[])) for i in range(6)]
sources_fourier = [ColumnDataSource(data=dict(x=[], y=[])) for i in range(3)]
sources_husid   = [ColumnDataSource(data=dict(x=[], y=[])) for i in range(3)]
sources_cav     = [ColumnDataSource(data=dict(x=[], y=[])) for i in range(3)]

#############
##  Lines  ##
#############
for i in range(3):
    plots_acc[i].line('x', 'y', source=sources_acc[i])
    plots_vel[i].line('x', 'y', source=sources_vel[i])
    plots_dis[i].line('x', 'y', source=sources_dis[i])

for i in range(len(axis_types)):
    for j in range(3):
        plots_fourier[i].line('x', 'y', source=sources_fourier[j], line_width=lwidth, color=Spectral6[j])
        plots_husid[i].line('x', 'y', source=sources_husid[j], line_width=lwidth, color=Spectral6[j])
        plots_cav[i].line('x', 'y', source=sources_cav[j], line_width=lwidth, color=Spectral6[j])

    for j in range(6):
        plots_sa[i].line('x', 'y', source=sources_sa[j], line_width=lwidth, color=Spectral6[j])
        plots_dva[i].line('x', 'y', source=sources_dva[j], line_width=lwidth, color=Spectral6[j])

####################
##  Plots inputs  ##
####################
options_ta = NumericInput(title='Initial period [s]', value=0., low=0., mode='float',
    sizing_mode='stretch_width', width=pwdith)

options_tb = NumericInput(title='Ending period [s]', value=1.5, low=0., mode='float',
    sizing_mode='stretch_width', width=pwdith)

options_xi = NumericInput(title='Damping ratio', value=0.05, low=0., mode='float',
    sizing_mode='stretch_width', width=pwdith)

options_ax = CheckboxGroup(labels=['Logarithmic X axis', 'Logarithmic Y axis'],
    sizing_mode='stretch_width', width=pwdith)

options_gr = CheckboxGroup(labels=['X axis minor grid', 'Y axis minor grid'],
    sizing_mode='stretch_width', width=pwdith)

inputs_plots = column([options_ta, options_tb, options_xi, options_ax, options_gr, div_plots],
    sizing_mode='stretch_width', width=pwdith)

#############
##  Table  ##
#############
columns = [TableColumn(field="fields", title="Parameter"), TableColumn(field="values", title="Value")]
table   = DataTable(columns=columns, sizing_mode='stretch_width', width=pwdith*4//3, reorderable=False)

###########
##  Map  ##
###########
plot_map = figure(x_axis_type="mercator", y_axis_type="mercator", sizing_mode='stretch_width', width=pwdith*4*2//3)
plot_map.add_tile("CartoDB Positron")

############
##  Tabs  ##
############
panel_records = [TabPanel(child=layout([[plots_acc[i]], [plots_vel[i]], [plots_dis[i]]],
    sizing_mode='stretch_width', width=pwdith*4)) for i in range(3)]

tabs_records  = Tabs(tabs=panel_records,
    sizing_mode='stretch_width', width=pwdith*4)

panel_plots   = [
    TabPanel(child=plots_sa[0]     , title='Sa Spectra'),
    TabPanel(child=plots_dva[0]    , title='D-V-A Spectrum'),
    TabPanel(child=plots_fourier[0], title='Fourier Spectrum'),
    TabPanel(child=plots_husid[0]  , title='Husid Plot'),
    TabPanel(child=plots_cav[0]    , title='Cumulative Absolute Velocity')
]

tabs_plots    = Tabs(tabs=panel_plots,
    sizing_mode='stretch_width', width=pwdith*3)

######################################
## Load and process data functions  ##
######################################
def load_event(event_name):
    stations_codes = []
    with np.load(os.path.join(dataPath, 'seismicDatabase', 'npz', event_name + '.npz'), allow_pickle=True) as f:
        event = {}
        for key, value in f.items():
            if not key.startswith('st'):
                continue

            station      = value.item()
            station_code = station['station_code']

            stations_codes.append(station_code)
            event[station_code] = station

    stations_codes = list(sorted(stations_codes))

    return event, stations_codes

def compute_vel_dis(station):
    velocities    = []
    displacements = []
    dt = station['dt']
    for i in range(3):
        n = len(station['acc_uncorrected_%i' %(i+1)])
        t = np.linspace(0., (n-1)*dt, n)
        acc = station['acc_uncorrected_%i' %(i+1)]
        vel = spi.cumtrapz(acc, x=t, initial=0.)
        dis = spi.cumtrapz(vel, x=t, initial=0.)
        velocities.append(vel)
        displacements.append(dis)

    return velocities, displacements

def compute_sa_spectra(station, tn, xi):
    spectra = seismic.SpectraRot(station['acc_uncorrected_1'], station['acc_uncorrected_2'], station['dt'], tn, xi, nTheta)

    spectrum_x       = spectra[nx]/g
    spectrum_y       = spectra[ny]/g
    spectrum_rotd50  = np.median(spectra, axis=0)/g
    spectrum_rotd0   = np.min(spectra, axis=0)/g
    spectrum_rotd100 = np.max(spectra, axis=0)/g
    spectrum_geom    = np.sqrt(spectrum_x*spectrum_y)

    return spectrum_x, spectrum_y, spectrum_rotd50, spectrum_rotd0, spectrum_rotd100, spectrum_geom

def compute_dva_spectra(sa_spectra, tn):
    dva_spectra = []
    for sa_spectrum in sa_spectra:
        dva_spectra.append(sa_spectrum*tn/(2.*np.pi))

    return dva_spectra

def compute_fourier_spectra(station):
    fourier_spectra = []
    dt = station['dt']
    Fs = 1./dt
    for i in range(3):
        L  = len(station['acc_uncorrected_%i' %(i+1)])
        freq = Fs*np.linspace(0., L//2, L//2)/L
        Fa = 2.*np.abs(spf.fft(station['acc_uncorrected_%i' %(i+1)])*dt)[:L//2]
        fourier_spectra.extend([freq, Fa])

    return fourier_spectra

def compute_husid_plot(station):
    husid_plots = []
    dt = station['dt']
    for i in range(3):
        n  = len(station['acc_uncorrected_%i' %(i+1)])
        t  = np.linspace(0., (n-1)*dt, n)
        ia = np.pi/(2.*g)*spi.cumtrapz(station['acc_uncorrected_%i' %(i+1)]**2, t, initial=0.)
        husid_plots.append(ia)

    return husid_plots

def compute_cav_plot(station):
    cav_plots = []
    dt = station['dt']
    for i in range(3):
        n   = len(station['acc_uncorrected_%i' %(i+1)])
        t   = np.linspace(0., (n-1)*dt, n)
        cav = spi.cumtrapz(np.abs(station['acc_uncorrected_%i' %(i+1)]), t, initial=0.)
        cav_plots.append(cav)

    return cav_plots

#######################
## Update functions  ##
#######################
def update_event(attrname, old, new):
    global event
    
    if old == new:
        return

    event, station_codes = load_event(new)

    select_station.options = station_codes
    if filter_sCode.value in station_codes:
        select_station.value = ''
        select_station.value = filter_sCode.value
    else:
        select_station.value = station_codes[0]
    
def update_station(attrname, old, new):
    global event

    station = event[select_station.value]

    if 0 in options_ax.active:
        tn = np.logspace(np.log10(max(options_ta.value, 0.001)), np.log10(options_tb.value), npoints)
    else:
        tn = np.linspace(options_ta.value, options_tb.value, npoints)

    vel, dis        = compute_vel_dis(station)
    sa_spectra      = compute_sa_spectra(station, tn, options_xi.value)
    dva_spectra     = compute_dva_spectra(sa_spectra, tn)
    fourier_spectra = compute_fourier_spectra(station)
    husid_plot      = compute_husid_plot(station)
    cav_plot        = compute_cav_plot(station)

    for i in range(3):
        n  = len(station['acc_uncorrected_%i' %(i+1)])
        dt = station['dt']
        t  = np.linspace(0., (n-1)*dt, n)
        sources_acc[i].data = dict(x=t, y=station['acc_uncorrected_%i' %(i+1)]/g)
        sources_vel[i].data = dict(x=t, y=vel[i])
        sources_dis[i].data = dict(x=t, y=dis[i])

        sources_fourier[i].data = dict(x=fourier_spectra[2*i], y=fourier_spectra[2*i+1])
        sources_husid[i].data   = dict(x=t, y=husid_plot[i])
        sources_cav[i].data     = dict(x=t, y=cav_plot[i])

    for i in range(6):
        sources_sa[i].data  = dict(x=tn, y=sa_spectra[i])
        sources_dva[i].data = dict(x=tn, y=dva_spectra[i])

def update_spectra_options(attrname, old, new):
    global event

    # assert type(options_ta.value) == float
    # assert type(options_tb.value) == float
    # assert type(options_xi.value) == float

    if options_ta.value >= options_tb.value:
        return

    station = event[select_station.value]

    if 0 in options_ax.active:
        tn = np.logspace(np.log10(max(options_ta.value, 0.001)), np.log10(options_tb.value), npoints)
    else:
        tn = np.linspace(options_ta.value, options_tb.value, npoints)

    sa_spectra = compute_sa_spectra(station, tn, options_xi.value)

####################################
##  Custom Javascript Functions   ##
####################################
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

########################
##  Filter functions  ##
########################
def filter_events():
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
        conditions &= (flatfile['Event type'] == eType.lower())

    if sCode != 'Any':
        conditions &= (flatfile['Station code'] == sCode)

    EarthquakeNames = flatfile[conditions]['Earthquake Name']
    EarthquakeNames = sorted(EarthquakeNames.drop_duplicates().values.tolist(), reverse=True)
    
    Alert.args = dict(n=len(EarthquakeNames))
    
    lettersandnumbers  = string.ascii_lowercase + string.digits
    button_filter.name = ''.join(random.choice(lettersandnumbers) for i in range(10))
    if len(EarthquakeNames) > 0:
        select_event.options = EarthquakeNames
        select_event.value   = EarthquakeNames[0]

########################
## Javascript events  ##
########################
button_filter.js_on_change('name', Alert)
button_filter.on_click(filter_events)
select_event.on_change('value', update_event)
select_station.on_change('value', update_station)
options_ta.on_change('value', update_spectra_options)
options_tb.on_change('value', update_spectra_options)
options_xi.on_change('value', update_spectra_options)
# checkbox_axis.on_change('active', update_axis)
# checkbox_grid.on_change('active', update_grid)

######################
##  Export website  ##
######################
_env = Environment(loader=FileSystemLoader('StrongMotionDatabase'))
FILE = _env.get_template("siberrisk_seismicdatabase.html")
curdoc().template = FILE

distribution = grid([[div_filter, None, None, None],
                     [filter_since, filter_minMw, filter_eType, button_filter],
                     [filter_until, filter_maxMw, filter_sCode, None],
                     [div_records],
                     [select_event, select_station, select_format, button_download],
                     [tabs_records],
                     [inputs_plots, tabs_plots],
                     [table, plot_map, None]], sizing_mode='stretch_width')

distribution.children[14] = (distribution.children[14][0], distribution.children[14][1], distribution.children[14][2]  , 1, 3)
distribution.children[15] = (distribution.children[15][0], distribution.children[15][1], distribution.children[14][2]+3, 1, 9)
distribution.children[16] = (distribution.children[16][0], distribution.children[16][1], distribution.children[16][2]  , 1, 4)
distribution.children[17] = (distribution.children[17][0], distribution.children[17][1], distribution.children[16][2]+4, 1, 8)

curdoc().add_root(distribution)
curdoc().title = 'Strong Motion Database ' + u'\u2013' + ' SIBER-RISK'
