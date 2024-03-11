import os
import numpy as np
import pandas as pd

###################
## Bokeh modules ##
###################
from bokeh.layouts import layout

from bokeh.models import DatePicker, Tabs, TabPanel
from bokeh.models.widgets import CheckboxGroup, DataTable, NumericInput, Select, TableColumn

from bokeh.plotting import figure

# Paths
currentDir = os.path.dirname(__file__)
libPath    = os.path.abspath(os.path.join(currentDir, '..', '..', 'lib'))
dataPath   = os.path.abspath(os.path.join(currentDir, '..', '..', '..', 'data'))

################
##  Settings  ##
################
g          = 9.81 # m/s**2
npoints    = 250
allowed    = False
pwdith     = 1
pheight    = 250
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
eTypes = ['All'] + list(map(lambda x: x.capitalize(), np.sort(flatfile['Event type'].unique()).tolist()))

# Magnitudes
min_mag = flatfile['Magnitude [Mw]'].min()
max_mag = flatfile['Magnitude [Mw]'].max()

# Station codes
sCodes = ['All'] + np.sort(flatfile['Station code'].unique()).tolist()

# Seismic events
seismic_events = list(reversed(flatfile['Earthquake Name'].unique().tolist()))
station_codes  = np.sort(flatfile[flatfile['Earthquake Name'] == seismic_events[0]]['Station code']).tolist()

#########################
##  Filter components  ##
#########################
filter_since = DatePicker(min_date=min_date, max_date=max_date, value=min_date,
    title='Show events since', sizing_mode='stretch_width', width=pwdith)

filter_until = DatePicker(min_date=min_date, max_date=max_date, value=max_date,
    title='Show events until', sizing_mode='stretch_width', width=pwdith)

filter_eType = Select(title='Event type', value=eTypes[0], options=eTypes,
    sizing_mode='stretch_width', width=pwdith*2)

filter_minMw = NumericInput(title='Show events larger or equal than', value=min_mag,
    sizing_mode='stretch_width', width=pwdith)

filter_maxMw = NumericInput(title='Show events smaller or equal than', value=max_mag,
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

####################
##  Plots inputs  ##
####################
options_ta = NumericInput(title='Initial period [s]', value=0.,
    sizing_mode='stretch_width', width=pwdith)

options_tb = NumericInput(title='Ending period [s]', value=1.5,
    sizing_mode='stretch_width', width=pwdith)

options_xi = NumericInput(title='Damping ratio', value=0.05,
    sizing_mode='stretch_width', width=pwdith)

options_ax = CheckboxGroup(labels=['Logarithmic X axis', 'Logarithmic Y axis'],
    sizing_mode='stretch_width', width=pwdith)

options_gr = CheckboxGroup(labels=['X axis minor grid', 'Y axis minor grid'],
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
#TODO: check if acceleration lengths are equal
def load_event(event_name):
    with np.load(os.path.join(dataPath, 'seismicDatabase', 'npz', event_name + '.npz'), allow_pickle=True) as f:
        event = {}
        for key, value in f.items():
            if not key.startswith('st'):
                continue

            station      = value.item()
            station_code = station['station_code']

            stations_codes.append(station_code)
            event[station_code] = station

    stations_codes = list(sorted(station_codes))

    return event, station_codes

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

    return spectrum_x, spectrum_y, spectrum_rotd50, spectrum_rod0, spectrum_rotd100, spectrum_geom

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
        husid_plots.extend([t, ia])

    return husid_plots

def compute_cav_plot(station):
    cav_plots = []
    dt = station['dt']
    for i in range(3):
        n   = len(station['acc_uncorrected_%i' %(i+1)])
        t   = np.linspace(0., (n-1)*dt, n)
        cav = spi.cumtrapz(np.abs(station['acc_uncorrected_%i' %(i+1)]), t, initial=0.)
        cav_plots.extend([t, cav])

    return cav_plots

#######################
## Update functions  ##
#######################

def update_event(attrname, old, new):
    global event
    
    if old == new:
        return

    event, station_codes = load_event(new)

    if filter_scode.value in station_codes:
        select_station.value = ''
        select_station.value = filter_scode.value
    else:
        select_station.value = station_codes[0]
    
def update_station(attrname, old, new):
    global event

    station = event[select_station.value]

    vel, dis        = compute_vel_dis(station)
    sa_spectra      = compute_sa_spectra(station, tn, xi)
    dva_spectra     = compute_dva_spectra(sa_spectra, tn)
    fourier_spectra = compute_fourier_spectra(station)
    husid_plot      = compute_husid_plot(station)
    cav_plot        = compute_cav_plot(station)

def update_spectra_options(attrname, old, new):
    global event

    if options_ta.value < 0. or options_ta.value >= options_tb.value:
        return

    if options_xi.value < 0.:
        return

    station = event[select_station.value]



    if option_ta.value <= 0.:
        Tn_lin = np.linspace(0., float(tb.value), npoints)
        Tn_log = np.hstack((0., np.logspace(np.log10(0.001), np.log10(float(tb.value)), npoints)))
    else:
        Tn_lin = np.linspace(float(ta.value), float(tb.value), npoints)
        Tn_log = np.logspace(np.log10(float(ta.value)), np.log10(float(tb.value)), npoints)
