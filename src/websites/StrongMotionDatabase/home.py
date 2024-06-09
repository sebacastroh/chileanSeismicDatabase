import os
import sys
import json
import random
import string
import datetime
import numpy as np
import pandas as pd
import scipy.fftpack as spf
import scipy.integrate as spi
from pyproj import Transformer

###################
## Bokeh modules ##
###################
from bokeh.events import ButtonClick

from bokeh.io import curdoc

from bokeh.layouts import column, grid, layout

from bokeh.models import ColumnDataSource, CustomJS, DatePicker, Tabs, TabPanel
from bokeh.models.widgets import Button, CheckboxGroup, DataTable, Div, NumericInput, Select, TableColumn

from bokeh.palettes import Spectral6

from bokeh.plotting import figure

from jinja2 import Environment, FileSystemLoader

# Paths
currentDir = os.path.dirname(__file__)
dataPath   = os.path.abspath(os.path.join(currentDir, '..', '..', '..', 'data'))
libPath    = os.path.abspath(os.path.join(currentDir, '..', '..', 'lib'))
srcPath    = os.path.abspath(os.path.join(currentDir, '..', '..'))

if not libPath in sys.path:
    sys.path.append(libPath)

import seismic

################
##  Settings  ##
################
g                  = 9.81 # m/s**2
npoints            = 250
allowed            = False
pwdith             = 1
pheight            = 250
lwidth             = 2
nTheta             = 91
nx                 = 0
ny                 = 45
formats            = ['Matlab (MAT-binary, *.mat)', 'Python (Numpy, *.npz)']
axis_types         = [('linear', 'linear'), ('linear', 'log'), ('log', 'linear'), ('log', 'log')]
spectrum_labels    = ['Component 1', 'Component 2', 'RotD50', 'RotD0', 'RotD100', 'Geometric mean']
transformer        = Transformer.from_crs('EPSG:4326', 'EPSG:3857')
station_attributes = [('starttime', 'Start time'),
    ('magnitude', 'Magnitude [Mw]'),
    ('hypocenter_lon', 'Longitude hypocenter'),
    ('hypocenter_lat', 'Latitude hypocenter'),
    ('event_type', 'Event type'),
    ('depth', 'Depth [km]'),
    ('station_name', 'Station name'),
    ('station_lon', 'Longitude station'),
    ('station_lat', 'Latitude station'),
    ('dt', 'Time interval (dt) [s]'),
    ('Rhypo', 'Hypocentral distance [km]'),
    ('Repi', 'Epicentral distance [km]'),
    ('Rrup', 'Rupture distance [km]'),
    ('Rjb', 'Joyner-Boore distance [km]'),
    ('vs30', 'Vs30 [m/s]'),
    ('hvsr', 'HVSR'),
    ('azimuth', 'Azimuth [o]'),
    ('last_update', 'Last update')]
seismic_component  = 'acc_filtered_'

########################
##  Global variables  ##
########################
event   = {}
ta      = None
tb      = None
xi      = None
plotted = False

############
##  Divs  ##
############
div_filter  = Div(text='<h2>Filter options</h2>', sizing_mode='stretch_width', width=pwdith)

div_records = Div(text='<h2>Strong Motion Database</h2>\nTo download several files at once please ' +
    'click <a href="downloadManager" target=' +
    '"_blank">here</a>. Examples to read the database files are available for <a href="examplePython.py"' +
    ' target="_blank">Python</a> and <a href="exampleMatlab.m" target="_blank">Matlab</a>.<br/>',
    sizing_mode='stretch_width', width=pwdith)

div_secPlot = Div(text='<h2>Plots</h2>\nFill the options below and then click the button to generate the plots for the selected station.<br/>',
    sizing_mode='stretch_width', width=pwdith, styles={'margin-bottom': '15px'})

div_plots   = Div(text='<b>TIP:</b> You can hide/show the curves by pressing their names at the legend box!',
    sizing_mode='stretch_width', width=pwdith)

div_secSumm = Div(text='<h2>Summary</h2>',
    sizing_mode='stretch_width', width=pwdith)

#########################
##  Filter components  ##
#########################
filter_since = DatePicker(title='Show events since', sizing_mode='stretch_width', width=pwdith)

filter_until = DatePicker(title='Show events until', sizing_mode='stretch_width', width=pwdith)

filter_eType = Select(title='Event type',sizing_mode='stretch_width', width=pwdith)

filter_minMw = NumericInput(title='Show events larger or equal than', mode='float', sizing_mode='stretch_width', width=pwdith)

filter_maxMw = NumericInput(title='Show events smaller or equal than', mode='float', sizing_mode='stretch_width', width=pwdith)

filter_sCode = Select(title='Recorded by station', sizing_mode='stretch_width', width=pwdith)

#########################
##  Select components  ##
#########################
select_event   = Select(title='Seismic events', sizing_mode='stretch_width', width=pwdith)

select_station = Select(title='Stations', sizing_mode='stretch_width', width=pwdith)

select_format = Select(title='File format', sizing_mode='stretch_width', width=pwdith)

###############
##  Buttons  ##
###############
button_filter   = Button(label='Apply filters', button_type='primary',
    sizing_mode='stretch_width', width=pwdith, align='end')

button_download = Button(label='Download', button_type='success',
    sizing_mode='stretch_width', width=pwdith, align='end')

button_plots    = Button(label='Generate plots', button_type='primary',
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
    tooltips=tooltips_a, height=pheight, sizing_mode='stretch_width', width=pwdith) for i in range(3)]

######################
##  Velocity plots  ##
######################
plots_vel = [figure(title='Velocity', x_axis_label='Time [s]', y_axis_label='Velocity [m/s]',
    tooltips=tooltips_v, height=pheight, sizing_mode='stretch_width', width=pwdith) for i in range(3)]

##########################
##  Displacement plots  ##
##########################
plots_dis = [figure(title='Displacement', x_axis_label='Time [s]', y_axis_label='Displacement [m]',
    tooltips=tooltips_d, height=pheight, sizing_mode='stretch_width', width=pwdith) for i in range(3)]

#########################
##  Sa Spectrum plots  ##
#########################
plots_sa = [figure(title='Pseudo-acceleration Response Spectrum',
    x_axis_label='Period [s]', y_axis_label='Spectral Acceleration [g]',
    tooltips=tooltips_sa, sizing_mode='stretch_width', width=pwdith,
    x_axis_type=axis_type[0], y_axis_type=axis_type[1]) for axis_type in axis_types]

##########################
##  DVA Spectrum plots  ##
##########################
plots_dva = [figure(title='Combined D-V-A Spectrum',
    x_axis_label='Period [s]', y_axis_label='Spectral Velocity [m/s]',
    tooltips=tooltips_sv, sizing_mode='stretch_width', width=pwdith,
    x_axis_type=axis_type[0], y_axis_type=axis_type[1]) for axis_type in axis_types]

##############################
##  Fourier Spectrum plots  ##
##############################
plots_fourier = [figure(title='Fourier Spectrum',
    x_axis_label='Frequency [Hz]', y_axis_label='Amplitude [m/s]',
    tooltips=tooltips_f, sizing_mode='stretch_width', width=pwdith,
    x_axis_type=axis_type[0], y_axis_type=axis_type[1]) for axis_type in axis_types]

###################
##  Husid plots  ##
###################
plots_husid = [figure(title='Husid Plot',
    x_axis_label='Time [s]', y_axis_label='Arias intensity [m/s]',
    tooltips=tooltips_ai, sizing_mode='stretch_width', width=pwdith,
    x_axis_type=axis_type[0], y_axis_type=axis_type[1]) for axis_type in axis_types]

#################
##  CAV plots  ##
#################
plots_cav = [figure(title='Cumulative Absolute Plot',
    x_axis_label='Time [s]', y_axis_label='CAV [m/s]',
    tooltips=tooltips_cav, sizing_mode='stretch_width', width=pwdith,
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
    plots_acc[i].line('x', 'y', source=sources_acc[i], legend_label='')
    plots_vel[i].line('x', 'y', source=sources_vel[i], legend_label='')
    plots_dis[i].line('x', 'y', source=sources_dis[i], legend_label='')

for i in range(len(axis_types)):
    for j in range(3):
        plots_fourier[i].line('x', 'y', source=sources_fourier[j], legend_label='Component %i' %(j+1), line_width=lwidth, color=Spectral6[j])
        plots_husid[i].line(  'x', 'y', source=sources_husid[j]  , legend_label='Component %i' %(j+1), line_width=lwidth, color=Spectral6[j])
        plots_cav[i].line(    'x', 'y', source=sources_cav[j]    , legend_label='Component %i' %(j+1), line_width=lwidth, color=Spectral6[j])

    for j in range(6):
        plots_sa[i].line( 'x', 'y', source=sources_sa[j] , legend_label=spectrum_labels[j], line_width=lwidth, color=Spectral6[j])
        plots_dva[i].line('x', 'y', source=sources_dva[j], legend_label=spectrum_labels[j], line_width=lwidth, color=Spectral6[j])

###################
##  Minor grids  ##
###################
for i in range(3):
    plots_acc[i].xgrid[0].minor_grid_line_color = '#e5e5e5'
    plots_acc[i].ygrid[0].minor_grid_line_color = '#e5e5e5'

    plots_vel[i].xgrid[0].minor_grid_line_color = '#e5e5e5'
    plots_vel[i].ygrid[0].minor_grid_line_color = '#e5e5e5'

    plots_dis[i].xgrid[0].minor_grid_line_color = '#e5e5e5'
    plots_dis[i].ygrid[0].minor_grid_line_color = '#e5e5e5'

###########################
##  Interactive legends  ##
###########################
for i in range(len(axis_types)):
    plots_sa[i].legend.click_policy      = 'hide'
    plots_dva[i].legend.click_policy     = 'hide'
    plots_fourier[i].legend.click_policy = 'hide'
    plots_husid[i].legend.click_policy   = 'hide'
    plots_cav[i].legend.click_policy     = 'hide'

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

inputs_plots = column([options_ta, options_tb, options_xi, button_plots, options_ax, options_gr, div_plots],
    sizing_mode='stretch_width', width=pwdith)

#############
##  Table  ##
#############
source_table  = ColumnDataSource(data=dict(fields=[], values=[]))
columns_table = [TableColumn(field='fields', title='Parameter'), TableColumn(field='values', title='Value')]
data_table    = DataTable(source=source_table, columns=columns_table,
    sizing_mode='stretch_width', width=pwdith, reorderable=False, height=2*pheight)

###########
##  Map  ##
###########
plot_map = figure(x_axis_type='mercator', y_axis_type='mercator', x_range=(0, 1), y_range=(0, 1),
    sizing_mode='stretch_width', width=pwdith)
plot_map.add_tile('CartoDB Positron')
source_hypo = ColumnDataSource(data=dict(lat=[], lon=[]))
source_sta  = ColumnDataSource(data=dict(lat=[], lon=[]))

plot_map.scatter(x='lon', y='lat', size=15, fill_color='blue'  , line_color='black', fill_alpha=0.8, source=source_sta , legend_label='Station')
plot_map.scatter(x='lon', y='lat', size=25, fill_color='yellow', line_color='black', fill_alpha=0.8, source=source_hypo, legend_label='Hypocenter', marker='star')
plot_map.legend.location = 'top_left'

############
##  Tabs  ##
############
panel_records = [TabPanel(child=layout([[plots_acc[i]], [plots_vel[i]], [plots_dis[i]]],
    sizing_mode='stretch_width', width=pwdith), title='Component %i' %(i+1)) for i in range(3)]

tabs_records  = Tabs(tabs=panel_records,
    sizing_mode='stretch_width', width=pwdith)

panel_plots = []
for i in range(len(axis_types)):
    panel_plots.append([
        TabPanel(child=plots_sa[i]     , title='Sa Spectra'),
        TabPanel(child=plots_dva[i]    , title='D-V-A Spectrum'),
        TabPanel(child=plots_fourier[i], title='Fourier Spectrum'),
        TabPanel(child=plots_husid[i]  , title='Husid Plot'),
        TabPanel(child=plots_cav[i]    , title='Cumulative Absolute Velocity')
    ])

tabs_plots    = Tabs(tabs=panel_plots[0],
    sizing_mode='stretch_width', width=pwdith)

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

            if not p_waves[select_event.value][station_code]['status']:
                continue

            stations_codes.append(station_code)
            event[station_code] = station

    stations_codes = list(sorted(stations_codes))

    return event, stations_codes

def compute_vel_dis(station):
    velocities    = []
    displacements = []
    dt            = station['dt']
    p_wave        = p_waves[select_event.value][select_station.value]['pos']
    for i in range(3):
        n = len(station[seismic_component + '%i' %(i+1)])
        if n == 0:
            velocities.append(np.empty(0))
            displacements.append(np.empty(0))
            continue

        t = np.linspace(0., (n-1)*dt, n)
        acc = station[seismic_component + '%i' %(i+1)]

        vel = spi.cumtrapz(acc, x=t, initial=0.)
        if p_wave > 0:
            vel -= vel[:p_wave].mean()

        dis = spi.cumtrapz(vel, x=t, initial=0.)
        if p_wave > 0:
            dis -= dis[:p_wave].mean()

        velocities.append(vel)
        displacements.append(dis)

    return velocities, displacements

def compute_sa_spectra(station, tn, xi):
    spectra = seismic.SpectraRot(station[seismic_component + '1'], station[seismic_component + '2'], station['dt'], tn, xi, nTheta)

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
        L  = len(station[seismic_component + '%i' %(i+1)])
        freq = Fs*np.linspace(0., L//2, L//2)/L
        Fa = 2.*np.abs(spf.fft(station[seismic_component + '%i' %(i+1)])*dt)[:L//2]
        fourier_spectra.extend([freq, Fa])

    return fourier_spectra

def compute_husid_plot(station):
    husid_plots = []
    dt = station['dt']
    for i in range(3):
        n  = len(station[seismic_component + '%i' %(i+1)])
        t  = np.linspace(0., (n-1)*dt, n)
        ia = np.pi/(2.*g)*spi.cumtrapz(station[seismic_component + '%i' %(i+1)]**2, t, initial=0.)
        husid_plots.append(ia)

    return husid_plots

def compute_cav_plot(station):
    cav_plots = []
    dt = station['dt']
    for i in range(3):
        n   = len(station[seismic_component + '%i' %(i+1)])
        t   = np.linspace(0., (n-1)*dt, n)
        cav = spi.cumtrapz(np.abs(station[seismic_component + '%i' %(i+1)]), t, initial=0.)
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
    global event, ta, tb, xi, plotted

    if select_station.value == '':
        return

    ta      = None
    tb      = None
    xi      = None
    plotted = False

    station  = event[select_station.value]
    vel, dis = compute_vel_dis(station)

    # Update source data
    for i in range(3):
        n  = len(station[seismic_component + '%i' %(i+1)])
        dt = station['dt']
        t  = np.linspace(0., (n-1)*dt, n)
        sources_acc[i].data = dict(x=t, y=station[seismic_component + '%i' %(i+1)]/g)
        sources_vel[i].data = dict(x=t, y=vel[i])
        sources_dis[i].data = dict(x=t, y=dis[i])

        sources_fourier[i].data = dict(x=[], y=[])
        sources_husid[i].data   = dict(x=[], y=[])
        sources_cav[i].data     = dict(x=[], y=[])

    for i in range(6):
        sources_sa[i].data  = dict(x=[], y=[])
        sources_dva[i].data = dict(x=[], y=[])

    # Update legend text
    for i in range(3):
        plots_acc[i].legend[0].items[0].update(label=dict(value = station['component_%i' %(i+1)]))
        plots_vel[i].legend[0].items[0].update(label=dict(value = station['component_%i' %(i+1)]))
        plots_dis[i].legend[0].items[0].update(label=dict(value = station['component_%i' %(i+1)]))

    for i in range(len(axis_types)):
        for j in range(3):
            plots_fourier[i].legend[0].items[j].update(label=dict(value = 'Component %i' %(j+1)))
            plots_husid[i].legend[0].items[j].update(  label=dict(value = 'Component %i' %(j+1)))
            plots_cav[i].legend[0].items[j].update(    label=dict(value = 'Component %i' %(j+1)))
    
        for j in range(2):
            plots_sa[i].legend[0].items[j].update( label=dict(value = 'Component %i' %(j+1)))
            plots_dva[i].legend[0].items[j].update(label=dict(value = 'Component %i' %(j+1)))

    # Update data table
    source_table.data = dict(fields=[attribute[1] for attribute in station_attributes],
        values=[station[attribute[0]] for attribute in station_attributes])

    # Update map plot
    hypo_x, hypo_y = transformer.transform(station['hypocenter_lat'], station['hypocenter_lon'])
    sta_x , sta_y  = transformer.transform(station['station_lat']   , station['station_lon'])

    xmean = 0.5*(hypo_x + sta_x)
    ymean = 0.5*(hypo_y + sta_y)

    dist = max([abs(hypo_x-xmean), abs(sta_x-xmean), abs(hypo_y-ymean), abs(sta_y-ymean)])

    xmin = xmean - dist - 30000
    xmax = xmean + dist + 30000

    ymin = ymean - dist - 30000
    ymax = ymean + dist + 30000

    plot_map.x_range.update(start=xmin, end=xmax)
    plot_map.y_range.update(start=ymin, end=ymax)

    source_hypo.data = dict(lat=[hypo_y], lon=[hypo_x])
    source_sta.data  = dict(lat=[sta_y] , lon=[sta_x])

def update_plots():
    global event, ta, tb, xi, plotted

    if (ta == options_ta.value) and (tb == options_tb.value) and (xi == options_xi.value):
        return

    ta = options_ta.value
    tb = options_tb.value
    xi = options_xi.value

    station = event[select_station.value]

    tn_log = np.logspace(np.log10(max(options_ta.value, 0.001)), np.log10(options_tb.value), int(npoints/2))
    tn_lin = np.linspace(options_ta.value, options_tb.value, int(npoints/2))
    tn     = np.unique(np.hstack((tn_log, tn_lin)))

    sa_spectra  = compute_sa_spectra(station, tn, options_xi.value)
    dva_spectra = compute_dva_spectra(sa_spectra, tn)

    if not plotted:
        fourier_spectra = compute_fourier_spectra(station)
        husid_plot      = compute_husid_plot(station)
        cav_plot        = compute_cav_plot(station)

    # Update source data
    if not plotted:
        for i in range(3):
            n  = len(station[seismic_component + '%i' %(i+1)])
            dt = station['dt']
            t  = np.linspace(0., (n-1)*dt, n)

            sources_fourier[i].data = dict(x=fourier_spectra[2*i], y=fourier_spectra[2*i+1])
            sources_husid[i].data   = dict(x=t, y=husid_plot[i])
            sources_cav[i].data     = dict(x=t, y=cav_plot[i])

    for i in range(6):
        sources_sa[i].data  = dict(x=tn, y=sa_spectra[i])
        sources_dva[i].data = dict(x=tn, y=dva_spectra[i])

    # Update legend text
    for i in range(len(axis_types)):
        if not plotted:
            for j in range(3):
                plots_fourier[i].legend[0].items[j].update(label=dict(value = station['component_%i' %(j+1)]))
                plots_husid[i].legend[0].items[j].update(  label=dict(value = station['component_%i' %(j+1)]))
                plots_cav[i].legend[0].items[j].update(    label=dict(value = station['component_%i' %(j+1)]))
    
        for j in range(2):
            plots_sa[i].legend[0].items[j].update( label=dict(value = station['component_%i' %(j+1)]))
            plots_dva[i].legend[0].items[j].update(label=dict(value = station['component_%i' %(j+1)]))

    plotted = True

def update_axis_options(attrname, old, new):
    if len(options_ax.active) == 0:
        ax = 0
    elif len(options_ax.active) == 2:
        ax = 3
    elif 0 in options_ax.active:
        ax = 2
    else:
        ax = 1

    tabs_plots.tabs = panel_plots[ax]

def update_grid_options(attrname, old, new):
    if 0 in options_gr.active:
        colorx = '#e5e5e5'
    else:
        colorx = None

    if 1 in options_gr.active:
        colory = '#e5e5e5'
    else:
        colory = None

    for i in range(len(axis_types)):
        plots_sa[i].xgrid[0].minor_grid_line_color = colorx
        plots_sa[i].ygrid[0].minor_grid_line_color = colory

        plots_dva[i].xgrid[0].minor_grid_line_color = colorx
        plots_dva[i].ygrid[0].minor_grid_line_color = colory

        plots_fourier[i].xgrid[0].minor_grid_line_color = colorx
        plots_fourier[i].ygrid[0].minor_grid_line_color = colory
        
        plots_husid[i].xgrid[0].minor_grid_line_color = colorx
        plots_husid[i].ygrid[0].minor_grid_line_color = colory

        plots_cav[i].xgrid[0].minor_grid_line_color = colorx
        plots_cav[i].ygrid[0].minor_grid_line_color = colory


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

Download = CustomJS(args=dict(source=select_event, extension=select_format),code="""
if (extension.value == 'Matlab (MAT-binary, *.mat)') {
    var filename = source.value.concat('.mat');
    var path = 'data/seismicDatabase/mat/';
} else {
    var filename = source.value.concat('.npz');
    var path = 'data/seismicDatabase/npz/';
}

var link = document.createElement('a');
link.setAttribute('download', filename);
link.href = '%s' + path + filename;
document.body.appendChild(link);
link.click();
link.remove();
""" %sys.argv[1])

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
button_plots.on_click(update_plots)
select_event.on_change('value', update_event)
select_station.on_change('value', update_station)
options_ax.on_change('active', update_axis_options)
options_gr.on_change('active', update_grid_options)
button_download.js_on_event(ButtonClick, Download)

################
##  Database  ##
################
with open(os.path.join(srcPath, 'data', 'p_waves.json')) as f:
    p_waves = json.load(f)

flatfile = pd.read_csv(os.path.join(dataPath, 'flatFile.csv'))
flatfile = flatfile[flatfile['Corrected records']]

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
select_event.update(value=seismic_events[0], options=seismic_events)
select_station.update(value=station_codes[0], options=station_codes)
select_format.update(value=formats[0], options=formats)

######################
##  Export website  ##
######################
_env = Environment(loader=FileSystemLoader('StrongMotionDatabase'))
FILE = _env.get_template('siberrisk_seismicdatabase.html')
curdoc().template = FILE

grid_filter = grid([[filter_since , filter_minMw  , filter_eType , button_filter  ] ,
                    [filter_until , filter_maxMw  , filter_sCode , None           ]], sizing_mode='stretch_width')
grid_select  = grid([[select_event, select_station, select_format, button_download]], sizing_mode='stretch_width')
grid_plots   = grid([[inputs_plots, tabs_plots]], sizing_mode='stretch_width')
grid_details = grid([[data_table  , plot_map]], sizing_mode='stretch_width')

grid_plots.children[0] = (grid_plots.children[0][0], 0, 0, 1, 1)
grid_plots.children[1] = (grid_plots.children[1][0], 0, 2, 1, 5)

curdoc().add_root(div_filter)
curdoc().add_root(grid_filter)
curdoc().add_root(div_records)
curdoc().add_root(grid_select)
curdoc().add_root(tabs_records)
curdoc().add_root(div_secPlot)
curdoc().add_root(grid_plots)
curdoc().add_root(div_secSumm)
curdoc().add_root(grid_details)

curdoc().title = 'Strong Motion Database ' + u'\u2013' + ' SIBER-RISK'
