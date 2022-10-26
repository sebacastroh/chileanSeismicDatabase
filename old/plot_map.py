# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 19:37:20 2020

@author: sebac
"""
import os
import shapefile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

import matplotlib.cm as cm
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

plt.close('all')

plt.ion()

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Times New Roman'
plt.rcParams['font.size'] = 14
plt.rcParams['axes.titlepad'] = 10
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{times}'

data_path = '.'#'D:\\Google Drive\\Trabajo\\SIBER-RISK Strong Motion Database\\data\\paper'
figure_path = '.'#'D:\\Google Drive\\Trabajo\\SIBER-RISK Strong Motion Database\paper\\figures'
df = pd.read_excel(os.path.join(data_path, 'flatFile_corrected_v2.xlsx'))

sf = shapefile.Reader(os.path.join(data_path, 'ne_10m_admin_0_countries'))
shp = sf.shape(2)

points = np.array(shp.points)

#resolution = 'small0512px'
resolution = 'large4096px'

savefig = True

arcgis_url1 = 'http://server.arcgisonline.com/arcgis/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}.jpg'
arcgis_url2 = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Physical_Map/MapServer/tile/{z}/{y}/{x}.jpg'

tiles1 = cimgt.GoogleTiles(url=arcgis_url1)
tiles2 = cimgt.GoogleTiles(url=arcgis_url2)

fig = plt.figure(figsize=(10.62, 9.82))#figsize=(15.84, 6.336))

a1 = plt.subplot2grid((3,2), (0,0), rowspan=2, fig=fig)#fig.add_subplot(1,2,1)

unique, index = np.unique(df['Earthquake Name'], return_index=True)
order = np.argsort(df['Magnitude [Mw]'][index])

a1.hist(df['Magnitude [Mw]'][index], bins=10, range=[4., 9.], edgecolor='k')
a1.set_xlabel('Magnitude [Mw]')
a1.set_ylabel('Number of events')
a1.legend(['Events processed until\nSeptember 1, 2020'])

#axins = zoomed_inset_axes(a1, 2., loc=7) # zoom-factor: 2.5, location: upper-left
axins = inset_axes(a1, 1.8, 1.8, loc=1, bbox_to_anchor=(0.46, 0.71), bbox_transform=a1.figure.transFigure)
axins.hist(df['Magnitude [Mw]'][index], bins=7, range=[5.5, 9.], edgecolor='k')
x1, x2, y1, y2 = 5.9, 9.1, 0., 80. # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2) # apply the y-limits
axins.xaxis.set_visible('False')
axins.yaxis.set_visible('False')
axins.set_xticks([6, 7, 8, 9])

mark_inset(a1, axins, loc1=3, loc2=4, fc="none", ec="0.5")

a2 = plt.subplot2grid((3,4), (0,2), rowspan=2, fig=fig, projection=ccrs.PlateCarree())#fig.add_subplot(1,4,3, projection=ccrs.PlateCarree())
a2.set_extent([-76., -64., -56, -14], crs=ccrs.PlateCarree())
# a2.background_img(name="natural-earth-1", resolution=resolution)
# a2.background_img(name="esri", resolution='large8192px')

a2.add_image(tiles1, 8)
a2.add_image(tiles2, 8, alpha=0.65)

print('Percentage of events with a moment magnitude lower than 6: %0.2f%%' %(100.*np.count_nonzero(df['Magnitude [Mw]'][index] < 6.)/len(df['Magnitude [Mw]'][index])))
print('Events equal or larger than 7: %i' %np.count_nonzero(df['Magnitude [Mw]'][index] >= 7.))

import matplotlib as mpl

boundaries=[4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9.]

cmaplist = [cm.jet((i-4.)/5.) for i in boundaries]

mag_cm = mpl.colors.LinearSegmentedColormap.from_list(
        'Magnitude cmap', cmaplist, N=len(boundaries))

for i in range(len(shp.parts)-1):
    a2.plot(points[shp.parts[i]:shp.parts[i+1],0], points[shp.parts[i]:shp.parts[i+1],1], 'k', lw=0.5, zorder=2)
a2.plot(points[shp.parts[i+1]:,0], points[shp.parts[i+1]:,1], 'k', lw=0.5, zorder=2)
a2.scatter(df['Hypocenter longitude'][index[order]], df['Hypocenter latitude'][index[order]], s=1.5, c=df['Magnitude [Mw]'][index[order]], cmap=mag_cm, vmin=4., vmax=9., zorder=3)
a2.set_aspect('equal')
clb2 = plt.colorbar(a2.collections[0], boundaries=boundaries)
clb2.ax.set_ylabel('Magnitude [Mw]', rotation=270, labelpad=15)

a2.set_xticks([-75, -70, -65], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
a2.xaxis.set_major_formatter(lon_formatter)

a2.set_yticks([-55, -50, -45, -40, -35, -30, -25, -20, -15], crs=ccrs.PlateCarree())
lat_formatter = LatitudeFormatter()
a2.yaxis.set_major_formatter(lat_formatter)

a3 = plt.subplot2grid((3,4), (0,3), rowspan=2, fig=fig, projection=ccrs.PlateCarree())#fig.add_subplot(1,4,4, projection=ccrs.PlateCarree())
a3.set_extent([-76., -64., -56, -14], crs=ccrs.PlateCarree())
# a3.background_img(name="natural-earth-1", resolution=resolution)
# a3.background_img(name="esri", resolution='large8192px')

a3.add_image(tiles1, 8)
a3.add_image(tiles2, 8, alpha=0.65)

unique, index, counts = np.unique(df['Station name'], return_index=True, return_counts=True)
order = np.argsort(counts)

boundaries= [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280]

cmaplist = [cm.jet(i/280) for i in boundaries]

events_cm = mpl.colors.LinearSegmentedColormap.from_list(
        'Events cmap', cmaplist, N=len(boundaries))

for i in range(len(shp.parts)-1):
    a3.plot(points[shp.parts[i]:shp.parts[i+1],0], points[shp.parts[i]:shp.parts[i+1],1], 'k', lw=0.5, zorder=2)
a3.plot(points[shp.parts[i+1]:,0], points[shp.parts[i+1]:,1], 'k', lw=0.5, zorder=2)
a3.scatter(df['Station longitude'][index[order]], df['Station latitude'][index[order]], s=1.5, c=counts[order], cmap=events_cm, vmin=counts.min(), vmax=counts.max(), zorder=3)
a3.set_aspect('equal')
clb3 = plt.colorbar(a3.collections[0], boundaries=boundaries)
clb3.ax.set_ylabel('Number of events by station', rotation=270, labelpad=15)

a3.set_xticks([-75, -70, -65], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
a3.xaxis.set_major_formatter(lon_formatter)

#a3.set_yticks([-55, -50, -45, -40, -35, -30, -25, -20, -15], crs=ccrs.PlateCarree())
#lat_formatter = LatitudeFormatter()
#a3.yaxis.set_major_formatter(lat_formatter)
a3.set_yticks([])

a1.set_title('(a)')
a2.set_title('(b)')
a3.set_title('(c)')

# hspace = 0.3*1.6
# fig.subplots_adjust(hspace)
fig.savefig(os.path.join(figure_path, 'histogram_database.pdf'), bbox_inches='tight', pad_inches=0, dpi=300)
