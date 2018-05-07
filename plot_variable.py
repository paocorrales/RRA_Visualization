#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  4 09:03:56 2018

@author: pao
"""
import time
import wrf
import numpy as np
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import glob
from mpl_toolkits.axes_grid1 import AxesGrid

start_time = time.clock()

path = "../2017-09-23_01_00_00/"
variable = "T2"

files = sorted(glob.glob(path + "anal00*"))
files = files[0:len(files)-1]

#Abro un .nc para obtener la latitud y la longitud
ncfile = Dataset(path + "anal00001")
lon = wrf.getvar(ncfile, "XLONG")
lat = wrf.getvar(ncfile, "XLAT")
llon = lon[0, 0]
llat = lat[0, 0]
rlon = lon[148, 98]
rlat = lat[148, 98]


#Elementos comunes a todas las figuras
levels = np.linspace(272, 310, 39)
ticks = np.linspace(272, 308, 10)

#Grafico
fig = plt.figure(figsize=(5, 13), dpi = 300)

grid = AxesGrid(fig, 111,
                nrows_ncols=(10, 6),
                axes_pad=0.05,
                cbar_mode='single',
                cbar_location='bottom',
                cbar_pad=0.05
                )

for f in range(len(files)):

    ncfile = Dataset(files[f])
    print(files[f])
    var = wrf.to_np(wrf.getvar(ncfile, variable, meta = False))

    m = Basemap(resolution = 'i', llcrnrlon = llon, llcrnrlat = llat, urcrnrlon = rlon, urcrnrlat = rlat,
                projection = 'lcc', lat_1 = -31.847992, lat_2 = -31.848, lat_0 = -31.848, lon_0 = -61.537,
                ax = grid[f])
    x, y = m(wrf.to_np(lon), wrf.to_np(lat))
    cf = grid[f].contourf(x, y, var, levels = levels, cmap = 'RdYlBu_r')
    grid[f].contour(x, y, var, levels = levels, colors = '#f4fbd2', linewidths = 0.1)
    m.drawcoastlines(linewidth = 0.5)
    m.drawcountries(linewidth = 0.5)
    cbar = grid.cbar_axes[0].colorbar(cf, ticks = ticks)
    cbar.ax.tick_params(labelsize = 8)
    for cax in grid.cbar_axes:
        cax.toggle_label(True)
        
plt.subplots_adjust(left  = 0.01, right = 0.99, bottom = 0.01, top = 0.99, wspace = 0.00001)    
fig.savefig(variable + "_wrfout.png")
plt.close()

print(time.clock() - start_time, "seconds")

#ncfile1 = Dataset("../2017-09-27_23_00_00/anal00001")
#ncfile2 = Dataset("../2017-09-27_23_00_00/anal00002")
#
#var1 = wrf.to_np(wrf.getvar(ncfile1, "PSFC", meta = False))
#var2 = wrf.to_np(wrf.getvar(ncfile2, "PSFC", meta = False))
#
#
#ncfile1 = Dataset("../fecha1/anal00001")
#t1 = wrf.getvar(ncfile1, "T")
#
#ncfile2 = Dataset("../fecha1/anal00002")
#t2 = wrf.getvar(ncfile2, "T")
#
#diff = t1-t2
#
#    plt.setp(ax.spines.values(), linewidth=0.6)
#
#    m = Basemap(resolution = 'i', llcrnrlon = llon, llcrnrlat = llat, urcrnrlon = rlon, urcrnrlat = rlat,
#                projection = 'lcc', lat_1 = -31.847992, lat_2 = -31.848, lat_0 = -31.848, lon_0 = -61.537)
#    x, y = m(wrf.to_np(lon), wrf.to_np(lat))
#    cf = m.contourf(x, y, diff[1,:,:], cmap = 'RdYlBu_r')
#   cbar = m.colorbar(cf, drawedges = True, location='bottom')
#    #cbar.ax.tick_params(labelsize = 3)
#    m.contour(x, y, diff[1,:,:], colors = '#f4fbd2', linewidths = 0.1)
#    m.drawcoastlines(linewidth = 0.5)
#    m.drawcountries(linewidth = 0.5)
#    plt.show()





























