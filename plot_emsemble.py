#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  4 09:03:56 2018

@author: pao
"""
#########################################################################
#      Plots        #
#########################################################################

import time
import wrf
import numpy as np
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import glob
from mpl_toolkits.axes_grid1 import AxesGrid

start_time = time.clock()

##########################################################################

path = "../2017-09-27_01_00_00/anal00*"
wd = "ANA" #ANA = analsis / GUESS = guess
variable = "T2"
figsize = (5, 4)
nrows_ncols = (2, 5) #Grilla de figuras
datetime = "2017-09-27_01_00_00"
###########################################################################

files = sorted(glob.glob(path))

#Abro el .nc que contiene la media del emsable
ncfile = Dataset(files[-1])
mean_ens = wrf.getvar(ncfile, variable, meta = False)

#Leo la latitud y la longitud
lon = wrf.getvar(ncfile, "XLONG")
lat = wrf.getvar(ncfile, "XLAT")
llon = lon[0, 0] #Si cambia el dominio hay que cambiar esto.
llat = lat[0, 0]
rlon = lon[148, 98]
rlat = lat[148, 98]

files = files[0:len(files)-1] #Saco el .nc con la media

#########################################################################
#      Grafica la variable elegida para cada miembro del emsable        #
#########################################################################

#Barra de  colores
ncfile = Dataset(files[0])
var = wrf.getvar(ncfile, variable, meta = False)
lev_max = np.max(var) + 0.01*np.max(var)
lev_min = np.min(var) + 0.01*np.min(var)
levels = np.around(np.linspace(lev_min, lev_max, 39), 2)
ticks = np.around(np.linspace(lev_min, lev_max, 10), 2)

#Grafico
fig = plt.figure(figsize = figsize, dpi = 300)
grid = AxesGrid(fig, 111,
                nrows_ncols = nrows_ncols,
                axes_pad = 0.05,
                cbar_mode = 'single',
                cbar_location = 'bottom',
                cbar_pad = 0.05
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
fig.savefig(wd + "_" + datetime + "_" + variable + "_miembros.pdf")
plt.close()

print("Plot varaible ready!")

#########################################################################
#   Grafica la diferencia entre cada miembro y la media del ensable     #
#########################################################################

#Barra de  colores
ncfile = Dataset(files[0])
var = wrf.getvar(ncfile, variable, meta = False) - mean_ens
scale = np.around((np.amax(np.abs(var))*1.1)/20, 2)
levels = np.linspace(-20, 20, 41)*scale
ticks = np.linspace(-20, 20, 11)*scale

#Grafico
fig = plt.figure(figsize = figsize, dpi = 300)
grid = AxesGrid(fig, 111,
                nrows_ncols = nrows_ncols,
                axes_pad = 0.05,
                cbar_mode = 'single',
                cbar_location = 'bottom',
                cbar_pad = 0.05
                )

for f in range(len(files)):

    ncfile = Dataset(files[f])
    print(files[f])
    var = wrf.to_np(wrf.getvar(ncfile, variable, meta = False))
    var = var - mean_ens

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
fig.savefig(wd + "_" + datetime + "_" +  variable + "_diff.pdf")
plt.close()
print("Plot_diff ready!")

print(time.clock() - start_time, "seconds")























