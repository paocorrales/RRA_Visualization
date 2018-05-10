#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  4 09:03:56 2018

@author: pao
"""
#########################################################################
#      Grafica la pp acum  elegida para cada miembro del emsable        #
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
variable = "RAINNC"
figsize = (5, 4)
nrows_ncols = (2, 5) #Grilla de figuras
datetime = "2017-09-27_01_00_00"
###########################################################################

files = sorted(glob.glob(path))

#Abro un .nc para obtener la latitud y la longitud
ncfile = Dataset(files[0])
lon = wrf.getvar(ncfile, "XLONG")
lat = wrf.getvar(ncfile, "XLAT")
llon = lon[0, 0] #Si cambia el dominio hay que cambiar esto.
llat = lat[0, 0]
rlon = lon[148, 98]
rlat = lat[148, 98]

#Leo los datos

files = files[0:len(files)-1]
alldata = np.empty((len(files), lon.shape[0], lon.shape[1]))

for f in range(len(files)):
    ncfile = Dataset(files[f])
    print(files[f])
    alldata[f,:,:] = wrf.to_np(wrf.getvar(ncfile, variable, meta = False)) + wrf.to_np(wrf.getvar(ncfile, "RAINC", meta = False))

#Barra de  colores
lev_max = np.max(alldata)
lev_min = 0
if lev_max == 0:
    levels = np.around(np.linspace(lev_min, 10, 39), 2)
    ticks = np.around(np.linspace(lev_min, 10, 10), 2)
else:
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

    m = Basemap(resolution = 'i', llcrnrlon = llon, llcrnrlat = llat, urcrnrlon = rlon, urcrnrlat = rlat,
                projection = 'lcc', lat_1 = -31.847992, lat_2 = -31.848, lat_0 = -31.848, lon_0 = -61.537,
                ax = grid[f])
    x, y = m(wrf.to_np(lon), wrf.to_np(lat))
    cf = grid[f].contourf(x, y, alldata[f,:,:], levels = levels, cmap = 'Blues')
    #grid[f].contour(x, y, var, levels = levels, colors = '#f4fbd2', linewidths = 0.1)
    m.drawcoastlines(linewidth = 0.5)
    m.drawcountries(linewidth = 0.5)
    cbar = grid.cbar_axes[0].colorbar(cf, ticks = ticks, boundaries=[-5] + levels + [5], extend = 'both')
    cbar.ax.tick_params(labelsize = 8)
    for cax in grid.cbar_axes:
        cax.toggle_label(True)
        
plt.subplots_adjust(left  = 0.01, right = 0.99, bottom = 0.01, top = 0.99, wspace = 0.00001)    
fig.savefig(wd + "_" + datetime + "_" + variable + "_miembros.pdf")

plt.close()
print("Plot variable ready!")
print(time.clock() - start_time, "seconds")























