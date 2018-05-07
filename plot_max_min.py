#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  7 14:16:07 2018

@author: pao
"""

import time
import wrf
import numpy as np
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import glob
#from mpl_toolkits.axes_grid1 import AxesGrid

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

#Abro la media del ensamble
ens_mean = Dataset(path + "anal00061")
var_mean = wrf.getvar(ens_mean, variable)

#Elementos comunes a todas las figuras
levels = np.linspace(272, 310, 39)
ticks = np.linspace(272, 308, 10)

#Inicializo variables 
               
all = np.empty((149, 99, 60))

for f in range(len(files)):
    ncfile = Dataset(files[f])
    print(files[f])
    all[:,:,f] = wrf.to_np(wrf.getvar(ncfile, variable, meta = False))
    
           
print(time.clock() - start_time, "seconds")    

min_var = np.amin(all, 2)
max_var = np.amax(all, 2)
spread_var = np.std(all, 2)

print(time.clock() - start_time, "seconds")  

#Grafico

fig = plt.figure(figsize=(5, 6), dpi = 300)

m = Basemap(resolution = 'i', llcrnrlon = llon, llcrnrlat = llat, urcrnrlon = rlon, urcrnrlat = rlat,
                projection = 'lcc', lat_1 = -31.847992, lat_2 = -31.848, lat_0 = -31.848, lon_0 = -61.537)
x, y = m(wrf.to_np(lon), wrf.to_np(lat))
        
ax1 = plt.subplot(221)
cf = m.contourf(x, y, min_var, cmap = 'RdYlBu_r')
m.contour(x, y, min_var, colors = '#f4fbd2', linewidths = 0.1)
m.drawcoastlines(linewidth = 0.5)
m.drawcountries(linewidth = 0.5)
cbar = m.colorbar(cf)
cbar.ax.tick_params(labelsize = 8)
ax1.set_title("Mínimo")

ax2 = plt.subplot(222)
cf = m.contourf(x, y, max_var, cmap = 'RdYlBu_r')
m.contour(x, y, max_var, colors = '#f4fbd2', linewidths = 0.1)
m.drawcoastlines(linewidth = 0.5)
m.drawcountries(linewidth = 0.5)
cbar = m.colorbar(cf)
cbar.ax.tick_params(labelsize = 8)
ax2.set_title("Máximo")

ax3 = plt.subplot(223)
cf = m.contourf(x, y, var_mean, cmap = 'RdYlBu_r')
m.contour(x, y, var_mean, colors = '#f4fbd2', linewidths = 0.1)
m.drawcoastlines(linewidth = 0.5)
m.drawcountries(linewidth = 0.5)
cbar = m.colorbar(cf)
cbar.ax.tick_params(labelsize = 8)
ax3.set_title("Media")

ax4 = plt.subplot(224)
cf = m.contourf(x, y, spread_var, cmap = 'RdYlBu_r')
m.contour(x, y, spread_var, colors = '#f4fbd2', linewidths = 0.1)
m.drawcoastlines(linewidth = 0.5)
m.drawcountries(linewidth = 0.5)
cbar = m.colorbar(cf)
cbar.ax.tick_params(labelsize = 8)
ax4.set_title("Spread")       

#plt.tight_layout(pad=0.4, w_pad=0.5,) 
plt.subplots_adjust(left  = 0.08, right = 0.92, bottom = 0.1, top = 0.9, wspace = 0.001)    
fig.savefig("../" + variable + "_wrfout_resumen.png")
       
        
        
###########################################################################        
        
#cf = m.contourf(x, y, max_mean, cmap = 'RdYlBu_r')
#m.contour(x, y, max_mean, colors = '#f4fbd2', linewidths = 0.1)
#m.drawcoastlines(linewidth = 0.5)
#m.drawcountries(linewidth = 0.5)
#cbar = m.colorbar(cf)
#cbar.ax.tick_params(labelsize = 8)
#plt.show()        
        