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
from matplotlib.colors import Normalize

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

start_time = time.clock()

##########################################################################

pathANA = "../ANA/2017-09-27_01_00_00/anal00*"
pathGUESS = "../GUESS/2017-09-27_01_00_00/gues*"
variable2D = np.array(("T2", "Q2"))
VariablesCalc = "M10" #Magnitud del viento a 10 metros
variableInterp = np.array(("M", "T")) #M = magnitud del viento
nivelInterp = np.array((850., 500.))
figsize = (4, 3)
datetime = "2017-09-27_01_00_00"
###########################################################################

filesA = sorted(glob.glob(pathANA))
a_mean = Dataset(filesA[-1])
filesG = sorted(glob.glob(pathGUESS))
g_mean = Dataset(filesG[-1])

T2_diff = wrf.getvar(g_mean, "T2", meta = False) - wrf.getvar(a_mean, "T2", meta = False)
Q2_diff = wrf.getvar(g_mean, "Q2", meta = False) - wrf.getvar(a_mean, "Q2", meta = False)

M10_diff = wrf.g_uvmet.get_uvmet10_wspd_wdir(g_mean)[0,:,:] - wrf.g_uvmet.get_uvmet10_wspd_wdir(a_mean)[0,:,:]

p_g = wrf.g_pressure.get_pressure_hpa(g_mean)
p_a = wrf.g_pressure.get_pressure_hpa(a_mean)

M_g = wrf.g_uvmet.get_uvmet_wspd_wdir(g_mean)[0,:,:,:]
M_a = wrf.g_uvmet.get_uvmet_wspd_wdir(a_mean)[0,:,:,:]
M850_diff = wrf.interplevel(M_g, p_g, 850.) - wrf.interplevel(M_a, p_a, 850.)

T_g = wrf.g_temp.get_tk(g_mean)
T_a = wrf.g_temp.get_tk(a_mean)
T500_diff = wrf.interplevel(T_g, p_g, 500.) - wrf.interplevel(T_a, p_a, 500.)

labelsize = 4
fontsize = 4

#Leo la latitud y la longitud
lon = wrf.getvar(g_mean, "XLONG")
lat = wrf.getvar(g_mean, "XLAT")
llon = lon[0, 0] #Si cambia el dominio hay que cambiar esto.
llat = lat[0, 0]
rlon = lon[148, 98]
rlat = lat[148, 98]

#########################################################################
#Grafica la diferencia entre la media del guess y la media del analisis #
#########################################################################

fig = plt.figure(figsize = figsize, dpi = 300)
    
m = Basemap(resolution = 'i', llcrnrlon = llon, llcrnrlat = llat, urcrnrlon = rlon, urcrnrlat = rlat,
                projection = 'lcc', lat_1 = -31.847992, lat_2 = -31.848, lat_0 = -31.848, lon_0 = -61.537)
x, y = m(wrf.to_np(lon), wrf.to_np(lat))
        
plt.subplot(231)
cf = m.contourf(x, y, T2_diff, cmap = 'seismic', 
                norm = MidpointNormalize(midpoint = 0, vmin = np.min(T2_diff), vmax = np.max(T2_diff)))
#m.contour(x, y, min_var, colors = '#f4fbd2', linewidths = 0.1)
m.drawcoastlines(linewidth = 0.5)
m.drawcountries(linewidth = 0.5)
cbar = m.colorbar(cf)
cbar.ax.tick_params(labelsize = labelsize)
plt.title("T2m (K)", fontsize = fontsize)

plt.subplot(232)
cf = m.contourf(x, y, Q2_diff*1000, cmap = 'seismic'    )
#m.contour(x, y, max_var, colors = '#f4fbd2', linewidths = 0.1)
m.drawcoastlines(linewidth = 0.5)
m.drawcountries(linewidth = 0.5)
cbar = m.colorbar(cf)
cbar.ax.tick_params(labelsize = labelsize)
plt.title("Q2m (g/kg)", fontsize = fontsize)

plt.subplot(233)
cf = m.contourf(x, y, M10_diff, cmap = 'seismic',
                norm = MidpointNormalize(midpoint = 0, vmin = np.min(M10_diff)-1, vmax = np.max(M10_diff)))
#m.contour(x, y, var_mean, colors = '#f4fbd2', linewidths = 0.1)
m.drawcoastlines(linewidth = 0.5)
m.drawcountries(linewidth = 0.5)
cbar = m.colorbar(cf)
cbar.ax.tick_params(labelsize = labelsize)
plt.title("Magnitud del viento \n a 10 m", fontsize = fontsize)

plt.subplot(234)
cf = m.contourf(x, y, M850_diff, cmap = 'seismic', 
                norm = MidpointNormalize(midpoint = 0, vmin = np.min(M850_diff), vmax = np.max(M850_diff)))
#m.contour(x, y, spread_var, colors = '#f4fbd2', linewidths = 0.1)
m.drawcoastlines(linewidth = 0.5)
m.drawcountries(linewidth = 0.5)
cbar = m.colorbar(cf)
cbar.ax.tick_params(labelsize = labelsize)
plt.title("Magnitud del viento \n en 850 hpa", fontsize = fontsize) 

plt.subplot(235)
cf = m.contourf(x, y, T500_diff, cmap = 'seismic', 
                norm = MidpointNormalize(midpoint = 0, vmin = np.min(T500_diff), vmax = np.max(T500_diff)))
#m.contour(x, y, spread_var, colors = '#f4fbd2', linewidths = 0.1)
m.drawcoastlines(linewidth = 0.5)
m.drawcountries(linewidth = 0.5)
cbar = m.colorbar(cf)
cbar.ax.tick_params(labelsize = labelsize)
plt.title("Temperatura (K) \n 500 hpa", fontsize = fontsize)       

#plt.tight_layout(pad=0.4, w_pad=0.5) 
plt.subplots_adjust(left  = 0.08, right = 0.92, bottom = 0.1, top = 0.9, wspace = 0.01)    
fig.savefig("Diff_" + datetime + "_guess_analisis.pdf")
plt.close()       

print("Difference guess-analysis ready!")
print(time.clock() - start_time, "seconds")  




















