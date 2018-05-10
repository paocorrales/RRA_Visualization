#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  7 14:16:07 2018

@author: pao
"""
###############################################################################
#   Calcula y grafica el maximo, minimo y spread del emsable mas su media     #
###############################################################################

import time
import wrf
import numpy as np
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import glob
from mpl_toolkits.axes_grid1 import AxesGrid

start_time = time.clock()

###############################################################################

wd = "ANA" #O GUESS
path = "../ANA/2017-09-27_01_00_00/anal00*"
variable2D = np.array(("T2", "Q2"))
variableInterp = np.array(("M", "T")) #M = magnitud del viento
nivelInterp = np.array((850., 500.))
figsize = (5, 4)
datetime = "2017-09-27_01_00_00"

###############################################################################
files = sorted(glob.glob(path))

#Abro un .nc para obtener la latitud y la longitud
ncfile = Dataset(files[0])
lon = wrf.getvar(ncfile, "XLONG")
lat = wrf.getvar(ncfile, "XLAT")
llon = lon[0, 0]
llat = lat[0, 0]
rlon = lon[148, 98]
rlat = lat[148, 98]

#                          Variables 2D
###############################################################################
for v in range(len(variable2D)):
    #Abro la media del ensamble
    ens_mean = Dataset(files[-1])
    var_mean = wrf.getvar(ens_mean, variable2D[v])

    #Inicializo variables 
    all_ens = np.empty((len(files)-1, lon.shape[0], lon.shape[1]))

    for f in range(len(files)-1):
        ncfile = Dataset(files[f])
        all_ens[f,:,:] = wrf.to_np(wrf.getvar(ncfile, variable2D[v], meta = False))
           
    #Calculo el maximo, minimo y spread del emsanble
    min_var = np.amin(all_ens, 0)
    max_var = np.amax(all_ens, 0)
    spread_var = np.std(all_ens, 0)  
    
    #Barra de  colores
    if variable2D[v] == "T2":
        levels = np.linspace(272, 310, 39)
        ticks = np.linspace(272, 308, 10)
    else:
        lev_max = np.max(max_var)
        lev_min = np.min(min_var)
        if lev_max == 0:
            levels = np.around(np.linspace(lev_min, 10, 39), 2)
            ticks = np.around(np.linspace(lev_min, 10, 10), 2)
        else:
            levels = np.around(np.linspace(lev_min, lev_max, 39), 2)
            ticks = np.around(np.linspace(lev_min, lev_max, 10), 2)
    
    #Grafico
    
    fig = plt.figure(figsize = figsize, dpi = 300)
    
    m = Basemap(resolution = 'i', llcrnrlon = llon, llcrnrlat = llat, urcrnrlon = rlon, urcrnrlat = rlat,
                    projection = 'lcc', lat_1 = -31.847992, lat_2 = -31.848, lat_0 = -31.848, lon_0 = -61.537)
    x, y = m(wrf.to_np(lon), wrf.to_np(lat))
            
    plt.subplot(221)
    cf = m.contourf(x, y, min_var, cmap = 'RdYlBu_r')
    #m.contour(x, y, min_var, colors = '#f4fbd2', linewidths = 0.1)
    m.drawcoastlines(linewidth = 0.5)
    m.drawcountries(linewidth = 0.5)
    cbar = m.colorbar(cf)
    cbar.ax.tick_params(labelsize = 8)
    plt.title("Mínimo")
    
    plt.subplot(222)
    cf = m.contourf(x, y, max_var, cmap = 'RdYlBu_r')
    #m.contour(x, y, max_var, colors = '#f4fbd2', linewidths = 0.1)
    m.drawcoastlines(linewidth = 0.5)
    m.drawcountries(linewidth = 0.5)
    cbar = m.colorbar(cf)
    cbar.ax.tick_params(labelsize = 8)
    plt.title("Máximo")
    
    plt.subplot(223)
    cf = m.contourf(x, y, var_mean, cmap = 'RdYlBu_r')
    #m.contour(x, y, var_mean, colors = '#f4fbd2', linewidths = 0.1)
    m.drawcoastlines(linewidth = 0.5)
    m.drawcountries(linewidth = 0.5)
    cbar = m.colorbar(cf)
    cbar.ax.tick_params(labelsize = 8)
    plt.title("Media")
    
    plt.subplot(224)
    cf = m.contourf(x, y, spread_var, cmap = 'RdYlBu_r')
    #m.contour(x, y, spread_var, colors = '#f4fbd2', linewidths = 0.1)
    m.drawcoastlines(linewidth = 0.5)
    m.drawcountries(linewidth = 0.5)
    cbar = m.colorbar(cf)
    cbar.ax.tick_params(labelsize = 8)
    plt.title("Spread")       
    
    #plt.tight_layout(pad=0.4, w_pad=0.5,) 
    plt.subplots_adjust(left  = 0.08, right = 0.92, bottom = 0.1, top = 0.9, wspace = 0.001)    
    fig.savefig(wd + "_" + datetime + "_" + variable2D[v] + "_max_min.pdf")
    plt.close()       
    print(variable2D[v] + " ready!")
    print(time.clock() - start_time, "seconds")  
        
###############################################################################        
          
#                      Variables Interpoladas a un nivel
###############################################################################
for v in range(len(variableInterp)):
    #Abro la media del ensamble
    ens_mean = Dataset(files[-1])
    p_mean = wrf.g_pressure.get_pressure_hpa(ens_mean)
    
    if variableInterp[v] == "M":
        var_mean = wrf.g_uvmet.get_uvmet_wspd_wdir(ens_mean)[0,:,:,:]
        var_mean = wrf.interplevel(var_mean, p_mean, nivelInterp[v])
    else:
        var_mean = wrf.getvar(ens_mean, variableInterp[v], meta = False)
        var_mean = wrf.interplevel(var_mean, p_mean, nivelInterp[v])

    #Inicializo variables 
    all_ens = np.empty((len(files)-1, lon.shape[0], lon.shape[1]))

    for f in range(len(files)-1):
        ncfile = Dataset(files[f])
        p =  wrf.g_pressure.get_pressure_hpa(ncfile)
        if variableInterp[v] == "M":
            tmp = wrf.g_uvmet.get_uvmet_wspd_wdir(ncfile)[0,:,:,:]
            all_ens[f,:,:] = wrf.interplevel(tmp, p, nivelInterp[v])
        elif variableInterp[v] == "T":
            tmp = wrf.g_temp.get_tk(ncfile)
            all_ens[f,:,:] = wrf.interplevel(tmp, p, nivelInterp[v])
        else:
            tmp = wrf.to_np(wrf.getvar(ncfile, variableInterp[v], meta = False))
            all_ens[f,:,:] = wrf.interplevel(tmp, p, nivelInterp[v])
           
    #Calculo el maximo, minimo y spread del emsanble
    min_var = np.amin(all_ens, 0)
    max_var = np.amax(all_ens, 0)
    spread_var = np.std(all_ens, 0)  
    
    #Barra de  colores
    if variable2D[v] == "T2":
        levels = np.linspace(272, 310, 39)
        ticks = np.linspace(272, 308, 10)
    else:
        lev_max = np.max(max_var)
        lev_min = np.min(min_var)
        if lev_max == 0:
            levels = np.around(np.linspace(lev_min, 10, 39), 2)
            ticks = np.around(np.linspace(lev_min, 10, 10), 2)
        else:
            levels = np.around(np.linspace(lev_min, lev_max, 39), 2)
            ticks = np.around(np.linspace(lev_min, lev_max, 10), 2)
    
    #Grafico
    
    fig = plt.figure(figsize = figsize, dpi = 300)
    
    m = Basemap(resolution = 'i', llcrnrlon = llon, llcrnrlat = llat, urcrnrlon = rlon, urcrnrlat = rlat,
                    projection = 'lcc', lat_1 = -31.847992, lat_2 = -31.848, lat_0 = -31.848, lon_0 = -61.537)
    x, y = m(wrf.to_np(lon), wrf.to_np(lat))
            
    plt.subplot(221)
    cf = m.contourf(x, y, min_var, cmap = 'RdYlBu_r')
    #m.contour(x, y, min_var, colors = '#f4fbd2', linewidths = 0.1)
    m.drawcoastlines(linewidth = 0.5)
    m.drawcountries(linewidth = 0.5)
    cbar = m.colorbar(cf)
    cbar.ax.tick_params(labelsize = 8)
    plt.title("Mínimo")
    
    plt.subplot(222)
    cf = m.contourf(x, y, max_var, cmap = 'RdYlBu_r')
    #m.contour(x, y, max_var, colors = '#f4fbd2', linewidths = 0.1)
    m.drawcoastlines(linewidth = 0.5)
    m.drawcountries(linewidth = 0.5)
    cbar = m.colorbar(cf)
    cbar.ax.tick_params(labelsize = 8)
    plt.title("Máximo")
    
    plt.subplot(223)
    cf = m.contourf(x, y, var_mean, cmap = 'RdYlBu_r')
    #m.contour(x, y, var_mean, colors = '#f4fbd2', linewidths = 0.1)
    m.drawcoastlines(linewidth = 0.5)
    m.drawcountries(linewidth = 0.5)
    cbar = m.colorbar(cf)
    cbar.ax.tick_params(labelsize = 8)
    plt.title("Media")
    
    plt.subplot(224)
    cf = m.contourf(x, y, spread_var, cmap = 'RdYlBu_r')
    #m.contour(x, y, spread_var, colors = '#f4fbd2', linewidths = 0.1)
    m.drawcoastlines(linewidth = 0.5)
    m.drawcountries(linewidth = 0.5)
    cbar = m.colorbar(cf)
    cbar.ax.tick_params(labelsize = 8)
    plt.title("Spread")       
    
    #plt.tight_layout(pad=0.4, w_pad=0.5,) 
    plt.subplots_adjust(left  = 0.08, right = 0.92, bottom = 0.1, top = 0.9, wspace = 0.001)    
    fig.savefig(wd + "_" + datetime + "_" + variableInterp[v] + "_" + str(nivelInterp[v]) + "_max_min.pdf")
    plt.close()       
    
    print(variable2D[v] + " in " + str(nivelInterp[v]) + " ready!")
    print(time.clock() - start_time, "seconds")  
        
###############################################################################

#                      Reflectividad (Colmax)
###############################################################################
#Abro la media del ensamble
ens_mean = Dataset(files[-1])
var_mean = wrf.g_dbz.get_max_dbz(ens_mean)

#Inicializo variables 
all_ens = np.empty((len(files)-1, lon.shape[0], lon.shape[1]))

for f in range(len(files)-1):
    ncfile = Dataset(files[f])
    all_ens[f,:,:] = wrf.g_dbz.get_max_dbz(ncfile)
       
#Calculo el maximo, minimo y spread del emsanble
min_var = np.amin(all_ens, 0)
max_var = np.amax(all_ens, 0)
spread_var = np.std(all_ens, 0)  

#Barra de  colores

lev_max = np.max(max_var)
lev_min = np.min(min_var)
if lev_max == 0:
    levels = np.around(np.linspace(lev_min, 10, 39), 2)
    ticks = np.around(np.linspace(lev_min, 10, 10), 2)
else:
    levels = np.around(np.linspace(lev_min, lev_max, 39), 2)
    ticks = np.around(np.linspace(lev_min, lev_max, 10), 2)
    
#Grafico

fig = plt.figure(figsize = figsize, dpi = 300)

m = Basemap(resolution = 'i', llcrnrlon = llon, llcrnrlat = llat, urcrnrlon = rlon, urcrnrlat = rlat,
                projection = 'lcc', lat_1 = -31.847992, lat_2 = -31.848, lat_0 = -31.848, lon_0 = -61.537)
x, y = m(wrf.to_np(lon), wrf.to_np(lat))
        
plt.subplot(221)
cf = m.contourf(x, y, min_var, cmap = 'RdYlBu_r')
#m.contour(x, y, min_var, colors = '#f4fbd2', linewidths = 0.1)
m.drawcoastlines(linewidth = 0.5)
m.drawcountries(linewidth = 0.5)
cbar = m.colorbar(cf)
cbar.ax.tick_params(labelsize = 8)
plt.title("Mínimo")

plt.subplot(222)
cf = m.contourf(x, y, max_var, cmap = 'RdYlBu_r')
#m.contour(x, y, max_var, colors = '#f4fbd2', linewidths = 0.1)
m.drawcoastlines(linewidth = 0.5)
m.drawcountries(linewidth = 0.5)
cbar = m.colorbar(cf)
cbar.ax.tick_params(labelsize = 8)
plt.title("Máximo")

plt.subplot(223)
cf = m.contourf(x, y, var_mean, cmap = 'RdYlBu_r')
#m.contour(x, y, var_mean, colors = '#f4fbd2', linewidths = 0.1)
m.drawcoastlines(linewidth = 0.5)
m.drawcountries(linewidth = 0.5)
cbar = m.colorbar(cf)
cbar.ax.tick_params(labelsize = 8)
plt.title("Media")

plt.subplot(224)
cf = m.contourf(x, y, spread_var, cmap = 'RdYlBu_r')
#m.contour(x, y, spread_var, colors = '#f4fbd2', linewidths = 0.1)
m.drawcoastlines(linewidth = 0.5)
m.drawcountries(linewidth = 0.5)
cbar = m.colorbar(cf)
cbar.ax.tick_params(labelsize = 8)
plt.title("Spread")       

#plt.tight_layout(pad=0.4, w_pad=0.5,) 
plt.subplots_adjust(left  = 0.08, right = 0.92, bottom = 0.1, top = 0.9, wspace = 0.001)    
fig.savefig(wd + "_" + datetime + "_reflectividad_COLMAX_max_min.pdf")
plt.close()       

print("Reflectivity ready!")
print(time.clock() - start_time, "seconds")  
        
###############################################################################                
          
        
        
