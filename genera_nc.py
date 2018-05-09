#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  9 10:44:02 2018

@author: pao
"""

import time
import wrf
import numpy as np
from netCDF4 import Dataset
import glob

path = "../2017-09-27_01_00_00/anal00*"

files = sorted(glob.glob(path))

ncfile = Dataset(files[0])

lon = wrf.getvar(ncfile, "XLONG")
lat = wrf.getvar(ncfile, "XLAT")

psfc = wrf.getvar(ncfile, "PSFC")

#Creat a list of ncfiles
nclist = []
for f in range(len(files)):
    nclist.append(Dataset(files[f])) 

p = wrf.g_pressure.get_pressure(nclist, method = 'join')/1000

start_time = time.clock()

path = "../2017-09-27_01_00_00/anal00*"
files = sorted(glob.glob(path))

interp_levels = [1000, 975, 925, 850, 700, 500, 300, 200]
tmp = np.empty((11,8,149,99))
for f in range(len(files)):
    ncfile =  Dataset(files[f])
    
    p = wrf.g_pressure.get_pressure(ncfile)/100
    t = wrf.g_temp.get_tk(ncfile)
    
    tmp[f,:,:,:] = wrf.vinterp(ncfile,
                   field = t,
                   vert_coord = "p",
                   interp_levels = interp_levels,
                   extrapolate = False,
                   field_type = "tk",
                   log_p = False)
    #tk[f,:,:,:] = tmp


file = 'test.nc'
newfile = Dataset("test.nc", "w", format="NETCDF4")

#create dimensions
newfile.createDimension('ens', len(files))
newfile.createDimension('level', len(interp_levels))
newfile.createDimension('lon', 99)
newfile.createDimension('lat', 149)

#define variables
lev = newfile.createVariable("LEV", 'f8', ('level',))
lat = newfile.createVariable('LAT', 'f8', ('lat','lon'))
lon = newfile.createVariable('LON', 'f8', ('lat','lon'))
tk = newfile.createVariable('TK','f8', ('ens','level','lat','lon'))

lev[:] = interp_levels
lon[:] = wrf.getvar(ncfile, "XLONG", meta = False)
lat[:] = wrf.getvar(ncfile, "XLAT", meta = False)
tk[:] = tmp

#close ncfile
newfile.close()
    
print(time.clock() - start_time, "seconds")