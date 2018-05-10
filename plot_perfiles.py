#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 11:33:21 2018

@author: pao
"""

#########################################################################
#                               Plots                                   #
#########################################################################

import time
import wrf
import numpy as np
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import glob

pathANA = "../ANA/2017-09-27_01_00_00/anal00*"
pathGUESS = "../GUESS/2017-09-27_01_00_00/gues*"
variable = np.array(("T", "Q", "U", "V", ))


filesA = sorted(glob.glob(pathANA))
filesG = sorted(glob.glob(pathGUESS))
forshape = Dataset(filesA[0])
forshape = wrf.getvar(forshape, "P", meta = False)

TA = TG = UA = UG = QA = QG = np.empty((len(filesA)-1, forshape.shape[0], forshape.shape[1], forshape.shape[2]))

for f in range(len(filesA)-1):
    g_ncfile = Dataset(filesG[f])
    a_ncfile = Dataset(filesA[f])
    
    TA[f,:,:,:] = wrf.g_temp.get_tk(a_ncfile)
    TG[f,:,:,:] = wrf.g_temp.get_tk(g_ncfile)
    TDiff = (TG - TA)**2
    
    UA[f,:,:,:] = wrf.g_uvmet.get_uvmet(a_ncfile)[0,:,:,:]
    UG[f,:,:,:] = wrf.g_uvmet.get_uvmet(g_ncfile)[0,:,:,:]
    UDiff = (UG - UA)**2
    
    QA[f,:,:,:] = wrf.getvar(a_ncfile, "QVAPOR", meta = False)
    QG[f,:,:,:] = wrf.getvar(g_ncfile, "QVAPOR", meta = False)
    QDiff = (QG - QA)**2
    
    
    
    
    #wrf.g_uvmet.get_uvmet(a_mean)[0,:,:,:]