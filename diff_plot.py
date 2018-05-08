import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from wrf import to_np, getvar, smooth2d, get_cartopy, cartopy_xlim, cartopy_ylim, latlon_coords

data_path1="/home/jruiz/share/RRA/TMP/ANA/2017-09-27_01_00_00/"
data_path2="/home/jruiz/share/RRA/TMP/ANA/2017-09-27_01_00_00/"

# Open the NetCDF files
ncfile1 = Dataset(data_path1 + "anal00009")
ncfile2 = Dataset(data_path2 + "anal00011")

# Get the sea level pressure
T2_1 = getvar(ncfile1, "T" ,  meta=False)
T2_2 = getvar(ncfile2, "T" ,  meta=False)

lon = getvar(ncfile1, "XLONG" ,  meta=False)
lat = getvar(ncfile1, "XLAT"  ,  meta=False)


diff=T2_2 - T2_1

plt.figure()

plt.pcolor(lon,lat,diff[10,:,:])

plt.colorbar()

plt.savefig('ana_pert1.png')

plt.show()

