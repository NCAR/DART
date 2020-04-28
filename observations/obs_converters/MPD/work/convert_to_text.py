# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

import numpy as np
from netCDF4 import Dataset

CCYY = 2019  ##date of data file
MM = 6
DD = 14
NUM_SITES = 5  ##number of MPD sites
OBS_ERR = 0.001 ###set observation error 1g/m3

for site in range(1, NUM_SITES + 1):
    f = Dataset('wv_mpd0{}.{:02d}{:02d}{:02d}.dafinal.Python.nc'.format(site, int(CCYY % 100), MM, DD))
    dat = np.array(f.variables['Absolute_Humidity'])
    dat = dat * 0.001  ##convert to kg/m3
    dat_mask = np.array(f.variables['Absolute_Humidity_mask'])

    nt, nz = dat.shape
    tt = np.array(f.variables['time_Absolute_Humidity'])
    zz = np.array(f.variables['range_Absolute_Humidity'])

    lat = np.array(f.variables['lidar_latitude'])
    lon = np.array(f.variables['lidar_longitude'])
    if(lon < 0):
        lon = 360.0 + lon
    elev = np.array(f.variables['lidar_elevation'])

    for t in range(nt-1):
        ii_sum = int(tt[t] / 60)
        hh = int(ii_sum / 60)
        ii = int(ii_sum % 60)
        for z in range(nz):
            if(dat_mask[t, z] == 0 and dat[t, z] >= 0):  ##QC mask and remove negative data
                textfile = "work/{:04d}{:02d}{:02d}{:02d}{:02d}".format(CCYY, MM, DD, hh, ii)
                f1 = open(textfile, "a")
                string = "{:9.5f} {:9.5f} {:8.1f} {:4d} {:02d} {:02d} {:02d} {:02d} {:02d} {:e} {:e}\n".format(lat, lon, elev + zz[z], CCYY, MM, DD, hh, ii, 0, dat[t, z], OBS_ERR)
                f1.write(string)
                f1.close()
