#!/usr/bin/env python3
import numpy as np
from netCDF4 import Dataset

ccyy = 2019  ##date of data file
mm = 6
dd = 14
otype = 1 ##defined in obs_kind_mod, index for MPD_ABSOLUTE_HUMIDITY
n_sites = 5  ##number of MPD sites
obs_err = 0.001 ###set observation error 1g/m3

for site in range(1, n_sites+1):
  f = Dataset('wv_mpd0{}.{:02d}{:02d}{:02d}.dafinal.Python.nc'.format(site, int(ccyy%100), mm, dd))
  dat = np.array(f.variables['Absolute_Humidity'])
  dat = dat * 0.001  ##convert to kg/m3
  dat_mask = np.array(f.variables['Absolute_Humidity_mask'])

  nt, nz = dat.shape
  tt = np.array(f.variables['time_Absolute_Humidity'])
  zz = np.array(f.variables['range_Absolute_Humidity'])

  lat = np.array(f.variables['lidar_latitude'])
  lon = np.array(f.variables['lidar_longitude'])
  if(lon<0):
    lon = 360.0 + lon
  elev = np.array(f.variables['lidar_elevation'])

  for t in range(nt-1):
    ii_sum = int(tt[t]/60)
    hh = int(ii_sum/60)
    ii = int(ii_sum%60)
    for z in range(nz):
      if(dat_mask[t,z]==0 and dat[t,z]>=0):  ##QC mask and remove negative data
        textfile = "work/{:04d}{:02d}{:02d}{:02d}{:02d}".format(ccyy, mm, dd, hh, ii)
        f1 = open(textfile, "a")
        string = "{:2d} {:9.5f} {:9.5f} {:8.1f} {:4d} {:02d} {:02d} {:02d} {:02d} {:02d} {:e} {:e}\n".format(otype, lat, lon, elev+zz[z], ccyy, mm, dd, hh, ii, 0, dat[t, z], obs_err)
        f1.write(string)
        f1.close()
