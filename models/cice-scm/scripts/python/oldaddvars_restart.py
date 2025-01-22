import xarray as xr

file = 'dart_restart.nc'
ds = xr.open_dataset(file)
data_vars=['aicen','vicen','vsnon']
ds = xr.open_dataset(file)
ds = ds[data_vars]
NCAT=5
for n in range(1,NCAT+1):
  ds['aice'+'{:02}'.format(n)]=ds.aicen.isel(ncat=n-1)
  ds['vice'+'{:02}'.format(n)]=ds.vicen.isel(ncat=n-1)
  ds['vsno'+'{:02}'.format(n)]=ds.vsnon.isel(ncat=n-1)
ds=ds.drop(data_vars)
encd ={'aice01': {'_FillValue': None},
       'aice02': {'_FillValue': None},
       'aice03': {'_FillValue': None},
       'aice04': {'_FillValue': None},
       'aice05': {'_FillValue': None},
       'vice01': {'_FillValue': None},
       'vice02': {'_FillValue': None},
       'vice03': {'_FillValue': None},
       'vice04': {'_FillValue': None},
       'vice05': {'_FillValue': None},
       'vsno01': {'_FillValue': None},
       'vsno02': {'_FillValue': None},
       'vsno03': {'_FillValue': None},
       'vsno04': {'_FillValue': None},
       'vsno05': {'_FillValue': None}}

ds.to_netcdf(file,mode='a',encoding= encd)


