import xarray as xr

ENS=80
ens=[]
for n in range(1,64):
  ens.append('mem'+'{:04}'.format(n))
for n in range(65,ENS+1):
  ens.append('mem'+'{:04}'.format(n))

data_vars=['aicen','vicen','vsnon']
NCAT=5

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

for enum,e in enumerate(ens):
  file = e+'/dart_restart.nc'
#  print(file)
  ds = xr.open_dataset(file)
  ds = ds[data_vars]
  for n in range(1,NCAT+1):
    ds['aice'+'{:02}'.format(n)]=ds.aicen.isel(ncat=n-1)
    ds['vice'+'{:02}'.format(n)]=ds.vicen.isel(ncat=n-1)
    ds['vsno'+'{:02}'.format(n)]=ds.vsnon.isel(ncat=n-1)
  ds=ds.drop(data_vars)
  ds.to_netcdf(file,mode='a',encoding=encd)


