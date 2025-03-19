import xarray as xr
import numpy as np

ENS=80
ens=[]
for n in range(1,64):
  ens.append('mem'+'{:04}'.format(n))
for n in range(65,ENS+1):
  ens.append('mem'+'{:04}'.format(n))

data_vars=['aicen','vicen','vsnon']

NCAT=5
expanded_vars=[]
for n in range(1,NCAT+1):
  expanded_vars.append('vice'+'{:02}'.format(n))
  expanded_vars.append('aice'+'{:02}'.format(n))
  expanded_vars.append('vsno'+'{:02}'.format(n))

#print(expanded_vars)

encd ={'aicen': {'_FillValue': None},
       'vicen': {'_FillValue': None},
       'vsnon': {'_FillValue': None}}

for enum,e in enumerate(ens):
  #print(e)
  
  file = e+'/restart_state.nc'
  ds = xr.open_dataset(file)
  aicen=np.zeros((5,4))
  vicen=np.zeros((5,4))
  vsnon=np.zeros((5,4))
  for n in range(1,NCAT+1):
    aicen[n-1,:]=ds['aice'+'{:02}'.format(n)].values
    vicen[n-1,:]=ds['vice'+'{:02}'.format(n)].values
    vsnon[n-1,:]=ds['vsno'+'{:02}'.format(n)].values
  ds['aicen']=xr.DataArray(data=aicen, dims=['ncat','ni'])
  ds['vicen']=xr.DataArray(data=vicen, dims=['ncat','ni'])
  ds['vsnon']=xr.DataArray(data=vsnon, dims=['ncat','ni'])
  ds=ds.drop(expanded_vars)
  ds.to_netcdf(file,encoding=encd)
  #print(ds)

#  comd = 'ncks -x -v'+expanded_vars+' -O '+e+'/dart_restart.nc '+e+'/dart_restart.nc'
#  os.system(comd)

