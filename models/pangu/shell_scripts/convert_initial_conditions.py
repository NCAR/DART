import xarray as xr
import numpy as np
from scipy.ndimage import zoom
%matplotlib inline



num_instances = 3
start_datetime = "2024-01-20-00"
data_dir = "era5_data/"

ds_sl = xr.open_dataset(f'{data_dir}/era5_ens_sfc.nc')
ds_pl = xr.open_dataset(f'{data_dir}/era5_ens_upper.nc')

# Calculate zoom factors, modify depend on your data resolution you downloaded
# (13 / data_lev_num, 721 / data_lat_num, 1440 / data_lon_num)
zoom_factor_sl = (721 / 361, 1440 / 720)
zoom_factor_pl = (1, 721 / 361, 1440 / 720)


for mem in range(num_instances):
    strmem = str(mem).zfill(4)
    savename_sfc = f"input_surface_{strmem}.forecast-{start_datetime}.npy"
    savename_upper = f"input_upper_{strmem}.forecast-{start_datetime}.npy"

    data_sfc = np.zeros((4,721,1440))
    
    data = ds_sl.msl[0,mem,:].values.astype(np.float32)
    resized_data = zoom(data, zoom_factor_sl)
    data_sfc[0,:] = resized_data
    
    data = ds_sl.u10[0,mem,:].values.astype(np.float32)
    resized_data = zoom(data, zoom_factor_sl)
    data_sfc[1,:] = resized_data
    
    data = ds_sl.v10[0,mem,:].values.astype(np.float32)
    resized_data = zoom(data, zoom_factor_sl)
    data_sfc[2,:] = resized_data
    
    data = ds_sl.t2m[0,mem,:].values.astype(np.float32)
    resized_data = zoom(data, zoom_factor_sl)
    data_sfc[3,:] = resized_data
    np.save(savename_sfc, data_sfc)
    
    data_upper = np.zeros((5,13,721,1440))
    
    # Calculate zoom factors
    zoom_factor_pl = (1, 721 / 361, 1440 / 720)
    
    data = ds_pl.z[0,mem,::-1,:,:].values.astype(np.float32)
    resized_data = zoom(data, zoom_factor_pl)
    data_upper[0,:] = resized_data
    
    data = ds_pl.q[0,mem,::-1,:,:].values.astype(np.float32)
    resized_data = zoom(data, zoom_factor_pl)
    data_upper[1,:] = resized_data
    
    data = ds_pl.t[0,mem,::-1,:,:].values.astype(np.float32)
    resized_data = zoom(data, zoom_factor_pl)
    data_upper[2,:] = resized_data
    
    data = ds_pl.u[0,mem,::-1,:,:].values.astype(np.float32)
    resized_data = zoom(data, zoom_factor_pl)
    data_upper[3,:] = resized_data
    
    data = ds_pl.v[0,mem,::-1,:,:].values.astype(np.float32)
    resized_data = zoom(data, zoom_factor_pl)
    data_upper[4,:] = resized_data
    
    np.save(savename_upper, data_upper)