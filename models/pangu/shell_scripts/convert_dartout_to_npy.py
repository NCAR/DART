import sys
import xarray as xr
import numpy as np
from datetime import datetime

# Function to calculate the average of horizontal neighbors
def average_horizontal_neighbors(array, var, lat, lon, level=None):
    neighbors = []
    lat_max, lon_max = array.shape[-2], array.shape[-1]
    
    # Collect the values of the 8 horizontal neighbors
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            if i == 0 and j == 0:
                continue
            lat_idx, lon_idx = lat + i, lon + j
            if 0 <= lat_idx < lat_max and 0 <= lon_idx < lon_max:
                if level is not None:
                    neighbor_value = array[var, level, lat_idx, lon_idx]
                else:
                    neighbor_value = array[var, lat_idx, lon_idx]
                if not np.isnan(neighbor_value):
                    neighbors.append(neighbor_value)
    
    # Calculate the average, if there are valid neighbors
    if neighbors:
        return np.mean(neighbors)
    else:
        # If no valid neighbors, return NaN (or any other default value)
        return np.nan

# Replace NaN values with the average of their horizontal neighbors
def replace_nan_with_neighbors(array):
    is_nan = np.isnan(array)
    nan_indices = np.argwhere(is_nan)
    
    for index in nan_indices:
        if array.ndim == 3:  # (var, lat, lon)
            var, lat, lon = index
            array[var, lat, lon] = average_horizontal_neighbors(array, var, lat, lon)
        elif array.ndim == 4:  # (var, level, lat, lon)
            var, level, lat, lon = index
            array[var, level, lat, lon] = average_horizontal_neighbors(array, var, lat, lon, level)
    
    # Handling any remaining NaNs (e.g., all neighbors were NaN)
    remaining_nans = np.isnan(array)
    array[remaining_nans] = np.nanmean(array)  # Example: fill with overall mean (excluding NaNs)
    
    return array
	
# list paths
out_path = sys.argv[3]
save_date = sys.argv[1]
ninst = sys.argv[2]
savename_sfc = out_path+'/input_surface'+ninst+'.postassim-'+save_date+'.npy'
savename_upper = out_path+'/input_upper'+ninst+'.postassim-'+save_date+'.npy'

# read data

path = 'pgoutput'+ninst+'.nc'
ds_dart = xr.open_dataset(path)
path = 'pginput'+ninst+'.nc'

path =  out_path+'/output_surface'+ninst+'.forecast-'+save_date+'.npy'  
sfc_template = np.load(path)
path =  out_path+'/output_upper'+ninst+'.forecast-'+save_date+'.npy'   
upper_template = np.load(path)
phb_upper_template = np.mean(upper_template[0, :, ::-1, :], axis=(-1, -2))

# !!!!!!!!!!
# era 5 data has the meridional coordinate in reverse!
# don't forget to change it back!

# save surface data
# surface states are not assimilated, use the forecast state as output
# maybe just copy the forecast.npy to output.npy?

sfc_arr = np.zeros_like(sfc_template)
sfc_arr[0,:] = sfc_template[0,:]
sfc_arr[1,:] = sfc_template[1,:]
sfc_arr[2,:] = sfc_template[2,:]
sfc_arr[3,:] = sfc_template[3,:]
sfc_arr = replace_nan_with_neighbors(sfc_arr)

np.save(savename_sfc, sfc_arr)

# save upper level data
# PH is not assimilated, use the forecast state as output
upper_arr = np.zeros_like(upper_template)
upper_arr[0,:] = upper_template[0, :]
upper_arr[1,:] = ds_dart.Q[0,:,::-1,:].values
upper_arr[2,:] = ds_dart.T[0,:,::-1,:].values
upper_arr[3,:] = ds_dart.U[0,:,::-1,:].values
upper_arr[4,:] = ds_dart.V[0,:,::-1,:].values
upper_arr = replace_nan_with_neighbors(upper_arr)

np.save(savename_upper, upper_arr)
