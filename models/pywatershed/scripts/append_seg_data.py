import os
import numpy  as np
import pandas as pd
from netCDF4 import Dataset

# Define file paths
# Change the file path. Ideally, these files should live in the 'work' directory
filepath = '/Users/gharamti/Documents/MATLAB/Derecho/HydroDART/pywatershed/'
table_f, source_f, dest_f = 'river_midpoints.csv', 'parameters_dis_seg.nc', 'parameters_dis_seg_app.nc'
NTWRK_FILE, PARAM_FILE, CSV_TABLE = map(lambda f: os.path.join(filepath, f), [dest_f, source_f, table_f])

if not os.path.isfile(PARAM_FILE) or not os.path.isfile(CSV_TABLE):
    raise FileNotFoundError(f"Missing file: {PARAM_FILE} or {CSV_TABLE}")

# Read CSV and sort by index
# Could read headers to check for the right columns
data = pd.read_csv(CSV_TABLE).iloc[:, [4, 5, 6]].values
data = data[np.argsort(data[:, 0])] # sorted array

lats, lons = data[:, 1], data[:, 2]

# Read NetCDF dimensions & variables
with Dataset(PARAM_FILE, 'r') as nc:
    n_seg, gages         = len(nc.dimensions['nsegment']), len(nc.dimensions['npoigages'])
    to_index, poi_gauges = nc.variables['tosegment'][:].astype(int), nc.variables['poi_gage_id'][:].astype(int)

# Compute upstream links
fromIndices, num_up_links, fromIndsStart, fromIndsEnd = [], np.zeros(n_seg, int), np.zeros(n_seg, int), np.zeros(n_seg, int)

k = 0
for i in range(n_seg):
    upstream = np.where(to_index == (i + 1))[0] + 1  # Convert to 1-based index
    num_up_links[i] = len(upstream)
    if upstream.size > 0:
        fromIndsStart[i], fromIndsEnd[i] = k + 1, k + len(upstream)
        fromIndices.extend(upstream)
        k += len(upstream)

# Create new NetCDF file
if os.path.isfile(NTWRK_FILE):
    os.remove(NTWRK_FILE)

with Dataset(NTWRK_FILE, 'w', format='NETCDF4') as nc_new, Dataset(PARAM_FILE, 'r') as nc_old:
    # Copy dimensions
    for dim_name, dim in nc_old.dimensions.items():
        nc_new.createDimension(dim_name, len(dim) if not dim.isunlimited() else None)
    nc_new.createDimension('index', len(fromIndices))  

    # Copy variables
    for var_name, var in nc_old.variables.items():
        new_var    = nc_new.createVariable(var_name, var.datatype, var.dimensions)
        new_var[:] = var[:]

    # Ignore attributes for now

    # Create new variables
    nc_new.createVariable('fromIndsStart', 'i8', ('nsegment',))[:]  = fromIndsStart
    nc_new.createVariable('fromIndsEnd'  , 'i8', ('nsegment',))[:]  = fromIndsEnd
    nc_new.createVariable('fromIndices'  , 'i8', ('index',))[:]     = np.array(fromIndices, int)
    nc_new.createVariable('num_up_links' , 'i8', ('nsegment',))[:]  = num_up_links
    nc_new.createVariable('seg_lats'     , 'f8', ('nsegment',))[:]  = lats
    nc_new.createVariable('seg_lons'     , 'f8', ('nsegment',))[:]  = lons
    nc_new.createVariable('poi_gauges'   , 'i8', ('npoigages',))[:] = poi_gauges

