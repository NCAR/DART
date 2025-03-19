# ------------------------------------------------------------------------- #
# IMPORT LIBRARIES                                                          #
# ------------------------------------------------------------------------- #
import xarray as xr
import os
import glob
import h5py
import numpy as np
import datetime

import geopandas as gpd
from shapely.geometry import Polygon
from shapely.geometry import Point
from scipy.interpolate import lagrange
from numpy.polynomial.polynomial import Polynomial

# ------------------------------------------------------------------------- #
# BEGIN FUNCTION DEFINITION                                                   #
# ------------------------------------------------------------------------- #

def determine_data_source(name, data_type):
    icesat2_names = ['ICESAT2', 'ICESAT-2', 'IS2', 'IS-2', 'Icesat2','Icesat-2','IceSat-2', 'icesat2','icesat-2','is2','is-2']

    if name in icesat2_names:
        if data_type in ['freeboard', 'fb']:
            print("You're converting freeboard data from NASA's ICESat-2 mission.")
        elif data_type in ['thickness','hi']:
            print("You're converting thickness data from NASA's ICESat-2 mission!")
        else:
            AssertionError('Data type not recognized.')
        read_is2_data(data_type)
    else:
        raise NotImplementedError
        

def read_is2_data(data_type, start_date = '2019.01.01', tracks = ['gt1r', 'gt1l', 'gt2r', 'gt2l', 'gt3r', 'gt3l'], 
                  input_path = '/glade/scratch/mollyw/external_data/ICESat-2/thickness/',
                  output_path = '/glade/scratch/mollyw/external_data/IS2_reduced/thickness/'):
    
    days_of_year = [31,29,31,30,31,30,31,31,30,31,30,31]
    year = start_date[0:4]
    dir_list = [name for name in sorted(os.listdir(input_path))]
    print(dir_list)
    
    if os.path.isdir(output_path) is False:
        print('Output directory presently does not exist! Generating directory now.')
        os.makedirs(output_path)
    
    if data_type in ['thickness', 'hi', 'sithick','SIT']:
        for dir_day in dir_list:
            files = sorted(glob.glob(input_path+dir_day+'/IS2SITDAT4_01*[0-9]*.nc'))
            if len(files) == 0: 
                FileNotFoundError('I did not find any IS2 thickness data in your system for this date. Please check if data is available in correct directories.')
            else:
                generate_thickness_netcdf(files, dir_day, output_path)

        check = os.listdir(output_path)
        if len(check)==0:
            FileNotFoundError('No files were converted! Please debug.')
        
    elif data_type in ['freeboard', 'fb']:
        for dir_day in dir_list:
            files = sorted(glob.glob(input_path+dir_day+'/ATL10-01*[0-9]*.h5'))
            if len(files) == 0:
                FileNotFoundError('I did not find any IS2 freeboard data in your system for this date. Please check if data is available in correct directories.')
            else:
                generate_freeboard_netcdf(files, tracks, dir_day, output_path)
                
        check = os.listdir(output_path)
        if len(check) == 0:
            FileNotFoundError('No files were converted! Please debug.')
            
    elif data_type in ['concentration', 'siconc','aice','SIC']:
        raise NotImplementedError('Not finished')
        
    elif data_type in ['motion']:
        raise NotImplementedError('Not currently supported')
        
    elif data_type in ['snow']:
        raise NotImplementedError('Not currently supported')
        
    elif data_type in ['ice_age', 'age']:
        raise NotImplementedError('Not currently supported')
        
    else:
        raise NotImplementedError("The type of data you wish to convert is not supported. Please contact DAReS if you wish to add support for your data type.")
        
        
def generate_freeboard_netcdf(files, tracks, dir_day, output_path):
    data_list = []
    for file in files:
        with h5py.File(file, mode='r') as f:
            for track in tracks:
                if len(f['%s/freeboard_beam_segment/' %track].keys()) < 2:
                    pass
                else:
                    # print(f['/%s/freeboard_beam_segment/latitude' %track])
                    # access the latitude
                    latvar = f['/%s/freeboard_beam_segment/latitude' %track]
                    # and rename, as a long list?
                    latitude = latvar[:]

                    # access the longitude
                    lonvar = f['/%s/freeboard_beam_segment/longitude' %track]
                    # and rename, as a long list?
                    longitude = lonvar[:]

                    # access the freeboard values 
                    dset_name = '/%s/freeboard_beam_segment/beam_fb_height' %track
                    datavar = f[dset_name]
                    # and rename as a long list
                    data = datavar[:]

                    # access the freeboard uncertainty
                    unc_var = f['/%s/freeboard_beam_segment/beam_fb_sigma' %track]
                    unc = unc_var[:]

                    # collect metadata
                    units = datavar.attrs['units']
                    long_name = datavar.attrs['long_name']
                    _FillValue = datavar.attrs['_FillValue']

                    # handle FillValue
                    data[data == _FillValue] = np.nan
                    unc[unc == _FillValue] = np.nan

                    # collect time information
                    timevar = f['/%s/freeboard_beam_segment/delta_time' %track]
                    time = timevar[:]

                    # make dataset, decode time to cftime object, and save in a list 
                    data_list.append(xr.decode_cf(xr.Dataset({"data": (["time"], data),
                                                              "error": (["time"], unc),
                                                              "lon": (["time"], longitude),
                                                              "lat": (["time"], latitude)
                                                             }, coords={"time": (['time'],
                                                                                 time,
                                                                                 {"units": "seconds since 2018-01-01"}
                                                                                )
                                                                       }
                                                            )
                                                 )
                                    )

    day_data = xr.concat(data_list, dim = 'time', compat = 'no_conflicts')
    day_data = day_data.dropna(dim = 'time', how = 'any')
    filename = 'IS2_fb_'+dir_day+'.nc'
    day_data.to_netcdf(output_path+filename)
    print(dir_day + ' done!')

def generate_thickness_netcdf(files, dir_day, output_path):
    data_list = []
    for file in files:
        ds = xr.open_dataset(file)

        year = int(dir_day[0:4])
        month = int(dir_day[5:7])
        day = int(dir_day[8:10])

        # access the latitude
        latitude = list(ds.latitude.values)

        # access the longitude
        longitude = list(ds.longitude.values)

        # access the thickness values 
        datavar = ds.ice_thickness
        data = list(datavar.values)

        # access the thickness uncertainty
        unc = list(ds.ice_thickness_unc.values)

        # collect time information
        time = [datetime.datetime(year,month,day) for i in range(0, len(data))]
                    
        # make dataset, decode time to cftime object, and save in a list 
        data_list.append(xr.Dataset({"data": (["time"], data),
                                     "error": (["time"], unc),
                                     "lon": (["time"], longitude),
                                     "lat": (["time"], latitude)}, 
                                     coords={"time": (['time'],time)}))

    day_data = xr.concat(data_list, dim = 'time', compat = 'no_conflicts')
    day_data = day_data.dropna(dim='time', how='any')
    filename = 'IS2_hi_'+dir_day+'.nc'
    day_data.to_netcdf(output_path+filename)
    print(dir_day + ' done!')

def bin_data(grid_file, data_file, reduce='mean', itd=False, write = False, output_path = None, thres=50):
    '''
    "bin_data" aggregates point data (satellite measurements, buoys, anything
    with a (lat,lon), into bins that correspond to locations on a grid; the
    function also reduces the data into one measurement per grid cell.

    Inputs:
    (1) grid_file: a .shp file containing the geometries of the grid you'd like to use
    (2) data_file: a netcdf file containing the (lat,lon) time, and data 
        values for each point measurement. There are no limits on the number
        of type of data included in the file (though strings may get tricky
        when using the reduction function (means, medians, etc). Pre-process
        as you dare.
    (3) reduce: the type of aggregate function used to reduce the data. Options
        include any possible "aggfunc" argument to geopandas.dissolve.
    (4) itd: flag for whether the "ITD" should be computed. ITD is here used
        to represent any histogram type distribution, it does not have to be
        related to ice thickness. Default is False. THIS IS NOT YET IMPLEMENTED.
    (5) thres: threshold for the number of points within a grid cell considered
        that is considered feasible for aggregation into a distribution.
        Default is 50.

    Outputs:
    (1) reduced_data: a GeoPandas DataFrame of containing the Polygon
        information for each gridcell, the reduced measurements included in
        data_file, and the locations and time of each measurement.

    Author: Molly M. Wieringa (mmw906@uw.edu, Univ. Washington ICG)
    '''
    # open the grid file
    grid = gpd.read_file(grid_file)
    crs='epsg:3411'

    # open the data file
    ds = xr.open_dataset(data_file)

    # convert the point data data into Points
    df = ds.to_dataframe().reset_index()

    # generate Points with an iterator
    points = [Point(x, y) for x, y in zip(df['lon'], df['lat'])]

    # convert to a GeoDataFrame
    point_data = gpd.GeoDataFrame(df, geometry=points, crs=crs)

    # combine the grid data and the point data into on GeoDF and drop NaNs
    binned_data = gpd.sjoin(grid, point_data, how='left',
                            predicate='contains').dropna()

    # pull and assign each bin a "cell_number". This is identical to the index,
    # but accessible to operate over for grouping functions
    index_list = list(binned_data.index.values)
    binned_data["cell_number"] = index_list

    # reduce the data with the a given function
    reduced_data = binned_data.dissolve(by='cell_number', aggfunc=reduce)

    if itd is True:
        raise NotImplementedError

        # now do the distribution: group by cell_number and aggregate in lists
        # temp = binned_data.groupby('cell_number').agg(lambda x: list(x))

        # determine the number of freeboard values per index
        # counts = [len(x) for x in temp["freeboard"]]

        # save counts as part of the GeoDataFrame
        # reduced_data["count"] = counts

        # filter out cells that do not have sufficent satellite returns
        # reduced_data = reduced_data.where(reduced_data["count"] > thres).dropna()

    if write is True:
        reduced_data.to_netcdf(output_path)
    else:
        print('No data files are being written. If you desire data files, set write=True')

    return reduced_data