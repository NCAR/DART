# ------------------------------------------------------------------------- #
# IMPORT LIBRARIES                                                          #
# ------------------------------------------------------------------------- #
import xarray as xr
import os
import glob
import h5py
import numpy as np

# ------------------------------------------------------------------------- #
# BEGIN FUNCTION DEFINITION                                                   #
# ------------------------------------------------------------------------- #


def determine_data_source(name, data_type):
    icesat2_names = ['ICESAT2', 'ICESAT-2', 'IS2', 'IS-2', 'Icesat2','Icesat-2','IceSat-2', 'icesat2','icesat-2','is2','is-2']

    if name in icesat2_names:
        if data_type in ['freeboard', 'fb']:
            print("You're converting freeboard data from NASA's ICESat-2 mission.")
        read_is2_data(data_type)
    else:
        raise NotImplementedError
        

def read_is2_data(data_type, start_date = '2019.01.01', tracks = ['gt1r', 'gt1l', 'gt2r', 'gt2l', 'gt3r', 'gt3l'], 
                  input_path = '/glade/scratch/mollyw/external_data/ICESat-2/',
                  output_path = '/glade/scratch/mollyw/external_data/IS2_reduced/freeboard/'):
    
    days_of_year = [31,29,31,30,31,30,31,31,30,31,30,31]
    year = start_date[0:4]
    dir_list = [name for name in sorted(os.listdir(input_path))]
    print(dir_list)
    
    if os.path.isdir(output_path) is False:
        print('Output directory presently does not exist! Generating directory now.')
        os.makedirs(output_path)
    
    if data_type in ['thickness', 'hi', 'sithick','SIT']:
        raise NotImplementedError('Not finished')
        
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
    filename = 'IS2_fb_'+dir_day+'.nc'
    day_data.to_netcdf(output_path+filename)
    print(dir_day + ' done!')