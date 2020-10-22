import argparse
import datetime
import pathlib
import xarray as xa

# python3 create_parameter_restart_file.py \
#    /Users/jamesmcc/Downloads/ \
#    --hydro_rst /Users/jamesmcc/Downloads/croton_NY/NWM/RESTART/HYDRO_RST.2011-08-26_00_00_DOMAIN1 \
#    --restart /Users/jamesmcc/Downloads/croton_NY/NWM/RESTART/RESTART.2011082600_DOMAIN1 \
#    --existing_variables qlink1 qlink1 SNEQV \
#    --new_variables mult_qBucket mult_qSfcLatRunoff some_dum_snow \
#    --values         1      5   15
#
# ncdump ~/Downloads/test_param_restart_file.nc 
#
# netcdf test_param_restart_file {
# dimensions:
# 	links = 249 ;
# 	Time = 1 ;
# 	south_north = 16 ;
# 	west_east = 15 ;
# variables:
# 	float mult_qBucket(links) ;
# 		mult_qBucket:_FillValue = NaNf ;
# 	float mult_qSfcLatRunoff(links) ;
# 		mult_qSfcLatRunoff:_FillValue = NaNf ;
# 	float some_dum_snow(Time, south_north, west_east) ;
# 		some_dum_snow:_FillValue = NaNf ;

# // global attributes:
# 		:his_out_counts = 9360 ;
# 		:Restart_Time = "2011-08-26_00:00:00" ;
# 		:Since_Date = "2010-08-01_00:00:00" ;
# 		:DTCT = 300.f ;
# 		:channel_only = 0 ;
# 		:channelBucket_only = 0 ;
# data:

#  mult_qBucket = 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
#     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
#     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
#     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
#     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
#     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
#     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
#     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
#     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
#     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
#     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ;

#  mult_qSfcLatRunoff = 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
#     5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
#     5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
#     5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
#     5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
#     5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
#     5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
#     5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
#     5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
#     5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
#     5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 ;

#  some_dum_snow =
#   15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
#   15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
#   15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
#   15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
#   15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
#   15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
#   15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
#   15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
#   15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
#   15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
#   15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
#   15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
#   15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
#   15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
#   15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
#   15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 ;
# }

# TODO(JLM): make either of the restart files optional, right now both are required even if not used.

def create_parameter_restart_file(
    out_path: str=None,
    out_mode: str='r',
    hydro_rst_file: str=None,
    restart_file: str=None,
    existing_variables: list=None,
    new_variables: list=None,
    values: list=None
):

    hydro_rst_time = None
    restart_time = None
    
    if hydro_rst_file is not None:
        hydro_rst = xa.open_dataset(hydro_rst_file)
        hydro_rst_time = datetime.datetime.strptime(hydro_rst.Restart_Time, '%Y-%m-%d_%H:%M:%S')

    if restart_file is not None:
        restart = xa.open_dataset(restart_file)
        restart_time_raw = restart.data_vars.get('Times').data[0].decode('utf-8')
        restart_time = datetime.datetime.strptime(restart_time_raw, '%Y-%m-%d_%H:%M:%S')

    if hydro_rst_time is not None and restart_time is not None:
        if hydro_rst_time != restart_time:
            raise ValueError('Times of the supplied restart files do not match.')

    the_date = hydro_rst_time if hydro_rst_time is not None else restart_time
    the_date_str = the_date.strftime('%Y-%m-%d_%H:%M')  # drop the seconds
    out_file = pathlib.Path(out_path) / ('parameter_restart.' + the_date_str + '.nc')

    if out_file.exists() and out_mode != 'w':
        raise FileExistsError(
            'The out_mode=' + out_mode +
            '  argument precludes overwriting/clobbering the existing file: ' +
            str(out_file)
        )
    
    new_dataset = None

    for new, exst, val in zip(new_variables, existing_variables, values):

        src = None
        if hydro_rst_file is not None and exst in hydro_rst.keys():
            src = hydro_rst
        if restart_file is not None and exst in restart.keys():
            src = restart

        if src is None:
            continue

        tmp_dataset = src[[exst]]
        tmp_dataset = tmp_dataset.rename({exst:new})
        tmp_dataset[new] = tmp_dataset[new]*float(0.0) + float(val)
        tmp_dataset[new].encoding.update({'_FillValue': -9999})

        if new_dataset is None:
            new_dataset = tmp_dataset
        else:
            new_dataset = new_dataset.merge(tmp_dataset)

        del tmp_dataset

    new_dataset.to_netcdf(out_file)

    return(out_file)
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create a WRF-Hydro parameter assimilation file.')

    parser.add_argument(
        'out_path',
        metavar='/path/to/',
        type=str,
        nargs=1,
        help='A path to write the parameter file to.'
    )

    parser.add_argument(
        '--hydro_rst',
        metavar='/path/to/HYDRO_RST.nc',
        type=str,
        nargs=1,
        help='A Hydro Restart file for hydro dimensions.'
    )

    parser.add_argument(
        '--restart',
        metavar='/path/to/RESTART.nc',
        type=str,
        nargs=1,
        help='A NoahMP Restart file for NoahMP dimensions.'
    )

    parser.add_argument(
        '--new_variables',
        nargs='*',
        help='A list of variables collated with values list.'
    )

    parser.add_argument(
        '--existing_variables',
        nargs='*',
        help='A list of variables collated with values list.'
    )

    parser.add_argument(
        '--values',
        nargs='*',
        help='A list of variables collated with variables list.'
    )

    args = parser.parse_args()
    out_path = pathlib.PosixPath(args.out_path[0])

    hydro_rst_file = args.hydro_rst
    if hydro_rst_file is not None:
        hydro_rst_file = pathlib.PosixPath(hydro_rst_file[0])

    restart_file = args.restart
    if restart_file is not None:
        restart_file = pathlib.PosixPath(restart_file[0])

    existing_variables = args.existing_variables
    new_variables = args.new_variables
    values = args.values

    # out_path = '/Users/jamesmcc/Downloads/'
    # hydro_rst_file = pathlib.PosixPath('/Users/jamesmcc/Downloads/croton_NY/NWM/RESTART/HYDRO_RST.2011-08-26_00_00_DOMAIN1')
    # restart_file = '/Users/jamesmcc/Downloads/croton_NY/NWM/RESTART/RESTART.2011082600_DOMAIN1'
    # new_variables = ['mult_qSfcLatRunoff', 'mult_qBucket']
    # existing_variables = ['qlink1', 'qlink2']
    # values = [1, 5]

    create_parameter_restart_file(
        out_path=out_path,
        hydro_rst_file=hydro_rst_file,
        restart_file=restart_file,
        existing_variables=existing_variables,
        new_variables=new_variables,
        values=values
    )
