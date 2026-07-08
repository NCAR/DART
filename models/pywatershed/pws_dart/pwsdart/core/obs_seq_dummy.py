import datetime
import pathlib
import re
import xarray
import datetime as dt
from .gregorian import gregorian

def obs_seq_dummy(
    hydro_rst: str,
    dir: pathlib.PosixPath    
):

    """
    Create observation sequence dummy files from hydro restart data.
    
    Args:
        hydro_rst: Path to the NetCDF restart file
        dir: Output directory for observation files
    """
    
    # Read the NetCDF file
    data = xarray.open_dataset(hydro_rst)
    
    # Get the time variable and convert to datetime
    time_var = data.time
    
    # Convert from days since 1979-01-01 to actual datetime
    # The time variable uses "days since 1979-01-01 00:00:00"
    time_datetime = time_var.values[0]  # Get the first (and only) time value
    
    # Convert numpy datetime64 to Python datetime
    if hasattr(time_datetime, 'astype'):
        time_datetime = time_datetime.astype('datetime64[s]').astype(datetime.datetime)
    
    # Extract date components
    year = time_datetime.year
    month = time_datetime.month
    day = time_datetime.day
    hour = time_datetime.hour
    minute = time_datetime.minute
    second = time_datetime.second
    
    # Create restart_time list in the format expected by gregorian function
    restart_time = [year, month, day, hour, minute, second]
    
    # Convert to gregorian
    gregorian_sec, gregorian_day = gregorian(restart_time)

# DO NOT EDIT THIS STRING
    obs_seq_dum = """ obs_sequence
obs_kind_definitions
           1
           1 STREAM_FLOW                    
  num_copies:            1  num_qc:            1
  num_obs:          1  max_num_obs:          1
observation                                                     
Data QC                                                         
  first:            1  last:          1
 OBS            1
  0.74756473302841187     
   1.0000000000000000     
          -1           -1          -1
obdef
loc3d
     4.949794212072877        0.7396891311859278         278.8099975585938      3
kind
           1
 gageID_linkID
 04233286           21983581
           1
     {0}     {1}
  1000000.0     
""".format(gregorian_sec, gregorian_day)

# Giagantic obs error variance ensures that this ob won't have an impact during 
# creation of the initial ensemble
# DO NOT EDIT THIS STRING
    obs_seq_identity_dum = """  obs_sequence
obs_kind_definitions
           0
  num_copies:            1  num_qc:            1
  num_obs:            1  max_num_obs:            1
observation                                                     
Data QC                                                         
  first:            1  last:            1
 OBS            1
  0.12345     
   1.0000000000000000     
          -1          -1          -1
obdef
loc3d
     4.949794212072877        0.7396891311859278         278.8099975585938      3
kind
        -130
  {0}     {1}
  1000000.0     
""".format(gregorian_sec, gregorian_day)

    
    # Create the hydro_file_list.txt and param_file_list.txt files
    obs_seq_dum_pathlib = dir / 'obs_seq.dum'
    obs_seq_out_pathlib = dir / 'obs_seq.out'
    with open(obs_seq_dum_pathlib, 'w') as opened_file:
        # Pending further upgrades which will help choose between these.
        #_ = opened_file.write(str(obs_seq_dum))
        _ = opened_file.write(str(obs_seq_identity_dum))


    obs_seq_out_pathlib.symlink_to(obs_seq_dum_pathlib)
