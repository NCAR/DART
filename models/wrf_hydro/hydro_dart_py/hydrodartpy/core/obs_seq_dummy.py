import datetime
import pathlib
import re
import xarray

from .gregorian import gregorian

def obs_seq_dummy(
    hydro_rst: str,
    dir: pathlib.PosixPath    
):

    data = xarray.open_dataset(hydro_rst)
    restart_time = data.attrs['Restart_Time']
    restart_time = re.split('[-_:]',restart_time)
    restart_time = [int(ii) for ii in restart_time]
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
