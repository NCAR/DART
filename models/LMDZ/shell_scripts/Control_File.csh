#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#****************************************************************************
# Set control variables for assimilation
# It will be used by job.csh and run_lmdz.csh 
# Tarkeshwar Singh Sep 2014
# ****************************************************************************
set num_proc       = 40
set host_file      = my-hosts2
set mpi_path       = /data/opt/mpi/openmpi-1.6.3
set Host_File_Path = /home/tk

# List of obs_seq.out files to be assimilated. 
set obs_seq_list   = olist

# store restart file at every $restart_store_freq day
set restart_store_freq = 2

# Define initial inflation parameters
set inf_initial    = 1.0
set inf_sd_initial = 0.6

# Logical parameter to control assimilation with or without sampling corrections
# 1 = with sampling corrections
# any others values means no sampling correction
set sampling_error_correction = 1

# Define the LMDZ exe and limit file name
set gcm_exe    = gcm_360x180x19_phylmd_seq.e
set ce0l_exe   = ce0l_360x180x19_phylmd_seq.e
set limit_file = limit.nc_360x180x19

# Path of model control files, data and excutables.
set DART_LMDZ5        = /home/tarkesh/DART/kodiak/models/LMDZ5
set DART_ics          = /home/tk/WORK/DART/monsoon_360x180x19/ENSEMBLES/clim_fnl_22May/ 
set LMDZ_DEF_PATH     = /home/tk/WORK/DART/monsoon_360x180x19/LMDZ_init
set OBS_PATH          = /home/tarkesh/DART/lanai/observations/NCEP/ascii_to_obs/obs_seq2010/ 
set DART_restarts_ics = /home/tk/WORK/DART/monsoon_360x180x19/FILTER_RUN/NCEP/EXP1/OUTPUT_20100928

exit 0


