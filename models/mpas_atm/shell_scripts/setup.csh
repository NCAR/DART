#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
#set echo
########################################################################################
# Set up parameters that are used for other scripts throughout the cycling period.
# Note: 1. Specify a full path name at each line.
#       2. Namelist options should be specified for all Namelist files separately. 
########################################################################################
# General configuration 
set EXPERIMENT_NAME = MPASg

set DATE_INI = 2017-01-11_00:00:00      # initial cycle of the entire period
set DATE_BEG = 2017-01-11_18:00:00      # start date to run this script
set DATE_END = 2017-01-12_00:00:00      # end date to run this script
set INTV_DAY = 0                        # cycling frequency - assimilation_period_days    in input.nml
set INTV_SEC = 21600                    # cycling frequency - assimilation_period_seconds in input.nml

# PBS setup
set   RUN_IN_PBS = yes          # Run on cheyenne? yes or no.
set  PROJ_NUMBER = Nxxxxxxx	# Your account key for cheyenne
set FILTER_NODES = 16           # Total no. of nodes for DART/filter (at least bigger than ensemble size)
set  MODEL_NODES = 8            # Total no. of nodes for MPAS/atmosphere_model 
set N_PROCS_ANAL = 36           # Number of mpi processors or tasks for filter; Reduce this for a larger memory
set N_PROCS_MPAS = 36           # Number of mpi processors for model (=> MODEL_NODES * N_PROCS for graph.info)
set        QUEUE = regular	# queue for filter and model runs
set  TIME_FILTER = 00:50:00	# wall clock time for mpi filter runs
set    TIME_MPAS = 00:10:00	# wall clock time for mpi model runs

# filter configuration
set     TIMESTEP = 30.          # config_dt [seconds]
set     LEN_DISP = 5000.        # Finest scale
set       CUTOFF = 0.10         # half-width localization radius
set         VLOC = 60000.       # half-width localization radius in height [meters] - will be mulplied by CUTOFF
set     ENS_SIZE = 100		# Ensemble size
set num_output_obs_members = 2  # output members in obs_seq.out
set num_output_state_members = 2  # output members in output states    
set binary_obs_seq = false	# binary or ascii obs_seq.final to produce
set DISTRIB_MEAN = false	# true for a large-memory job; false otherwise
set CONVERT_OBS  = true		# convert_all_obs_verticals_first = .true. in &assim_tools_nml
set CONVERT_STAT = false	# convert_all_state_verticals_first = .true. in &assim_tools_nml
set ADAPTIVE_INF = true         # adaptive_inflation - If true, this script only supports
                                # spatially-varying state space prior inflation.
                                # And you also need to edit inf_sd_initial, inf_damping,
                                # inf_lower_bound, and inf_sd_lower_bound in &filter_nml.
set     INFL_OUT = output_priorinf
set      INFL_IN = input_priorinf

# SST
set   SST_UPDATE = false
set   SST_FNAME  = sfc_update.nc

# First Guess (for cold-start runs and initial ensemble)
set       VTABLE = Vtable.GFS

# init_template_filename in &model_nml in input.nml
set F_TEMPLATE   = /path_to_your_mpas/init.nc
set LBC_DIR      = /path_to_your_lbc_files                              # only for regional MPAS

# Directories
set ROOT_DIR     = /path_to_your_directory/MPAS_DART
set OUTPUT_DIR   = ${ROOT_DIR}/${EXPERIMENT_NAME}			# output directory
set RUN_DIR      = ${ROOT_DIR}/${EXPERIMENT_NAME} 			# Run MPAS/DART cycling
set OBS_DIR      = ${ROOT_DIR}/OBS_SEQ            			# obs_seq.out
set GRID_DIR     = /path_to_your_directory/GRID                  	# Grid info
set SST_DIR      = /path_to_your_directory/SST          			# sfc_update.nc
set ENS_DIR      = member

set DART_DIR     = /path_to_your_directory/Manhattan/models/mpas_atm/work	# DART executables
set NML_DIR      = /path_to_your_directory/Manhattan/models/mpas_atm/data 	# namelist files
set MPAS_DIR     = /path_to_your_mpas/MPASV7                         		# MPAS executables
set INIT_DIR     = ${ROOT_DIR}/INIT_ENS                            		# MPAS init.nc
set WPS_DIR      = /path_to_your_directory/WPSV3                                # WPS executables

# Namelist files - template files can be found in $NML_DIR
set NML_INIT     = namelist.init_atmosphere      # Namelist for init_atmosphere_model
set NML_MPAS     = namelist.atmosphere.global	 # or namelist.atmosphere.regional for regional MPAS
set NML_WPS      = namelist.wps			 # Namelist for WPS
set NML_DART     = input.nml			 # Namelist for DART
set STREAM_INIT  = streams.init_atmosphere	 # I/O list for init_atmosphere_model
set STREAM_ATM   = streams.atmosphere		 # I/O list for atmosphere_model

# Commands (do not need modification unless moving to new system)
set REMOVE = '/bin/rm -rf'
set   COPY = 'cp -pf'
set   MOVE = 'mv -f'
set   LINK = 'ln -sf'
