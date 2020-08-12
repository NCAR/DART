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
set EXPERIMENT_NAME = TEST	        # What do you call this experiment?
set DATE_INI = 2016-08-19_00:00:00      # initial cycle of the entire period
set DATE_BEG = 2016-08-19_00:00:00      # start date to run this script
set DATE_END = 2016-08-19_06:00:00      # end date to run this script
set INTV_DAY = 0                        # cycling frequency - assimilation_period_days    in input.nml
set INTV_SEC = 21600                    # cycling frequency - assimilation_period_seconds in input.nml

#  Settings specific to initializing the ensemble from grib files
set INIT_GRIB_FILE_LIST  = /glade/p/work/USER/MPAS-DART/rundir/griblist  #  Name of variable with grib files (full path needed)
set TIME_INIT            = 00:50:00
set INIT_FORECAST_LENGTH = 12           #  Number of hours to integrate initial ensemble

# PBS setup
set   RUN_IN_PBS = yes          # Run on cheyenne? yes or no.
set  PROJ_NUMBER = xxxxxxxx	# Your account key for cheyenne
set FILTER_NODES = 35           # Total no. of nodes for DART/filter (at least bigger than ensemble size)
set  MODEL_NODES = 4            # Total no. of nodes for MPAS/atmosphere_model 
set       N_CPUS = 36		# Number of cpus per node (default = 36)
set      N_PROCS = 32		# Number of mpi processors (=> MODEL_NODES * N_PROCS for graph.info)
set        QUEUE = economy	# queue for filter and model runs
set  TIME_FILTER = 00:30:00	# wall clock time for mpi filter runs
set    TIME_MPAS = 00:10:00	# wall clock time for mpi model runs

set    HPSS_SAVE = no           # Backup in HPSS? yes or no. If yes, edit below.

# Ensemble configuration
set    MPAS_GRID = x1.40962   	# All grid parameters will be changed based on this.
set     ENS_SIZE = 30 	 	# Ensemble size
set ADAPTIVE_INF = true         # adaptive_inflation - If true, this script only supports
                                # spatially-varying state space prior inflation.
                                # And you also need to edit inf_sd_initial, inf_damping,
                                # inf_lower_bound, and inf_sd_lower_bound in &filter_nml.
set     INFL_OUT = output_priorinf
set      INFL_IN = input_priorinf

set   SST_UPDATE = false
set   SST_FNAME  = ${MPAS_GRID}.sfc_update.nc

# First Guess (for cold-start runs and initial ensemble)
set       VTABLE = Vtable.GFS

# Directories
set ENS_DIR      = member
set ROOT_DIR     = /glade/scratch/USER
set OUTPUT_DIR   = ${ROOT_DIR}/MPAS-DART/output                         # output directory
set CSH_DIR      = ${ROOT_DIR}/MPAS-DART/TEST/Scripts	                # shell scripts
set RUN_DIR      = ${ROOT_DIR}/MPAS-DART/rundir          		# Run MPAS/DART cycling
set SST_DIR      = ${ROOT_DIR}/MPAS-DART/SST				# sfc_update.nc
set OBS_DIR      = ${ROOT_DIR}/MPAS-DART/obs            		# obs_seq.out
set INIT_DIR     = ${ROOT_DIR}/MPAS-DART/INIT				# initial ensemble
set FG_DIR       = ${ROOT_DIR}/FG  					# First Guess
set GRID_DIR     = ${ROOT_DIR}						# Grid info
set HPSS_DIR     = MPAS_DART/${EXPERIMENT_NAME}      			# hpss archives

set MPAS_DIR     = /glade/u/home/USER/MPAS_v5.0  				# MPAS executables
set DART_DIR     = /glade/u/home/USER/DART-Man/models/mpas_atm/work		# DART executables
set WPS_DIR      = /glade/u/home/wrfhelp/PRE_COMPILED/WPSV3.9_intel_serial  	#  WPS executables

# Namelist files
set NML_INIT     = namelist.init_atmosphere      # Namelist for init_atmosphere_model
set NML_MPAS     = namelist.atmosphere		 # Namelist for atmosphere_model
set NML_WPS      = namelist.wps			 # Namelist for WPS
set NML_DART     = input.nml			 # Namelist for DART
set STREAM_ATM   = streams.atmosphere		 # I/O list for atmosphere_model
set STREAM_INIT  = streams.init_atmosphere	 # I/O list for init_atmosphere_model

# Commands (do not need modification unless moving to new system)
set HSICMD = 'hsi put -P'
set REMOVE = '/bin/rm -rf'
set   COPY = 'cp -pf'
set   MOVE = 'mv -f'
set   LINK = 'ln -sf'
unalias cd
unalias ls


