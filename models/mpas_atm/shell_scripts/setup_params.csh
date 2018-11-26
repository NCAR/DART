#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
#set echo
########################################################################################
# Set up parameters that are used for other scripts throughout the cycling period.
# Note: 1. Specify a full path name at each line.
#       2. Namelist options should be specified for all Namelist files separately. 
########################################################################################
# General configuration 
set EXPERIMENT_NAME = x4.30k  		# What do you call this experiment?
set DATE_INI = 2017-05-15_00:00:00      # initial cycle of the entire period
set DATE_BEG = 2017-06-02_00:00:00      # start date to run this script
set DATE_END = 2017-06-15_18:00:00      # end date to run this script
set INTV_DAY = 0                        # cycling frequency - assimilation_period_days    in input.nml
set INTV_SEC = 21600                    # cycling frequency - assimilation_period_seconds in input.nml

#  Settings specific to initializing the ensemble from grib files
set INIT_GRIB_FILE_LIST  = /glade/p/mmm/syha/FNL/UNGRIB/griblist  #  Name of variable with grib files (full path needed)
set TIME_INIT            = 00:10:00
set INIT_FORECAST_LENGTH = 240          #  Number of hours to integrate initial ensemble

# PBS setup
set   RUN_IN_PBS = yes          # Run on cheyenne? yes or no.
set  PROJ_NUMBER = NMMM0040	# Your account key for cheyenne
set FILTER_NODES = 18          # Total no. of nodes for DART/filter (at least bigger than ensemble size)
set  MODEL_NODES = 4            # Total no. of nodes for MPAS/atmosphere_model 
set N_PROCS_ANAL = 36           # Number of mpi processors or tasks for filter; Reduce this for a larger memory
set N_PROCS_MPAS = 36           # Number of mpi processors for model (=> MODEL_NODES * N_PROCS for graph.info)
set        QUEUE = economy	# queue for filter and model runs
set  TIME_FILTER = 00:40:00	# wall clock time for mpi filter runs
set    TIME_MPAS = 00:10:00	# wall clock time for mpi model runs

set    HPSS_SAVE = yes          # Backup in HPSS? yes or no. If yes, edit below.

# Ensemble configuration
set    MPAS_GRID = x4.133890  	# This global grid can provide LBCs for a limited-area MPAS.
set     TIMESTEP = 180.         # config_dt
set     LEN_DISP = 30000.       # Finest scale
set     ENS_SIZE = 100 		# Ensemble size
set ADAPTIVE_INF = true         # adaptive_inflation - If true, this script only supports
                                # spatially-varying state space prior inflation.
                                # And you also need to edit inf_sd_initial, inf_damping,
                                # inf_lower_bound, and inf_sd_lower_bound in &filter_nml.
set     INFL_OUT = output_priorinf
set      INFL_IN = input_priorinf

set   SST_UPDATE = true 
set   SST_FNAME  = x4.133890.sfc_update.spring2017.nc

# First Guess (for cold-start runs and initial ensemble)
set       VTABLE = Vtable.GFS

# Directories
set ROOT_DIR     = /glade2/scratch2/syha
set OUTPUT_DIR   = ${ROOT_DIR}/MPAS_DART/Cycle                          # output directory
set RUN_DIR      = ${ROOT_DIR}/MPAS_DART/Cycle           		# Run MPAS/DART cycling
set OBS_DIR      = ${ROOT_DIR}/MPAS_DART/OBS/OBS_SEQ			# obs_seq.out
set INIT_DIR     = ${ROOT_DIR}/MPAS_DART/INIT_ENS			# initial ensemble
set GRID_DIR     = /glade/p/work/syha/MPAS_INIT/GRID/x4.133890.gwd 	# Grid info
set SST_DIR      = /glade/p/mmm/syha/MPAS_DART/SST			# sfc_update.nc
set HPSS_DIR     = MPAS_DART/${EXPERIMENT_NAME}      			# hpss archives
set ENS_DIR      = member

set CSH_DIR      = /glade/p/mmm/syha/Manhattan/models/mpas_atm/shell_scripts
set DART_DIR     = /glade/p/mmm/syha/Manhattan/models/mpas_atm/work	# DART executables
set NML_DIR      = /glade/p/mmm/syha/Manhattan/models/mpas_atm/data	# DART executables
set MPAS_DIR     = /glade/p/mmm/syha/MPASV52_tend_vnTiedtke		# MPAS executables
set WPS_DIR      = /glade/p/mmm/syha/WPSV39                             #  WPS executables

# Namelist files
set NML_INIT     = namelist.init_atmosphere      # Namelist for init_atmosphere_model
set NML_MPAS     = namelist.atmosphere		 # Namelist for atmosphere_model
set NML_WPS      = namelist.wps			 # Namelist for WPS
set NML_DART     = input.nml			 # Namelist for DART
set STREAM_INIT  = streams.init_atmosphere	 # I/O list for init_atmosphere_model
set STREAM_ATM   = streams.atmosphere		 # I/O list for atmosphere_model
set STREAM_TEND  = streams.atmosphere.tend	 # I/O list for atmosphere_model
set STREAM_ENS1  = streams.atmosphere.ens1	 # I/O list for atmosphere_model

# Commands (do not need modification unless moving to new system)
set HSICMD = 'hsi put -P'
set REMOVE = '/bin/rm -rf'
set   COPY = 'cp -pf'
set   MOVE = 'mv -f'
set   LINK = 'ln -sf'
unalias cd
unalias ls

# Create LBCs for a limited-area version of MPAS => create_lbc.csh
set fcmd = /glade/p/mmm/syha/MPAS-Tools/limited_area/mpas_to_mpas/a.out.bdy_only.uncoupled.cheyenne
set lopt = "--use-reconstruct-winds"		# Wind option for LBCs
set rgrd = ${RUN_DIR}/CONUS.x20.835586.init.nc	# Regional IC
# Ex. ${fcmd} ${lopt} init.nc ${rgrd} global_lbc.2017-05-15_06.nc	# init.nc for MPAS_GRID

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

