#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# TJH   ADAPTIVE_INFLATION is disconnected from input.nml
# TJH   ASSIM_INT_HOURS  is implicit in (ALL) the scripts except assim_advance.csh
#                        ASSIM_INT_MINUTES support needs to be added to param.csh,
#                        it is referenced in assim_advance.csh but not declared in param.csh

# Set up environment. Current settings are for NCAR's Derecho
module load nco          # set this appropriately #%%%#
module load ncl/6.6.2    # set this appropriately #%%%#

#  Set the assimilation parameters
set NUM_ENS            = 50
set ASSIM_INT_MINUTES  = 0   # 0 means use ASSIM_INT_HOURS
set ASSIM_INT_HOURS    = 6   # ignored if ASSIM_INT_MINUTES > 0
set IC_PERT_SCALE      = 0.25
set ADAPTIVE_INFLATION = 1   # set to 1 if using adaptive inflation to tell the scripts to look for the files
set NUM_DOMAINS        = 1

#  Directories where things are run
#  IMPORTANT : Scripts provided rely on this directory structure and names relative to BASE_DIR.
#              Do not change, otherwise tutorial will fail.    
set BASE_DIR         = /glade/derecho/scratch/USER/WORK_DIR     # set this appropriately #%%%#
set RUN_DIR          = ${BASE_DIR}/rundir
set TEMPLATE_DIR     = ${BASE_DIR}/template
set OBSPROC_DIR      = ${BASE_DIR}/obsproc
set OUTPUT_DIR       = ${BASE_DIR}/output
set ICBC_DIR         = ${BASE_DIR}/icbc
set POST_STAGE_DIR   = ${BASE_DIR}/post
set OBS_DIAG_DIR     = ${BASE_DIR}/obs_diag
set PERTS_DIR        = ${BASE_DIR}/perts

#  Assign path to DART, WRF, WPS and WRFDA build
set SHELL_SCRIPTS_DIR = ${BASE_DIR}/scripts
set DART_DIR          = /glade/work/USER/DART                     # set this appropriately #%%%#
set WRF_DM_SRC_DIR    = /glade/work/USER/WRFV3                    # set this appropriately #%%%#
set WPS_SRC_DIR       = /glade/work/USER/WPS                      # set this appropriately #%%%#
set VAR_SRC_DIR       = /glade/work/USER/WRFDA                    # set this appropriately #%%%#

# for generating wrf template files
set GEO_FILES_DIR     = /glade/u/home/wrfhelp/WPS_GEOG            # set this appropriately #%%%#
set GRIB_DATA_DIR     = ${ICBC_DIR}/grib_data                     # set this appropriately #%%%#
set GRIB_SRC          = 'GFS'                                     # set this appropriately #%%%#

# list of variables for extraction and cycling
set extract_vars_a   = ( U V PH THM MU QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP QNICE QNRAIN \
                         U10 V10 T2 Q2 PSFC TSLB SMOIS TSK RAINC RAINNC GRAUPELNC )
set extract_vars_b   = ( U V W PH THM MU QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP QNICE QNRAIN \
                         U10 V10 T2 Q2 PSFC TSLB SMOIS TSK RAINC RAINNC GRAUPELNC \
                         REFL_10CM VT_DBZ_WT )
set cycle_vars_a     =   ( U V PH THM MU QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP QNICE QNRAIN \
                           U10 V10 T2 Q2 PSFC TSLB SMOIS TSK )
set increment_vars_a = ( U V PH THM MU QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP QNICE QNRAIN U10 V10 T2 Q2 PSFC )

#  Diagnostic parameters
set OBS_VERIF_DAYS      = 7

#  Generic queuing system parameters
set SUPER_PLATFORM          = derecho
set COMPUTER_CHARGE_ACCOUNT = YOUR_ACCT                  # set this appropriately #%%%#
set EMAIL                   = YOUR_EMAIL                 # set this appropriately #%%%#

if ( $SUPER_PLATFORM == 'derecho') then
   # Derecho values (uses 'PBS' queueing system) 
   # Set these appropriately for your PBS system  #%%%#  
   set FILTER_QUEUE       = main
   set FILTER_PRIORITY    = premium
   set FILTER_TIME        = 0:35:00
   set FILTER_NODES       = 2
   set FILTER_PROCS       = 128
   set FILTER_MPI         = 128

   set ADVANCE_QUEUE      = main
   set ADVANCE_PRIORITY   = premium
   set ADVANCE_TIME       = 0:20:00
   set ADVANCE_NODES      = 1
   set ADVANCE_PROCS      = 128
   set ADVANCE_MPI        = 128
else
   # 'LSF' queueing system example
   # Set these appropriately for your LSF or Slurm system #%%%# 
   set FILTER_QUEUE        = regular
   set FILTER_TIME         = 0:25
   set FILTER_CORES        = 512
   set ADVANCE_QUEUE       = regular
   set ADVANCE_TIME        = 0:18
   set ADVANCE_CORES       = 64
   set FILTER_PTILE        = 16
   set ADVANCE_PTILE       = 16
endif

#  System specific commands
setenv   REMOVE 'rm -rf'
setenv   COPY 'cp -pfr'
setenv   MOVE 'mv -f'
setenv   LINK 'ln -fs'
setenv   WGET /usr/bin/wget
setenv   LIST 'ls'

exit 0

