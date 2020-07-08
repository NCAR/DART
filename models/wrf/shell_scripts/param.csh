#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# TJH   ADAPTIVE_INFLATION is disconnected from input.nml
# TJH   ASSIM_INT_HOURS  is implicit in (ALL) the scripts except assim_advance.csh
#                        ASSIM_INT_MINUTES support needs to be added to param.csh,
#                        it is referenced in assim_advance.csh but not declared in param.csh

# Set up environment. Current settings are for NCAR's Cheyenne
module load mpt          # set this appropriately #%%%#
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
set BASE_DIR         = /glade2/scratch2/USER/WORK_DIR     # set this appropriately #%%%#
set RUN_DIR          = ${BASE_DIR}/rundir
set TEMPLATE_DIR     = ${BASE_DIR}/template
set OBSPROC_DIR      = ${BASE_DIR}/obsproc
set OUTPUT_DIR       = ${BASE_DIR}/output
set ICBC_DIR         = ${BASE_DIR}/icbc
set POST_STAGE_DIR   = ${BASE_DIR}/post
set OBS_DIAG_DIR     = ${BASE_DIR}/obs_diag
set PERTS_DIR        = ${BASE_DIR}/perts

#  Directories that can be used by many things
set SHELL_SCRIPTS_DIR = ${BASE_DIR}/scripts
set DART_DIR          = /glade/p/work/USER/DART_manhattan           # set this appropriately #%%%#
set WRF_DM_SRC_DIR    = /glade/p/work/USER/WRFV3_dmpar              # set this appropriately #%%%#
set WPS_SRC_DIR       = /glade/p/work/USER/WPS                      # set this appropriately #%%%#
set VAR_SRC_DIR       = /glade/p/work/USER/WRFDA                    # set this appropriately #%%%#

# for generating wrf template files
set GEO_FILES_DIR     = /glade/p/work/USER/WPS        # set this appropriately #%%%#
set GRIB_DATA_DIR     = /glade/p/work/USER/WPS/GRIB   # set this appropriately #%%%#
set GRIB_SRC          = 'GFS'                         # set this appropriately #%%%#

# list of variables for extraction and cycling
set extract_vars_a   = ( U V PH T MU QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP QNICE QNRAIN \
                         U10 V10 T2 Q2 PSFC TSLB SMOIS TSK RAINC RAINNC GRAUPELNC )
set extract_vars_b   = ( U V W PH T MU QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP QNICE QNRAIN \
                         U10 V10 T2 Q2 PSFC TSLB SMOIS TSK RAINC RAINNC GRAUPELNC \
                         REFL_10CM VT_DBZ_WT )
set cycle_vars_a     =   ( U V PH T MU QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP QNICE QNRAIN \
                           U10 V10 T2 Q2 PSFC TSLB SMOIS TSK )
set increment_vars_a = ( U V PH T MU QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP QNICE QNRAIN U10 V10 T2 Q2 PSFC )

#  Diagnostic parameters
set OBS_VERIF_DAYS      = 7

#  Generic queuing system parameters
set SUPER_PLATFORM      = cheyenne

# TJH consistent way of checking the SUPER_PLATFORM and injecting that
#     header information into the scripts ... rather than have scripts
#     that have redundant blocks in them ...
#
set COMPUTER_CHARGE_ACCOUNT = YOUR_ACCT                  # set this appropriately #%%%#
set EMAIL                   = YOUR_EMAIL@SOMEPLACE.COM   # set this appropriately #%%%#

if ( $SUPER_PLATFORM == 'cheyenne') then
   # cheyenne values (uses 'PBS' queueing system) 
   # set this appropriately #%%%#  ... ALL OF THESE if using PBS
   set FILTER_QUEUE       = regular
   set FILTER_TIME        = 0:35:00
   set FILTER_NODES       = 10
   set FILTER_PROCS       = 36
   set FILTER_MPI         = 36
   set ADVANCE_QUEUE       = regular
   set ADVANCE_TIME        = 0:20:00
   set ADVANCE_NODES       = 3
   set ADVANCE_PROCS       = 36
   set ADVANCE_MPI         = 36
else
   # yellowstone (uses 'LSF' queueing system)
   # set this appropriately #%%%#  ... ALL OF THESE if using LSF
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
# TJH ... The LINK command probably should not have the force option.
# TJH ... and if the LINK fails, should it die right there?
setenv   REMOVE 'rm -rf'
setenv   COPY 'cp -pfr'
setenv   MOVE 'mv -f'
setenv   LINK 'ln -fs'
setenv   WGET /usr/bin/wget
setenv   LIST 'ls'

exit 0

