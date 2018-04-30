#!/bin/csh

   #  Set the assimilation parameters
   set NUM_ENS            = 50
   set ASSIM_INT_HOURS    = 6   # ignored if assim_int_minutes > 0
   set IC_PERT_SCALE      = 0.25
   set ADAPTIVE_INFLATION = 1   # set to 1 if using adaptive inflaton to tell the scripts to look for the files
   set NUM_DOMAINS        = 1

   #  Directories where things are run
   set BASE_DIR         = /glade2/scratch2/USER/WORK_DIR     # set this appropriately #%%%#
   set RUN_DIR          = ${BASE_DIR}/rundir
   set TEMPLATE_DIR     = ${BASE_DIR}/template
   set OBSPROC_DIR      = ${BASE_DIR}/obsproc
   set OUTPUT_DIR       = ${BASE_DIR}/output
   set ICBC_DIR         = ${BASE_DIR}/icbc
   set POST_STAGE_DIR   = ${BASE_DIR}/post
   set OBS_DIAG_DIR     = ${BASE_DIR}/obs_diag

   #  Directories that can be used by many things
   set SHELL_SCRIPTS_DIR = ${BASE_DIR}/scripts
   set DART_DIR          = /glade/p/work/USER/DART_manhattan           # set this appropriately #%%%#
   set WRF_SRC_DIR       = /glade/p/work/USER/WRFV3                    # set this appropriately #%%%#
   set WPS_SRC_DIR       = /glade/p/work/USER/WPS                      # set this appropriately #%%%#
   set VAR_SRC_DIR       = /glade/p/work/USER/WRFDA                    # set this appropriately #%%%#

   # for generating wrf template files
   set GEO_FILES_DIR     =  /glade/p/work/USER/WPS       # set this appropriately #%%%#
   set GRIB_DATA_DIR     = /glade/p/work/USER/WPS/GRIB   # set this appropriately #%%%#

   # list of variables for extraction and cycling
   ################################################################
   # IMPORTANT - SHOULD MATCH those set in new_advance_model.csh  #
   ################################################################
   set extract_vars_a = ( U V PH T MU QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP QNICE QNRAIN \
                       U10 V10 T2 Q2 PSFC TSLB SMOIS TSK RAINC RAINNC GRAUPELNC )
   set cycle_vars_a =   ( U V PH T MU QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP QNICE QNRAIN \
                     U10 V10 T2 Q2 PSFC TSLB SMOIS TSK )
   set increment_vars_a = ( U V PH T MU QVAPOR QCLOUD QRAIN QICE QSNOW QGRAUP QNICE QNRAIN H_DIABATIC U10 V10 T2 Q2 PSFC )

   #  Diagnostic parameters
   set OBS_VERIF_DAYS      = 7

   #  Generic queuing system parameters
#   set SUPER_PLATFORM      = yellowstone
   set SUPER_PLATFORM      = cheyenne

   #  CHARGE ACCOUNTS
   set NCAR_GAU_ACCOUNT    = YOUR_ACCT   # set this appropriately #%%%#
   set CNCAR_GAU_ACCOUNT   = YOUR_ACCT   # set this appropriately #%%%#

# yellowstone parameters
   set FILTER_QUEUE        = regular
   set FILTER_TIME         = 0:25
   set FILTER_CORES        = 512
   set ADVANCE_QUEUE       = regular
   set ADVANCE_TIME        = 0:18
   set ADVANCE_CORES       = 64
   set NCAR_FILTER_PTILE   = 16         
   set NCAR_ADVANCE_PTILE  = 16
# cheyenne parameters
   set CFILTER_QUEUE       = regular
   set CFILTER_TIME        = 0:35:00
   set CEMAIL              = YOUR_EMAIL@SOMEPLACE.COM      # set this appropriately #%%%#
   set CFILTER_NODES       = 30
   set CFILTER_PROCS       = 30
   set CFILTER_MPI         = 30
   set CADVANCE_QUEUE       = regular
   set CADVANCE_TIME        = 00:16:00
   set CDADVANCE_TIME       = 00:35:00
   set CADVANCE_NODES       = 3
   set CADVANCE_PROCS       = 36
   set CADVANCE_MPI         = 36

   #  System specific commands
   setenv   REMOVE 'rm -rf'
   setenv   COPY 'cp -pfr'
   setenv   MOVE 'mv -f'
   setenv   LINK 'ln -fs'
   setenv   WGET /usr/bin/wget
