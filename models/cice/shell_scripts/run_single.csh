#!/bin/csh -f
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#----------------------------------------------------------------------
# This was written by Yongfei Zhang to explore running a single
# instance of the development version of CESM ... cesm1_5_beta06c
#----------------------------------------------------------------------

set myname = "yfzhang"
set casedir = "$WORK/cesmcases/"
set rundir  = "$SCRATCH/cesmruns/"

setenv MYCASENAME  "Dtest_single"
setenv RES         T62_g16
setenv MACH        yellowstone 
setenv COMPSET     DTEST


cd $casedir

if ( -d $casedir/$MYCASENAME) then
    echo "case existed"
    echo "to delete the old case, use \
          rm -rf $casedir/$MYCASENAME" ;exit
endif

set scriptdir = "/glade/scratch/bitz/darttest/cesm1_5_beta06c/cime/scripts/"

${scriptdir}/create_newcase  -res  $RES \
                             -mach $MACH \
                             -compset $COMPSET \
                             -case $MYCASENAME \
                             -project P93300065

cd $MYCASENAME

set stream_year_align = 2000
set stream_year_first = 2000
set stream_year_last  = 2000

./xmlchange -file env_build.xml   -id CESMSCRATCHROOT   -val "$rundir/"
./xmlchange -file env_run.xml     -id DIN_LOC_ROOT      -val "$SCRATCH/inputdata_cam"
./xmlchange -file env_run.xml     -id RUNDIR            -val "$rundir/$MYCASENAME/run"
./xmlchange -file env_run.xml     -id STOP_OPTION       -val ndays
./xmlchange -file env_run.xml     -id STOP_N            -val 1
./xmlchange -file env_run.xml     -id RESUBMIT          -val 1

#./xmlchange -file env_run.xml     -id DATM_MODE         -val CPLHIST3HrWx
#./xmlchange -file env_run.xml     -id DATM_CPLHIST_YR_START -val $stream_year_first
#./xmlchange -file env_run.xml     -id DATM_CPLHIST_YR_END   -val $stream_year_last
#./xmlchange -file env_run.xml     -id DATM_CPLHIST_YR_ALIGN -val $stream_year_align
./xmlchange -file env_run.xml     -id RUN_STARTDATE     -val ${stream_year_first}-01-01

./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val 64
./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val 64
./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val 64
./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val 64
./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val 64
./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val 64
./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val 64
./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val 64


./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val 1
./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val 1
./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val 1
./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val 1
./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val 1
./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val 1
./xmlchange -file env_mach_pes.xml -id NTHRDS_ROF -val 1
./xmlchange -file env_mach_pes.xml -id NTHRDS_WAV -val 1

./xmlchange -file env_mach_pes.xml -id TOTALPES            -val 64
./xmlchange -file env_mach_pes.xml -id MAX_TASKS_PER_NODE  -val 16

./case.setup
echo "case setup finished"

./case.build

exit 0


