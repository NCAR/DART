#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# This file is just to see if CESM can be bullt on whatever architecture.
# Single instance, no DART.
# You are expected to read this and understand it before executing it.

# You will have to modify this to point to your CESM repository
set CESMCLONE = /glade/work/thoar/CESM/my_cesm_sandbox

# You will have to modify this to suit your own purpose.
set CASE = clm5bgc-crop_simple
set CASEDIR = /glade/work/thoar/cases
set COMPSET = 2000_DATM%GSWP3v1_CLM50%BGC-CROP_SICE_SOCN_MOSART_SGLC_SWAV
set CASEROOT = ${CASEDIR}/${CASE}

# Just echo what compsets are available
${CESMCLONE}/cime/scripts/query_config  --compsets clm

# Prevent accidental removal of an experiment.
if ( -e ${CASEROOT} ) then
   echo "WARNING" 
   echo "WARNING : The CASE ${CASEROOT} exists. Stopping."
   echo "WARNING : If you want to perform the experiment,"
   echo "WARNING : you have to manually remove the CASEROOT, EXEDIR, RUNDIR"
   echo "WARNING : or change the CASE name to something else."
   exit
endif


echo " Starting create_newcase ..."

${CESMCLONE}/cime/scripts/create_newcase \
    --case ${CASEROOT} \
    --compset ${COMPSET} \
    --mach cheyenne \
    --res f09_f09_mg17 \
    --project P86850054 \
    --run-unsupported || exit 1

cd ${CASEROOT}

echo " Finished create_newcase ..."
echo " Starting case.setup     ..."

foreach FILE ( *xml )
   cp -v ${FILE} ${FILE}.original
end

./xmlchange RUN_TYPE=startup
./xmlchange RUN_STARTDATE=2000-01-01
./xmlchange DOUT_S=FALSE
./xmlchange CALENDAR=GREGORIAN

@ nthreads = 1
set nodes_per_instance = 2

@ atm_tasks = -1 * ${nodes_per_instance}
@ cpl_tasks = -1 * ${nodes_per_instance}
@ ocn_tasks = -1 * ${nodes_per_instance}
@ wav_tasks = -1 * ${nodes_per_instance}
@ glc_tasks = -1 * ${nodes_per_instance}
@ ice_tasks = -1 * ${nodes_per_instance}
@ rof_tasks = -1 * ${nodes_per_instance}
@ lnd_tasks = -1 * ${nodes_per_instance}
@ esp_tasks = -1 * ${nodes_per_instance}

./xmlchange ROOTPE_ATM=0,NTHRDS_ATM=$nthreads,NTASKS_ATM=$atm_tasks
./xmlchange ROOTPE_CPL=0,NTHRDS_CPL=$nthreads,NTASKS_CPL=$cpl_tasks
./xmlchange ROOTPE_OCN=0,NTHRDS_OCN=$nthreads,NTASKS_OCN=$ocn_tasks
./xmlchange ROOTPE_WAV=0,NTHRDS_WAV=$nthreads,NTASKS_WAV=$wav_tasks
./xmlchange ROOTPE_GLC=0,NTHRDS_GLC=$nthreads,NTASKS_GLC=$glc_tasks
./xmlchange ROOTPE_ICE=0,NTHRDS_ICE=$nthreads,NTASKS_ICE=$ice_tasks
./xmlchange ROOTPE_ROF=0,NTHRDS_ROF=$nthreads,NTASKS_ROF=$rof_tasks
./xmlchange ROOTPE_LND=0,NTHRDS_LND=$nthreads,NTASKS_LND=$lnd_tasks
./xmlchange ROOTPE_ESP=0,NTHRDS_ESP=$nthreads,NTASKS_ESP=$esp_tasks

# At one point these had no effect when run before case.setup
# With release-clm5.0.04,  these are best _before_ case.setup
./xmlchange --subgroup case.run --id JOB_QUEUE          --val premium
./xmlchange --subgroup case.run --id JOB_WALLCLOCK_TIME --val 0:20

./case.setup || exit 2

# create a couple extra history files 

set fname = "user_nl_clm"
setenv stop_option         nhours
setenv stop_n              24
@ clm_dtime = 1800
@ h1nsteps = $stop_n * 3600 / $clm_dtime
echo "hist_empty_htapes = .false."                                     >! ${fname}
echo "hist_fincl1 = 'NEP','TOTECOSYSC','TOTVEGC','TLAI'"               >> ${fname}
echo "hist_fincl2 = 'NEP','FSH','EFLX_LH_TOT_R'"                       >> ${fname}
echo "hist_fincl3 = 'TV','TLAI','PBOT','TBOT','TSA','RH2M_R','SNOWDP'" >> ${fname}
echo "hist_nhtfrq = -$stop_n,1,-$stop_n"                               >> ${fname}
echo "hist_mfilt  = 1,$h1nsteps,1"                                     >> ${fname}
echo "hist_avgflag_pertape = 'A','A','I'"                              >> ${fname}
echo "hist_dov2xy = .true.,.true.,.false."                             >> ${fname}
echo "hist_type1d_pertape = ' ',' ',' '"                               >> ${fname}

echo " Finished case.setup     ..."
echo " Starting case.build     ..."

./case.build || exit 3

echo " Finished case.build     ..."
echo ""
echo "  cd  ${CASEROOT}"
echo "  ./case.submit"
echo ""
echo "  if that works, try a do_nothing assimilation cycle."

# These modify env_run.xml so they can be executed AFTER a successful clm cycle.
if ( 1 == 2 ) then
   ./xmlchange DATA_ASSIMILATION_LND=TRUE
   ./xmlchange DATA_ASSIMILATION_CYCLES=1
   ./xmlchange DATA_ASSIMILATION_SCRIPT=do_nothing.csh
   echo "echo Hello from `hostname`" >! do_nothing.csh
   chmod 755 do_nothing.csh
endif

exit 0

