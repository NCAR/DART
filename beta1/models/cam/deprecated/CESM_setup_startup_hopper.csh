#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# ---------------------
# Purpose
# ---------------------
#
# This script is designed to configure and build a multi-instance CESM model
# that has CAM, CLM, and CICE as active components over a single data ocean,
# and will use DART to assimilate observations at regular intervals.
# This script does not build DART.

# ---------------------
# How to set up the script
# ---------------------
# -- Either edit and run this script in the $DART/models/cam/shell_scripts
#    directory where it first appears,
#    or copy it to somewhere that it will be preserved and run it there.
#    It will create a 'case' directory, where the model will be built,
#    and an execution directory, where each forecast and assimilation will
#    take place.  The short term archiver will use a third directory for
#    storage of model output until it can be moved to long term storage (HPSS)
# -- Examin the whole script to identify things to change for your experiments.
# -- Provide any initial files needed by your run:
#       inflation
#       sampling error correction
#       CAM/CLM/CICE initial ensemble
#       ...
# -- Run this script.
# -- Edit the DART input.nml that appears in the $caseroot directory.
# -- Submit the job using $caseroot/${case}.${mach}.submit
#      ($mach may not be needed for cesm releases after cesm1_1_beta04)
#
# ---------------------
# Important features
# ---------------------
#
# If you want to change something in your case other than the runtime
# settings, you need to delete everything and start the run from scratch.
#
# ./${CASENAME}.*.clean_build
# ./configure -cleanall
#
# ====================================================================
# === IMPORTANT modifications to the distribution code before ANYTHING
# ====================================================================

# The cesm1_1_beta04 lt_archive script did not create parent dirs
# if they did not already exist.  To fix the script, edit:
#   ${cesmroot}/scripts/ccsm_utils/Tools/lt_archive.csh 
# and change 'mkdir' to 'mkdir -p'.  This is fixed in more recent
# versions of the code.

# ====================================================================
# ====  Set case options
# ====================================================================

# case will be used many ways;
#    directory and file names, both locally and on HPSS, and
#    script names; so consider it's length and information content.
# num_instances:  Number of ensemble members
# reuse_existing_case:
#    false; Remove $caseroot and $exeroot and rebuild
#    true;  configure -cleannamelist

setenv case                 Exp1
setenv compset              F_2000
setenv resolution           f09_f09
setenv cesmtag              cesm1_1_beta04
setenv num_instances        4
setenv reuse_existing_case  false

# ====================================================================
# define machines and directories
#
# mach            Computer name
# cesm_datadir    Root path of the public CESM data files
# cesmroot        Location of the cesm code base
#                 For cesm1_1_beta04 on hopper, MUST use 'nscollin' value provided.
# DARTroot        Location of DART code tree.
#                    Executables, scripts and input in $DARTroot/models/cam/...
# caseroot        Your (future) cesm case directory, where this CESM+DART will be built.
#                    Preferably not a frequently scrubbed location.
#                 caseroot will be deleted if reuse_existing_case is false (below)
#                 So this script, and other useful things should be kept elsewhere.
# exeroot         (Future) Run-time directory; scrubbable, large amount of space needed.
# archdir         (Future) Short-term archive directory
#                    until the long-term archiver moves it to permanent storage.
# note:
#         hopper organizes the scratch space into 26 subdirs based on the
#         first letter of the username, and then below that has the actual
#         user's scratch dir.  $firstlet below is the first letter of the user name.
# ====================================================================

setenv mach         hopp2
setenv cesm_datadir /project/projectdirs/ccsm1/inputdata
setenv cesmroot     /global/u1/n/nscollin/${cesmtag}

set firstlet = `echo $USER | cut -c1`

setenv DARTroot     /global/u1/${firstlet}/${USER}/DART/development

setenv caseroot     /global/u1/${firstlet}/${USER}/cases/${case}
setenv exeroot      /scratch/scratchdirs/${USER}/case}
setenv archdir      /scratch/scratchdirs/${USER}/archive/${case}

# ======================
# configure settings
# ====================================================================

# yyyy-mm-dd
setenv run_startdate 2008-10-31

setenv sst_dataset ${cesm_datadir}/atm/cam/sst/sst_HadOIBl_bc_0.9x1.25_1850_2011_c110307.nc
setenv year_start  1850
setenv year_end    2010

# ====================================================================
# runtime settings --  How many assimilation steps will be done after this one
#
# stop_option   Units for determining the forecast length between assimilations
#               Changing stop_option requires changes to user_nl_cam below.
# stop_n        Number of time units in the forecast
# ====================================================================

setenv resubmit      0
setenv stop_option   nhours
setenv stop_n        6

# ====================================================================
# job settings
#
# timewall    can be changed during a series by changing the ${case}.${mach}.run
# queue   can be changed during a series by changing the ${case}.${mach}.run
#         lrg_ queues are used in order to fit more instances on each node.
#         FV 1-degree can comfortably fit 4 instances on 1 lrg_ node (~60 gbyte)
#         On bluefire the regular queue (or higher) is probably necessary,
#         because it appears that there are not separate queues for the lrg memory
#         and the regular memory nodes.  So economy jobs requesting smaller numbers
#         of processors seem to prevent this lrg_economy job (20 nodes for 1-degree)
#         from running for long periods.
# ====================================================================

setenv timewall     1:00
setenv queue        regular

# ====================================================================
# set these standard commands based on the machine you are running on.
# ====================================================================

# NERSC "hopper" 
set   MOVE = '/bin/mv -fv'
set   COPY = '/bin/cp -fv --preserve=timestamps'
set   LINK = '/bin/ln -fvs'
set REMOVE = '/bin/rm -fr'


# ====================================================================
# Create the case.
#
# For list of the pre-defined cases: ./create_newcase -list
# To create a variant case, see the CESM documentation and carefully
# incorporate any needed changes into this script.
# ====================================================================

if ("${reuse_existing_case}" == "false") then
   echo "removing old files from ${caseroot} and ${exeroot}"
   ${REMOVE} ${caseroot}
   ${REMOVE} ${exeroot}
   ${cesmroot}/scripts/create_newcase -case ${caseroot} -mach ${mach} \
                   -res ${resolution} -compset ${compset} -skip_rundb

   if ( $status != 0 ) then
      echo "ERROR: Case could not be created."
      exit 1
   endif
else
   cd ${caseroot}
   ./configure  -cleannamelist
endif

# ====================================================================
# Configure the case.
# ====================================================================

cd ${caseroot}

./xmlchange -file env_build.xml    -id EXEROOT        -val ${exeroot}
./xmlchange -file env_build.xml    -id USE_ESMF_LIB   -val TRUE
./xmlchange -file env_build.xml    -id ESMF_LIBDIR    -val /scratch/scratchdirs/nscollin/esmf-mpi

# num_tasks_per_instance = #tasks_node / #instances_node
# hopper: #tasks_node = 4 with threading of 6/task = 24 processors
#         #instances_node = 1-degree: 2 on standard memory nodes 
#                           (using lrg mem nodes requires asking for ALL of them)
set num_tasks_per_node = 4
set num_tasks_per_instance = 2
set num_threads = 6
# This is hard-wiring for the current (1/17/2011) multi-instance CESM restriction
# that all instances must be advanced simultaneously.  Work is underway to relax that.
@ total_nt = $num_instances * $num_tasks_per_instance
echo "total MPI tasks requested = $total_nt"
echo "num_threads = $num_threads"

# Atm gets all the nodes and runs.
# Lnd gets all the nodes and runs.
# The other components divide them up.
# ? ? ?  What about sglc?  It's a stub, and doesn't matter what pes are assigned to it.
# This algorithm figures out whether there are enough processors requested
# to run each component on whole nodes, or the components need to share some nodes.
# The hard-coded numbers (and ratios between them) are estimates; change them if you
# know better.
@ atm_pes  = $total_nt
@ large_small = $total_nt / (8 * $num_tasks_per_node)
if ($large_small > 0) then
   # Large_small > 0 means there are at least 8 nodes requested.
   # Allot whole nodes to the major components.
   # hopp2 doesn't use many 'tasks', so it needs a different distribution than on bluefire.
   # In particular, we need to keep num_instances > num_tasks for atm, lnd, and ice.
   @ docn_pes = $num_tasks_per_node
   @ cpl_pes  = ( $total_nt / (3 * $num_tasks_per_node) ) * $num_tasks_per_node
   @ cice_pes = $total_nt - ($cpl_pes + $docn_pes)
   @ lnd_pes  = $total_nt 
else
   # 40% cpl,  40% cice, 20% docn, 100% lnd.  These may occupy fractions of nodes.
   @ cpl_pes  = (2 * $total_nt) / 5
   @ cice_pes = (2 * $total_nt) / 5
   @ docn_pes = $total_nt - ($cpl_pes + $cice_pes)
   @ lnd_pes  = $total_nt 
endif

echo "task layout"
echo "[0 ......................... ATM ............................. $atm_pes]"

# first pe of each component (counted from 0)
@ atm_rootpe  = 0
@ lnd_rootpe  = 0
@ cpl_rootpe  = 0
@ cice_rootpe = $cpl_rootpe  + $cpl_pes
@ docn_rootpe = $cice_rootpe + $cice_pes
echo "[$cpl_rootpe ... CPL ... $cice_rootpe ... ICE ... $docn_rootpe ... OCN ... $total_nt]"
echo "[$lnd_rootpe ... LND ..................................................... $total_nt]"
echo ""
echo "ATM gets $atm_pes"
echo "ICE gets $cice_pes"
echo "LND gets $lnd_pes"
echo "CPL gets $cpl_pes"
echo "OCN gets $docn_pes"
echo ""

./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $atm_pes
./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val $num_threads
./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val $atm_rootpe
./xmlchange -file env_mach_pes.xml -id  NINST_ATM -val $num_instances

./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val $cpl_pes
./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val $num_threads
./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL -val $cpl_rootpe

./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $cice_pes
./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val $num_threads
./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val $cice_rootpe
./xmlchange -file env_mach_pes.xml -id  NINST_ICE -val $num_instances

./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $docn_pes
./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val $num_threads
./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN -val $docn_rootpe

./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $lnd_pes
./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val $num_threads
./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val $lnd_rootpe
./xmlchange -file env_mach_pes.xml -id  NINST_LND -val $num_instances

./xmlchange -file env_conf.xml -id RUN_TYPE                -val startup
./xmlchange -file env_conf.xml -id RUN_STARTDATE           -val $run_startdate
./xmlchange -file env_conf.xml -id DOCN_SSTDATA_FILENAME   -val $sst_dataset
./xmlchange -file env_conf.xml -id DOCN_SSTDATA_YEAR_START -val $year_start
./xmlchange -file env_conf.xml -id DOCN_SSTDATA_YEAR_END   -val $year_end
./xmlchange -file env_conf.xml -id CLM_CONFIG_OPTS         -val '-rtm off'
# The river transport model ON is useful only when using an active ocean or
# land surface diagnostics.

# Substantial archiving changes exist in the Tools/st_archive.sh script.
# Do not change the CALENDAR or the CONTINUE_RUN
# DOUT_S     is to turn on/off the short-term archiving
# DOUT_L_MS  is to store to the HPSS (formerly "MSS")

./xmlchange -file env_run.xml -id CONTINUE_RUN               -val FALSE
./xmlchange -file env_run.xml -id RESUBMIT                   -val $resubmit
./xmlchange -file env_run.xml -id STOP_OPTION                -val $stop_option
./xmlchange -file env_run.xml -id STOP_N                     -val $stop_n
./xmlchange -file env_run.xml -id CALENDAR                   -val GREGORIAN
./xmlchange -file env_run.xml -id DOUT_S_ROOT                -val ${archdir}
./xmlchange -file env_run.xml -id DOUT_S                     -val TRUE
./xmlchange -file env_run.xml -id DOUT_S_SAVE_INT_REST_FILES -val TRUE
./xmlchange -file env_run.xml -id DOUT_L_MS                  -val TRUE
./xmlchange -file env_run.xml -id DOUT_L_MSROOT              -val "csm/${case}"
./xmlchange -file env_run.xml -id DOUT_L_HTAR                -val FALSE

# ====================================================================
# Create namelist template: user_nl_cam, user_nl_clm
# ====================================================================

cd ${caseroot}

cat <<EOF >! user_nl_cam
&camexp
 inithist             = 'ENDOFRUN'
 div24del2flag        = 4
 empty_htapes         = .true.
 fincl1               = 'PHIS:I'
 nhtfrq               = -$stop_n
 iradae               = -$stop_n
 aerodep_flx_datapath = '${cesm_datadir}/atm/cam/chem/trop_mozart_aero/aero'
 aerodep_flx_file     = 'aerosoldep_monthly_1849-2006_1.9x2.5_c090803.nc'
 aerodep_flx_cycle_yr = 2000
 aerodep_flx_type     = 'CYCLICAL'
/
EOF

cat <<EOF >! user_nl_clm
&clmexp
  fatmgrid = '${cesm_datadir}/lnd/clm2/griddata/griddata_0.9x1.25_070212.nc'
  faerdep  = '/scratch/scratchdirs/nscollin/cesm_datafiles/aerosoldep_rcp4.5_monthly_1849-2104_0.9x1.25_c100407.nc'
  outnc_large_files = .true.
  hist_empty_htapes = .true.
/
EOF
#  faerdep  = '${cesm_datadir}/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_rcp4.5_monthly_1849-2104_0.9x1.25_c100407.nc'
#  This resolution doesn't exist on hopper yet.

# ====================================================================
# Update source files if need be
#    Ideally, using DART will not require any source mods.
#    Until then, this script accesses source mods from a hard-wired location below.
#    Those may eventually be packaged into the DART repository.
#    If you have additional source mods, they will need to be merged into any DART
#    mods and put in the SourceMods subdirectory found in the 'case' directory.
# ====================================================================
cp /global/u1/n/nscollin/cesm_mods/seq*F90           $caseroot/SourceMods/src.drv
cp /global/u1/n/nscollin/cesm_mods/hist*F90          $caseroot/SourceMods/src.clm
cp /global/u1/n/nscollin/cesm_mods/ccsm_comp_mod.F90 $caseroot/SourceMods/src.drv

# this one needs a recursive copy to get all the files in the subdirs
#${COPY} -r  ~thoar/${cesmtag}/SourceMods/* ${caseroot}/SourceMods/
#if ( $status == 0) then
#   echo "FYI - Local Source Modifications used for this case:"
#   ls -lr ${caseroot}/SourceMods/*
#else
#   echo "FYI - No SourceMods for this case"
#endif
#
# ====================================================================
# Configure
# ====================================================================

cd ${caseroot}

./configure -case

if ( $status != 0 ) then
   echo "ERROR: Case could not be configured."
   exit 2
endif

# ====================================================================
# Stage the required parts of DART in the caseroot directory.
# ====================================================================

cd ${caseroot}

if ("${reuse_existing_case}" == "false") then
   ${MOVE} Tools/st_archive.sh Tools/st_archive.sh.orig
endif
${COPY} ${DARTroot}/models/cam/shell_scripts/st_archive.sh Tools/

# The cesm1_1_beta04 release had an error in that it did not
# provide the lt_archive.csh script, and the one in the repos
# did not have the -p flag, which is a good idea. So, for now ...
${COPY} ${cesmroot}/scripts/ccsm_utils/Tools/lt_archive.csh Tools/

${COPY} ${DARTroot}/models/cam/shell_scripts/assimilate.csh  .
${COPY} ${DARTroot}/models/cam/work/input.nml                .

# Ensure that the input.nml ensemble size matches the number of instances.
# WARNING: the output files contain ALL ensemble members ==> BIG

ex input.nml <<ex_end
g;ens_size ;s;= .*;= $num_instances;
g;num_output_state_members ;s;= .*;= $num_instances;
g;num_output_obs_members ;s;= .*;= $num_instances;
wq
ex_end

# ====================================================================
# Update the scripts that build the namelists.
# The active components scripts need to support the multi-instance naming.
# ====================================================================

echo ''
echo 'Editing the Buildconf/{cam,cice,clm}.buildnml.csh files'
echo ''

cd ${caseroot}/Buildconf

${COPY}  cam.buildnml.csh  cam.buildnml.csh.orig
${COPY} cice.buildnml.csh cice.buildnml.csh.orig
${COPY}  clm.buildnml.csh  clm.buildnml.csh.orig

# The CAM buildnml script only needs changing in one place.

ex cam.buildnml.csh <<ex_end
/cam_inparm/
/ncdata/
s;= '.*';= "cam_initial_\${atm_inst_counter}.nc";
wq
ex_end

# The CICE buildnml script only needs changing in one place.

ex cice.buildnml.csh <<ex_end
/setup_nml/
/ice_ic/
s;= '.*';= "ice_restart_\${ice_inst_counter}.nc";
wq
ex_end

# The CLM buildnml script needs changing in MULTIPLE places.

@ n = 1
while ($n <= $num_instances)
   set inst = `printf "%04d" $n`
   ex clm.buildnml.csh <<ex_end
/lnd_in_$inst/
/finidat/
s;= '.*';= "clm_restart_${n}.nc";
wq
ex_end
   @ n++
end

chmod 0755 clm.buildnml.csh

# ====================================================================
# The *.run script must be modified to call the DART assimilate script.
# The modifications are contained in a "here" document that MUST NOT
# expand the wildcards etc., before it is run. This is achieved by
# double-quoting the characters used to delineate the start/stop of
# the "here" document. No kidding. It has to be "EndOfText",
# not 'EndOfText' or EndOfText.
# ====================================================================

cd ${caseroot}

echo ''
echo 'Adding the call to assimilate.csh to the *.run script.'
echo ''

cat << "EndOfText" >! add_to_run.txt

# -------------------------------------------------------------------------
# START OF DART: if CESM finishes correctly (pirated from ccsm_postrun.csh);
# perform an assimilation with DART.
# -------------------------------------------------------------------------

set CplLogFile = `ls -1t cpl.log* | head -n 1`
if ($CplLogFile == "") then
 echo 'ERROR: Model did not complete - no cpl.log file present - exiting'
 echo 'ERROR: Assimilation will not be attempted.'
 exit -4
endif

grep 'SUCCESSFUL TERMINATION' $CplLogFile
if ( $status == 0 ) then
  ${CASEROOT}/assimilate.csh

  if ( $status == 0 ) then
     echo "`date` -- DART HAS FINISHED"
  else
     echo "`date` -- DART FILTER ERROR - ABANDON HOPE"
     exit -5
  endif
endif

# END OF DART BLOCK
# -------------------------------------------------------------------------

"EndOfText"

# Now that the "here" document is created,
# determine WHERE to insert it -- ONLY IF it is not already there.

grep "ABANDON HOPE" ${case}.${mach}.run
set STATUSCHECK = $status

if ( ${STATUSCHECK} == 0 ) then
   echo "DART block already present in ${case}.${mach}.run"
else if ( ${STATUSCHECK} == 1 ) then

   set MYSTRING = `grep --line-number "CSM EXECUTION HAS FINISHED" ${case}.${mach}.run`
   set MYSTRING = `echo $MYSTRING | sed -e "s#:# #g"`

   @ origlen = `cat ${case}.${mach}.run | wc -l`
   @ keep = $MYSTRING[1]
   @ lastlines = $origlen - $keep

   mv ${case}.${mach}.run ${case}.${mach}.run.orig

   head -n $keep      ${case}.${mach}.run.orig >! ${case}.${mach}.run
   cat                add_to_run.txt           >> ${case}.${mach}.run
   tail -n $lastlines ${case}.${mach}.run.orig >> ${case}.${mach}.run

else
   echo "ERROR in grep of ${case}.${mach}.run: aborting"
   echo "status was ${STATUSCHECK}"
   exit 10
endif

# ====================================================================
# Edit the run script to reflect project, queue, and wallclock
# ====================================================================

echo ''
echo 'Updating the run script to set the project number, wallclock time'
echo 'and queue name.'
echo ''

if ($?proj) then
  set PROJ=`grep '^#PBS ' $case.$mach.run | grep -e '-P' `
  sed -e s/$PROJ[3]/$proj/ < $case.$mach.run >! temp
  /bin/mv temp  $case.$mach.run
endif

if ($?timewall) then
  set TIMEWALL=`grep '^#PBS ' $case.$mach.run | grep walltime `
  sed -e /"${TIMEWALL}"/s/=.\$/=$timewall/ < $case.$mach.run >! temp
  /bin/mv temp  $case.$mach.run
endif

if ($?queue) then
  set QUEUE=`grep '^#PBS ' $case.$mach.run | grep -e '-q' `
  sed -e s/$QUEUE[3]/$queue/ < $case.$mach.run >! temp
  /bin/mv temp  $case.$mach.run
endif

chmod 0755 ${case}.${mach}.run

# ====================================================================
# IMPORTANT: All resubmits must be type 'startup'.
# Change Tools/ccsm_postrun.csh line 83 to CONTINUE_RUN -val FALSE'
# ====================================================================

cd ${caseroot}/Tools

echo ''
echo 'Changing Tools/ccsm_postrun.csh such that all the resubmits are "startup",'
echo 'which means CONTINUE_RUN should be FALSE in ccsm_postrun.csh'
echo ''

ex ccsm_postrun.csh <<ex_end
/use COMP_RUN_BARRIERS as surrogate for timing run logical/
/CONTINUE_RUN/
s;TRUE;FALSE;
wq
ex_end

# ====================================================================
# build
# ====================================================================

cd ${caseroot}

echo ''
echo 'Building the case'
echo ''

./$case.$mach.build

if ( $status != 0 ) then
   echo "ERROR: Case could not be built."
   exit 3
endif

# ====================================================================
# Stage the required parts of DART in the execution root directory,
# now that EXEROOT exists.
# ====================================================================

foreach FILE ( filter cam_to_dart dart_to_cam )
   ${COPY} ${DARTroot}/models/cam/work/${FILE} ${exeroot}/
   if ( $status != 0 ) then
      echo "ERROR: ${DARTroot}/models/cam/work/${FILE} not copied to ${exeroot}"
      echo "ERROR: ${DARTroot}/models/cam/work/${FILE} not copied to ${exeroot}"
      exit 3
   endif
end

# ====================================================================
# Stage the restarts now that the run directory exists
# ====================================================================

# Perfect model observations assimilation
#set stagedir = /scratch/scratchdirs/nscollin/ned_datafiles
# Real observations assimilation
set stagedir = /scratch/scratchdirs/nscollin/tim_datafiles

echo ''
echo "Copying the restart files from ${stagedir}"
echo 'into the CESM run directory.'
echo ''

@ n = 1
while ($n <= $num_instances)
  echo "Staging restarts for instance $n of $num_instances"
  ${COPY} ${stagedir}/CAM/caminput_${n}.nc ${exeroot}/run/cam_initial_${n}.nc
  ${COPY} ${stagedir}/CLM/clminput_${n}.nc ${exeroot}/run/clm_restart_${n}.nc
  ${COPY} ${stagedir}/ICE/iceinput_${n}.nc ${exeroot}/run/ice_restart_${n}.nc
  @ n++
end

# This script will copy an existing prior_inflate_restart to the run dir if found.
if (-f ${stagedir}/DART/prior_inflate_restart) then
   ${COPY} ${stagedir}/DART/prior_inflate_restart ${exeroot}/run/prior_inflate_restart.${run_startdate}-00000
   echo ''
   echo "${stagedir}/DART/prior_inflate_restart has been copied to "
   echo "${exeroot}/run/prior_inflate_restart.${run_startdate}-00000"
   echo 'If that has the wrong state vector, you will need to replace it before running.'
   echo ''
else
   echo ''
   echo 'If using inflation in DART you may need to copy an inflation restart file'
   echo "to ${exeroot}/run/prior_inflate_restart.${run_startdate}-00000"
   echo 'before running.  It must include the exact fields as your DART state vector.'
   echo "You can make one with ${DARTroot}/models/cam/work/fill_inflation_restart"
   echo ''
endif


# only warn people if a precomputed final_full for this number of instances 
# does not already exist.
if (! -f ${DARTroot}/system_simulation/final_full_precomputed_tables/final_full.${num_instances}) then
   echo ''
   echo 'If you are using the DART sampling error correction feature'
   echo 'the assimilate.csh script will expect to copy this file:'
   echo "${DARTroot}/system_simulation/final_full_precomputed_tables/final_full.${num_instances}"
   echo 'and it does not exist for your number of ensemble members.'
   echo "Generate one by building and running ${DARTroot}/system_simulation/work/full_error"
   echo 'with the namelist set to your ensemble size.'
   echo ''
endif

# ====================================================================
# What to do next
# ====================================================================

echo ''
echo 'Time to check the case.'
echo "cd into ${caseroot}"
echo 'Modify what you like in input.nml, make sure the observation directory'
echo 'names set in assimilate.csh match those on your system, and submit'
echo 'the CESM job by running:'
echo "./$case.$mach.submit"
echo ''

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

