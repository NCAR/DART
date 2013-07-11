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
# This version of the setup script is designed for a B compset where
# CAM, POP, and CLM are all active and going to be assimilated separately.
# The assimilation happens once a day for all components at 0Z.
#
# It depends on a couple patched utility files in ~nancy/cesm1_1_1
# Copy them over to your own $HOME/cesm1_1_1 dir before running this script.
#
# This script is designed to configure and build a multi-instance CESM model.
# It will use DART to assimilate observations at regular intervals.
# This script will build the DART executables first if they are not found.
# 
# Until the binary POP restart files have been converted from big-endian
# to little-endian, you MUST COMPILE DART with intel and -convert big-endian
# contact dart@ucar.edu if you want to use another compiler.
#
# This script relies heavily on the information in:
# http://www.cesm.ucar.edu/models/cesm1.1/cesm/doc/usersguide/book1.html
#
# ---------------------
# How to set up the script
# ---------------------
# -- Either edit and run this script in the $DART/models/CESM/shell_scripts
#    directory where it first appears,
#    or copy it to somewhere that it will be preserved and run it there.
#    It will create a 'case' directory, where the model will be built,
#    and an execution directory, where each forecast and assimilation will
#    take place.  The short term archiver will use a third directory for
#    storage of model output until it can be moved to long term storage (HPSS)
# -- Examine the whole script to identify things to change for your experiments.
# -- Provide any initial files needed by your run:
#       inflation
#       sampling error correction
#       CESM initial ensemble
#       ...
# -- Run this script.
# -- Edit the DART input.nml that appears in the $CASEROOT directory.
# -- Submit the job using $CASEROOT/${case}.submit
#
# ---------------------
# Important features
# ---------------------
#
# If you want to change something in your case other than the runtime
# settings, it is safest to delete everything and start the run from scratch.
# For the brave, read
#
# http://www.cesm.ucar.edu/models/cesm1.1/cesm/doc/usersguide/x1142.html
#
# and you may be able to salvage something with
# ./cesm_setup -clean
# ./cesm_setup
# ./${case}.clean_build
# ./${case}.build
#
# ==============================================================================
# ====  Set case options
# ==============================================================================

# case will be used many ways;
#    directory and file names, both locally and on HPSS, and
#    script names; so consider it's length and information content.
# num_instances:  Number of ensemble members

setenv case                 cesm_testme
setenv compset              B_2000_CAM5
setenv resolution           0.9x1.25_gx1v6
setenv cesmtag              cesm1_1_1
setenv num_instances        3

# ==============================================================================
# define machines and directories
#
# mach            Computer name
# cesm_datadir    Root path of the public CESM data files
# cesmroot        Location of the cesm code base
#                 For cesm1_1 on yellowstone
# DARTroot        Location of DART code tree.
#                    Executables, scripts and input in $DARTroot/models/dev/...
# caseroot        Your (future) cesm case directory, where this CESM+DART will be built.
#                    Preferably not a frequently scrubbed location.
#                    This script will delete any existing caseroot, so this script,
#                    and other useful things should be kept elsewhere.
# rundir          (Future) Run-time directory; scrubbable, large amount of space needed.
# exeroot         (Future) directory for executables - scrubbable, large amount of space needed.
# archdir         (Future) Short-term archive directory
#                    until the long-term archiver moves it to permanent storage.
# ==============================================================================

setenv mach         yellowstone
setenv cesm_datadir /glade/p/cesm/cseg/inputdata
setenv cesmroot     /glade/p/cesm/cseg/collections/$cesmtag
setenv caseroot     /glade/p/work/${USER}/cases/${case}
setenv exeroot      /glade/scratch/${USER}/${case}/bld
setenv rundir       /glade/scratch/${USER}/${case}/run
setenv archdir      /glade/scratch/${USER}/archive/${case}

setenv DARTroot     /glade/u/home/${USER}/DART

set RTM_stagedir = /glade/scratch/thoar/DART_POP_RESTARTS/2004-01-01-00000
set CLM_stagedir = /glade/scratch/thoar/DART_POP_RESTARTS/CLM_2004-01-01-00000/cesm_test
set CAM_stagedir = /glade/p/cesm/cseg/inputdata/atm/cam/inic/fv
set POP_stagedir = /glade/p/work/aliciak/DART_IC/CCSM4_ensembles/rest/2004-01-01-00000

# ==============================================================================
# configure settings ... run_startdate format is yyyy-mm-dd
# ==============================================================================

setenv refyear     2004
setenv refmon      01
setenv refday      01
setenv run_reftod  00000
setenv run_refdate $refyear-$refmon-$refday

# ==============================================================================
# runtime settings --  How many assimilation steps will be done after this one
#                      plus archiving options
#
# resubmit      How many job steps to run on continue runs (will be 0 initially)
# stop_option   Units for determining the forecast length between assimilations
# stop_n        Number of time units in the first forecast
# assim_n       Number of time units between assimilations
# ==============================================================================

setenv resubmit      10
setenv stop_option   nhours
setenv stop_n        72
setenv assim_n       24
setenv short_term_archiver on
setenv long_term_archiver  on

# ==============================================================================
# job settings
#
# timewall   can be changed during a series by changing the ${case}.run
# queue      can be changed during a series by changing the ${case}.run
#
# TJH: How many T62_gx1v6 CESM instances can fit on 1 node?
# ==============================================================================

setenv ACCOUNT      P8685nnnn
setenv timewall     0:50
setenv queue        economy
setenv ptile        15

# ==============================================================================
# set these standard commands based on the machine you are running on.
# ==============================================================================

set nonomatch       # suppress "rm" warnings if wildcard does not match anything

# The FORCE options are not optional.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following
switch ("`hostname`")
   case be*:
      # NCAR "bluefire"
      set   MOVE = '/usr/local/bin/mv -fv'
      set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
      set   LINK = '/usr/local/bin/ln -fvs'
      set REMOVE = '/usr/local/bin/rm -fr'

   breaksw
   default:
      # NERSC "hopper", NWSC "yellowstone"
      set   MOVE = '/bin/mv -fv'
      set   COPY = '/bin/cp -fv --preserve=timestamps'
      set   LINK = '/bin/ln -fvs'
      set REMOVE = '/bin/rm -fr'

   breaksw
endsw

# ==============================================================================
# some simple error checking before diving into the work
# ==============================================================================

# fatal idea to make caseroot the same dir as where this setup script is
# since the build process removes all files in the caseroot dir before
# populating it.  try to prevent shooting yourself in the foot.
if ( $caseroot == `dirname $0` ) then
   echo "ERROR: the setup script should not be located in the caseroot"
   echo "directory, because all files in the caseroot dir will be removed"
   echo "before creating the new case.  move the script to a safer place."
   exit -1
endif

# make sure these directories exist
set musthavedirs = "cesm_datadir cesmroot DARTroot"
foreach VAR ( $musthavedirs )
   # VAR is the shell variable name, DIR is the value
   set DIR = `eval echo \${$VAR}`
   if ( ! -d $DIR ) then
      echo "ERROR: directory '$DIR' not found"
      echo " In the setup script check the setting of: $VAR "
      exit -1
   endif
end

# make sure there is a filter in these dirs.   try to build
# them if we can't find filter already built for each model.
set musthavefiles = "cam POP clm CESM"
foreach MODEL ( $musthavefiles )
   set targetdir = $DARTroot/models/$MODEL/work
   if ( ! -x $targetdir/filter ) then
      echo "WARNING: executable file 'filter' not found"
      echo " Looking for: $targetdir/filter "
      echo " Trying to rebuild all model files now."
      (cd $targetdir; ./quickbuild.csh -mpi)
      if ( ! -x $targetdir/filter ) then
         echo "ERROR: executable file 'filter' not found"
         echo " Unsuccessfully tried to rebuild: $targetdir/filter "
         echo " Required DART assimilation executables are not found "
         exit -1
      endif
   endif
end

# ==============================================================================
# Create the case.
#
# For list of the pre-defined cases: ./create_newcase -list
# To create a variant case, see the CESM documentation and carefully
# incorporate any needed changes into this script.
# ==============================================================================

   echo "removing old files from ${caseroot}"
   echo "removing old files from ${exeroot}"
   echo "removing old files from ${rundir}"
   ${REMOVE} ${caseroot}
   ${REMOVE} ${exeroot}
   ${REMOVE} ${rundir}
   ${cesmroot}/scripts/create_newcase -case ${caseroot} -mach ${mach} \
                   -res ${resolution} -compset ${compset}

   if ( $status != 0 ) then
      echo "ERROR: Case could not be created."
      exit -1
   endif

# ==============================================================================
# Configure the case - this creates the CASEROOT directory.
# ==============================================================================

cd ${caseroot}

# Save a copy for debug purposes
foreach FILE ( *xml )
   if ( ~ -e        ${FILE}.original ) then
      ${COPY} $FILE ${FILE}.original
   endif
end

if ( $num_instances < 10) then

   # This is only for the purpose of debugging the code.
   # A more efficient layout must be found
   @ atm_pes = $ptile * $num_instances * 4
   @ ocn_pes = $ptile * $num_instances * 4
   @ lnd_pes = $ptile * $num_instances * 4
   @ ice_pes = $ptile * $num_instances * 1
   @ glc_pes = $ptile * $num_instances
   @ rof_pes = $ptile * $num_instances
   @ cpl_pes = $ptile * 4

else

   # This is only for the purpose of debugging the code.
   # A more efficient layout must be found
   #
   @ atm_pes = $ptile * $num_instances * 2
   @ ocn_pes = $ptile * $num_instances * 2
   @ lnd_pes = $ptile * $num_instances * 2
   @ ice_pes = $ptile * $num_instances
   @ glc_pes = $ptile * $num_instances
   @ rof_pes = $ptile * $num_instances
   @ cpl_pes = $ptile * $num_instances

endif

#echo "task partitioning ... atm+ocn // lnd+ice+glc+rof"
echo ""
echo "ATM  gets $atm_pes"
echo "CPL  gets $cpl_pes"
echo "ICE  gets $ice_pes"
echo "LND  gets $lnd_pes"
echo "GLC  gets $glc_pes"
echo "DROF gets $rof_pes"
echo "OCN  gets $ocn_pes"
echo ""

./xmlchange NTHRDS_CPL=1,NTASKS_CPL=$cpl_pes
./xmlchange NTHRDS_GLC=1,NTASKS_GLC=$glc_pes,NINST_GLC=1
./xmlchange NTHRDS_ATM=1,NTASKS_ATM=$atm_pes,NINST_ATM=$num_instances
./xmlchange NTHRDS_LND=1,NTASKS_LND=$lnd_pes,NINST_LND=$num_instances
./xmlchange NTHRDS_ICE=1,NTASKS_ICE=$ice_pes,NINST_ICE=$num_instances
./xmlchange NTHRDS_ROF=1,NTASKS_ROF=$rof_pes,NINST_ROF=$num_instances
./xmlchange NTHRDS_OCN=1,NTASKS_OCN=$ocn_pes,NINST_OCN=$num_instances
./xmlchange ROOTPE_ATM=0
./xmlchange ROOTPE_OCN=0
./xmlchange ROOTPE_CPL=0
./xmlchange ROOTPE_LND=0
./xmlchange ROOTPE_ICE=0
./xmlchange ROOTPE_GLC=0
./xmlchange ROOTPE_ROF=0

# http://www.cesm.ucar.edu/models/cesm1.1/cesm/doc/usersguide/c1158.html#run_start_stop
# "A hybrid run indicates that CESM is initialized more like a startup, but uses
# initialization datasets from a previous case. This is somewhat analogous to a
# branch run with relaxed restart constraints. A hybrid run allows users to bring
# together combinations of initial/restart files from a previous case (specified
# by $RUN_REFCASE) at a given model output date (specified by $RUN_REFDATE).
# Unlike a branch run, the starting date of a hybrid run (specified by $RUN_STARTDATE)
# can be modified relative to the reference case. In a hybrid run, the model does not
# continue in a bit-for-bit fashion with respect to the reference case. The resulting
# climate, however, should be continuous provided that no model source code or
# namelists are changed in the hybrid run. In a hybrid initialization, the ocean
# model does not start until the second ocean coupling (normally the second day),
# and the coupler does a "cold start" without a restart file.

./xmlchange RUN_TYPE=hybrid
./xmlchange RUN_STARTDATE=$run_refdate
./xmlchange START_TOD=$run_reftod
./xmlchange RUN_REFDATE=$run_refdate
./xmlchange RUN_REFTOD=$run_reftod
./xmlchange GET_REFCASE=FALSE
./xmlchange EXEROOT=${exeroot}

./xmlchange CALENDAR=GREGORIAN

./xmlchange STOP_OPTION=$stop_option
./xmlchange STOP_N=$stop_n
./xmlchange CONTINUE_RUN=FALSE
./xmlchange RESUBMIT=0
./xmlchange PIO_TYPENAME=pnetcdf

./xmlchange CLM_CONFIG_OPTS='-bgc cn'

if ($short_term_archiver == 'off') then
   ./xmlchange DOUT_S=FALSE
else
   ./xmlchange DOUT_S=TRUE
   ./xmlchange DOUT_S_ROOT=${archdir}
   ./xmlchange DOUT_S_SAVE_INT_REST_FILES=FALSE
endif
if ($long_term_archiver == 'off') then
   ./xmlchange DOUT_L_MS=FALSE
else
   ./xmlchange DOUT_L_MS=TRUE
   ./xmlchange DOUT_L_MSROOT="csm/${case}"
   ./xmlchange DOUT_L_HTAR=FALSE
endif

# level of debug output, 0=minimum, 1=normal, 2=more, 3=too much, valid values: 0,1,2,3 (integer)

./xmlchange DEBUG=FALSE
./xmlchange INFO_DBUG=0


# ==============================================================================
# Set up the case.
# This creates the EXEROOT and RUNDIR directories.
# ==============================================================================

./cesm_setup

if ( $status != 0 ) then
   echo "ERROR: Case could not be set up."
   exit -2
endif

# this should be removed when the files are fixed:
echo "PATCHING ADDITIONAL FILES HERE - SHOULD BE REMOVED WHEN FIXED"
echo caseroot is ${caseroot}
if (    -d     ~/${cesmtag} ) then
   ${COPY} ~/${cesmtag}/*buildnml.csh ${caseroot}/Buildconf/
   ls -l ~/${cesmtag}/*buildnml.csh ${caseroot}/Buildconf/*buildnml.csh
   ${COPY} ~/${cesmtag}/preview_namelists ${caseroot}
   ls -l ~/${cesmtag}/preview_namelists ${caseroot}/preview_namelists
endif

# ==============================================================================
# Modify namelist templates for each instance.
# ==============================================================================

@ inst = 1
while ($inst <= $num_instances)

   # instance now includes the leading underscore
   set instance  = `printf _%04d $inst`

   # special for some files: no leading underscore, 2 digits
   set instance2  = `printf %02d $inst`

   # ===========================================================================
   set fname = "user_nl_cam${instance}"
   # ===========================================================================
   # For a HOP TEST ... empty_htapes = .false.
   # For a HOP TEST ... use a default fincl1

   echo " inithist      = 'DAILY'"                      >> ${fname}
   echo " ncdata        = 'cam_initial${instance}.nc'"  >> $fname
   echo " empty_htapes  = .true. "                      >> ${fname}
   echo " fincl1        = 'PHIS:I' "                    >> ${fname}
   echo " nhtfrq        = -$assim_n "                   >> ${fname}
   echo " mfilt         = 1 "                           >> ${fname}

   # ===========================================================================
   set fname = "user_nl_pop2${instance}"
   # ===========================================================================

   # POP Namelists
   # init_ts_suboption = 'data_assim'   for non bit-for-bit restarting (assimilation mode)
   # init_ts_suboption = 'rest'         for
   # init_ts_suboption = 'spunup'       for
   # init_ts_suboption = 'null'         for
   # For a HOP TEST (untested)... tavg_file_freq_opt = 'nmonth' 'nday' 'once'"
   # For a HOP TEST ... cool to have restart files every day, not just for end.

   echo "init_ts_suboption = 'null'" >> $fname

   # ===========================================================================
   set fname = "user_nl_cice${instance}"
   # ===========================================================================
   # CICE Namelists
   
   echo "ice_ic = 'b40.20th.005_ens${instance2}.cice.r.2004-01-01-00000.nc'" >> $fname

   # ===========================================================================
   set fname = "user_nl_clm${instance}"
   # ===========================================================================

   # Customize the land namelists
   # The history tapes are a work in progress. If you write out the instantaneous
   # flux variables every 30 minutes to the .h1. file, the forward observation
   # operators for these fluxes should just read them from the .h1. file rather
   # than trying to create them from the (incomplete DART) CLM state.
   # For a HOP TEST ... hist_empty_htapes = .false.
   # For a HOP TEST ... use a default hist_fincl1
   #
#  set CLM_stagedir = /glade/scratch/afox/bptmp/MD_40_PME/run/MD_40_PME

   @ thirtymin = $assim_n * 2

   echo "finidat = '${CLM_stagedir}.clm2${instance}.r.${run_refdate}-${run_reftod}.nc'" >> $fname
   echo "hist_empty_htapes = .true."                 >> $fname
   echo "hist_fincl1 = 'NEP'"                        >> $fname
   echo "hist_fincl2 = 'NEP','FSH','EFLX_LH_TOT_R'"  >> $fname
   echo "hist_nhtfrq = -$assim_n,1,"                 >> $fname
   echo "hist_mfilt  = 1,$thirtymin"                 >> $fname
   echo "hist_avgflag_pertape = 'A','A'"             >> $fname

   # ===========================================================================
   set fname = "user_nl_rtm${instance}"
   # ===========================================================================
   # RIVER RUNOFF CAN START FROM AN OLD CLM RESTART FILE

   echo "finidat_rtm = '${RTM_stagedir}/b40.20th.005_ens${instance2}.clm2.r.${run_refdate}-${run_reftod}.nc'" >> $fname

   @ inst ++
end

# This is expected to fail at the CLM stage (until clm.buildnml.csh is fixed)
# "build-namelist ERROR:: Can NOT set both -clm_startfile option AND finidat on namelist"
# problem is ... I need to specify finidat.

./preview_namelists

# ==============================================================================
# Update source files if need be
#    Ideally, using DART will not require any source mods.
#    Until then, this script accesses source mods from a hard-wired location below.
#    Those may eventually be packaged into the DART repository.
#    If you have additional source mods, they will need to be merged into any DART
#    mods and put in the SourceMods subdirectory found in the 'case' directory.
# ==============================================================================

if (    -d     ~/${cesmtag}/SourceMods ) then
   ${COPY} -r  ~/${cesmtag}/SourceMods/* ${caseroot}/SourceMods/
else
   echo "ERROR - No SourceMods for this case."
   echo "ERROR - No SourceMods for this case."
   echo "DART requires modifications to several src files."
   echo "These files can be downloaded from:"
   echo "http://www.image.ucar.edu/pub/DART/CESM/DART_SourceMods_cesm1_1_1.tar"
   echo "untar these into your HOME directory - they will create a"
   echo "~/cesm_1_1_1  directory with the appropriate SourceMods structure."
   exit -4
endif

# ==============================================================================
# build
# ==============================================================================

echo ''
echo 'Building the case'
echo ''

./${case}.build

if ( $status != 0 ) then
   echo "ERROR: Case could not be built."
   exit -5
endif

# ==============================================================================
# Stage the restarts now that the run directory exists
# I ran this compset once without setting user_nl_atm_00nn to see where the
# initial files come from.
# ==============================================================================

cd ${rundir}

echo ''
echo "Copying the restart files from the staging directories"
echo 'into the CESM run directory and creating the pointer files.'
echo ''

# TJH FIXME ... simply put full path in the pointer file? traceability?
# TJH FIXME ... what do we do with the *.hdr file
# After a run completes, the 'normal' pointer files created are:
# rpointer.drv
# rpointer.atm_0001           cesm_2000.cam_0001.r.2004-01-04-00000.nc
# rpointer.lnd_0001           cesm_2000.clm2_0001.r.2004-01-04-00000.nc
# rpointer.ice_0001           cesm_2000.cice_0001.r.2004-01-04-00000.nc
# rpointer.rof_0001           cesm_2000.rtm_0001.r.2004-01-04-00000.nc
# rpointer.ocn_0001.ovf       cesm_2000.pop_0001.ro.2004-01-04-00000
# rpointer.ocn_0001.restart   cesm_2000.pop_0001.r.2004-01-04-00000.nc       RESTART_FMT=nc

@ inst = 1
while ($inst <= $num_instances)
   set n4 = `printf %04d $inst`
   set n2 = `printf %02d $inst`

   echo ''
   echo "Staging restarts for instance $inst of $num_instances"

   ${LINK} ${CAM_stagedir}/cami-mam3_0000-01-01_0.9x1.25_L30_c100618.nc      cam_initial_${n4}.nc
   ${LINK} ${POP_stagedir}/b40.20th.005_ens${n2}.pop.r.2004-01-01-00000      .
   ${LINK} ${POP_stagedir}/b40.20th.005_ens${n2}.pop.r.2004-01-01-00000.hdr  .
   ${LINK} ${POP_stagedir}/b40.20th.005_ens${n2}.pop.ro.2004-01-01-00000     .
   ${LINK} ${POP_stagedir}/b40.20th.005_ens${n2}.cice.r.2004-01-01-00000.nc  .

#  ${LINK} ${RTM_stagedir}/b40.20th.005_ens${n2}.clm2.r.2004-01-01-00000.nc  .
#  ${LINK} ${CLM_stagedir}.clm2_$instance.r.${run_refdate}-${run_reftod}.nc  .

   echo "cam_initial_${n4}.nc"                                         >! rpointer.atm_${n4}
#  echo "${stagedir}.clm2_$instance.r.${run_refdate}-${run_reftod}.nc" >! rpointer.lnd_${n4}
   echo "b40.20th.005_ens${n2}.cice.r.2004-01-01-00000.nc"             >! rpointer.ice_${n4}
#  echo "b40.20th.005_ens${n2}.clm2.r.2004-01-01-00000.nc"             >! rpointer.rof_${n4}

   echo "b40.20th.005_ens${n2}.pop.ro.2004-01-01-00000"                >! rpointer.ocn_${n4}.ovf
   echo "b40.20th.005_ens${n2}.pop.r.2004-01-01-00000"                 >! rpointer.ocn_${n4}.restart
   echo "RESTART_FMT=bin"                                              >> rpointer.ocn_${n4}.restart

   @ inst ++
end

# ==============================================================================
# Edit the run script to reflect project, queue, and wallclock
# ==============================================================================

cd ${caseroot}

echo ''
echo 'Updating the run script to set wallclock and queue.'
echo ''

if ( ~ -e  ${case}.run.original ) then
   ${COPY} ${case}.run ${case}.run.original
endif

source Tools/ccsm_getenv
set BATCH = `echo $BATCHSUBMIT | sed 's/ .*$//'`
switch ( $BATCH )
   case bsub*:
      # NCAR "bluefire", "yellowstone"
      set TIMEWALL=`grep BSUB ${case}.run | grep -e '-W' `
      set    QUEUE=`grep BSUB ${case}.run | grep -e '-q' `
      sed -e "s/ptile=[0-9][0-9]*/ptile=$ptile/" \
          -e "s/$TIMEWALL[3]/$timewall/" \
          -e "s/$QUEUE[3]/$queue/" < ${case}.run >! temp.$$
      ${MOVE} temp.$$  ${case}.run
   breaksw

   default:

   breaksw
endsw

# ==============================================================================
# The *.run script must be modified to call the DART assimilate script.
# The modifications are contained in a "here" document that MUST NOT
# expand the wildcards etc., before it is run. This is achieved by
# double-quoting the characters used to delineate the start/stop of
# the "here" document. No kidding. It has to be "EndOfText",
# not 'EndOfText' or EndOfText.
# ==============================================================================

echo ''
echo 'Adding the call to assimilate.csh to the *.run script.'
echo ''

cat << "EndOfText" >! add_to_run.txt

# -------------------------------------------------------------------------
# START OF DART: if CESM finishes correctly (pirated from ccsm_postrun.csh);
# perform an assimilation with DART.

set CplLogFile = `ls -1t cpl.log* | head -n 1`
if ($CplLogFile == "") then
   echo 'ERROR: Model did not complete - no cpl.log file present - exiting.'
   echo 'ERROR: Assimilation will not be attempted.'
   setenv LSB_PJL_TASK_GEOMETRY "{(0)}"
   setenv EXITCODE -1
   mpirun.lsf ${CASEROOT}/shell_exit.sh
   exit -1
endif

grep 'SUCCESSFUL TERMINATION' $CplLogFile
if ( $status == 0 ) then
   ${CASEROOT}/assimilate.csh

   if ( $status == 0 ) then
      echo "`date` -- DART HAS FINISHED"
   else
      echo "`date` -- DART FILTER ERROR - ABANDON HOPE"
      setenv LSB_PJL_TASK_GEOMETRY "{(0)}"
      setenv EXITCODE -3
      mpirun.lsf ${CASEROOT}/shell_exit.sh
      exit -3
   endif
else
   echo 'ERROR: Model did not complete successfully - exiting.'
   echo 'ERROR: Assimilation will not be attempted.'
   setenv LSB_PJL_TASK_GEOMETRY "{(0)}"
   setenv EXITCODE -2
   mpirun.lsf ${CASEROOT}/shell_exit.sh
   exit -2
endif

# END OF DART BLOCK
# -------------------------------------------------------------------------
"EndOfText"

# Now that the "here" document is created,
# determine WHERE to insert it -- ONLY IF it is not already there.

grep "ABANDON HOPE" ${case}.run
set STATUSCHECK = $status

if ( ${STATUSCHECK} == 0 ) then
   echo "DART block already present in ${case}.run"
else if ( ${STATUSCHECK} == 1 ) then

   set MYSTRING = `grep --line-number "CSM EXECUTION HAS FINISHED" ${case}.run`
   set MYSTRING = `echo $MYSTRING | sed -e "s#:# #g"`

   @ origlen = `cat ${case}.run | wc -l`
   @ keep = $MYSTRING[1]
   @ lastlines = $origlen - $keep

   head -n $keep      ${case}.run    >! temp.$$
   cat                add_to_run.txt >> temp.$$
   tail -n $lastlines ${case}.run    >> temp.$$

   ${MOVE} temp.$$ ${case}.run
   ${REMOVE} add_to_run.txt

else
   echo "ERROR in grep of ${case}.run: aborting"
   echo "status was ${STATUSCHECK}"
   exit -6
endif

chmod 0744 ${case}.run

# ==============================================================================
# Stage the required parts of DART in the CASEROOT directory.
# ==============================================================================

# The standard CESM short-term archiving script may need to be altered
# to archive addtional or subsets of things, or to reduce the amount of
# data that is sent to the long-term archive.  Put a version of st_archive.sh
# in  ${DARTroot}/models/CESM/shell_scripts when/if necessary

if ( ~ -e  Tools/st_archive.sh.original ) then
   ${COPY} Tools/st_archive.sh Tools/st_archive.sh.original
else
   echo "a Tools/st_archive.sh backup copy already exists"
endif

${COPY} ${DARTroot}/models/CESM/shell_scripts/st_archive.sh           Tools/
${COPY} ${DARTroot}/models/CESM/shell_scripts/assimilate.csh          .
${COPY} ${DARTroot}/shell_scripts/shell_exit.sh                       .
${COPY} ${DARTroot}/models/CESM/shell_scripts/cam_assimilate.csh      .
${COPY} ${DARTroot}/models/CESM/shell_scripts/pop_assimilate.csh      .
${COPY} ${DARTroot}/models/CESM/shell_scripts/clm_assimilate.csh      .


# ==============================================================================
# Construct a set of scripts for restarting broken runs without having to
# run this entire script again and rebuild CESM.  This section creates 5 scripts
# which do the following tasks:
#  1) resetting the xml files to run (or rerun) the first step of the experiment
#  2) the xml changes you need to make between steps 1 and 2
#  3) restage initial case files to start over
#  4) restore the files from the last successful cesm advance to restart
#      in the middle of a run
#  5) restage the dart executables from the dartroot directory to the
#      CESM bld directory
# ==============================================================================


#-----

cat << EndOfText >! xml_changes_for_step1.sh
#!/bin/sh

# this script changes the env_run options that are needed for
# the first job step.
# this script was autogenerated by $0
# using the variables set in that script

./xmlchange STOP_OPTION=$stop_option
./xmlchange STOP_N=$stop_n
./xmlchange CONTINUE_RUN=FALSE
./xmlchange RESUBMIT=0
exit 0

EndOfText
chmod 0775 xml_changes_for_step1.sh

#-----

cat << EndOfText >! xml_changes_for_stepN.sh
#!/bin/sh

# this script changes the env_run options that are needed for
# any jobs step after the first one.
# this script was autogenerated by $0
# using the variables set in that script

./xmlchange STOP_OPTION=$stop_option
./xmlchange STOP_N=$assim_n
./xmlchange CONTINUE_RUN=TRUE
./xmlchange RESUBMIT=$resubmit
exit 0

EndOfText
chmod 0775 xml_changes_for_stepN.sh

#-----

cat << EndOfText >! reset_step1.sh
#!/bin/sh

# this script removes the current contents of the run directory and
# replaces the initial staged files needed to start the experiment over.
# this script was autogenerated by $0
# using the variables set in that script

# before removing everything, be sure we make it to the run dir
cd ${rundir}
if [[ "\`pwd\`" != ${rundir} ]]; then
   echo did not change to run directory successfully.  exiting.
   exit -1
fi
${REMOVE} ${rundir}/*

let inst=1
while [[ \$inst -le $num_instances ]]
do
   # the instance numbers
   echo staging files for instance \$inst

   n4=\`printf %04d \$inst\`
   n2=\`printf %02d \$inst\`

   ${LINK} ${CAM_stagedir}/cami-mam3_0000-01-01_0.9x1.25_L30_c100618.nc       cam_initial_\${n4}.nc
   ${LINK} ${POP_stagedir}/b40.20th.005_ens\${n2}.pop.r.2004-01-01-00000      .
   ${LINK} ${POP_stagedir}/b40.20th.005_ens\${n2}.pop.r.2004-01-01-00000.hdr  .
   ${LINK} ${POP_stagedir}/b40.20th.005_ens\${n2}.pop.ro.2004-01-01-00000     .
   ${LINK} ${POP_stagedir}/b40.20th.005_ens\${n2}.cice.r.2004-01-01-00000.nc  .

   echo "cam_initial_\${n4}.nc"                                > rpointer.atm_\${n4}
   echo "b40.20th.005_ens\${n2}.cice.r.2004-01-01-00000.nc"    > rpointer.ice_\${n4}

   echo "b40.20th.005_ens\${n2}.pop.ro.2004-01-01-00000"       > rpointer.ocn_\${n4}.ovf
   echo "b40.20th.005_ens\${n2}.pop.r.2004-01-01-00000"        > rpointer.ocn_\${n4}.restart
   echo "RESTART_FMT=bin"                                     >> rpointer.ocn_\${n4}.restart

   let inst=inst+1
done

date > ${rundir}/make_cam_inflation_cookie

cd ${caseroot}

# reset the env_run options to start a new run (or start over)
./xml_changes_for_step1.sh

exit 0

EndOfText
chmod 0775 reset_step1.sh

#-----

cat << EndOfText >! reset_last_successful_step.sh
#!/bin/sh

# this script removes the current contents of the run directory and
# restores the files from the last successfully archived directory.
# this script was autogenerated by $0
# using the variables set in that script

lastarchivedir=\`ls -1dt ${archdir}/.sta2/* | head -n 1\`
if [[ ! -d \$lastarchivedir ]]; then
  lastarchivedir=\`ls -1dt ${archdir}/rest/* | head -n 1\`
  if [[ ! -d \$lastarchivedir ]]; then
    echo cannot find last archive directory in ${archdir}/.sta2 
    echo or in ${archdir}/rest.  exiting.
    exit -1
  fi
fi

# before removing everything, be sure we make it to the run dir
cd ${rundir}
if [[ "\`pwd\`" != ${rundir} ]]; then
   echo did not change to run directory successfully.  exiting.
   exit -1
fi
${REMOVE} ${rundir}/*

${COPY} \${lastarchivedir}/*  .

let inst=1
while [[ \$inst -le $num_instances ]]
do
   # instance string includes the leading underscore
   inst_string=\`printf _%04d \$inst\`
   timetag=\`basename \${lastarchivedir}\`

   ${LINK} ${case}.cam\${inst_string}.i.\${timetag}.nc cam_initial\${inst_string}.nc

   let inst=inst+1
done
exit 0

# reset the env_run options for a continue run
./xml_changes_for_stepN.sh

# you may want to reset the RESUBMIT value, depending on how many
# runs worked before failing.
# ./xmlchange RESUBMIT=$resubmit

EndOfText
chmod 0775 reset_last_successful_step.sh

#-----

# ==============================================================================
# Stage the DART executables in the CESM execution root directory: EXEROOT
# ==============================================================================

cat << EndOfText >! refresh_dart_files.sh
#!/bin/sh

# this script copies over the dart executables and namelists to the
# proper directory.  if you have to update any dart code or namelists,
# do it in the $DARTroot directory and then rerun refresh_dart_files.sh
# this script was autogenerated by $0
# using the variables set in that script

${COPY} ${DARTroot}/models/cam/work/cam_to_dart   ${exeroot}/.
${COPY} ${DARTroot}/models/cam/work/dart_to_cam   ${exeroot}/.
${COPY} ${DARTroot}/models/cam/work/filter        ${exeroot}/filter_cam
${COPY} ${DARTroot}/models/cam/work/input.nml                cam_input.nml

${COPY} ${DARTroot}/models/clm/work/clm_to_dart   ${exeroot}/.
${COPY} ${DARTroot}/models/clm/work/dart_to_clm   ${exeroot}/.
${COPY} ${DARTroot}/models/clm/work/filter        ${exeroot}/filter_clm
${COPY} ${DARTroot}/models/clm/work/input.nml                clm_input.nml

${COPY} ${DARTroot}/models/POP/work/pop_to_dart   ${exeroot}/.
${COPY} ${DARTroot}/models/POP/work/dart_to_pop   ${exeroot}/.
${COPY} ${DARTroot}/models/POP/work/filter        ${exeroot}/filter_pop
${COPY} ${DARTroot}/models/POP/work/input.nml                pop_input.nml

${COPY} ${DARTroot}/models/CESM/work/cesm_to_dart ${exeroot}/.
${COPY} ${DARTroot}/models/CESM/work/dart_to_cesm ${exeroot}/.
${COPY} ${DARTroot}/models/CESM/work/filter       ${exeroot}/filter_cesm
${COPY} ${DARTroot}/models/CESM/work/input.nml               input.nml

if [[ -x update_namelists.sh ]]; then
  ./update_namelists.sh
fi

exit 0

EndOfText
chmod 0775 refresh_dart_files.sh

./refresh_dart_files.sh


# ==============================================================================
# fix up the namelists to be sure they are consistent with the
# ensemble size, suggest settings for num members in the output
# diagnostics files, etc.  the user is free to update these after
# setup and before running.
# ==============================================================================

cat << EndOfText >! update_namelists.sh
#!/bin/sh

# this script makes certain namelist settings consistent with the number
# of ensemble members built by the setup script.
# this script was autogenerated by $0
# using the variables set in that script

# Ensure that the input.nml ensemble size matches the number of instances.
# WARNING: the output files contain ALL ensemble members ==> BIG

ex cam_input.nml <<ex_end
g;ens_size ;s;= .*;= $num_instances;
g;num_output_state_members ;s;= .*;= $num_instances;
g;num_output_obs_members ;s;= .*;= $num_instances;
wq
ex_end

ex clm_input.nml <<ex_end
g;ens_size ;s;= .*;= $num_instances;
g;num_output_state_members ;s;= .*;= $num_instances;
g;num_output_obs_members ;s;= .*;= $num_instances;
g;casename ;s;= .*;= "../$case",;
wq
ex_end

# num_output_state_members intentionally not set for POP.
ex pop_input.nml <<ex_end
g;ens_size ;s;= .*;= $num_instances;
g;num_output_obs_members ;s;= .*;= $num_instances;
wq
ex_end

ex input.nml <<ex_end
g;ens_size ;s;= .*;= $num_instances;
g;num_output_state_members ;s;= .*;= $num_instances;
g;num_output_obs_members ;s;= .*;= $num_instances;
wq
ex_end

exit 0

EndOfText
chmod 0775 update_namelists.sh

./update_namelists.sh

#=========================================================================
# Stage the files needed for SAMPLING ERROR CORRECTION
#
# The sampling error correction is a lookup table.
# Each ensemble size has its own (static) file.
# It is only needed if any
# input.nml:&assim_tools_nml:sampling_error_correction = .true.,
#=========================================================================

set nmls = "cam_input.nml pop_input.nml clm_input.nml input.nml"
foreach N ( $nmls )
  set  MYSTRING = `grep sampling_error_correction $N`
  set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
  set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
  set SECSTRING = `echo $MYSTRING[2] | tr '[:upper:]' '[:lower:]'`
  
  if ( $SECSTRING == true ) then
     set SAMP_ERR_FILE = ${DARTroot}/system_simulation/final_full_precomputed_tables/final_full.${ensemble_size}
     if (  -e   ${SAMP_ERR_FILE} ) then
      ${COPY} ${SAMP_ERR_FILE} .
      break   # we only need to copy it once if anyone has SEC on.
     else
        echo "ERROR: no sampling error correction file for this ensemble size."
        echo "ERROR: looking for ${SAMP_ERR_FILE} in"
        echo "ERROR: ${DARTroot}/system_simulation/final_full_precomputed_tables"
        echo "ERROR: one can be generated for any ensemble size; see docs"
        exit -3
     endif
  endif
end
  

# ==============================================================================
# Initial setup for the default inflation scenario.
# ==============================================================================
# CAM usually uses adaptive state-space prior inflation. The initial settings
# are in the filter_nml and ... during an assimilation experiment, the output
# from one assimilation is the input for the next. To facilitate this operationally,
# it is useful to specify an initial file of inflation values for the first
# assimilation step. However, I can think of no general way to do this. The
# utility that creates the initial inflation values (fill_inflation_restart)
# needs the model size from model_mod. To get that, CAM needs a 'cam_phis.nc'
# file which we generally don't have at this stage of the game (it exists after
# a model advance). So ... until I think of something better ... I am making a
# cookie file that indicates this is the very first assimilation. If this
# cookie file exists, the assimilate.csh script will make the inflation restart
# file before it performs the assimilation. After the first assimilation takes
# place, the cookie file must be 'eaten' so that subsequent assimilations do not
# overwrite whatever _should_ be there.

date >! ${rundir}/make_cam_inflation_cookie

# ==============================================================================
# What to do next
# ==============================================================================

echo ""
echo "Time to check the case."
echo ""
echo "cd into ${caseroot}"
echo "1) edit ${caseroot}/Buildconf/clm.buildnml.csh ... remove 'hybrid' portion of line 86"
echo "2) edit ${caseroot}/Buildconf/rtm.buildnml.csh ... comment out line 36"
echo "3) If the ${case}.run script still contains:"
echo '     #BSUB -R "select[scratch_ok > 0]" '
echo "   around line 9, delete it."
echo ""
echo "Check the streams listed in the streams text files.  If more or different"
echo 'dates need to be added, then do this in the $CASEROOT/user_*files*'
echo "then invoke 'preview_namelists' so you can check the information in the"
echo "CaseDocs or ${rundir} directories."
echo ""
echo "Modify what you like in xxx_input.nml, make sure the observation directory"
echo "names set in xxx_assimilate.csh match those on your system, and submit"
echo "the CESM job by running:"
echo "./${case}.submit"
echo ""
echo "For continued submissions after the initial (hybrid) startup,"
echo " run the xml_changes_for_stepN.sh script, which makes the following"
echo " changes to the env_run variables:"
echo ""
echo "  ./xmlchange CONTINUE_RUN=TRUE"
echo "  ./xmlchange RESUBMIT=<number_of_cycles_to_run>"
echo "  ./xmlchange STOP_N=$assim_n"
echo ""
echo "If you have to start over, 'reset_step1.sh' will reset the original files"
echo " for the first timestep. 'reset_last_successful_step.sh' will set up the"
echo " files from the last fully successful model advance/assimilation."
echo " If you need to recompile any part of the DART system, 'refresh_dart_files.sh'"
echo " will copy them into the correct places."
echo ""

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

