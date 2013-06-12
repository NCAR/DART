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
# that has POP as active components
# and will use DART to assimilate observations at regular intervals.
# This script does not build DART. It works best if the appropriate DART
# executables have been built, however.
#
# This script relies heavily on the information in:
# http://www.cesm.ucar.edu/models/cesm1.1/cesm/doc/usersguide/book1.html
#
# ---------------------
# How to set up the script
# ---------------------
# -- Either edit and run this script in the $DART/models/POP/shell_scripts
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
#       POP initial ensemble
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

# ==============================================================================
# ====  Set case options
# ==============================================================================

# case will be used many ways;
#    directory and file names, both locally and on HPSS, and
#    script names; so consider it's length and information content.
# num_instances:  Number of ensemble members

setenv case                 pop_test
setenv compset              GIAF
setenv resolution           T62_gx1v6
setenv cesmtag              cesm1_1_1
setenv num_instances        4

# ==============================================================================
# define machines and directories
#
# mach            Computer name
# cesm_datadir    Root path of the public CESM data files
# cesmroot        Location of the cesm code base
#                    i.e. cesm1_1_1 on yellowstone
# DARTroot        Location of DART code tree.
#                    Executables, scripts and input in $DARTroot/models/POP/...
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

setenv DARTroot     /glade/u/home/${USER}/svn/DART/dev

# ==============================================================================
# configure settings ... run_startdate format is yyyy-mm-dd
# ==============================================================================

setenv stream_year_first 2004
setenv stream_year_last  2004
setenv stream_year_align 2004
setenv refyear     2004
setenv refmon      01
setenv refday      01
setenv run_reftod  00000
setenv run_refdate $refyear-$refmon-$refday

# ==============================================================================
# runtime settings --  How many assimilation steps will be done after this one
#
# stop_option   Units for determining the forecast length between assimilations
# stop_n        Number of time units in the forecast
# ==============================================================================

setenv resubmit      0
setenv stop_option   ndays
setenv stop_n        2

# ==============================================================================
# job settings
#
# timewall   can be changed during a series by changing the ${case}.run
# queue      can be changed during a series by changing the ${case}.run
#
# TJH: How many T62_gx1v6 POP instances can fit on 1 node?
# ==============================================================================

setenv ACCOUNT      P8685nnnn
setenv timewall     0:20
setenv queue        regular
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

# Try to lay out the tasks in some sort of reasonable way.
# The ocean root pe is the only one that is not zero.

if ($num_instances == 4) then

   # This is only for the purpose of debugging the code.
   # A more efficient layout must be done when running a full assimilation.

   @ cpl_pes = $ptile
   @ atm_pes = $ptile * $num_instances
   @ ice_pes = $ptile * $num_instances
   @ lnd_pes = $ptile
   @ rof_pes = $ptile
   @ glc_pes = $ptile
   @ ocn_pes = $ptile * $num_instances

else

   # layout for 30 members.
   # Made in conjunction with Mike Levy, David Bailey and Jim Edwards.

   @ cpl_pes = 450
   @ atm_pes = 450
   @ ice_pes = 450
   @ lnd_pes = 450
   @ rof_pes = 450
   @ glc_pes = 450
   @ ocn_pes = 1800

endif

echo "task partitioning ... ice and ocean run at the same time."
echo ""
echo "CPL gets $cpl_pes"
echo "ATM gets $atm_pes"
echo "ICE gets $ice_pes"
echo "LND gets $lnd_pes"
echo "ROF gets $rof_pes"
echo "GLC gets $glc_pes"
echo "OCN gets $ocn_pes"
echo ""

./xmlchange ROOTPE_CPL=0,NTHRDS_CPL=1,NTASKS_CPL=$cpl_pes
./xmlchange ROOTPE_ATM=0,NTHRDS_ATM=1,NTASKS_ATM=$atm_pes,NINST_ATM=$num_instances
./xmlchange ROOTPE_ICE=0,NTHRDS_ICE=1,NTASKS_ICE=$ice_pes,NINST_ICE=$num_instances
./xmlchange ROOTPE_LND=0,NTHRDS_LND=1,NTASKS_LND=$lnd_pes,NINST_LND=1
./xmlchange ROOTPE_ROF=0,NTHRDS_ROF=1,NTASKS_ROF=$rof_pes,NINST_ROF=1
./xmlchange ROOTPE_GLC=0,NTHRDS_GLC=1,NTASKS_GLC=$glc_pes,NINST_GLC=1
./xmlchange              NTHRDS_OCN=1,NTASKS_OCN=$ocn_pes,NINST_OCN=$num_instances
./xmlchange ROOTPE_OCN=$ice_pes

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
./xmlchange RUN_REFCASE=b40.20th.005
./xmlchange RUN_REFDATE=$run_refdate
./xmlchange RUN_REFTOD=$run_reftod
./xmlchange BRNCH_RETAIN_CASENAME=FALSE
./xmlchange GET_REFCASE=FALSE
./xmlchange EXEROOT=${exeroot}
./xmlchange PIO_TYPENAME=pnetcdf

# DOUT_S     is to turn on/off the short-term archiving
# DOUT_L_MS  is to store to the HPSS (formerly "MSS")

./xmlchange DOUT_S_ROOT=${archdir}
./xmlchange DOUT_S=TRUE
./xmlchange DOUT_S_SAVE_INT_REST_FILES=FALSE
./xmlchange DOUT_L_MS=FALSE
./xmlchange DOUT_L_MSROOT="csm/${case}"
./xmlchange DOUT_L_HTAR=FALSE

./xmlchange DATM_MODE=CPLHIST3HrWx
./xmlchange DATM_CPLHIST_CASE=$case
./xmlchange DATM_CPLHIST_YR_ALIGN=$refyear
./xmlchange DATM_CPLHIST_YR_START=$refyear
./xmlchange DATM_CPLHIST_YR_END=$refyear

# The streams files were generated with a NO_LEAP calendar in mind.
# We need to test these with a GREGORIAN calendar.

./xmlchange CALENDAR=GREGORIAN
./xmlchange STOP_OPTION=$stop_option
./xmlchange STOP_N=$stop_n
./xmlchange CONTINUE_RUN=FALSE
./xmlchange RESUBMIT=$resubmit

# The river transport model ON is useful only when using an active ocean or
# land surface diagnostics. Setting ROF_GRID to 'null' turns off the RTM.
# TJH - guidance needed from Alicia

# ./xmlchange ROF_GRID='null'
# ./xmlchange ROF_GRID='r05'

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

# ==============================================================================
# Modify namelist templates for each instance.
#
# hist_empty_htapes = .true.     suppresses the creation of all history files
# hist_fincl1 = 'TG',            except the first one, which will have one variable
# hist_nhtfrq = -$stop_n,        create one every $stop_n HOURS
# hist_mfilt  =  1,              with precisely one day in it
# hist_avgflag_pertape = 'I'     use instantaneous values - no average
#
# ==============================================================================

@ inst = 1
while ($inst <= $num_instances)

   set instance  = `printf %04d $inst`
   set instance2 = `printf %02d $inst`

   # ===========================================================================
   set fname = "user_nl_datm_$instance"
   # ===========================================================================

   echo "streams = 'datm.streams.txt.CPLHIST3HrWx.Solar_$instance             $stream_year_align $stream_year_first $stream_year_last'," >> $fname
   echo "          'datm.streams.txt.CPLHIST3HrWx.nonSolarNonPrecip_$instance $stream_year_align $stream_year_first $stream_year_last'"  >> $fname
   echo "dtlimit = 1.5, 1.5"               >> $fname
   echo "fillalgo = 'nn', 'nn'"            >> $fname
   echo "fillmask = 'nomask','nomask'"     >> $fname
   echo "mapalgo = 'bilinear','bilinear'"  >> $fname
   echo "mapmask = 'nomask','nomask'"      >> $fname
   echo "taxmode = 'cycle','cycle'"        >> $fname
   echo "tintalgo = 'linear','linear'"     >> $fname
   echo "restfils = 'unset'"               >> $fname
   echo "restfilm = 'unset'"               >> $fname

   # ===========================================================================
   set fname = "user_nl_cice_$instance"
   # ===========================================================================
   # CICE Namelists
   # this is only used for a hybrid start, else rpointers are used.

   echo "ice_ic = 'b40.20th.005_ens"$instance2".cice.r.2004-01-01-00000.nc'" >> $fname

   # ===========================================================================
   set fname = "user_nl_pop2_$instance"
   # ===========================================================================

   # POP Namelists
   # init_ts_suboption = 'data_assim'   for non bit-for-bit restarting (assimilation mode)
   # init_ts_suboption = 'null'         for 'perfect' restarting/forecasting

   echo "init_ts_suboption = 'data_assim'" >> $fname

   @ inst ++
end

# DLND

echo "streams = 'drof.streams.txt.rof.diatren_iaf_rx1" 1 1948 2009"'" >> user_nl_drof

# ==============================================================================
# to create custom streamfiles ...
# "To modify the contents of a stream txt file, first use preview_namelists to
#  obtain the contents of the stream txt files in CaseDocs, and then place a copy
#  of the modified stream txt file in $CASEROOT with the string user_ prepended."
#
# -or-
#
# we copy a template stream txt file from the
# $DARTroot/models/POP/shell_scripts directory and modify one for each instance.
#
# ==============================================================================

./preview_namelists

# This gives us a stream txt file for each instance that we can
# modify for our own purpose.

foreach FILE (CaseDocs/*streams*)
   set FNAME = $FILE:t

   switch ( ${FNAME} )
      case *presaero*:
         echo "Using default prescribed aerosol stream.txt file ${FNAME}"
         breaksw
      case *diatren*:
         echo "Using default runoff stream.txt file ${FNAME}"
         breaksw
      case *\.Precip_*:
         echo "Precipitation in nonSolarNonPrecip stream.txt file - not ${FNAME}"
         breaksw
      default:
         ${COPY} $FILE user_$FNAME
         chmod   644   user_$FNAME
         breaksw
   endsw

end

# Replace each default stream txt file with one that uses the CAM DATM
# conditions for a default year and modify the instance number.

foreach FNAME (user*streams*)
   set name_parse = `echo $FNAME | sed 's/\_/ /g'`
   @ filename_index = $#name_parse - 1
   @ instance_index = $#name_parse
   set streamname = $name_parse[$filename_index]
   set    n = `echo $name_parse[$instance_index] | sed 's/0*//'`

   if (-e $DARTroot/models/POP/shell_scripts/user_$streamname*template) then

      echo "Copying DART template for $FNAME and changing instances."

      ${COPY} $DARTroot/models/POP/shell_scripts/user_$streamname*template $FNAME

      sed s/NINST/$n/g $FNAME >! out.$$
      ${MOVE} out.$$ $FNAME

   else
      echo "DIED Looking for a DART stream txt template for $FNAME"
      echo "DIED Looking for a DART stream txt template for $FNAME"
      exit -3
   endif

end

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
   echo "DART requires modifications to several src.pop2/ files."
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
#
# set stagedir = /glade/scratch/thoar/DART_POP_RESTARTS/2004-01-01-00000
# ==============================================================================

set stagedir = /glade/p/work/aliciak/DART_IC/CCSM4_ensembles/rest/2004-01-01-00000

echo ''
echo "Copying the restart files from ${stagedir}"
echo 'into the CESM run directory.'
echo ''

@ inst = 1
while ($inst <= $num_instances)
   set n4 = `printf %04d $inst`
   set n2 = `printf %02d $inst`

   echo ''
   echo "Staging restarts for instance $inst of $num_instances"

   #AK: Note that the pop ocean must have the .hdr file to describe the
   #    metadata for binary restarts.
   ${COPY} ${stagedir}/rpointer.ocn_${n4}.restart                       ${rundir}
   ${COPY} ${stagedir}/rpointer.ocn_${n4}.ovf                           ${rundir}
   ${COPY} ${stagedir}/rpointer.ice_${n4}                               ${rundir}
   ${COPY} ${stagedir}/b40.20th.005_ens${n2}.pop.r.2004-01-01-00000     ${rundir}
   ${COPY} ${stagedir}/b40.20th.005_ens${n2}.pop.r.2004-01-01-00000.hdr ${rundir}
   ${COPY} ${stagedir}/b40.20th.005_ens${n2}.pop.ro.2004-01-01-00000    ${rundir}
   ${COPY} ${stagedir}/b40.20th.005_ens${n2}.cice.r.2004-01-01-00000.nc ${rundir}

   @ inst ++
end

if (  -e   ${stagedir}/prior_inflate_restart* ) then
   ${COPY} ${stagedir}/prior_inflate_restart* ${rundir}/.
endif

# ==============================================================================
# Edit the run script to reflect project, queue, and wallclock
# ==============================================================================

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
   exit -1
endif

grep 'SUCCESSFUL TERMINATION' $CplLogFile
if ( $status == 0 ) then
   ${CASEROOT}/assimilate.csh

   if ( $status == 0 ) then
      echo "`date` -- DART HAS FINISHED"
   else
      echo "`date` -- DART FILTER ERROR - ABANDON HOPE"
      exit -3
   endif
else
   echo 'ERROR: Model did not complete successfully - exiting.'
   echo 'ERROR: Assimilation will not be attempted.'
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

   ${MOVE}   temp.$$ ${case}.run
   ${REMOVE} add_to_run.txt

else
   echo "ERROR in grep of ${case}.run: aborting"
   echo "status was ${STATUSCHECK}"
   exit -6
endif

chmod 0755 ${case}.run

# ==============================================================================
# Stage the required parts of DART in the CASEROOT directory.
# ==============================================================================

# The standard CESM short-term archiving script may need to be altered
# to archive addtional or subsets of things, or to reduce the amount of
# data that is sent to the long-term archive.

if ( ~ -e  Tools/st_archive.sh.original ) then
   ${COPY} Tools/st_archive.sh Tools/st_archive.sh.original
endif

# NOTE: the assimilate.csh script and input.nml must be modified for your
#       situation. The script has variables that point to the location of
#       the observations sequence files and the DART working directory
#       and may be customized for a more efficient PE layout for DART.
#       If you are running this, you know what to do with input.nml.
#       If you don't, you should give up now. Really.

${COPY} ${DARTroot}/models/POP/shell_scripts/st_archive.sh   Tools/
${COPY} ${DARTroot}/models/POP/shell_scripts/assimilate.csh  assimilate.csh
${COPY} ${DARTroot}/models/POP/work/input.nml                input.nml

# ==============================================================================
# Stage the DART executables in the CESM execution root directory: EXEROOT
# ==============================================================================

foreach FILE ( filter pop_to_dart dart_to_pop )
   ${COPY} ${DARTroot}/models/POP/work/${FILE} ${exeroot}/
   if ( $status != 0 ) then
      echo "ERROR: ${DARTroot}/models/POP/work/${FILE} not copied to ${exeroot}"
      echo "ERROR: ${DARTroot}/models/POP/work/${FILE} not copied to ${exeroot}"
      exit -7
   endif
end

# ==============================================================================
# What to do next
# ==============================================================================

echo ''
echo "Time to check the case."
echo ''
echo "cd into ${caseroot}"
echo "Modify what you like in input.nml, make sure the observation directory"
echo "names set in assimilate.csh match those on your system, and submit"
echo "the CESM job by running:"
echo "./${case}.submit"
echo ''
echo "For continued submissions after the initial (hybrid) startup,"
echo "make the following changes to the env_run variables:"
echo ''
echo "  ./xmlchange -file env_run.xml -id STOP_N        -val 1"
echo "  ./xmlchange -file env_run.xml -id CONTINUE_RUN  -val TRUE"
echo "  ./xmlchange -file env_run.xml -id RESUBMIT      -val <your_favorite_number>"
echo ''
echo "Check the streams listed in the streams text files.  If more or different"
echo 'dates need to be added, then do this in the $CASEROOT/user_*files*'
echo "then invoke 'preview_namelists' so you can check the information in the"
echo "CaseDocs or ${rundir} directories."
echo ''

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

