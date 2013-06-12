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
# that has CLM running under a set of data atmospheres (from a DART/CAM assim).
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
# -- Either edit and run this script in the $DART/models/clm/shell_scripts
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
#       CLM initial ensemble
#       ...
# -- Check to make sure the stream txt files point to the proper stream files for
#       each ensemble member.
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

setenv case                 clm_test
setenv compset              I_2000_CN
setenv resolution           f19_f19
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
#                    Executables, scripts and input in $DARTroot/models/clm/...
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
# The initial ensemble can be set by specifying the 'finidat' variable in the
# user_nl_clm_${instance}. A FULL pathname to the file is required. This is nice
# for two reasons - one is that you don't need to copy the files and rename them
# (tedious), the second is that the full pathname provides a means of tracking 
# the origin of the initial ensemble.
# ==============================================================================

set stagedir = /glade/scratch/afox/bptmp/MD_40_PME/run/MD_40_PME

# ==============================================================================
# configure settings ... run_startdate format is yyyy-mm-dd
# ==============================================================================

setenv stream_year_first 2000
setenv stream_year_last  2000
setenv stream_year_align 2000
setenv refyear     2000
setenv refmon      01
setenv refday      31
setenv run_reftod  00000
setenv run_refdate $refyear-$refmon-$refday

# ==============================================================================
# runtime settings --  How many assimilation steps will be done after this one
#                      plus archiving options
#
# resubmit      How many job steps to run on continue runs (will be 0 initially)
# stop_option   Units for determining the forecast length between assimilations
# stop_n        Number of time units in the model_advance
# assim_n       Number of time units between assimilations
#
# ==============================================================================

setenv resubmit            0
setenv stop_option         nhours
setenv stop_n              24
setenv assim_n             24
setenv short_term_archiver on
setenv long_term_archiver  off

# ==============================================================================
# job settings
#
# timewall   can be changed during a series by changing the ${case}.run
# queue      can be changed during a series by changing the ${case}.run
#
# TJH How many f19_f19 CLM instances can fit on 1 'regular' node?
# ==============================================================================

setenv ACCOUNT      P8685nnnn
setenv timewall     0:20
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

# make sure these variables point to valid directories
foreach VAR ( cesm_datadir cesmroot DARTroot )
   # VAR is the shell variable name, DIR is the value
   set DIR = `eval echo \${$VAR}`
   if ( ! -d $DIR ) then
      echo "ERROR: directory '$DIR' not found"
      echo " In the setup script check the setting of: $VAR "
      exit -1
   endif
end

# make sure there is a filter in these dirs
foreach MODEL ( clm )
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
   # A more efficient layout must be done when running a full assimilation.

   @ cpl_pes  = $ptile
   @ atm_pes  = $ptile * $num_instances * 2
   @ glc_pes  = $ptile
   @ rof_pes  = $ptile * $num_instances * 2
   @ lnd_pes  = $ptile * $num_instances * 2
   @ ice_pes  = $ptile
   @ ocn_pes  = $ptile

else

   # This is only for the purpose of debugging the code.
   # Every component runs sequentially on all pes.
   # A more efficient layout must be found

   @ cpl_pes = $ptile * $num_instances
   @ atm_pes = $ptile * $num_instances
   @ glc_pes = $ptile * $num_instances
   @ rof_pes = $ptile * $num_instances
   @ lnd_pes = $ptile * $num_instances
   @ ice_pes = $ptile * $num_instances
   @ ocn_pes = $ptile * $num_instances 
   
endif

echo "task layout"
echo ""
echo "CPL gets $cpl_pes"
echo "ATM gets $atm_pes"
echo "GLC gets $glc_pes"
echo "ROF gets $rof_pes"
echo "LND gets $lnd_pes"
echo "ICE gets $ice_pes"
echo "OCN gets $ocn_pes"
echo ""

./xmlchange NTHRDS_CPL=1,NTASKS_CPL=$cpl_pes
./xmlchange NTHRDS_ATM=1,NTASKS_ATM=$atm_pes,NINST_ATM=$num_instances
./xmlchange NTHRDS_GLC=1,NTASKS_GLC=$glc_pes,NINST_GLC=1
./xmlchange NTHRDS_ROF=1,NTASKS_ROF=$rof_pes,NINST_ROF=$num_instances
./xmlchange NTHRDS_LND=1,NTASKS_LND=$lnd_pes,NINST_LND=$num_instances
./xmlchange NTHRDS_ICE=1,NTASKS_ICE=$ice_pes,NINST_ICE=1
./xmlchange NTHRDS_OCN=1,NTASKS_OCN=$ocn_pes,NINST_OCN=1
./xmlchange ROOTPE_CPL=0
./xmlchange ROOTPE_ATM=0
./xmlchange ROOTPE_GLC=0
./xmlchange ROOTPE_ROF=0
./xmlchange ROOTPE_LND=0
./xmlchange ROOTPE_ICE=0
./xmlchange ROOTPE_OCN=0

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

./xmlchange RUN_TYPE=startup
./xmlchange RUN_STARTDATE=$run_refdate
./xmlchange RUN_REFCASE=${case}
./xmlchange RUN_REFDATE=$run_refdate
./xmlchange RUN_REFTOD=$run_reftod
./xmlchange BRNCH_RETAIN_CASENAME=TRUE
./xmlchange GET_REFCASE=FALSE
./xmlchange EXEROOT=${exeroot}
./xmlchange PIO_TYPENAME=pnetcdf


./xmlchange DATM_MODE=CPLHIST3HrWx
./xmlchange DATM_CPLHIST_CASE=$case
./xmlchange DATM_CPLHIST_YR_ALIGN=$refyear
./xmlchange DATM_CPLHIST_YR_START=$refyear
./xmlchange DATM_CPLHIST_YR_END=$refyear

# Do not change the CALENDAR or the CONTINUE_RUN

./xmlchange CALENDAR=GREGORIAN
./xmlchange STOP_OPTION=$stop_option
./xmlchange STOP_N=$stop_n
./xmlchange CONTINUE_RUN=FALSE
./xmlchange RESUBMIT=$resubmit

# The river transport model ON is useful only when using an active ocean or
# land surface diagnostics. Setting ROF_GRID to 'null' turns off the RTM.

./xmlchange ROF_GRID='null'
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

# ==============================================================================
# Create namelist template: user_nl_clm
# Example user_nl_clm namelist adding and removing fields on primary history file
# hist_fincl1 = 'COSZEN', 'DECL'
# hist_fexcl1 = 'TG', 'TV', 'TSOI', 'H2OSOI'
# DART needs the lon,lat,levgrnd,lonatm,latatm,lonrof,latrof DIMENSION
# information from the .h0. history file - nothing else.
#
# hist_empty_htapes = .true.     suppresses the creation of all history files
# hist_fincl1 = 'TG',            except the first one, which will have one variable
# hist_nhtfrq = -$stop_n,        create one every $stop_n HOURS
# hist_mfilt  =  1,              with precisely one day in it
# hist_avgflag_pertape = 'I'     use instantaneous values - no average
#
# The fincl2 history tape has the half-hourly flux tower observations.
# The observation operators in obs_def_tower_mod.f90
# are going to read from the .h1. history file for these values.
# ==============================================================================

@ inst = 1
while ($inst <= $num_instances)

   set instance  = `printf %04d $inst`

   # ===========================================================================
   set fname = "user_nl_datm_$instance"
   # ===========================================================================

   echo "dtlimit = 1.5, 1.5, 1.5"                    >> $fname
   echo "fillalgo = 'nn', 'nn', 'nn'"                >> $fname
   echo "fillmask = 'nomask','nomask','nomask'"      >> $fname
   echo "mapalgo = 'bilinear','bilinear','bilinear'" >> $fname
   echo "mapmask = 'nomask','nomask','nomask'"       >> $fname
   echo "streams = 'datm.streams.txt.CPLHIST3HrWx.Solar_$instance             $stream_year_align $stream_year_first $stream_year_last'," >> $fname
   echo "          'datm.streams.txt.CPLHIST3HrWx.Precip_$instance            $stream_year_align $stream_year_first $stream_year_last'," >> $fname
   echo "          'datm.streams.txt.CPLHIST3HrWx.nonSolarNonPrecip_$instance $stream_year_align $stream_year_first $stream_year_last'"  >> $fname
   echo "taxmode = 'cycle','cycle','cycle'"          >> $fname
   echo "tintalgo = 'coszen','nearest','linear'"     >> $fname
   echo "restfils = 'unset'"                         >> $fname
   echo "restfilm = 'unset'"                         >> $fname

   # ===========================================================================
   set fname = "user_nl_clm_$instance"
   # ===========================================================================

   # Customize the land namelists
   # The initial ensemble can be set by specifying the 'finidat' variable in the
   # user_nl_clm_${instance}. A FULL pathname to the file is required. This is nice
   # for two reasons - one is that you don't need to copy the files and rename them
   # (tedious), the second is that the full pathname provides a means of tracking 
   # the origin of the initial ensemble.
 
   @ thirtymin = $assim_n * 2

   echo "finidat = '${stagedir}.clm2_$instance.r.${run_refdate}-${run_reftod}.nc'" >> $fname
   echo "hist_empty_htapes = .false."                >> $fname
   echo "hist_fincl1 = 'NEP'"                        >> $fname
   echo "hist_fincl2 = 'NEP','FSH','EFLX_LH_TOT_R'"  >> $fname
   echo "hist_nhtfrq = -$stop_n,1,"                  >> $fname
   echo "hist_mfilt  = 1,$thirtymin"                 >> $fname
   echo "hist_avgflag_pertape = 'A','A'"             >> $fname

   @ inst ++
end

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
      default:
         ${COPY} $FILE user_$FNAME
         chmod   644   user_$FNAME
         breaksw
   endsw

end

# Replace each default stream txt file with one that uses the CLM DATM
# conditions for a default year and modify the instance number.

foreach FNAME (user*streams*)
   set name_parse = `echo $FNAME | sed 's/\_/ /g'`
   @ instance_index = $#name_parse
   @ filename_index = $#name_parse - 1
   set streamname = $name_parse[$filename_index]
   set instance   = $name_parse[$instance_index]

   if (-e $DARTroot/models/clm/shell_scripts/user_$streamname*template) then

      echo "Copying DART template for $FNAME and changing instances, refyear"

      ${COPY} $DARTroot/models/clm/shell_scripts/user_$streamname*template $FNAME

      sed s/NINST/$instance/g   $FNAME >! out.$$
      sed s/REFYEAR/$refyear/g  out.$$ >! $FNAME
      \rm -f out.$$

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
# IF NEED BE : Stage the restarts now that the run directory exists
# The initial ensemble can be set by specifying the 'finidat' variable in the
# user_nl_clm_${instance}. A FULL pathname to the file is required. This is nice
# for two reasons - one is that you don't need to copy the files and rename them
# (tedious), the second is that the full pathname provides a means of tracking 
# the origin of the initial ensemble. If you choose to copy and rename them,
# here is some code that may help.  Use at your own risk.
# ==============================================================================

# set stagedir = /glade/scratch/afox/bptmp/MD_40_PME/run

# echo ''
# echo "Copying the restart files from ${stagedir}"
# echo ''
#
# @ inst = 1
# while ($inst <= $num_instances)
#
#    echo "Staging restarts for instance $inst of $num_instances"
#
#    set LANDFILE = `printf ${stagedir}/MD_40_PME.clm2_%04d.r.2000-01-31-00000.nc $inst`
#    set LND_RESTART_FILENAME = `printf "${case}.clm2_%04d.r.%04d-%02d-%02d-%05d.nc" $inst $refyear $refmon $refday $run_reftod`
#
#    ${COPY} ${LANDFILE} ${rundir}/${LND_RESTART_FILENAME}
#
#    @ inst ++
# end

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
   setenv LSB_PJL_TASK_GEOMETRY "{(0)}"
   setenv EXITCODE -1
   mpirun.lsf ${CASEROOT}/shell_exit.sh
   exit $EXITCODE
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
      exit $EXITCODE
   endif
else
   echo 'ERROR: Model did not complete successfully - exiting.'
   echo 'ERROR: Assimilation will not be attempted.'
   setenv LSB_PJL_TASK_GEOMETRY "{(0)}"
   setenv EXITCODE -2
   mpirun.lsf ${CASEROOT}/shell_exit.sh
   exit $EXITCODE
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
else
   echo "a Tools/st_archive.sh backup copy already exists"
endif

# NOTE: the assimilate.csh script and input.nml must be modified for your
#       situation. The script has variables that point to the location of
#       the observations sequence files and the DART working directory
#       and may be customized for a more efficient PE layout for DART.
#       If you are running this, you know what to do with input.nml.
#       If you don't, you should give up now. Really.

${COPY} ${DARTroot}/models/clm/shell_scripts/st_archive.sh   Tools/
${COPY} ${DARTroot}/models/clm/shell_scripts/assimilate.csh  .
${COPY} ${DARTroot}/models/clm/work/input.nml                .
${COPY} ${DARTroot}/shell_scripts/shell_exit.sh              .

# ==============================================================================
# Stage the DART executables in the CESM execution root directory: EXEROOT
# ==============================================================================

foreach FILE ( perfect_model_obs filter clm_to_dart dart_to_clm )
   ${COPY} ${DARTroot}/models/clm/work/${FILE} ${exeroot}/
   if ( $status != 0 ) then
      echo "ERROR: ${DARTroot}/models/clm/work/${FILE} not copied to ${exeroot}"
      echo "ERROR: ${DARTroot}/models/clm/work/${FILE} not copied to ${exeroot}"
      exit -7
   endif
end

# ==============================================================================
# What to do next
# ==============================================================================

echo ""
echo "Time to check the case."
echo ""
echo "cd into ${caseroot}"
echo ""
echo "Modify what you like in input.nml, make sure the observation directory"
echo "names set in assimilate.csh match those on your system, and submit"
echo "the CESM job by running:"
echo "./${case}.submit"
echo ""
echo "For continued submissions after the initial startup,"
echo "make the following changes to the env_run variables:"
echo ""
echo "  ./xmlchange -file env_run.xml -id CONTINUE_RUN  -val TRUE"
echo "  ./xmlchange -file env_run.xml -id RESUBMIT      -val <your_favorite_number>"
echo ""
echo "Check the streams listed in the streams text files.  If more or different"
echo 'dates need to be added, then do this in the $CASEROOT/user_*files*'
echo "then invoke 'preview_namelists' so you can check the information in the"
echo "CaseDocs or ${rundir} directories."
echo ""

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

