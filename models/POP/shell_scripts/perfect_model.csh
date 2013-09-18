#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.

echo "`date` -- BEGIN GENERATE TRUE STATE"

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

      set BASEOBSDIR = /glade/proj3/image/Observations/WOD09
   breaksw

   case ys*:
      # NCAR "yellowstone"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'

      set BASEOBSDIR = /glade/p/image/Observations/WOD09
   breaksw

   default:
      # NERSC "hopper"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'

      set BASEOBSDIR = /scratch/scratchdirs/nscollin/ACARS
   breaksw
endsw

#-------------------------------------------------------------------------
# Determine time of model state ... from file name of first member
# of the form "./${MYCASE}.pop.r.2000-01-06-00000.nc"
#
# Piping stuff through 'bc' strips off any preceeding zeros.
#-------------------------------------------------------------------------

set FILE = `head -n 1 rpointer.ocn.restart`
set FILE = $FILE:t
set FILE = $FILE:r
set MYCASE = `echo $FILE | sed -e "s#\..*##"`
set OCN_DATE_EXT = `echo $FILE:e`
set OCN_DATE     = `echo $FILE:e | sed -e "s#-# #g"`
set OCN_YEAR     = `echo $OCN_DATE[1] | bc`
set OCN_MONTH    = `echo $OCN_DATE[2] | bc`
set OCN_DAY      = `echo $OCN_DATE[3] | bc`
set OCN_SECONDS  = `echo $OCN_DATE[4] | bc`
set OCN_HOUR     = `echo $OCN_DATE[4] / 3600 | bc`

echo "valid time of model is $OCN_YEAR $OCN_MONTH $OCN_DAY $OCN_SECONDS (seconds)"
echo "valid time of model is $OCN_YEAR $OCN_MONTH $OCN_DAY $OCN_HOUR (hours)"

#-------------------------------------------------------------------------
# Create temporary working directory for the assimilation and go there
#-------------------------------------------------------------------------

set temp_dir = assimilate_pop
echo "temp_dir is $temp_dir"

if ( -d $temp_dir ) then
   ${REMOVE} $temp_dir/*
else
   mkdir -p $temp_dir
endif
cd $temp_dir

#-----------------------------------------------------------------------------
# Get observation sequence file ... or die right away.
# The observation file names have a time that matches the stopping time of POP.
#-----------------------------------------------------------------------------

set YYYYMM   = `printf %04d%02d                ${OCN_YEAR} ${OCN_MONTH}`
set OBSFNAME = `printf obs_seq.0Z.%04d%02d%02d ${OCN_YEAR} ${OCN_MONTH} ${OCN_DAY}`
set OBS_FILE = ${BASEOBSDIR}/${YYYYMM}/${OBSFNAME}

if (  -e   ${OBS_FILE} ) then
   ${COPY} ${OBS_FILE} obs_seq.in
else
   echo "ERROR ... no observation file $OBS_FILE"
   echo "ERROR ... no observation file $OBS_FILE"
   exit -1
endif

#=========================================================================
# Block 1: Populate a run-time directory with the input needed to run DART.
#=========================================================================

echo "`date` -- BEGIN COPY BLOCK"

if (  -e   ${CASEROOT}/input.nml ) then
   ${COPY} ${CASEROOT}/input.nml .
else
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   exit -2
endif

echo "`date` -- END COPY BLOCK"

# Eat the cookie regardless
${REMOVE} ../make_pop_inflation_cookie

#=========================================================================
# Block 2: convert 1 pop restart file to a DART initial conditions file.
# At the end of the block, we have a DART restart file  perfect_ics
# that came from the pointer file ../rpointer.ocn.restart
#
# DART namelist settings appropriate/required:
# &perfect_model_obs_nml:  restart_in_file_name    = 'perfect_ics'
# &pop_to_dart_nml:        pop_to_dart_output_file = 'dart_ics'
#=========================================================================

echo "`date` -- BEGIN POP-TO-DART"

set member = 1

   set MYTEMPDIR = member_${member}
   mkdir -p $MYTEMPDIR
   cd $MYTEMPDIR

   # make sure there are no old output logs hanging around
   $REMOVE output.${member}.pop_to_dart

   set OCN_RESTART_FILENAME = ${MYCASE}.pop.r.${OCN_DATE_EXT}.nc
   set     OCN_NML_FILENAME = pop2_in
   set     DART_IC_FILENAME = perfect_ics
   set    DART_RESTART_FILE = perfect_restart

   sed -e "s#dart_ics#../${DART_IC_FILENAME}#" \
       -e "s#dart_restart#../${DART_RESTART_FILE}#" < ../input.nml >! input.nml

   ${LINK} ../../$OCN_RESTART_FILENAME pop.r.nc
   ${LINK} ../../$OCN_NML_FILENAME     pop_in

   ${EXEROOT}/pop_to_dart >! output.${member}.pop_to_dart

   if ($status != 0) then
      echo "ERROR ... DART died in 'pop_to_dart' ... ERROR"
      echo "ERROR ... DART died in 'pop_to_dart' ... ERROR"
      exit -3
   endif

   cd ..

echo "`date` -- END POP-TO-DART"

#=========================================================================
# Block 3: Advance the model and harvest the synthetic observations.
# Will result in a single file : 'perfect_restart' which we don't need
# for a perfect model experiment with CESM.
#
# DART namelist settings required:
# &perfect_model_obs_nml:           async                  = 0,
# &perfect_model_obs_nml:           adv_ens_command        = "./no_model_advance.csh",
# &perfect_model_obs_nml:           restart_in_file_name   = 'perfect_ics'
# &perfect_model_obs_nml:           restart_out_file_name  = 'perfect_restart'
# &perfect_model_obs_nml:           obs_sequence_in_name   = 'obs_seq.in'
# &perfect_model_obs_nml:           obs_sequence_out_name  = 'obs_seq.out'
# &perfect_model_obs_nml:           init_time_days         = -1,
# &perfect_model_obs_nml:           init_time_seconds      = -1,
# &perfect_model_obs_nml:           first_obs_days         = -1,
# &perfect_model_obs_nml:           first_obs_seconds      = -1,
# &perfect_model_obs_nml:           last_obs_days          = -1,
# &perfect_model_obs_nml:           last_obs_seconds       = -1,
#
#=========================================================================

# POP always needs a pop_in and a pop.r.nc to start.

set OCN_RESTART_FILENAME = ${MYCASE}.pop.r.${OCN_DATE_EXT}.nc

${LINK} ../$OCN_RESTART_FILENAME pop.r.nc
${LINK} ../pop2_in               pop_in

echo "`date` -- BEGIN PERFECT_MODEL_OBS"
${EXEROOT}/perfect_model_obs || exit -4
echo "`date` -- END PERFECT_MODEL_OBS"

${MOVE} True_State.nc    ../pop_True_State.${OCN_DATE_EXT}.nc
${MOVE} obs_seq.out      ../obs_seq.${OCN_DATE_EXT}.out
${MOVE} dart_log.out     ../pop_dart_log.${OCN_DATE_EXT}.out

#=========================================================================
# Block 4: Update the pop restart files.
#=========================================================================

# not needed ... perfect_model_obs does not update the model state.

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

echo "`date` -- END   GENERATE TRUE STATE"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

