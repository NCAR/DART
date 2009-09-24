#!/usr/local/bin/tcsh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2009, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next few lines under version control, do not edit>
# $URL$
# $Id$
# $Revision$
# $Date$

#-----------------------------------------------------------------------------
# assimilate.csh
#-----------------------------------------------------------------------------

set ensemble_size = ${NINST_OCN}

# Create temporary working directory for the assimilation
set temp_dir = assimilate_dir.$$
echo "temp_dir is $temp_dir"

# Create a clean temporary directory and go there
\rm -rf  $temp_dir
mkdir -p $temp_dir
cd       $temp_dir

#-------------------------------------------------------------------------
# Determine time of model state ... from file name of first member
# of the form "./${CASE}.pop.$ensemble_member.r.2000-01-06-00000.nc"
#-------------------------------------------------------------------------

set FILE = `head -1 ../rpointer.ocn.1.restart`
set FILE = $FILE:t
set FILE = $FILE:r
set OCN_DATE_EXT = `echo $FILE:e`
set OCN_DATE_STR = `echo $FILE:e | sed -e "s#-# #g"`
set OCN_DATE = `echo $OCN_DATE_STR`
@ OCN_YEAR    = $OCN_DATE[1]
@ OCN_MONTH   = $OCN_DATE[2]
@ OCN_DAY     = $OCN_DATE[3]
@ OCN_SECONDS = $OCN_DATE[4]

echo "valid time of model is $OCN_YEAR $OCN_MONTH $OCN_DAY $OCN_SECONDS"

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
#-----------------------------------------------------------------------------

set DARTDIR = ${HOME}/DART/models/POP/work

set DART_OBS_DIR = `printf %04d%02d ${OCN_YEAR} ${OCN_MONTH}`
set  OBSDIR = /ptmp/dart/Obs_sets/GTSPP/${DART_OBS_DIR}

#-------------------------------------------------------------------------
# Populate a run-time directory with the bits needed to run DART 
# Either get them from the CCSM 'run' directory or some stock repository
# The grid files are absolute paths ... so they need not move.
#-------------------------------------------------------------------------

foreach FILE ( input.nml filter pop_to_dart dart_to_pop ) 

   if ( -e   ../${FILE} ) then
      /usr/local/bin/cp -pv ../${FILE} .
   else if ( -e ${DARTDIR}/${FILE} ) then
      /usr/local/bin/cp -pv    ${DARTDIR}/${FILE} .
   else
      echo "DART required file $FILE not found ... ERROR"
      stop
   endif

end

# move inflation info (if they exist) ...

#-------------------------------------------------------------------------
# Block 1: convert N POP restart files to DART initial conditions file(s)
# pop_to_dart is serial code, we can do all of these at the same time
# and just wait for them to finish IFF it were not for the fact we'd have
# to have unique namelists for all of them.
#
# At the end of the block, we have DART restart files  filter_ics.[1-N]
# that came from pointer files ../rpointer.ocn.[1-N].restart
#
# DART namelist settings appropriate/required:
# &filter_nml:           restart_in_file_name    = 'filter_ics'
# &ensemble_manager_nml: single_restart_file_in  = '.false.'
# &pop_to_dart_nml:      pop_to_dart_output_file = 'dart.ud',
#-------------------------------------------------------------------------

set member = 1
while ( $member <= $ensemble_size )

   set DART_IC_FILE = `printf filter_ics.%04d $member`
   set OCN_RESTART_FILENAME = `head -1 ../rpointer.ocn.$member.restart`
   ln -sf ../$OCN_RESTART_FILENAME pop.r.nc
   ln -sf ../pop2_in.$member       pop_in

   ./pop_to_dart || exit 1

   mv dart.ud $DART_IC_FILE
   @ member++
end

#-------------------------------------------------------------------------
# Block 2: Actually run the assimilation.
# Will result in a set of files : 'filter_restart.xxxx'
#
# DART namelist settings required:
# &filter_nml:           async                  = 0,
# &filter_nml:           adv_ens_command        = "./no_model_advance.csh",
# &filter_nml:           restart_in_file_name   = 'filter_ics'
# &filter_nml:           restart_out_file_name  = 'filter_restart'
# &filter_nml:           obs_sequence_in_name   = 'obs_seq.out'
# &filter_nml:           obs_sequence_out_name  = 'obs_seq.final'
# &filter_nml:           init_time_days         = -1,
# &filter_nml:           init_time_seconds      = -1,
# &filter_nml:           first_obs_days         = -1,
# &filter_nml:           first_obs_seconds      = -1,
# &filter_nml:           last_obs_days          = -1,
# &filter_nml:           last_obs_seconds       = -1,
# &ensemble_manager_nml: single_restart_file_in = '.false.'
#
#-------------------------------------------------------------------------

# Determine proper observation sequence file.

set OBSFNAME = `printf obs_seq.0Z.%04d%02d%02d ${OCN_YEAR} ${OCN_MONTH} ${OCN_DAY}`
# set OBSFNAME = obs_seq.0Z.20000106.1obs
set OBS_FILE = ${OBSDIR}/${OBSFNAME} 

ln -sfv ${OBS_FILE} obs_seq.out

mpirun.lsf ./filter || exit 2

mv Prior_Diag.nc      ../Prior_Diag.${OCN_DATE_EXT}.nc
mv Posterior_Diag.nc  ../Posterior_Diag.${OCN_DATE_EXT}.nc
mv obs_seq.final      ../obs_seq.${OCN_DATE_EXT}.final
mv dart_log.out       ../dart_log.${OCN_DATE_EXT}.out

#  move inflation info (if they exist) ...

#-------------------------------------------------------------------------
# Block 3: Update the POP restart files ... sequentially (sigh) ...
#
# DART namelist settings required:
# &filter_nml:           restart_out_file_name  = 'filter_restart'
# &ensemble_manager_nml: single_restart_file_in = '.false.'
# &dart_to_pop_nml:      dart_to_pop_input_file = 'dart.ic',
# &dart_to_pop_nml:      advance_time_present   = .false.
#-------------------------------------------------------------------------

set member = 1
while ( $member <= $ensemble_size )

   set DART_RESTART_FILE = `printf filter_restart.%04d $member`
   set OCN_RESTART_FILENAME = `head -1 ../rpointer.ocn.$member.restart`
   ln -sf ../$OCN_RESTART_FILENAME pop.r.nc
   ln -sf ../pop2_in.$member       pop_in

   ln -sf $DART_RESTART_FILE dart.ic

   ./dart_to_pop || exit 3

   @ member++
end

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

ls -lrt

