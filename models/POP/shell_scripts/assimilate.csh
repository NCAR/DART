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

set   MOVE = '/usr/local/bin/mv -fv'
set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
set   LINK = '/usr/local/bin/ln -fvs'
set REMOVE = '/usr/local/bin/rm -fr'

set ensemble_size = ${NINST_OCN}

# Create temporary working directory for the assimilation
set temp_dir = assimilate_dir
echo "temp_dir is $temp_dir"

# Create a clean temporary directory and go there
if ( -d $temp_dir ) then
   ${REMOVE} $temp_dir/*
else
   mkdir -p $temp_dir
endif
cd $temp_dir

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
set  OBSDIR = /ptmp/dart/Obs_sets/WOD/${DART_OBS_DIR}

#-------------------------------------------------------------------------
# DART COPY BLOCK
# Populate a run-time directory with the bits needed to run DART 
# Either get them from the CCSM 'run' directory or some stock repository
# The grid files are absolute paths ... so they need not move.
#-------------------------------------------------------------------------

foreach FILE ( input.nml filter pop_to_dart dart_to_pop )

   if ( -e    ../${FILE} ) then
      ${COPY} ../${FILE} .
   else if ( -e ${DARTDIR}/${FILE} ) then
      ${COPY}   ${DARTDIR}/${FILE} .
   else
      echo "DART required file $FILE not found ... ERROR"
      stop
   endif

end

#-------------------------------------------------------------------------
# INFLATION COPY BLOCK
# This file is only relevant if 'inflation' is turned on -
# i.e. if inf_flavor(1) /= 0 - AND we are in a 'restart' mode.
#
# filter_nml
# inf_flavor                  = 2,                       0,
# inf_initial_from_restart    = .true.,                  .false.,
# inf_in_file_name            = 'prior_inflate_ics',     'post_inflate_ics',
#
# This is a 'test' configuration for this script. We are simply
# assuming that the namelist values are set such that we need this file,
# and that it is called 'prior_inflate_ics'. Since the inflation file is
# essentially a duplicate of the model state ... it is slaved to a specific 
# geometry. I created the file offline for the gx1v6 geometry on bluefire. 
# The inflation values are all unity.
#
# The strategy is to use an inflation file from CENTRALDIR if one exists - 
# if not, grab one from the same place we get the DART bits - 
# if not, grab the default one and hope for the best.
#
# The thought being: the first time through, CENTRALDIR doesn't have one
# so it will either be copied from the DART tree or use a vanilla value.
# After an assimilation, the output file will be copied back to CENTRALDIR
# to be used for subsequent assimilations.
#-------------------------------------------------------------------------

foreach FILE ( prior_inflate_ics ) 

   if ( -e    ../${FILE} ) then
      ${LINK} ../${FILE} .
   else if ( -e ${DARTDIR}/${FILE} ) then
      ${LINK}   ${DARTDIR}/${FILE} .
   else
      echo "DART auxiliary file $FILE not found ... starting from 1.0"
      set INFLATEFILE =  /ptmp/thoar/gx1v6_prior_inflate_restart.be
      ${COPY} ${INFLATEFILE} prior_inflate_ics
   endif

end

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
   ${LINK} ../$OCN_RESTART_FILENAME pop.r.nc
   ${LINK} ../pop2_in.$member       pop_in

   ./pop_to_dart || exit 1

   ${MOVE} dart.ud $DART_IC_FILE
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
set OBS_FILE = ${OBSDIR}/${OBSFNAME} 

${LINK} ${OBS_FILE} obs_seq.out

# FIXME: special for trying out non-monotonic task layouts.
setenv ORG_PATH "${PATH}"
setenv LSF_BINDIR /contrib/lsf/tgmpatch
setenv PATH ${LSF_BINDIR}:${PATH}
setenv ORG_TASK_GEOMETRY "${LSB_PJL_TASK_GEOMETRY}"
setenv NANCY_GEOMETRY_126 \
	"{(0,29,30,31,1,32,33,34,2,35,36,37,3,38,39,40,4,41,42,43,5)(44,45,46,6,47,48,49,7,50,51,52,8,53,54,55,9,56,57,58,10,59)(60,61,11,62,63,64,12,65,66,67,13,68,69,70,14,71,72,73,15,74,75)(76,16,77,78,79,17,80,81,82,18,83,84,85,19,86,87,88,20,89,90,91)(21,92,93,94,22,95,96,97,23,98,99,100,24,101,102,103,25,104,105,106,26)(107,108,109,27,110,111,112,28,113,114,115,116,117,118,119,120,121,122,123,124,125)}"
setenv LSB_PJL_TASK_GEOMETRY "${NANCY_GEOMETRY_126}"

which mpirun.lsf

mpirun.lsf ./filter || exit 2

${MOVE} Prior_Diag.nc      ../Prior_Diag.${OCN_DATE_EXT}.nc
${MOVE} Posterior_Diag.nc  ../Posterior_Diag.${OCN_DATE_EXT}.nc
${MOVE} obs_seq.final      ../obs_seq.${OCN_DATE_EXT}.final
${MOVE} dart_log.out       ../dart_log.${OCN_DATE_EXT}.out

# FIXME: should add something for posterior_inflate_restarts ...
if ( -e prior_inflate_restart ) then
   # 1) rename file to reflect current date
   # 2) must link generic name to specific date so next day
   # we can just grab the generic name without having to
   # parse the 'day before the current day'
   # 3) move both to CENTRALDIR so the DART INFLATION BLOCK works next time
   ${MOVE} prior_inflate_restart prior_inflate.${OCN_DATE_EXT}.restart.be 
   ${LINK} prior_inflate.${OCN_DATE_EXT}.restart.be prior_inflate_ics
   ${MOVE} prior_inflate.${OCN_DATE_EXT}.restart.be prior_inflate_ics ..
else
   echo "No prior_inflate_restart for ${OCN_DATE_EXT}"
endif

if ( -e prior_inflate_diag ) then
   ${MOVE} prior_inflate_diag ../prior_inflate.${OCN_DATE_EXT}.diag
else
   echo "No prior_inflate_diag for ${OCN_DATE_EXT}"
endif

# FIXME: special for trying out non-monotonic task layouts.
setenv PATH "${ORG_PATH}"
setenv LSB_PJL_TASK_GEOMETRY "${ORG_TASK_GEOMETRY}"

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
   ${LINK} ../$OCN_RESTART_FILENAME pop.r.nc
   ${LINK} ../pop2_in.$member       pop_in

   ${LINK} $DART_RESTART_FILE dart.ic

   ./dart_to_pop || exit 3

   @ member++
end

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

ls -lrt

