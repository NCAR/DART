#!/usr/local/bin/tcsh
#
# DART software - Copyright © 2004 - 2010 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$

# The FORCE options are not optional. 
# the VERBOSE options are useful for debugging.
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

   if ( -e ${CASEROOT}/${FILE} ) then
      ${COPY}   ${CASEROOT}/${FILE} .
   else if ( -e    ../${FILE} ) then
      ${COPY} ../${FILE} .
   else if ( -e ${DARTDIR}/${FILE} ) then
      ${COPY}   ${DARTDIR}/${FILE} .
   else
      echo "DART required file $FILE not found ... ERROR"
      exit 1
   endif

end

#-------------------------------------------------------------------------
# This is the file for the sampling error correction.
# Each ensemble size has its own file.
# It is static - it does not need to be archived, etc.
# It is only needed if 
# input.nml:&assim_tools_nml:sampling_error_correction = .true.,
#-------------------------------------------------------------------------

set SAMP_ERR_FILE = ${DARTDIR}/system_simulation/final_full.${ensemble_size}

if ( -e ${SAMP_ERR_FILE}/ ) then
   ${COPY} ${SAMP_ERR_FILE} .
else
   echo "WARNING: no sampling error correction file for this ensemble size."
   echo "warning: looking for system_simulation/final_full.${ensemble_size}"
endif

#-------------------------------------------------------------------------
# DART INFLATION BLOCK
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
# The strategy is to use the LATEST inflation file from CENTRALDIR if one exists - 
#
# After an assimilation, the output file will be copied back to CENTRALDIR
# to be used for subsequent assimilations.
#-------------------------------------------------------------------------

foreach FILE ( prior post ) 

   # These files may or may not exist. This causes some complexity.
   # So - we look for the 'newest' and use it. And Pray.

   (ls -rt1 ../${FILE}_inflate.*.restart.* | tail -1 >! latestfile) > & /dev/null
   set nfiles = `cat latestfile | wc -l`

   if ( $nfiles > 0 ) then
      set latest = `cat latestfile`
      ${LINK} $latest ${FILE}_inflate_ics
   else
      # MUST HAVE inf_initial_from_restart = .false.
      echo "WARNING: no incoming ${FILE}_inflate.YYYY-MM-DD-00000.restart.endiansuffix"
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

   # Each member will do its job in its own directory.
   # That way, we can do N of them simultaneously -
   # they all read their OWN 'input.nml' ... the output
   # filenames must inserted into the appropriate input.nml

   set MYTEMPDIR = member_${member}
   mkdir -p $MYTEMPDIR
   cd $MYTEMPDIR

   set OCN_RESTART_FILENAME = `head -1 ../../rpointer.ocn.$member.restart`
   ${LINK} ../../$OCN_RESTART_FILENAME pop.r.nc
   ${LINK} ../../pop2_in.$member       pop_in

   # the slash in the filename screws up 'sed' ... unless
   set DART_IC_FILE = `printf ..\\/filter_ics.%04d $member`

   sed -e "s/dart.ud/${DART_IC_FILE}/" < ../input.nml >! input.nml

   ../pop_to_dart &

   cd ..

   @ member++
end

wait

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

# POP always needs a pop_in and a pop.r.nc to start.

set OCN_RESTART_FILENAME = `head -1 ../rpointer.ocn.1.restart`
${LINK} ../$OCN_RESTART_FILENAME pop.r.nc
${LINK} ../pop2_in.1  pop_in

# Determine proper observation sequence file.

set OBSFNAME = `printf obs_seq.0Z.%04d%02d%02d ${OCN_YEAR} ${OCN_MONTH} ${OCN_DAY}`
set OBS_FILE = ${OBSDIR}/${OBSFNAME} 

${LINK} ${OBS_FILE}   obs_seq.out

# FIXME: special for trying out non-monotonic task layouts.
setenv ORG_PATH "${PATH}"
setenv LSF_BINDIR /contrib/lsf/tgmpatch
setenv PATH ${LSF_BINDIR}:${PATH}
setenv ORG_TASK_GEOMETRY "${LSB_PJL_TASK_GEOMETRY}"

# layout 1: rr by node
setenv NANCY_GEOMETRY_126_2NODES_RR \
	"{(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124)(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83,85,87,89,91,93,95,97,99,101,103,105,107,109,111,113,115,117,119,121,123,125)}";
	
# layout 2: flat
setenv NANCY_GEOMETRY_126_2NODES_FL \
	"{(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62)(63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125)}";

# layout 3: rr sub block by stride
setenv NANCY_GEOMETRY_126_2NODES_STR \
	"{(0,64,2,66,4,68,6,70,8,72,10,74,12,76,14,78,16,80,18,82,20,84,22,86,24,88,26,90,28,92,30,94,32,96,34,98,36,100,38,102,40,104,42,106,44,108,46,110,48,112,50,114,52,116,54,118,56,120,58,122,60,124,62)(1,65,3,67,5,69,7,71,9,73,11,75,13,77,15,79,17,81,19,83,21,85,23,87,25,89,27,91,29,93,31,95,33,97,35,99,37,101,39,103,41,105,43,107,45,109,47,111,49,113,51,115,53,117,55,119,57,121,59,123,61,125,63)}";

setenv NANCY_GEOMETRY_126_3NODES \
	"{(0,29,30,31,1,32,33,34,2,35,36,37,3,38,39,40,4,41,42,43,5,44,45,46,6,47,48,49,7,50,51,52,8,53,54,55,9,56,57,58,10,59)(60,61,11,62,63,64,12,65,66,67,13,68,69,70,14,71,72,73,15,74,75,76,16,77,78,79,17,80,81,82,18,83,84,85,19,86,87,88,20,89,90,91)(21,92,93,94,22,95,96,97,23,98,99,100,24,101,102,103,25,104,105,106,26,107,108,109,27,110,111,112,28,113,114,115,116,117,118,119,120,121,122,123,124,125)}"

setenv NANCY_GEOMETRY_126_6NODES \
	"{(0,29,30,31,1,32,33,34,2,35,36,37,3,38,39,40,4,41,42,43,5)(44,45,46,6,47,48,49,7,50,51,52,8,53,54,55,9,56,57,58,10,59)(60,61,11,62,63,64,12,65,66,67,13,68,69,70,14,71,72,73,15,74,75)(76,16,77,78,79,17,80,81,82,18,83,84,85,19,86,87,88,20,89,90,91)(21,92,93,94,22,95,96,97,23,98,99,100,24,101,102,103,25,104,105,106,26)(107,108,109,27,110,111,112,28,113,114,115,116,117,118,119,120,121,122,123,124,125)}"

setenv NANCY_GEOMETRY_126_7NODES \
         "{(0,7,14,21,28,35,42,49,56,63,70,77,84,91,98,105,112,119)(1,8,15,22,29,36,43,50,57,64,71,78,85,92,99,106,113,120)(2,9,16,23,30,37,44,51,58,65,72,79,86,93,100,107,114,121)(3,10,17,24,31,38,45,52,59,66,73,80,87,94,101,108,115,122)(4,11,18,25,32,39,46,53,60,67,74,81,88,95,102,109,116,123)(5,12,19,26,33,40,47,54,61,68,75,82,89,96,103,110,117,124)(6,13,20,27,34,41,48,55,62,69,76,83,90,97,104,111,118,125)}"

# layout: flat
setenv NANCY_GEOMETRY_54_1NODE \
	"{(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53)}";

setenv LSB_PJL_TASK_GEOMETRY "${NANCY_GEOMETRY_54_1NODE}"

which mpirun.lsf

mpirun.lsf ./filter || exit 2

${MOVE} Prior_Diag.nc      ../Prior_Diag.${OCN_DATE_EXT}.nc
${MOVE} Posterior_Diag.nc  ../Posterior_Diag.${OCN_DATE_EXT}.nc
${MOVE} obs_seq.final      ../obs_seq.${OCN_DATE_EXT}.final
${MOVE} dart_log.out       ../dart_log.${OCN_DATE_EXT}.out

# Accomodate any possible inflation files 

foreach INFLATION ( prior post )

   if ( -e ${INFLATION}_inflate_restart ) then
      # 1) rename file to reflect current date
      # 2) move to CENTRALDIR so the DART INFLATION BLOCK works next time
   
      ${MOVE} ${INFLATION}_inflate_restart ../${INFLATION}_inflate.${OCN_DATE_EXT}.restart.be 

   else
      echo "No ${INFLATION}_inflate_restart for ${OCN_DATE_EXT}"
   endif

   if ( -e ${INFLATION}_inflate_diag ) then
      ${MOVE} ${INFLATION}_inflate_diag ../${INFLATION}_inflate.${OCN_DATE_EXT}.diag
   else
      echo "No ${INFLATION}_inflate_diag for ${OCN_DATE_EXT}"
   endif

end

# FIXME: special for trying out non-monotonic task layouts.
setenv PATH "${ORG_PATH}"
setenv LSB_PJL_TASK_GEOMETRY "${ORG_TASK_GEOMETRY}"

#-------------------------------------------------------------------------
# Block 3: Update the POP restart files ... simultaneously ...
#
# DART namelist settings required:
# &filter_nml:           restart_out_file_name  = 'filter_restart'
# &ensemble_manager_nml: single_restart_file_in = '.false.'
# &dart_to_pop_nml:      dart_to_pop_input_file = 'dart.ic',
# &dart_to_pop_nml:      advance_time_present   = .false.
#-------------------------------------------------------------------------

set member = 1
while ( $member <= $ensemble_size )

   # Each member will do its job in its own directory.
   # That way, we can do N of them simultaneously -
   # they all read their OWN 'input.nml' ... the output
   # filenames must inserted into the appropriate input.nml

   set MYTEMPDIR = member_${member}
   mkdir -p $MYTEMPDIR
   cd $MYTEMPDIR

   set DART_RESTART_FILE = `printf filter_restart.%04d $member`
   ${LINK} ../$DART_RESTART_FILE dart.ic

   set OCN_RESTART_FILENAME = `head -1 ../../rpointer.ocn.$member.restart`
   ${LINK} ../../$OCN_RESTART_FILENAME pop.r.nc
   ${LINK} ../../pop2_in.$member       pop_in

   cp -f ../input.nml .

   ../dart_to_pop &

   cd ..

   @ member++
end

wait

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

ls -lrt

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

