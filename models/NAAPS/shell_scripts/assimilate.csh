#!/bin/tcsh 
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# The FORCE options are not optional. 
# the VERBOSE options are useful for debugging.
set   MOVE = '/bin/mv -fv'
set   COPY = '/bin/cp -fv --preserve=timestamps'
set   LINK = '/bin/ln -fvs'
set REMOVE = '/bin/rm -fr'

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

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
#-----------------------------------------------------------------------------
#-------------------------------------------------------------------------
# DART COPY BLOCK
# Populate a run-time directory with the bits needed to run DART 
# Either get them from the CCSM 'run' directory or some stock repository
# The grid files are absolute paths ... so they need not move.
#-------------------------------------------------------------------------

set DARTDIR = /shared/aerosol_heck2/users/sessions/_DART/kodiak/models/NAAPS/work

foreach FILE ( input.nml.template filter naaps_to_dart dart_to_naaps )

   if ( -e ${DARTDIR}/${FILE} ) then
      ${COPY}   ${DARTDIR}/${FILE} .
   else
      echo "DART required file $FILE not found ... ERROR"
      exit 1
   endif

end

#-------------------------------------------------------------------------
# Determine time of model state ... input namelist
#-------------------------------------------------------------------------

set FILE = `grep dtg input.nml.template`
set MODEL_DATE = `echo $FILE | sed -e "s#[',]##g"`
set DTG         = $MODEL_DATE[3]
set MODEL_YEAR    = `echo $DTG | cut -c1-4`
set MODEL_MONTH   = `echo $DTG | cut -c5-6`
set MODEL_DAY     = `echo $DTG | cut -c7-8`
set MODEL_HOURS   = `echo $DTG | cut -c9-10`

echo "assimilate.csh: valid time of model is $MODEL_YEAR $MODEL_MONTH $MODEL_DAY $MODEL_HOURS"
set YYYYMM = "${MODEL_YEAR}${MODEL_MONTH}"
set OBSDIR = /shared/aerosol_jack/users/sessions/products/DART/OBS_SEQ/${YYYYMM}

set FILE = `grep ens_size input.nml.template | sed -e "s#[',]##g"`
set ensemble_size = $FILE[3]

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

   (ls -rt1 ../${FILE}_inflate.*.restart.* | tail -n 1 >! latestfile) > & /dev/null
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
# Block 1: convert N NAAPS restart files to DART initial conditions file(s)
# naaps_to_dart is serial code, we can do all of these at the same time
# and just wait for them to finish IFF it were not for the fact we'd have
# to have unique namelists for all of them.
#
# At the end of the block, we have DART restart files  filter_ics.[1-N]
#
# DART namelist settings appropriate/required:
# &ensemble_manager_nml: single_restart_file_in    = '.false.'
# &filter_nml:           restart_in_file_name      = 'filter_ics'
# &naaps_to_dart_nml:    naaps_to_dart_output_file = 'dart_ics',
#-------------------------------------------------------------------------

set member = 1
while ( $member <= $ensemble_size )

   # Each member will do its job in its own directory.
   # That way, we can do N of them simultaneously -
   # they all read their OWN 'input.nml' ... the output
   # filenames must inserted into the appropriate input.nml

   set MYTEMPDIR = `printf member_%04d $member`
   mkdir -p $MYTEMPDIR
   cd $MYTEMPDIR

   # the slash in the filename screws up 'sed' ... unless
   set DART_IC_FILE = `printf ../filter_ics.%04d $member`

   sed -e "s#dart_ics#${DART_IC_FILE}#" < ../input.nml.template >! tmp.nml
   sed -e "s#MEMBER_NUMBER#$member#" < tmp.nml >! input.nml
   $REMOVE tmp.nml

   ../naaps_to_dart &

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

# Determine proper observation sequence file.

set OBSFNAME = "${YYYYMM}/obs_seq.out.${DTG}"
set OBS_FILE = /shared/aerosol_jack/users/sessions/products/DART/OBS_SEQ/${OBSFNAME}

${LINK} ${OBS_FILE} obs_seq.out

sed -e "s/MEMBER_NUMBER/1/" < input.nml.template >! input.nml

./filter || exit 2

${MOVE} preassim.nc      ../preassim.${DTG}.nc
${MOVE} analysis.nc      ../analysis.${DTG}.nc
${MOVE} obs_seq.final    ../obs_seq.${DTG}.final
${MOVE} dart_log.out     ../dart_log.${DTG}.out

# Accomodate any possible inflation files 

foreach INFLATION ( prior post )

   if ( -e ${INFLATION}_inflate_restart ) then
      # 1) rename file to reflect current date
      # 2) move to CENTRALDIR so the DART INFLATION BLOCK works next time
   
      ${MOVE} ${INFLATION}_inflate_restart ../${INFLATION}_inflate.${DTG}.restart.be 

   else
      echo "No ${INFLATION}_inflate_restart for ${DTG}"
   endif

   if ( -e ${INFLATION}_inflate_diag ) then
      ${MOVE} ${INFLATION}_inflate_diag ../${INFLATION}_inflate.${DTG}.diag
   else
      echo "No ${INFLATION}_inflate_diag for ${DTG}"
   endif

end

#-------------------------------------------------------------------------
# Block 3: Update the naaps restart files ... simultaneously ...
#
# DART namelist settings required:
# &ensemble_manager_nml: single_restart_file_in   = '.false.'
# &filter_nml:           restart_out_file_name    = 'filter_restart'
# &dart_to_naaps_nml:    dart_to_naaps_input_file = 'dart_restart',
# &dart_to_naaps_nml:    advance_time_present     = .false.
#-------------------------------------------------------------------------

set member = 1
while ( $member <= $ensemble_size )

   # Each member will do its job in its own directory.
   # That way, we can do N of them simultaneously -
   # they all read their OWN 'input.nml' ... the output
   # filenames must inserted into the appropriate input.nml

   set MYTEMPDIR = `printf member_%04d $member`
   mkdir -p $MYTEMPDIR
   cd $MYTEMPDIR

   set DART_RESTART_FILE = `printf ../filter_restart.%04d $member`
   sed -e "s#dart_restart#${DART_RESTART_FILE}#" < ../input.nml.template >! tmp.nml
   sed -e "s#MEMBER_NUMBER#$member#" < tmp.nml >! input.nml
   $REMOVE tmp.nml

   ../dart_to_naaps &

   cd ..

   @ member++
end

wait

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

exit 0


