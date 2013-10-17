#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# The FORCE options are not optional.
# the VERBOSE options are useful for debugging.
set   MOVE = '/usr/local/bin/mv -fv'
set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
set   LINK = '/usr/local/bin/ln -fvs'
set REMOVE = '/usr/local/bin/rm -fr'

set ensemble_size = ${NINST_ATM}

# Create temporary working directory for the assimilation
set temp_dir = assimilate_dir
echo "temp_dir is $temp_dir"
mkdir -p $temp_dir
cd $temp_dir

#-------------------------------------------------------------------------
# Determine time of model state ... from file name of first member
# of the form "./${CASE}.cam_${ensemble_member}.r.2000-01-06-00000.nc"
#-------------------------------------------------------------------------

set FILE = `head -n 1 ../rpointer.atm_0001`
set FILE = $FILE:t
set FILE = $FILE:r
set MYCASE = `echo $FILE | sed -e "s#\..*##"`
set MODEL_DATE_EXT = `echo $FILE:e`
set MODEL_DATE     = `echo $FILE:e | sed -e "s#-# #g"`
set MODEL_YEAR     = $MODEL_DATE[1]
set MODEL_MONTH    = $MODEL_DATE[2]
set MODEL_DAY      = $MODEL_DATE[3]
set MODEL_SECONDS  = $MODEL_DATE[4]
set MODEL_HOUR     = `echo $MODEL_DATE[4] / 3600 | bc`

echo "valid time of model is $MODEL_YEAR $MODEL_MONTH $MODEL_DAY $MODEL_SECONDS (seconds)"
echo "valid time of model is $MODEL_YEAR $MODEL_MONTH $MODEL_DAY $MODEL_HOUR (hours)"

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
#-----------------------------------------------------------------------------

set DARTDIR = ${HOME}/svn/DART/dev
set  OBSDIR = /glade/proj3/image/Observations/Synthetic/UVT_Set2_12H

#-------------------------------------------------------------------------
# DART COPY BLOCK
# Populate a run-time directory with the bits needed to run DART
#-------------------------------------------------------------------------

foreach FILE ( input.nml.zagar filter cam_to_dart dart_to_cam )
   if (  -e   ${DARTDIR}/models/cam/work/${FILE} ) then
      ${COPY} ${DARTDIR}/models/cam/work/${FILE} .
   else
      echo "DART required file ${DARTDIR}/${FILE} not found ... ERROR"
      exit 1
   endif
end

${MOVE} input.nml.zagar input.nml

${COPY} /glade/proj3/DART/raeder/FV1deg_4.0/cam_phis.nc .

#-------------------------------------------------------------------------
# DART SAMPLING ERROR CORRECTION BLOCK
# This stages the files needed for the sampling error correction.
# Each ensemble size has its own (static) file which does not need to be archived.
# It is only needed if
# input.nml:&assim_tools_nml:sampling_error_correction = .true.,
#-------------------------------------------------------------------------

set  MYSTRING = `grep sampling_error_correction input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set SECSTRING = `echo $MYSTRING[2] | tr 'A-Z' 'a-z'`

if ( $SECSTRING == true ) then
   set SAMP_ERR_FILE = ${DARTDIR}/system_simulation/final_full_precomputed_tables/final_full.${ensemble_size}
   if (  -e   ${SAMP_ERR_FILE} ) then
      ${COPY} ${SAMP_ERR_FILE} .
   else
      echo "ERROR: no sampling error correction file for this ensemble size."
      echo "ERROR: looking for ${SAMP_ERR_FILE}"
      exit 2
   endif
else
   echo "Sampling Error Correction not requested for this assimilation."
endif

#-------------------------------------------------------------------------
# DART INFLATION BLOCK
# This stages the files that contain the inflation values.
# The inflation values change through time and should be archived.
#
# This file is only relevant if 'inflation' is turned on -
# i.e. if inf_flavor(1) /= 0 AND inf_initial_from_restart = .TRUE.
#
# filter_nml
# inf_flavor                  = 2,                       0,
# inf_initial_from_restart    = .true.,                  .false.,
# inf_in_file_name            = 'prior_inflate_ics',     'post_inflate_ics',
# inf_out_file_name           = 'prior_inflate_restart', 'post_inflate_restart',
# inf_diag_file_name          = 'prior_inflate_diag',    'post_inflate_diag',
#
# NOTICE: the archiving scripts more or less require the names of these
# files to be as listed above. When being archived, the filenames get a
# unique extension (describing the assimilation time) appended to them.
#
# The inflation file is essentially a duplicate of the model state ...
# it is slaved to a specific geometry. The initial files are created
# offline with values of unity. For the purpose of this script, they are
# thought to be the output of a previous assimilation, so they should be
# named something like prior_inflate_restart.YYYY-MM-DD-SSSSS
#
# The first inflation file can be created with 'fill_inflation_restart'
# which can be built in the usual DART manner.
#
# The strategy is to use the LATEST inflation file from the CESM 'rundir'.
# After an assimilation, the new inflation values/files will be moved to
# the CESM rundir to be used for subsequent assimilations. If the short-term
# archiver has worked correctly, only the LATEST files will available. Of
# course, it is not required to have short-term archiving turned on, so ...
#-------------------------------------------------------------------------

set  MYSTRING = `grep inf_flavor input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,]# #g"`
set  PRIOR_INF = $MYSTRING[2]
set  POSTE_INF = $MYSTRING[3]

set  MYSTRING = `grep inf_initial_from_restart input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,]# #g"`
set  PRIOR_TF = `echo $MYSTRING[2] | tr [A-Z] [a-z]`
set  POSTE_TF = `echo $MYSTRING[3] | tr [A-Z] [a-z]`

# its a little tricky to remove both styles of quotes from the string.

set  MYSTRING = `grep inf_in_file_name input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set  PRIOR_INF_IFNAME = $MYSTRING[2]
set  POSTE_INF_IFNAME = $MYSTRING[3]

set  MYSTRING = `grep inf_out_file_name input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set  PRIOR_INF_OFNAME = $MYSTRING[2]
set  POSTE_INF_OFNAME = $MYSTRING[3]

set  MYSTRING = `grep inf_diag_file_name input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set  PRIOR_INF_DIAG = $MYSTRING[2]
set  POSTE_INF_DIAG = $MYSTRING[3]

# IFF we want PRIOR inflation:

if ( $PRIOR_INF > 0 ) then

   if ($PRIOR_TF == ".false.") then
      echo "ERROR: inf_flavor(1) = $PRIOR_INF, yet inf_initial_from_restart = $PRIOR_TF"
      echo "ERROR: fix input.nml to reflect whether you want prior inflation or not."
      exit 3
   endif

   # Look for the output from the previous assimilation
   (ls -rt1 ../${PRIOR_INF_OFNAME}.* | tail -n 1 >! latestfile) > & /dev/null
   set nfiles = `cat latestfile | wc -l`

   # If one exists, use it as input for this assimilation
   if ( $nfiles > 0 ) then
      set latest = `cat latestfile`
      ${LINK} $latest ${PRIOR_INF_IFNAME}
   else
      echo "ERROR: Requested prior inflation but specified no incoming prior inflation file."
      echo "ERROR: expected something like ../${PRIOR_INF_OFNAME}.YYYY-MM-DD-SSSSS"
      exit 4
   endif
else
   echo "Prior Inflation not requested for this assimilation."
endif

# POSTERIOR: We look for the 'newest' and use it - IFF we need it.

if ( $POSTE_INF > 0 ) then

   if ($POSTE_TF == ".false.") then
      echo "ERROR: inf_flavor(2) = $POSTE_INF, yet inf_initial_from_restart = $POSTE_TF"
      echo "ERROR: fix input.nml to reflect whether you want posterior inflation or not."
      exit 5
   endif

   # Look for the output from the previous assimilation
   (ls -rt1 ../${POSTE_INF_OFNAME}.* | tail -n 1 >! latestfile) > & /dev/null
   set nfiles = `cat latestfile | wc -l`

   # If one exists, use it as input for this assimilation
   if ( $nfiles > 0 ) then
      set latest = `cat latestfile`
      ${LINK} $latest ${POSTE_INF_IFNAME}
   else
      echo "ERROR: Requested POSTERIOR inflation but specified no incoming POSTERIOR inflation file."
      echo "ERROR: expected something like ../${POSTE_INF_OFNAME}.YYYY-MM-DD-SSSSS"
      exit 6
   endif
else
   echo "Posterior Inflation not requested for this assimilation."
endif

#-------------------------------------------------------------------------
# Block 1: convert N cam restart files to DART initial conditions file(s)
# cam_to_dart is serial code, we can do all of these at the same time
# and just wait for them to finish IFF it were not for the fact we'd have
# to have unique namelists for all of them.
#
# At the end of the block, we have DART restart files  filter_ic_old.[1-N]
# that came from pointer files ../rpointer.atm_[1-N]
#
# DART namelist settings appropriate/required:
# &filter_nml:           restart_in_file_name    = 'filter_ic_old'
# &ensemble_manager_nml: single_restart_file_in  = '.false.'
# &cam_to_dart_nml:      cam_to_dart_output_file = 'dart_ics',
#-------------------------------------------------------------------------

set member = 1
while ( ${member} <= ${ensemble_size} )

   # Each member will do its job in its own directory.
   # That way, we can do N of them simultaneously -
   # they all read their OWN 'input.nml' ... the output
   # filenames must inserted into the appropriate input.nml

   set MYTEMPDIR = member_${member}
   mkdir -p $MYTEMPDIR
   cd $MYTEMPDIR

   set POINTER_FILENAME = `printf rpointer.atm_%04d ${member}`
   set MODEL_RESTART_FILENAME = `head -n 1 ../../${POINTER_FILENAME}`
   set MODEL_INITIAL_FILENAME = `echo ${MODEL_RESTART_FILENAME} | sed "s#\.r\.#\.i\.#"`
   ${LINK} ../../$MODEL_INITIAL_FILENAME caminput.nc
   ${LINK} ../cam_phis.nc .

   # TJH can we use a .h0. file instead of some arbitrary cam_phis.nc

   set DART_IC_FILE = `printf ../filter_ic_old.%04d ${member}`

   sed -e "s#dart_ics#${DART_IC_FILE}#" < ../input.nml >! input.nml

   echo "starting cam_to_dart for member ${member} at "`date`
   ../cam_to_dart >! output.${member}.cam_to_dart &
   echo "finished cam_to_dart for member ${member} at "`date`

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
# &filter_nml:           restart_in_file_name   = 'filter_ic_old'
# &filter_nml:           restart_out_file_name  = 'filter_ic_new'
# &filter_nml:           obs_sequence_in_name   = 'obs_seq.out'
# &filter_nml:           obs_sequence_out_name  = 'obs_seq.final'
# &filter_nml:           init_time_days         = -1,
# &filter_nml:           init_time_seconds      = -1,
# &filter_nml:           first_obs_days         = -1,
# &filter_nml:           first_obs_seconds      = -1,
# &filter_nml:           last_obs_days          = -1,
# &filter_nml:           last_obs_seconds       = -1,
# &ensemble_manager_nml: single_restart_file_in = .false.
#
#-------------------------------------------------------------------------

# cam always needs a cam_initial.nc and a cam_history.nc to start.

set MODEL_RESTART_FILENAME = `head -n 1 ../rpointer.atm_0001`
set MODEL_INITIAL_FILENAME = `echo ${MODEL_RESTART_FILENAME} | sed "s#\.r\.#\.i\.#"`
set MODEL_HISTORY_FILENAME = `echo ${MODEL_RESTART_FILENAME} | sed "s#\.r\.#\.h0\.#"`

${LINK} ../$MODEL_INITIAL_FILENAME caminput.nc
#${LINK} ../$MODEL_RESTART_FILENAME cam_restart.nc
#${LINK} ../$MODEL_HISTORY_FILENAME cam_history.nc

# stage the proper observation sequence file.

set OBSFNAME = `printf syn_obs_seq${MODEL_YEAR}${MODEL_MONTH}${MODEL_DAY}%02d ${MODEL_HOUR}`
set OBS_FILE = ${OBSDIR}/${OBSFNAME}

${LINK} ${OBS_FILE} obs_seq.out

mpirun.lsf ./filter || exit 7

${MOVE} Prior_Diag.nc      ../Prior_Diag.${MODEL_DATE_EXT}.nc
${MOVE} Posterior_Diag.nc  ../Posterior_Diag.${MODEL_DATE_EXT}.nc
${MOVE} obs_seq.final      ../obs_seq.${MODEL_DATE_EXT}.final
${MOVE} dart_log.out       ../dart_log.${MODEL_DATE_EXT}.out

# Accomodate any possible inflation files
# 1) rename file to reflect current date
# 2) move to CENTRALDIR so the DART INFLATION BLOCK works next time and
#    that they can get archived.

foreach FILE ( ${PRIOR_INF_OFNAME} ${POSTE_INF_OFNAME} ${PRIOR_INF_DIAG} ${POSTE_INF_DIAG} )
   if ( -e ${FILE} ) then
      ${MOVE} ${FILE} ../${FILE}.${MODEL_DATE_EXT}
   else
      echo "No ${FILE} for ${MODEL_DATE_EXT}"
   endif
end

#-------------------------------------------------------------------------
# Block 3: Update the cam restart files ... simultaneously ...
#
# DART namelist settings required:
# &filter_nml:           restart_out_file_name  = 'filter_ic_new'
# &ensemble_manager_nml: single_restart_file_in = '.false.'
# &dart_to_cam_nml:      dart_to_cam_input_file = 'dart_restart',
# &dart_to_cam_nml:      advance_time_present   = .false.
# &atm_in_xxxx:ncdata = 'cam_initial_x.nc'
#-------------------------------------------------------------------------

set member = 1
while ( ${member} <= ${ensemble_size} )

   # Each member will do its job in its own directory.
   # After they are all done, we can move them.

   set MYTEMPDIR = member_${member}
   mkdir -p $MYTEMPDIR
   cd $MYTEMPDIR

   set DART_RESTART_FILE = `printf filter_ic_new.%04d ${member}`
   ${LINK} ../$DART_RESTART_FILE dart_restart

   set ATM_POINTER_FILENAME = `printf rpointer.atm_%04d ${member}`
   set LND_POINTER_FILENAME = `printf rpointer.lnd_%04d ${member}`
   set ICE_POINTER_FILENAME = `printf rpointer.ice_%04d ${member}`

   set ATM_RESTART_FILENAME = `head -n 1 ../../${ATM_POINTER_FILENAME}`
   set LND_RESTART_FILENAME = `echo ${ATM_RESTART_FILENAME} | sed "s#\.cam_#\.clm2_#"`
   set ICE_RESTART_FILENAME = `echo ${ATM_RESTART_FILENAME} | sed "s#\.cam_#\.cice_#"`

   set ATM_INITIAL_FILENAME = `echo ${ATM_RESTART_FILENAME} | sed "s#\.r\.#\.i\.#"`

   ${LINK} ../../$ATM_INITIAL_FILENAME caminput.nc

   echo "starting dart_to_cam for member ${member} at "`date`
   ../dart_to_cam >! output.${member}.dart_to_cam &
   echo "finished dart_to_cam for member ${member} at "`date`

   cd ..

   @ member++
end

wait

#-------------------------------------------------------------------------
# Block 4: The cam files have now been updated, move them into position.
#-------------------------------------------------------------------------

set member = 1
while ( ${member} <= ${ensemble_size} )

   cd member_${member}

   set ATM_POINTER_FILENAME = `printf rpointer.atm_%04d ${member}`
   set LND_POINTER_FILENAME = `printf rpointer.lnd_%04d ${member}`
   set ICE_POINTER_FILENAME = `printf rpointer.ice_%04d ${member}`

   set ATM_RESTART_FILENAME = `head -n 1 ../../${ATM_POINTER_FILENAME}`
   set LND_RESTART_FILENAME = `echo ${ATM_RESTART_FILENAME} | sed "s#\.cam_#\.clm2_#"`
   set ICE_RESTART_FILENAME = `echo ${ATM_RESTART_FILENAME} | sed "s#\.cam_#\.cice_#"`

   set ATM_INITIAL_FILENAME = `echo ${ATM_RESTART_FILENAME} | sed "s#\.r\.#\.i\.#"`

   ${COPY} ../../$ATM_INITIAL_FILENAME ../../cam_initial_${member}.nc
   ${COPY} ../../$LND_RESTART_FILENAME ../../clm_restart_${member}.nc
   ${COPY} ../../$ICE_RESTART_FILENAME ../../ice_restart_${member}.nc

   cd ..

   @ member++
end

#-------------------------------------------------------------------------
# Now that everything is staged, we have to communicate the current
# model time to the drv_in&seq_timemgr_inparm namelist
# which is built from CASEROOT/user_nl_drv by the *.run script
#-------------------------------------------------------------------------

ex ${CASEROOT}/Buildconf/cpl.buildnml.csh << ex_end
g; start_ymd;s;=[ ]*.*;= ${MODEL_YEAR}${MODEL_MONTH}${MODEL_DAY};
g; start_tod;s;=[ ]*.*;= $MODEL_SECONDS;
wq
ex_end

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

