#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# A shell script to run the MPAS-A(tmostphere) model from DART input.
#
# This script is called by advance_model.template in driver_mpas_dart.csh
# after the analysis step is done at each cycle.
#
# This script performs the following:
# 1.  Creates a temporary directory to run an MPAS-A realization (see options)
# 2.  Gets all the files necessary for the model run in the member directory.
# 3.  Updates an MPAS namelist from a template with new dates.
# 4.  Runs the MPAS-A model in a restart mode until the target time is reached.
# 5.  Regional MPAS is also supported, if chosen.
# 6.  Checks for incomplete runs.
# 7.  Saves horizontal winds if save_wind = true (for extended forecasts later).
#
# Note: 1. This script supports MPAS V7+ and the Manhattan release of DART. 
#          It is NOT backward compatible for older versions.
#       2. MPAS is run in a restart mode during the cycles, which means
#          both input and output of the model run are restart files.
#          This also means that one should not delete the member directory as
#          one needs to keep the restart file for each member during the cycle.
#       3. The input analysis file should be provided through update_mpas_states 
#          (which is supposed to be run before running this script).
#       4. For the required data to run this script, check the section of 'dependencies'.
#       5. Anything specific to the experiment is supposed to be provided in ${CENTRALDIR}/.
#       6. Regional MPAS is not supported in DART yet. 
#          For more info, please contact Soyoung Ha (syha@ucar.edu).
#
# Arguments for this script (created by 'filter' or 'perfect_model_obs') are:
# 1) ensemble member number
# 2) maximum ensemble member number
# 
# This script can loop over all the ensemble members unless the two input 
# arguments are identical for a specific ensemble member number.
#----------------------------------------------------------------------
set ensemble_member = $1
set ensemble_max    = $2

# Do you want to save horizontal winds in the analysis file?
#-------------------------------------------------------------------
set save_wind = false

# mpi command
#-------------------------------------------------------------------
#set mpicmd = "mpi -n 4"			# Mac OS
set mpicmd = "mpiexec_mpt"			# dplace -s 1"	# Cheyenne

# Other commands
#-------------------------------------------------------------------
set  REMOVE = 'rm -rf'
set    COPY = 'cp -pf'
set    MOVE = 'mv -f'
set    LINK = 'ln -sf'
unalias cd
unalias ls

# The run-time directory for the entire experiment is called CENTRALDIR;
#-------------------------------------------------------------------
set CENTRALDIR = `pwd`

# Copy necessary input files/executables/files common
# to all model advances to a clean, temporary directory.
#-------------------------------------------------------------------
if ( ! -r ${CENTRALDIR}/input.nml ) then
     echo ABORT\: advance_model.csh could not find required readable dependency ${CENTRALDIR}/input.nml
     exit 1
endif

foreach f ( namelist.atmosphere streams.atmosphere )
if ( ! -r ${CENTRALDIR}/$f ) then
     echo ABORT\: advance_model.csh could not find required readable dependency ${CENTRALDIR}/$f
     echo The file is assumed to be edited for your own configuration.
     exit 1
endif
end

if ( ! -x ${CENTRALDIR}/advance_time ) then
     echo ABORT\: advance_model.csh could not find required executable dependency ${CENTRALDIR}/advance_time
     exit 1
endif

if ( ! -d ${CENTRALDIR}/MPAS_RUN ) then
      echo ABORT\: advance_model.csh could not find required data directory ${CENTRALDIR}/MPAS_RUN, 
      echo         which contains all the default input files for running MPAS/atmosphere_model.
      exit 1
endif

if ( ! -x ${CENTRALDIR}/MPAS_RUN/atmosphere_model ) then
     echo ABORT\: advance_model.csh could not find required executable dependency 
     echo         ${CENTRALDIR}/MPAS_RUN/atmosphere_model
     exit 1
endif

# A list of input analysis file names 
#-------------------------------------------------------------------
set inlist  = `grep update_output_file_list ${CENTRALDIR}/input.nml | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g" | sed -e 's/"//g'`

set sample = `head -1 $inlist`
set dhead  = `dirname $sample | cut -c1-6`
if($dhead != "member") then
   echo "Check temp_dir below. The directory name cannot start with 'member'."
   echo "Input ensemble directories are named as $dhead instead."
   exit
endif

# Common input files based on the model configuration
#-------------------------------------------------------------------

   # Get the grid info files - now for PIO
   set fs_grid = `grep config_block_decomp_file_prefix ${CENTRALDIR}/namelist.atmosphere | awk '{print $3}' | sed -e "s/'//g"`

   # Surface update
   set if_sfc_update = `grep config_sst_update ${CENTRALDIR}/namelist.atmosphere | awk '{print $3}'`
   set fsfc = `sed -n '/<stream name=\"surface\"/,/\/>/{/Scree/{p;n};/##/{q};p}' ${CENTRALDIR}/streams.atmosphere | \
               grep filename_template | awk -F= '{print $2}' | awk -F$ '{print $1}' | sed -e 's/"//g'`

   # Sanity check - A switch for cycling
   set if_DAcycling = `grep config_do_DAcycling ${CENTRALDIR}/namelist.atmosphere | wc -l`
   if($if_DAcycling == 0) then
      echo "Please add config_do_DAcycling = .true. in \&restart"
      echo "in ${CENTRALDIR}/namelist.atmosphere."
      exit -1
   endif

#----------------------------------------------------------------------
# A main section for the model integration
#----------------------------------------------------------------------
while( $ensemble_member <= $ensemble_max )

   # Create a new temp directory for each member unless requested to keep and it exists already.
   set temp_dir = 'member'${ensemble_member}

   if(! -d $temp_dir) mkdir -p $temp_dir  || exit 1
   cd $temp_dir                           || exit 1

   # Get the program and necessary auxiliary files for the model
   ${LINK} ${CENTRALDIR}/MPAS_RUN/atmosphere_model     .         || exit 1
   ${LINK} ${CENTRALDIR}/MPAS_RUN/*BL                  .         || exit 1
   ${LINK} ${CENTRALDIR}/MPAS_RUN/*DATA                .         || exit 1
   ${LINK} ${CENTRALDIR}/MPAS_RUN/stream_list.atmosphere.* .     || exit 1
   ${LINK} ${CENTRALDIR}/advance_time                  .         || exit 1

   # Get the files specific for this experiment
   ${COPY} ${CENTRALDIR}/input.nml                     .         || exit 1
   ${LINK} ${CENTRALDIR}/streams.atmosphere            .         || exit 1
   ${COPY} ${CENTRALDIR}/namelist.atmosphere           .         || exit 1
   ${LINK} ${CENTRALDIR}/${fs_grid}*                   .	 || exit 1
   if( $if_sfc_update == .true. || $if_sfc_update == true ) then
       ${LINK} ${CENTRALDIR}/${fsfc} .
       ls -lL $fsfc						 || exit 1
   endif

   # Input analysis file
   set input_file = `head -n $ensemble_member ${CENTRALDIR}/${inlist}  | tail -1`
   set input_file = `basename $input_file`

   # Analysis time
   set anal_utc = `ncdump -v xtime $input_file | tail -2 | head -1 | cut -d";" -f1 | sed -e 's/"//g'`
   set    tanal = `echo ${anal_utc} 0 | ./advance_time`

   # Target forecast time (= next analysis time)
   set assim_days = `grep assimilation_period_days    input.nml | awk '{print $3}' | cut -d ',' -f1`
   set assim_secs = `grep assimilation_period_seconds input.nml | awk '{print $3}' | cut -d ',' -f1`
   set   targ_utc = `echo ${anal_utc} ${assim_days}d${assim_secs}s -w | ./advance_time`
   set      tfcst = `echo ${targ_utc} 0 | ./advance_time`

   @ addday     = $assim_secs / 86400
   @ modsec     = $assim_secs % 86400
   @ assim_days = $assim_days + $addday

   @ assim_hour = $assim_secs / 3600
   @ hoursec    = $assim_hour * 3600
   @ assim_secs = $assim_secs - $hoursec

   @ assim_min  = $assim_secs / 60
   @ minsec     = $assim_min * 60
   @ assim_secs = $assim_secs - $minsec

   # Forecast length
   set intv_utc = `echo $assim_days + 100 | bc | cut -b2-3`_`echo $assim_hour + 100 | bc | cut -b2-3`:`echo $assim_min + 100 | bc | cut -b2-3`:`echo $assim_secs + 100 | bc | cut -b2-3`

   # Prefix of each file name (for restart and diagnostics)
   set fhead = `sed -n '/<immutable_stream name=\"restart\"/,/\/>/{/Scree/{p;n};/##/{q};p}' streams.atmosphere | \
                    grep filename_template | awk -F= '{print $2}' | awk -F$ '{print $1}' | sed -e 's/"//g'`

   set fdiag = `sed -n '/<stream name=\"diagnostics\"/,/<\/stream>/{/Scree/{p;n};/##/{q};p}' streams.atmosphere | \
                    grep filename_template | awk -F= '{print $2}' | awk -F$ '{print $1}' | sed -e 's/"//g'`

   cat >! script.sed << EOF
   /config_start_time/c\
    config_start_time   = '${anal_utc}'
   /config_run_duration/c\
    config_run_duration = '${intv_utc}'
EOF


   echo  ${input_file}
   ls -l ${input_file}							|| exit
   if ( -e namelist.atmosphere )  ${REMOVE} namelist.atmosphere
   sed -f script.sed ${CENTRALDIR}/namelist.atmosphere >! namelist.atmosphere

   set is_it_regional = `grep config_apply_lbcs namelist.atmosphere | awk '{print $3}'`
   if ( $is_it_regional == true ) then
        echo Regional MPAS is not supported yet. Exit.
        exit
   endif

   # clean out any old log files
   if ( -e log.0000.out ) ${REMOVE} log.*

   # Run the model
   $mpicmd ./atmosphere_model

   # Check the output status
   ls -lrt >! list.${tanal}.txt
  
   # Model output at the target time
   set output_file = ${fhead}`echo ${targ_utc} | sed -e 's/:/\./g'`.nc
   set   diag_file = ${fdiag}`echo ${targ_utc} | sed -e 's/:/\./g'`.nc
   set date_utc = `ncdump -v xtime ${output_file} | tail -2 | head -1 | cut -d";" -f1 | sed -e 's/"//g'`

   # Check if the model was succefully completed.
   if($date_utc != $targ_utc) then
      echo $ensemble_member >>! ${CENTRALDIR}/blown.${tanal}_${tfcst}.out
      echo "Model failure! Check file " ${CENTRALDIR}/blown.${tanal}_${tfcst}.out
      exit 1
   endif

   if($save_wind == true) then
 
   set if_u_used = `grep use_u_for_wind input.nml | awk '{print $3}' | cut -d ',' -f1`
   if ( $if_u_used == .false. ) then
        ncks -O -v xtime,u ${input_file} analysis.uedge.${tanal}.nc
        ls -l analysis.uedge.${tanal}.nc
   else
        ncks -O -v xtime,uReconstructZonal,uReconstructMeridional ${input_file} analysis.uv.${tanal}.nc
        ls -l analysis.uv.${tanal}.nc
   endif

   endif	#($save_wind == true) then

   # Change back to the top directory.
   #-------------------------------------------------------------------
   cd $CENTRALDIR
   ls -l ${temp_dir}/${output_file}		|| exit
   echo ${temp_dir}/${output_file} >> list.${tfcst}.txt

   echo "Ensemble Member $ensemble_member completed"

   # Now repeat the entire process for other ensemble members
   #-------------------------------------------------------------------
   @ ensemble_member = $ensemble_member + 1

end

exit 0
