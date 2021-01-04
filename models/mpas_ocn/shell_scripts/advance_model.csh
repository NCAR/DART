#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# Shell script to run the MPAS-A(tmostphere) model from DART input.
#
# This script performs the following:
# 1.  Creates a temporary directory to run a MPAS-A realization (see options)
# 2.  Copies or links the necessary files into the temporary directory
# 3.  Converts DART state vectors to mpas input
# 4.  Updates an MPAS namelist from a template withe new dates
# 5.  Runs MPAS-A model until target time
# 7.  Checks for incomplete runs
# 8.  Converts mpas output to DART state vectors
#
# Arguments are (created by 'filter' or 'perfect_model_obs' and include):
# 1) the process number of caller,
# 2) the number of ensemble members/state copies belonging to that process, and 
# 3) the name of the control_file for that process.
# 
# Note: For the required data to run this script, 
#       check the section of 'dependencies'.
#
# FIXME: As of September 30, 2011, it has not been decided if we will run 
#        MPAS-A as a restart mode or not. So this script was written
#        assuming that we take the model output (instead of the restart file)
#        at target time to use it as a background for the next filter. 
#        If we switch to use a restart file later, we should edit this script to have
#        config_restart_interval instead of config_output_interval, and 
#        check config_restart_name instead of config_output_name, and
#        set up config_do_restart = .true. for namelist.input.
#        We may also need to change the dart_to_model converter to write 
#        coupled variables (with dry density) in restart file.
#
# If this script finishes and the 'control_file' still exists, it is
# an ERROR CONDITION and means one or more of the ensemble members did
# not advance properly. Despite our best attempts to trap on this
# condition, some MPI installations simply hang, some properly terminate.
#
# This script loops over all the entries in the control_file to advance 
# any/all of the ensemble members.  The number of trips through the 
# loop is the second argument to this script. The control_file contains 
# the information about which ensemble members are to be advanced by THIS TASK.
# Sometimes it may be just one ensemble member, sometimes all of them.
# Read DART/doc/html/filter_async_modes.html and the mpi_intro.html
# for an overview.
#

set      process = $1
set   num_states = $2
set control_file = $3

#----------------------------------------------------------------------
# Block 1: copy necessary input files/executables/files common
#          to all model advances to a clean, temporary directory.
#          These will be used by ALL of the ensemble
#          members being advanced by this script.
#----------------------------------------------------------------------

# The run-time directory for the entire experiment is called CENTRALDIR;
set CENTRALDIR = `pwd`

# Create a clean temporary directory and go there
set delete_temp_dir = false

# set this to true if you want to maintain complete individual input/output 
# for each member (to carry through non-updated fields)
set individual_members = true

# next line ensures that the last cycle leaves everything in the temp dirs
if ( $individual_members == true ) set delete_temp_dir = false

set  REMOVE = '/bin/rm -rf'
set    COPY = '/bin/cp -p'
set    MOVE = '/bin/mv -f'
set    LINK = '/bin/ln -sf'
unalias cd
unalias ls

# if process 0 go ahead and check for dependencies here
if ( $process == 0 ) then

   foreach fn ( advance_time dart_to_model model_to_dart )
   if ( ! -x ${CENTRALDIR}/$fn ) then
     echo ABORT\: advance_model.csh could not find required executable dependency ${CENTRALDIR}/$fn
     exit 1
   endif
   end

   if ( ! -d ${CENTRALDIR}/MPAS_RUN ) then
      echo ABORT\: advance_model.csh could not find required data directory ${CENTRALDIR}/MPAS_RUN, 
      echo         which contains all the MPAS run-time input files
      exit 1
   endif

   if ( ! -x ${CENTRALDIR}/MPAS_RUN/nhyd_atmos_model.exe ) then
     echo ABORT\: advance_model.csh could not find required executable dependency 
     echo         ${CENTRALDIR}/MPAS_RUN/nhyd_atmos_model.exe
     exit 1
   endif

   if ( ! -r ${CENTRALDIR}/input.nml ) then
     echo ABORT\: advance_model.csh could not find required readable dependency ${CENTRALDIR}/input.nml
     exit 1
   endif

   if ( ! -r ${CENTRALDIR}/namelist.input ) then
     echo ABORT\: advance_model.csh could not find required readable dependency ${CENTRALDIR}/namelist.input
     exit 1
   endif

endif # process 0 dependency checking


# Loop through each state
set state_copy = 1
set ensemble_member_line = 1
set      input_file_line = 2
set     output_file_line = 3

while($state_copy <= $num_states)
   
   set ensemble_member = `head -$ensemble_member_line ${CENTRALDIR}/$control_file | tail -1`
   set input_file      = `head -$input_file_line      ${CENTRALDIR}/$control_file | tail -1`
   set output_file     = `head -$output_file_line     ${CENTRALDIR}/$control_file | tail -1`
   
   # Create a new temp directory for each member unless requested to keep and it exists already.
   set temp_dir = 'advance_temp'${ensemble_member}

   if ( $delete_temp_dir == "true" ) then
        if( -d $temp_dir ) ${REMOVE} $temp_dir || exit 1
   endif

   if(! -d $temp_dir) mkdir -p $temp_dir  || exit 1
   cd $temp_dir  || exit 1

   # Get the program and necessary files for the model
   ${LINK} ${CENTRALDIR}/MPAS_RUN/nhyd_atmos_model.exe .         || exit 1
   ${LINK} ${CENTRALDIR}/MPAS_RUN/*BL                  .         || exit 1
   ${LINK} ${CENTRALDIR}/advance_time                  .         || exit 1

   # Get the namelists
   ${COPY} ${CENTRALDIR}/input.nml      .                        || exit 1
   ${COPY} ${CENTRALDIR}/namelist.input namelist.input.template  || exit 1

   # Get the grid info files
   set fs_grid = `grep config_decomp_file_prefix namelist.input.template | awk '{print $3}' | sed -e "s/'//g"`
   set  f_grid = `basename $fs_grid .part.`
   ${LINK} ${CENTRALDIR}/MPAS_RUN/${fs_grid}* .
   ${LINK} ${CENTRALDIR}/MPAS_RUN/${f_grid} .

   # Get the in/out file names for converters and the model
   set f1 = `grep  dart_to_model_input_file input.nml | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g"`
   set f2 = `grep model_to_dart_output_file input.nml | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g"`
   set f3 = `grep   model_analysis_filename input.nml | awk '{print $3}' | cut -d ',' -f1 | sed -e "s/'//g"`

   # model output file name 
   set fhead = `grep config_output_name namelist.input.template | awk '{print $3}' | cut -d . -f1 | sed -e "s/'//g"`

   #----------------------------------------------------------------------
   # Block 2: move/convert the DART state vector to the model netcdf file.
   #----------------------------------------------------------------------
   if(! -e ${CENTRALDIR}/$f3) then
      echo ABORT\: ${CENTRALDIR}/$f3 does not exist.
      exit
   endif

   ${COPY} ${CENTRALDIR}/$f3 .           || exit 2
   ${MOVE} ${CENTRALDIR}/$input_file $f1 || exit 2
   ${CENTRALDIR}/dart_to_model >&! out.dart_to_model 

   # The program dart_to_model has created an ascii file named mpas_time.
   # Time information is extracted from the file.
   set curr_utc = `head -1 mpas_time | tail -1`
   set targ_utc = `head -2 mpas_time | tail -1`
   set intv_utc = `head -3 mpas_time | tail -1`

   #----------------------------------------------------------------------
   # Block 3: advance the model
   #          Make sure the file name is consistent in the namelist.input.
   #----------------------------------------------------------------------
   cat >! script.sed << EOF
   /config_start_time/c\
   config_start_time = '$curr_utc'
   /config_run_duration/c\
   config_run_duration ='$intv_utc'
   /config_output_interval/c\
   config_output_interval = '$intv_utc'
   /config_input_name/c\
   config_input_name = '$f3'
   /config_frames_per_outfile/c\
   config_frames_per_outfile = 1
EOF
# The EOF on the line above MUST REMAIN in column 1.

   sed -f script.sed namelist.input.template >! namelist.input

   # clean out any old rsl files
   if ( -e log.0000.out ) ${REMOVE} log.*

   # run MPAS here
   # mpi run on bluefire
   mpirun.lsf /usr/local/bin/launch ./nhyd_atmos_model.exe || exit 3
   # mpi run on Mac OS
   #mpiexec -n 4 ./nhyd_atmos_model.exe || exit 3
  
   # Model output at the target time
   set     fout = `ls -1 ${fhead}.*.nc | tail -1`
   set date_utc = `echo $fout | awk -F. '{print $(NF-1)}'`
   set targ_grg = `echo $date_utc 0 -g | advance_time`
   set targ_day = $targ_grg[1]
   set targ_sec = $targ_grg[2]

   # Check if the model was succefully completed.
   if($date_utc != $targ_utc) then
      echo $ensemble_member >>! ${CENTRALDIR}/blown_${targ_day}_${targ_sec}.out
      echo "Model failure! Check file " ${CENTRALDIR}/blown_${targ_day}_${targ_sec}.out
      exit 1
   endif

   #-------------------------------------------------------------------
   # Block 4: Convert your model output to a DART format ics file,
   #          then move it back to CENTRALDIR
   #-------------------------------------------------------------------
   ${MOVE} input.nml input.nml.template
   cat >! script1.sed << EOF
   /model_analysis_filename/c\
   model_analysis_filename      = '$fout'
EOF
   sed -f script1.sed input.nml.template >! input.nml
   ${CENTRALDIR}/model_to_dart >&! out.model_to_dart 
   ${MOVE} $f2 ${CENTRALDIR}/$output_file || exit 4

   # Update the template file for the next cycle 
   # only for $delete_temp_dir == false
   ${MOVE} $fout $f3 || exit 5

   # FIXME: Do we want to clean up the directory before moving up?
   ${REMOVE} ${fhead}.*.nc mpas_time log.* 

   # Change back to original directory.
   cd $CENTRALDIR

   # delete the temp directory for each member if desired
   # If all goes well, there should be no need to keep this directory.
   # If you are debugging, you may want to keep this directory. 

   if ( $delete_temp_dir == true )  ${REMOVE} $temp_dir
   echo "Ensemble Member $ensemble_member completed"

   # and now repeat the entire process for any other ensemble member that
   # needs to be advanced by this task.
   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line = $input_file_line + 3
   @ output_file_line = $output_file_line + 3

end

# MANDATORY - Remove the control_file to signal completion. If it still
# exists in CENTRALDIR after all the ensemble members have been advanced,
# it means one or more of the advances failed and is an ERROR CONDITION.

${REMOVE} $control_file

exit 0


