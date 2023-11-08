#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download


# Email build results
TO_EMAIL_ADDRESSES='anorcio@ucar.edu'
# Log results of all quickbuilds and compilers
LOGFILE="all_quickbuilds_results.log"
# Log results of failed builds only
FAILURE_LOGFILE="all_quickbuilds_failures.log"

# Store compiler arguments in an array
FCS=( $@ )

# Check if the script is running in a batch job or interactive session
if [[ -z $PBS_ENVIRONMENT ]]; then
  echo "ERROR: You must run this in a batch job"
  echo "       qsub submit_me.sh"
  echo "       or an interactive session"
  exit 2
fi

# If no argument is provided, run with all compilers as default
if [[ "${#FCS[@]}" == 0 ]]; then
  FCS=( intel gnu cce nvhpc )
fi

# Include current date/time to log files
current_date="$(printf '\n%s\n' "$(date)")"
echo "$current_date" > "$LOGFILE"
echo "$current_date" > "$FAILURE_LOGFILE"

# Function to run quickbuild.sh for specific compilers
function compiler_operation {
  FC="$1"
  # Set up DART environment variable and do a cleanup in case DART already exists
  export DART="/glade/derecho/scratch/$USER/tmp/$FC/DART"
  if [[ -n $FC ]]; then
    rm -rf /glade/derecho/scratch/$USER/tmp/$FC
  fi
  # Create new tmp directory for a specific compiler and clone DART
  mkdir /glade/derecho/scratch/$USER/tmp/"$FC"
  git clone 'https://github.com/NCAR/DART.git' "$DART"

  # Check if DART directory exists. If not, return an error
  if [[ ! -d $DART ]]; then 
    echo "No DART directory: $DART"
    return 1
  fi

  # Specify the mkmf template for each compiler
  if [[ $FC == "intel" ]]; then
    template_name="mkmf.template.intel.linux"
  elif [[ $FC == "gnu" ]]; then
    template_name="mkmf.template.gfortran"
  elif [[ $FC == "cce" ]]; then
    template_name="mkmf.template.cce"
  elif [[ $FC == "nvhpc" ]]; then
    template_name="mkmf.template.nvhpc"
  else
    # Return an error for invalid arguments
    echo "$FC is not a valid argument. Compiler is not supported on Derecho."
    return 1
  fi

  # Log the compiler being processed
  printf '\nProcessing %s\n' "$FC"

  # Copy appropriate mkmf template to build_templates directory in DART
  cp "$DART/build_templates/$template_name" $DART/build_templates/mkmf.template

  # cd into DART. Note that fixsystem, preprocess, and modify input.nml were taken from the current version of build_everything/run_all_quickbuilds.sh  
  cd $DART

  # Run fixsystem to avoid all make commands from altering mpi_*_utilities_mod.f90 simultaneously
  cd assimilation_code/modules/utilities; ./fixsystem $FC
  cd -

  # Build preprocess once for the given compiler
  pp_dir=$DART/assimilation_code/programs/preprocess
  cd $pp_dir
  $DART/build_templates/mkmf -x -p $pp_dir/preprocess \
        -a $DART $pp_dir/path_names_preprocess
  cd -

  # Modify input.nml files to use local versions of obs_def_mod.f90 and obs_kind_mod.f90
  find . -name input.nml -exec sed -i -e "/^[[:space:]]*#/! s|.*output_obs_def_mod_file.*|output_obs_def_mod_file = './obs_def_mod.f90'|g" \
        -e "/^[[:space:]]*#/! s|.*output_obs_qty_mod_file.*|output_obs_qty_mod_file = './obs_kind_mod.f90'|g" \
        -e "/^[[:space:]]*#/! s|.*output_obs_kind_mod_file.*|output_obs_qty_mod_file = './obs_kind_mod.f90'|g" {} \;

  # Store the current directory and initialize arrays to hold process IDs, directories, and status codes 
  my_dir=$(pwd)
  pids=()
  dirs=()
  status=()

  # Find all quickbuild.sh executable files and remove './' and 'quickbuild.sh'
  files_to_process=( $(find $DART -executable -type f -name quickbuild.sh | sed -E 's#(\./|quickbuild\.sh)##g') )

  # Iterate over each file to and run quickbuild.sh
  for f in "${files_to_process[@]}"; do
    # cd to */quickbuild.sh file and run quickbuild.sh in the background
    cd $f; ./quickbuild.sh & 
    # Record the PID and directory of the each process then cd back to starting directory
    pids+=( "$!" )
    dirs+=( "$f" )
    cd $my_dir
  done

  # Wait for all background processes to finish and record their exit statuses
  for pid in ${pids[@]}; do
    wait ${pid}
    status+=( "$?" )
  done

  # Check the status of each build process and log results
  OVERALL_EXIT=0
  i=0
  for st in ${status[@]}; do
    # Display failed vs. passed processes
    if [[ ${st} -ne 0 ]]; then
      log_msg=$(printf '%s\n' "$FC RESULT: $i ${dirs[$i]} FAILED")
      echo "$log_msg"
      # Indicate at least one failure
      OVERALL_EXIT=1
    else
      echo "$FC RESULT: $i ${dirs[$i]} PASSED"
    fi
    ((i+=1))
  done

  return "$OVERALL_EXIT"
}

# Run compiler_operation for each compiler in the background and record each PID
background_processes=( )
for FC in "${FCS[@]}"; do
  compiler_operation "$FC" &
  background_processes+=( "$!" )
done &>"$LOGFILE"

# Wait for all compiler_operation to finish and record exit statuses
exit_statuses=( )
for process in "${background_processes[@]}"; do
  wait "$process"
  exit_statuses+=( "$?" )
done &>"$LOGFILE"

# Filter and append all failed builds to all_quickbuilds_failures.log
grep -aE '[[:blank:]]+FAILED[[:blank:]]*$' "$LOGFILE" >> "$FAILURE_LOGFILE"

# Email build results with 'all_quickbuilds_results.log' and 'all_quickbuilds_failures.log' as file attachments
mail_subject=$(printf 'TEST: Run all quickbuilds results: %s' "$(date)")
echo "$mail_subject" | mail -s "$mail_subject" -a "$LOGFILE" -a "$FAILURE_LOGFILE" $TO_EMAIL_ADDRESSES

# Teardown: Remove DART directory
for FC in "${FCS[@]}"; do
  if [[ -n $FC ]]; then
    rm -rf /glade/derecho/scratch/$USER/tmp/"$FC"
  fi
done

# For testing purposes, indicate when script is done
echo "Done"

# Check for any failed operations
for exit_status in "${exit_statuses[@]}"; do
  if [[ $exit_status -ne 0 ]]; then
    exit 1
  fi
done


exit 0
