#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# Log results of all quickbuilds and compilers
LOGFILE="all_quickbuilds_results.log"
# Log results of failed builds only
FAILURE_LOGFILE="all_quickbuilds_failures.log"

if [ $# -ne 1 ]; then
   echo "ERROR: expecting 1 argument"
   exit 1
fi

# Check if the script is running in a batch job or interactive session
if [[ -z $PBS_ENVIRONMENT ]]; then
  echo "ERROR: You must run this in a batch job"
  echo "       qsub submit_me.sh"
  echo "       or an interactive session"
  exit 2
fi

# Include current date/time to log files
current_date="$(printf '\n%s\n' "$(date)")"
echo "$current_date" > "$LOGFILE"
echo "$current_date" > "$FAILURE_LOGFILE"

# Specify the mkmf template for each compiler
if [[ $FC == "intel" ]]; then
  template_name="mkmf.template.intel.linux"
elif [[ $FC == "gcc" ]]; then
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

# Log the compiler being processed
printf '\nProcessing %s\n' "$FC"

module load $FC

cd $DART
cp build_templates/$template_name build_templates/mkmf.template

# Run fixsystem to avoid all make commands from altering mpi_*_utilities_mod.f90 simultaneously
cd assimilation_code/modules/utilities; ./fixsystem $FC
cd -

# Build preprocess once
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
    if [[ $FC == "gcc" ]]; then  # HK this should be needed for nvhpc too
        cd $f; ./quickbuild.sh mpif08 &
    else
        cd $f; ./quickbuild.sh &
    fi

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
