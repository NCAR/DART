#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# Usage: run_all_quickbuilds.sh compiler [gcc intel nvhpc]

if [ $# -ne 1 ]; then
   echo "ERROR: expecting 1 argument"
   exit 1
fi

compiler=$1

# Check if the script is running in a batch job or interactive session
if [[ -z $PBS_ENVIRONMENT ]]; then
  echo "ERROR: You must run this in a batch job"
  echo "       qsub submit_me.sh"
  echo "       or an interactive session"
  exit 2
fi

# Specify the mkmf template for each compiler
if [[ $compiler == "intel" ]]; then
  mkmf_template="mkmf.template.intel.linux"
elif [[ $compiler == "gcc" ]]; then
  mkmf_template="mkmf.template.gfortran"
elif [[ $compiler == "cce" ]]; then
  mkmf_template="mkmf.template.cce"
elif [[ $compiler == "nvhpc" ]]; then
  mkmf_template="mkmf.template.nvhpc"
else
  echo "$compiler is not a valid argument"
  exit 1
fi
  

test_dir="/glade/derecho/scratch/$USER/build_everything/$compiler"
if [[ -d $test_dir ]]; then
  echo "Directory exists: $test_dir"
  exit 2
fi

mkdir -p $test_dir
cd $test_dir
git clone 'https://github.com/NCAR/DART.git'
cd DART
export DART=$(git rev-parse --show-toplevel)

# mkmf for chosen compiler
module load $compiler
cp build_templates/$mkmf_template build_templates/mkmf.template

# Run fixsystem to avoid all make commands from altering mpi_*_utilities_mod.f90 simultaneously
cd assimilation_code/modules/utilities; ./fixsystem $compiler
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
    if [[ $compiler == "gcc" ]]; then
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
    log_msg=$(printf '%s\n' "$compiler RESULT: $i ${dirs[$i]} FAILED")
    echo "$log_msg"
    # Indicate at least one failure
    OVERALL_EXIT=1
  else
    echo "$compiler RESULT: $i ${dirs[$i]} PASSED"
  fi
  ((i+=1))
done

mv $test_dir  $test_dir.$(date +"%FT%H%M")
