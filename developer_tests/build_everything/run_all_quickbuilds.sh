#!/usr/bin/env bash

# options to set
#   1. DART directory
#   2. compiler

export DART=/glade/scratch/hkershaw/DART/Compile_all/DART
FC=intel


if [[ ! -d $DART ]] ; then 
  echo "No DART directory: " $DART
  exit 1
fi 

cp mkmf.template.compile_all $DART/build_templates/mkmf.template
cd $DART

# run fixsystem 
# running this once at the beginning otherwise all make commands will
# try and alter the mpi_*_utilities_mod.f90 simultaneously
cd assimilation_code/modules/utilities; ./fixsystem $FC
cd -


# build preprocess once
pp_dir=$DART/assimilation_code/programs/preprocess
cd $pp_dir
$DART/build_templates/mkmf -x -p $pp_dir/preprocess \
      -a $DART $pp_dir/path_names_preprocess
cd -

# local versions of obs_def_mod.f90 and obs_kind_mod.f90
find . -name input.nml -exec sed -i '' -e "/^[[:space:]]*#/! s|.*output_obs_def_mod_file.*|output_obs_def_mod_file = './obs_def_mod.f90'|g" \
      -e "/^[[:space:]]*#/! s|.*output_obs_qty_mod_file.*|output_obs_qty_mod_file = './obs_kind_mod.f90'|g" \
      -e "/^[[:space:]]*#/! s|.*output_obs_kind_mod_file.*|output_obs_qty_mod_file = './obs_kind_mod.f90'|g" {} \;  

my_dir=$(pwd)
pids=()
dirs=()
status=()

while read f; do

cd $f; ./quickbuild.sh & 
pids+=( "$!" )
dirs+=( "$f" )
cd $my_dir

done < $DART/developer_tests/build_everything/all_quickbuilds

for pid in ${pids[@]}; do
  #echo "${pid}"
  wait ${pid}
  status+=( "$?" )
done

# looping through the status arr to check exit code for each
i=0
for st in ${status[@]}; do
    if [[ ${st} -ne 0 ]]; then
        echo "$i ${dirs[$i]} failed"
        OVERALL_EXIT=1
    else
        echo "$i  ${dirs[$i]} finished"
    fi
    ((i+=1))
done

