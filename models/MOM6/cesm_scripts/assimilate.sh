#!/usr/bin/env bash 

main() {

comp_name=OCN

get_cesm_info

get_obs_sequence

set_dart_input_nml

run_filter

cleanup

}



#-------------------------------
# Functions
#-------------------------------

#-------------------------------
# info needed from CESM
# This info is avaiable in the python case object
#  it would be good to have the case passed to _do_assimilate
#  in case_run.py
#-------------------------------
get_cesm_info() {
   
exeroot=$(./xmlquery EXEROOT --value) 
rundir=$(./xmlquery RUNDIR --value)
assimilate=$(./xmlquery DATA_ASSIMILATION_$comp_name --value)
ensemble_size=$(./xmlquery NINST_$comp_name --value)
   
}

#-------------------------------
# set filter input.nml options
#-------------------------------
set_dart_input_nml() {

echo "setting up input.nml for DART" 

}

#-------------------------------
# grab the observation sequence file
#-------------------------------
get_obs_sequence() {

echo "grab obs_seq.out"

}

#-------------------------------
# run filter
#-------------------------------
run_filter() {

if [ "$assimilate" = TRUE ]; then
   mpirun "$exeroot"/filter
fi

}

#-------------------------------
# run filter
#-------------------------------
cleanup() {

echo "stashing restart files for the next cycle"

}

#-------------------------------

main "@"
