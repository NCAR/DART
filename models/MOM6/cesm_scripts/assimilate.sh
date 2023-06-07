#!/usr/bin/env bash 

# CESM caseroot is the first argument to the script

main() {

set -ex

caseroot=$1

dart_build_dir=/glade/scratch/hkershaw/DART/MOM6/DART/models/MOM6/work
comp_name=OCN
obs_dir=/glade/scratch/hkershaw/DART/MOM6/Observations/201502

echo "DART dart_build_dir" $dart_build_dir

cd $caseroot
get_cesm_info

cd $rundir
get_obs_sequence

setup_dart

setup_dart_input_nml

setup_template_files

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
case=$(./xmlquery CASE --value)
   
}


#-------------------------------
# Setup dart executables
#  This step should go in case.setup/case.build
#-------------------------------
setup_dart() {

# Should these checks just be for CONTINUE_RUN?
[ ! -f "$exeroot"/filter ] || cp $dart_build_dir/filter  $exeroot


}

#-------------------------------
# set filter input.nml options
#-------------------------------
setup_dart_input_nml() {

echo "setting up input.nml for DART" 

# list restart files
ls $case.mom6.r* > filter_input_list.txt
cp filter_input_list.txt  filter_output_list.txt

# What about other input.nml in the rundir?
cp $dart_build_dir/input.nml $rundir/input.nml
}

#-------------------------------
# set template files for filter
# restart, static filenames depend on the case
# ocean_geometry.nc is always called ocean_geometry.nc
#-------------------------------
setup_template_files() {

ln -sf $(head -1 filter_input_list.txt) mom6.r.nc

ln -sf $(ls $case.mom6.static* | head -1) mom6.static.nc
}

#-------------------------------
# grab the observation sequence file
#-------------------------------
get_obs_sequence() {

echo "grab obs_seq.out"

ln -sf $obs_dir/obs_seq.0Z.20150206 obs_seq.out
}

#-------------------------------
# run filter
#-------------------------------
run_filter() {

echo "running filter"
if [ "$assimilate" = TRUE ]; then
   mpirun -n 512 "$exeroot"/filter
fi

}

#-------------------------------
# cleanup
#-------------------------------
cleanup() {

echo "stashing restart files for the next cycle"

}

#-------------------------------

main "$@"
