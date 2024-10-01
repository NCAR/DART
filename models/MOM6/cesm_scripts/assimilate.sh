#!/usr/bin/env bash 

# CESM caseroot is the first argument to the script

main() {

set -ex

caseroot=$1

dart_build_dir=/glade/derecho/scratch/hkershaw/MOM6/DART/models/MOM6/work
comp_name=OCN
obs_dir=/glade/p/cisl/dares/Observations/WOD13
ntasks=1152  # should pull this from CIME

echo "DART dart_build_dir" $dart_build_dir

cd $caseroot
get_cesm_info

cd $rundir

get_model_time_from_filename

get_obs_sequence

setup_dart

setup_dart_input_nml

setup_template_files

run_filter

date_stamp_output

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
# get the model time from the filename
#  used to get the correct obs_seq.out
#-------------------------------
get_model_time_from_filename() {

pointer=rpointer.ocn_0001
filename=$(head -n 1 $pointer)

date1=${filename##*.r.}  # cut off case.mom6.r.
date="${date1%._*}"      # cut off _ensemble member.nc

array=(${date//-/ })

year=${array[0]}
month=${array[1]}
day=${array[2]}
seconds=${array[3]}

echo "model time is $year $month $day $seconds"

YYYYMMDD=${year}${month}${day}
YYYYMM=${year}${month}
obs_filename=obs_seq.0Z.${YYYYMMDD}
obs_file=${obs_dir}/${YYYYMM}/${obs_filename}

}


#-------------------------------
# Setup dart executables
#  This step should go in case.setup/case.build
#-------------------------------
setup_dart() {

# Should this check just be for CONTINUE_RUN?
if [ ! -f "$exeroot"/filter ]
then
  cp $dart_build_dir/filter  $exeroot
fi

}

#-------------------------------
# set filter input.nml options
#-------------------------------
setup_dart_input_nml() {

echo "setting up input.nml for DART" 

# list restart files

cat rpointer.ocn_???? > filter_input_list.txt
cp filter_input_list.txt  filter_output_list.txt

# Store MOM6 nml so it is not overwritten by dart input.nml
mv input.nml input.nml.mom6
cp $dart_build_dir/input.nml $rundir/input.nml
}

#-------------------------------
# set template files for filter
# restart, static filenames depend on the case
# ocean_geometry.nc is always called ocean_geometry.nc
#-------------------------------
setup_template_files() {

ln -sf $(head -1 filter_input_list.txt) mom6.r.nc

ln -sf $(ls $case.mom6.h.static* | head -1) mom6.static.nc
}

#-------------------------------
# grab the observation sequence file
#-------------------------------
get_obs_sequence() {

echo "grab obs_seq.out" $obs_file

ln -sf $obs_file obs_seq.out
}

#-------------------------------
# run filter
#-------------------------------
run_filter() {

echo "running filter"
if [ "$assimilate" = TRUE ]; then
   mpirun -n $ntasks "$exeroot"/filter
fi

}

#-------------------------------
# append timestamp to filter output
#-------------------------------
date_stamp_output() {

mv obs_seq.final "$case".obs_seq.final.${YYYYMMDD}
mv dart_log.out "$case".dart_log.out.${YYYYMMDD}
mv dart_log.nml "$case".dart_log.nml.${YYYYMMDD}

# possible dart output files:
netcdf=(\
"preassim_mean.nc"           "preassim_sd.nc" \
"preassim_priorinf_mean.nc"  "preassim_priorinf_sd.nc" \
"preassim_postinf_mean.nc"   "preassim_postinf_sd.nc" \
"postassim_mean.nc"          "postassim_sd.nc" \
"postassim_priorinf_mean.nc" "postassim_priorinf_sd.nc" \
"postassim_postinf_mean.nc"  "postassim_postinf_sd.nc" \
"output_mean.nc"             "output_sd.nc" \
"output_priorinf_mean.nc"    "output_priorinf_sd.nc" \
"output_postinf_mean.nc"     "output_postinf_sd.nc")

for file in ${netcdf[@]}; do
  if [ -f $file ]; then
     mv "$file" "$case"."$file".${YYYYMMDD}
  fi
done

}


#-------------------------------
# cleanup
# restore mom6 input.nml for next cycle
#-------------------------------
cleanup() {

echo "TODO: stashing restart files for the next cycle"
mv input.nml.mom6 input.nml

}

#-------------------------------

main "$@"
