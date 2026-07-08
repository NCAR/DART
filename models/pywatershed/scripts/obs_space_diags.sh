#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id: $
#
#-------------------------------------------------------------------------------
# DESCRIPTION:
#
# This script is run after submitting the output summary script: "DART_cleanup.csh" or "collect_output_sched.sh" 
# These scripts summarize the DART diagnostic files in output directery into large netcdf files (run_dir/all_*)
# ...........
# The following script first creates a diagnostics directory in your scratch. The output from 
# "DART_cleanup.csh" are moved. 
# DART executables "obs_diag" and "obs_seq_to_netcdf" from the work directory are linked to the diagnostics directory. 
# The input.nml is edited such that the obs_seq.final files are listed in a file for 
# the intended period. 

# A matlab script "run_HydroDARTdiags.m" can then be run to display the results.     
# 
#--------------------------------------------------------------------------------
#

if [[ $# -lt 3 ]]; then 
   echo "Please enter the following:"  
   echo "(1) RUN Direcotry" 
   echo "(2) Diag Directory"
   echo "(3) First date for obs diag: e.g., 2018-08-15_00" 
   echo "(4) Last date for obs_diag: e.g., 2018-09-20_00"
   echo "Usage: ./obs_space_diags.sh RUN_DIR DIAG_DIR 2018-09-07_00 2018-10-07_23"
   exit 0
fi

RUN_DIR=$1
DIAG_DIR=$2
home_dir=/PATH/TO/DART/models/pywatershed/work/
first_year=`echo  $3 | cut -c 1-4`
first_month=`echo $3 | cut -c 6-7`
first_day=`echo   $3 | cut -c 9-10`
first_hour=`echo  $3 | cut -c 12-13`

last_year=`echo   $4 | cut -c 1-4`
last_month=`echo  $4 | cut -c 6-7`
last_day=`echo    $4 | cut -c 9-10`
last_hour=`echo   $4 | cut -c 12-13`

mkdir $DIAG_DIR 

if [[ -e ${DIAG_DIR}/obs_diag_output.nc ]]; then 
   echo "An obs_diag_output.nc file already exists in ${scratch_dir}/diagnostics/${which_dir} .."
   read -p "Do you wish to continue and overwrite the current file? " -n 1 -r
   echo 

   if [[ ! $REPLY =~ ^[Yy]$ ]]; then
      echo "Exiting .."
      exit 1
   fi
fi

ln -sf ${home_dir}/streamflow_obs_diag          ${DIAG_DIR}
ln -sf ${home_dir}/obs_seq_to_netcdf ${DIAG_DIR}

rm -f ${DIAG_DIR}/dart_log.out   || true
rm -f ${DIAG_DIR}/dart_log.nml   || true 
rm -f ${DIAG_DIR}/LargeInnov.txt || true
rm -f ${DIAG_DIR}/list_2_diag    || true

mv ${RUN_DIR}/all_*              ${DIAG_DIR} || echo "Error moving the summarized DART fils exiting"

# Add the list of obs-seq-final files to the input.nml
sed -i "s#.*obs_sequence_list.*#    obs_sequence_list = ''#" ${RUN_DIR}/input.nml
ln -sf ${RUN_DIR}/input.nml ${DIAG_DIR}/input.nml

cd ${DIAG_DIR}

for data_dir in `ls ${RUN_DIR}/output`; do
    echo 'Adding obs_seq file at time: ' ${data_dir}
    ls -d -1 ${RUN_DIR}/output/${data_dir}/obs_seq.final* >> list_2_diag

    sed -i "s#.*obs_sequence_name.*#    obs_sequence_name = '${RUN_DIR}/output/${data_dir}/obs_seq.final.${data_dir}'#" input.nml
done

sed -i "s#.*obs_sequence_name.*#    obs_sequence_name = ''#" input.nml
sed -i "s#.*obs_sequence_list.*#    obs_sequence_list = '${DIAG_DIR}/list_2_diag'#" input.nml

# Adjust the beginning and end dates for diagnosing
sed -i "s#.*first_bin_center.*#    first_bin_center = $first_year, $first_month, $first_day, $first_hour, 0, 0#" input.nml
sed -i "s#.*last_bin_center.*#    last_bin_center = $last_year, $last_month, $last_day, $last_hour, 0, 0#" input.nml

# Regions to show
sed -i "s#.*nregions.*#    nregions = 1#" input.nml
sed -i "s#.*reg_names.*#    reg_names = '${exp_domain}'#" input.nml

# hlevel edges
sed -i "s#.*hlevel.*#    hlevel_edges = -100, 10000000#" input.nml

# Make sure we add all obs to the epoch file
sed -i "s#.*append_to_netcdf.*#   append_to_netcdf  = .true.#" input.nml

# Run obs_diag and obs_seq_to_netcdf ..
./streamflow_obs_diag 
./obs_seq_to_netcdf

rm dart_log.out dart_log.nml LargeInnov.txt || true

# Lunch MATLAB and display the results 
echo "You may now launch matlab 'matlab -nodesktop' and change your directory to: ${home_dir}/../matlab"
echo "In order to view the results, you need to edit/run: run_HydroDARTdiags.m"

# <next few lines under version control, do not edit>
# $URL: $
# $Revision: $
# $Date: $

