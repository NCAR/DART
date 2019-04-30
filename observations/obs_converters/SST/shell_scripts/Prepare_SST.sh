#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# This script is based on code donated to DART by Romain Escudier
# who was at Rutgers at the time. Thanks Romain!
#
# DART $Id$

. ./parameters_SST
. functions.sh

# module purge
# module load netcdf/4.3.0-gcc4.4.7
# module load gsl/2.1-gcc4.4.7 nco/4.6.1

echo " "
echo "generating netcdf files in : ${DIR_OUT}"
echo "generating dart   files in : ${DIR_OUT_DART}"
echo "representation error is ${ERROR_REP}"
echo " "

# Get name of mask
mask_file=../masks/Mask_${REGION}${FILE_SUFF}.nc

# Create directory for outputs
mkdir -p ${DIR_OUT_DART}

# Get the number of days to prepare
n_days=$(get_timediff_dates ${STARTDATE} ${ENDDATE})

# Loop on days
for i_day in $( seq 0 $(($n_days)) ) ; do

    # Get date of this loop iteration
    my_date=$(get_date_from_cycle ${i_day} $STARTDATE 1)
    my_year=${my_date:0:4} # Year
    echo "Preparing date : ${my_date}"
    mkdir -p ${DIR_OUT}/${my_year}

    # Region selection
    file_in=${DIR_IN}/${my_year}/${FILE_PREF}${my_date}${FILE_SUFF}.nc
    file_out=${DIR_OUT}/${my_year}/${FILE_PREF}${my_date}${FILE_SUFF}_${REGION}.nc

    if [ -f ${file_out} ]; then
       echo "WARNING: Region Selection output file already exists."
       echo "WARNING: ${file_out}"
       echo "WARNING: not subsetting again."
    else
       echo "input     netcdf file in : ${file_in}"
       echo "subsetted netcdf file in : ${file_out}"
       echo "using mask               : ${mask_file}"

       ./prepare_SST_file_NWA.sh ${file_in} ${file_out} ${mask_file} || exit 1

       echo "Finished prepare_SST_file_NWA.sh"
    fi

    echo " "

    # Transform to DART format
    file_dart_out=${DIR_OUT_DART}/${FILE_DART_PREF}${my_date}${FILE_DART_SUFF}
    if [ -f ${file_dart_out} ]; then
       echo "WARNING: DART observation sequence file already exists."
       echo "WARNING: ${file_dart_out}"
       echo "WARNING: Not creating it again."
    else
       echo "creating ${file_dart_out}"

       sed -e "s;<SST_IN_FILE>;${file_out};g" \
           -e "s;<SST_OUT_FILE>;${file_dart_out};g" \
           -e "s;<SST_TO_OBS_LOG>;dart_sst2seq_err${ERROR_REP}_${my_date}.out;g" \
           -e "s;<SST_TO_OBS_NML>;temp.nml;g" \
           -e "s;<ERROR_REP>;${ERROR_REP};g" input.nml.template > input.nml

       ../work/sst_to_obs || exit 2
       rm -f temp.nml
    fi
    echo " "
done

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

