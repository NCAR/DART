#!/bin/bash
#
# This code was donated to DART by Romain Escudier, who was at Rutgers at the time.
# It is not protected by the DART copyright agreement. Thanks Romain!
#
# DART $Id$

#---------------------------------------------------------------------------------------------#
# Get parameter value from netCDF file

get_param_from_nc() { FILENAME=$1 ; VARNAME=$2
   VARVALUE=$(ncdump -v ${VARNAME} ${FILENAME} | awk '{FS="data:"; RS="ceci1est8une9valeur5impossible"; print $2}');
   VARVALUE=${VARVALUE#*=}; VARVALUE=${VARVALUE%;*\}}
   echo ${VARVALUE}
}

#---------------------------------------------------------------------------------------------#
# Main program

set -e

SST_file_in=$1
SST_file_out=$2
mask_NWA=$3

if [ ! -f ${SST_file_out} ]; then

   #module load gsl/2.1-gcc4.4.7 nco/4.6.1 2> /dev/null
   # Create output file
   cp ${SST_file_in} ${SST_file_out}

   # Get time from infile (will be erased by append in ncks)
   my_time=$(get_param_from_nc ${SST_file_out} time)

   # Append the ROMS mask to output file
   ncks -A -v mask_roms ${mask_NWA} ${SST_file_out}

   ## Compute the composite mask
   ncap2 -O -s "mask=mask_roms*mask" ${SST_file_out} ${SST_file_out}

   # Select area
   ncks -O -x -v mask_roms -d lon,-102.66,-43.81 -d lat,7.78,54.83 ${SST_file_out} ${SST_file_out}

   # Put the correct time
   ncap2 -O -s "time[time]=time-time+${my_time}" ${SST_file_out} ${SST_file_out}

fi

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

