#!/bin/sh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#
# take a set of year-long WOD files and convert 
# the contents into a single year-long obs_seq file.
#
# another script will split it into daily files.
#

# set the first and last years.
let start_year=1998
let end_year=2005

EXEDIR='../work'
DATDIR='../data'

# end of things you should have to set in this script

let totalyears=${end_year}-${start_year}+1


# loop over each year
let y=1
let year=$start_year

while (( y <= totalyears )) ; do

  # status/debug - comment in or out as desired.
  echo starting processing for ${year} 

  ls -1 ${DATDIR}/*/*/*${year} > ${year}list

  sed -e "s/YYYY/$year/" input.nml.template > input.nml

  ${EXEDIR}/wod_to_obs

  # advance the loop counter
  let y=y+1
  let year=year+1
 
done

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

