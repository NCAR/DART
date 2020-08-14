#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# this does not actually generate wrf initial condition files; but once
# those exist, they have to be converted from wrf netcdf files to dart
# restart files.  that's what this script does - loop over all 50 wrf
# input netcdf files and convert them to dart state vector files, with
# the name 'filter_ic_old.00nn'.  it also copies one of the input files
# to the generic name 'wrfinput_d01', which the model_mod needs when
# it initializes - it doesn't use data from this file, but it needs to
# read in the grid information and number of domains from this file.

set i=1
while ($i <= 50)
 rm -f wrfinput_d01
 ln -s wrfinput_d01_148403_0_${i} wrfinput_d01
 ./wrf_to_dart
 set n=`printf "%04d"  $i`
 mv dart_wrf_vector filter_ic_old.$n
 @ i++
end

cp wrfinput_d01_148403_0_1 wrfinput_d01

exit 0


