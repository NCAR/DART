#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# A set of examples on how to simply tweak existing netCDF files.

set fname = input_file.nc

# The ON2_UNCERTAINTY variable in the netcdf files have IEEE NaN values,
# but none of the required metadata to interpret them correctly.
# These 2 lines will add the required attributes so that NaNs are replaced with
# a fill value that can be queried and checked for.
# Since the ON2_UNCERTAINTY is a standard deviation, it is enough to make it negative
   
ncatted -a _FillValue,ON2_UNCERTAINTY,o,f,NaN        input_file.nc
ncatted -a _FillValue,ON2_UNCERTAINTY,m,f,-1.0       input_file.nc

exit 0


