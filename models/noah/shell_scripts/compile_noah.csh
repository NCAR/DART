#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
 

gfortran  -ffree-form -ffree-line-length-none -g -O2 -O0 -fconvert=big-endian \
          -ffpe-trap=invalid,zero,overflow,underflow \
          -I/Users/thoar/gnu_4.5.0/include -D__GFORTRAN__ -D_SIMPLE_DRIVER_ \
          -L/Users/thoar/gnu_4.5.0/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lm \
          -o noah_1d_driver \
          ../src/module_model_constants.F90 \
          ../src/module_sfcdif_wrf.F90  \
          ../src/module_sf_noahlsm.F90 \
          ../src/module_Noahlsm_utility.F90 \
          ../src/module_netcdf_io.F90 \
          ../src/kwm_date_utilities.F90 \
          ../src/module_ascii_io.F90 \
          ../src/module_io.F90 \
          ../src/simple_driver.F90 

exit 0


