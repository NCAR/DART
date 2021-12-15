#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

set HDFDIR = /Users/thoar/intel_16.0.0
set INCLUDES = "-I${HDFDIR}/include"
set LIBRARIES = "-L${HDFDIR}/lib -lhdf5 -lhdf5_fortran"
set FFLAGS = '-O0'
set FFLAGS = '-g -C -check noarg_temp_created -fpe0 -fp-model precise \
              -ftrapuv -traceback -warn declarations,uncalled,unused'

\rm readhdf5 readhdf5lt

# readhdf5 does not work on my OSX software stack 
ifort readhdf5.f90 -o readhdf5 ${FFLAGS} ${INCLUDES} ${LIBRARIES}


set LIBRARIES = "-L${HDFDIR}/lib -lhdf5_hl -lhdf5 -lhdf5hl_fortran -lhdf5_fortran"

# Use for testing conflicts with the netcdf libs - more realistic
# set LIBRARIES = "-L${HDFDIR}/lib -lhdf5_hl -lhdf5 -lhdf5hl_fortran -lhdf5_fortran -lnetcdff -lnetcdf"

# readhdf5lt WORKS on my OSX software stack 
ifort readhdf5lt.f90 -o readhdf5lt ${FFLAGS} ${INCLUDES} ${LIBRARIES}

exit $status


