#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# Test whether the default MPI (MPICH2) is necessary for filter to hang.
#
# setenv MP_MPILIB pempi

# ./preprocess || exit 99
# 
# csh ./mkmf_cam_to_dart
# make
# csh ./mkmf_dart_to_cam
# make
csh ./mkmf_filter -mpi
make

# set echo verbose
echo '=============================================='
if ($#argv == 1) then
   if (! -d $1) mkdir $1
   cp filter $1  || exit 19

   cd $1:h
   if (-f $1:h/filter) then
      rm filter 
   endif
   ln -s $1:t/* .
else
   echo "Not copying executables anywhere"
endif

exit 0


