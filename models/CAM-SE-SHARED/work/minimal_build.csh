#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# Test whether the default MPI (MPICH2) is necessary for filter to hang.
#
# setenv MP_MPILIB pempi

# \rm -f ../../../obs_def/obs_def_mod.f90
# \rm -f ../../../obs_kind/obs_kind_mod.f90
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

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

