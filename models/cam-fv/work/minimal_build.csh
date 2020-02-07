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

# Set modules to be consistent with CESM2_1's.
# module purge
# module restore cesm2_1_0

# \rm -f ../../../obs_def/obs_def_mod.f90
# \rm -f ../../../obs_kind/obs_kind_mod.f90
# ./preprocess || exit 99

# csh ./mkmf_obs_diag
# make
#  >> Need to fix the mpi aspect of building mpi and nonmpi programs.
# csh ./mkmf_fill_inflation_restart
# Trying to follow Tim's 2018-12-31 advice to build everything with mpif90:
# csh ./mkmf_fill_inflation_restart -mpi
# make

csh ./mkmf_filter -mpi
make

# set echo verbose
echo '=============================================='
if ($#argv == 1) then
   if (! -d $1) mkdir $1
   cp filter $1  || exit 19
#    cp fill_inflation_restart $1  || exit 20

   cd $1:h
   if (-f $1:h/filter) then
      rm filter 
   endif
#    if (-f $1:h/fill_inflation_restart) then
#       rm fill_inflation_restart 
#    endif
   ln -s $1:t/* .
else
   echo "Not copying executables anywhere"
endif

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

