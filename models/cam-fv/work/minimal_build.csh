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

set cdebug = 0

if ($#argv == 0) then
   echo "Usage: ./minimal_build.csh output_dir"
   echo '       if output_dir = multinode (/glade/scratch/$user/...) filter will be copied '
   echo '       to output_dir and to a local dir = $output_dir:t.'
   echo '       Otherwise, the only copy of filter will be in $output_dir'
   exit

else if ($#argv == 1) then
   set local_dir = $1:t
   if (-d $local_dir) then
      echo "Directory $local_dir already exists.  Exiting"
      exit
   else
      mkdir $local_dir
      set out_file = ${local_dir}/minimal_build.out
   endif
else
   set out_file = minimal_build.out
endif

echo "Build date:"              >! $out_file
date                            >> $out_file

# module swap mpt/2.21            >>& $out_file
module list                     >>& $out_file

echo "PATH = "                  >> $out_file
echo $PATH | sed -e "s#:#\n#g " >> $out_file

if (-f filter) then
   mv filter filter.$$
   echo "\n Moved existing filter to filter.$$\n"
endif

csh ./mkmf_filter -mpi          >>& $out_file
make                            >>& $out_file

if ( $cdebug ) then 
   echo 'preserving .o and .mod files for debugging' >> $out_file
   mv *.o *.mod Makefile .cppdefs $local_dir
else
   \rm -f *.o *.mod Makefile .cppdefs
endif

echo "md5sum filter:" >> $out_file
md5sum filter >> $out_file

# set echo verbose
echo '==============================================' >> $out_file
if ($#argv == 1) then
   mv filter $local_dir
   ln -s $local_dir/filter .
   ls -l $local_dir >> $out_file

   # Test for multi node directory name, and copy filter there.
   if ($1:h != $1:t) then
      if (! -d $1) mkdir $1
      cp ${local_dir}/filter $1  || exit 19
#       cp fill_inflation_restart $1  || exit 20

      cd $1:h
      if (-f filter) then
         rm filter 
      endif
#       if (-f $1:h/fill_inflation_restart) then
#          rm fill_inflation_restart 
#       endif
      ln -s ${local_dir}/filter .
   endif
else
   echo "Not copying executables anywhere" >> $out_file
endif

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

