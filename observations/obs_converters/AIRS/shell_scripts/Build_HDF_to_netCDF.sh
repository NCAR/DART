#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# This file is intended to provide guidance on how to compile the 
# "HDF4 CF CONVERSION TOOLKIT" from the HDF-EOS Tools and Information Center
#
# The URL of the HDF-EOS Tools and Information Center is:
# http://hdfeos.org/software/h4cflib.php
# 
# The URL of the "HDF4 CF CONVERSION TOOLKIT" is:
# http://hdfeos.org/software/h4cflib/h4cflib_1.3.tar.gz
#
# This is not a substitute for the README and INSTALL contained in the tar file.
#
# My habit is to install software for my personal use in my $HOME/local directory.

export MYHOME=/glade/work/thoar

./configure \
            --prefix=${MYHOME}/local/h4cf_1.3 \
         --with-hdf4=${MYHOME}/local/hdf-eos \
         --with-jpeg=${MYHOME}/local/hdf-eos \
         --with-zlib=${MYHOME}/local/hdf-eos \
      --with-hdfeos2=${MYHOME}/local/hdf-eos \
       --with-netcdf=$NETCDF \
        --with-szlib=/usr/local/szip \
            CPPFLAGS=-I${MYHOME}/local/hdf-eos/include  \
             LDFLAGS=-L${MYHOME}/local/hdf-eos/lib || exit 1

make          || exit 2
make check    || exit 3
make install  || exit 4

exit

#  The best way to get the most current configure options is to use configure:
#  ./configure --help
#
#  Here is a recap of just the environment variables, there are many more options
#
#  Some influential environment variables:
#    CC          C compiler command
#    CFLAGS      C compiler flags
#    LDFLAGS     linker flags, e.g. -L<lib dir> if you have libraries in a
#                nonstandard directory <lib dir>
#    LIBS        libraries to pass to the linker, e.g. -l<library>
#    CPPFLAGS    (Objective) C/C++ preprocessor flags, e.g. -I<include dir> if
#                you have headers in a nonstandard directory <include dir>
#    LT_SYS_LIBRARY_PATH
#                User-defined run-time library search path.
#    CPP         C preprocessor
#    CXX         C++ compiler command
#    CXXFLAGS    C++ compiler flags
#    CXXCPP      C++ preprocessor
