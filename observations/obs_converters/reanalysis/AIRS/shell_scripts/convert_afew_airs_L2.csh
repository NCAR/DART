#!/bin/csh 
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id: convert_afew_gpsro.csh 12571 2018-05-03 15:20:58Z nancy@ucar.edu $
#
# example script for converting only a few days of AIRS temp/moisture retrieval
# observations to DART observation sequence format.
#
# for longer time periods see the convert_many_airs_L2.csh script.
#
#
# this script calls the convert_airs_L2 script with 5 args:
#
#  - the date in YYYYMMDD format
#  - the processing directory location, relative to the 'work' dir.
#  - whether to download the data automatically from the nasa web site.
#     'yes' will do the download, 'no' assumes the data is already downloaded
#     and on local disk.  (downloading data requires signing up for a 
#     username/password to access the site, and then setting the username 
#     and password in the convert_airs_L2 script before running it.)
#  - whether to convert the data.  set to 'yes' to make obs_seq files (the
#     usual use of this script). 'no' if just downloading or just cleaning up.
#  - whether to delete the data automatically from the local disk after the
#     conversion is done.  valid values are 'yes' or 'no'.
#

# examples of common use follow.  

## download only:
#./convert_airs_L2.csh 20100701 ../output_240 yes no no
#./convert_airs_L2.csh 20100702 ../output_240 yes no no
#./convert_airs_L2.csh 20100703 ../output_240 yes no no

## convert only.  assume all data already downloaded:
./convert_airs_L2.csh 20100701 ../output_240 no yes no
./convert_airs_L2.csh 20100702 ../output_240 no yes no
./convert_airs_L2.csh 20100703 ../output_240 no yes no
./convert_airs_L2.csh 20100704 ../output_240 no yes no
./convert_airs_L2.csh 20100705 ../output_240 no yes no

## download and convert, not removing files:
#./convert_airs_L2.csh 20071001 ../output_240 yes yes no
#./convert_airs_L2.csh 20071002 ../output_240 yes yes no
#./convert_airs_L2.csh 20071003 ../output_240 yes yes no

## clean up only after verifying conversion worked:
#./convert_airs_L2.csh 20071001 ../output_240 no no yes
#./convert_airs_L2.csh 20071002 ../output_240 no no yes
#./convert_airs_L2.csh 20071003 ../output_240 no no yes

## download, convert, and clean up all in one go:
#./convert_airs_L2.csh 20071001 ../output_240 yes yes yes
#./convert_airs_L2.csh 20071002 ../output_240 yes yes yes
#./convert_airs_L2.csh 20071003 ../output_240 yes yes yes


exit 0

# <next few lines under version control, do not edit>
# $URL: https://svn-dares-dart.cgd.ucar.edu/DART/branches/rma_trunk/observations/obs_converters/gps/shell_scripts/convert_afew_gpsro.csh $
# $Revision: 12571 $
# $Date: 2018-05-03 09:20:58 -0600 (Thu, 03 May 2018) $

