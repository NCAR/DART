#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# example script for converting only a few days of GPS Radiosonde Occultation
# observations to DART observation sequence format.
#
# for longer time periods see the convert_many_gpsro.csh script.
#
#
# this script calls the gpsro_to_obsseq script with 6 args:
#
#  - the date in YYYYMMDD format
#  - the processing directory location, relative to the 'work' dir.
#  - whether to download the data automatically from the cosmic web site.
#     'yes' will do the download, 'no' assumes the data is already downloaded
#     and on local disk.  (downloading data requires signing up for a 
#     username/password to access the site, and then setting the username 
#     and password in the gpsro_to_obsseq script before running it.)
#  - whether to convert the data.  set to 'yes' to make obs_seq files (the
#     usual use of this script). 'no' if just downloading or just cleaning up.
#  - whether to delete the data automatically from the local disk after the
#     conversion is done.  valid values are 'yes' or 'no'.
#  - the name of the text file containing the satellite names to convert
#    data from. 
#

# examples of common use follow.  

# download only:
./gpsro_to_obsseq.csh 20071001 ../cosmic yes no no satlist
./gpsro_to_obsseq.csh 20071002 ../cosmic yes no no satlist
./gpsro_to_obsseq.csh 20071003 ../cosmic yes no no satlist

# convert only.  assume all data already downloaded:
./gpsro_to_obsseq.csh 20071001 ../cosmic no yes no satlist
./gpsro_to_obsseq.csh 20071002 ../cosmic no yes no satlist
./gpsro_to_obsseq.csh 20071003 ../cosmic no yes no satlist

# download and convert, not removing files:
./gpsro_to_obsseq.csh 20071001 ../cosmic yes yes no satlist
./gpsro_to_obsseq.csh 20071002 ../cosmic yes yes no satlist
./gpsro_to_obsseq.csh 20071003 ../cosmic yes yes no satlist

# clean up only after verifying conversion worked:
./gpsro_to_obsseq.csh 20071001 ../cosmic no no yes satlist
./gpsro_to_obsseq.csh 20071002 ../cosmic no no yes satlist
./gpsro_to_obsseq.csh 20071003 ../cosmic no no yes satlist

# download, convert, and clean up all in one go:
./gpsro_to_obsseq.csh 20071001 ../cosmic yes yes yes satlist
./gpsro_to_obsseq.csh 20071002 ../cosmic yes yes yes satlist
./gpsro_to_obsseq.csh 20071003 ../cosmic yes yes yes satlist

exit 0


