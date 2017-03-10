#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Main script:
#    generate multiple days of gps observations
#
# calls the cosmic_to_obsseq script with 5 args:
#
#  - the date in YYYYMMDD format
#  - the processing directory location, relative to the 'work' dir.
#  - whether to download the data automatically from the cosmic web site.
#     'yes' will do the download, 'no' assumes the data is already downloaded
#     and on local disk.  (downloading data requires signing up for a 
#     username/password to access the site, and then setting the username 
#     and password in the cosmic_to_obsseq script before running it.)
#  - whether to convert the data.  set to 'yes' to make obs_seq files (the
#     usual use of this script). 'no' if just downloading or just cleaning up.
#  - whether to delete the data automatically from the local disk after the
#     conversion is done.  valid values are 'yes' or 'no'.
#

# examples of common use follow.  

# download only:
./cosmic_to_obsseq.csh 20071001 ../cosmic yes no no
./cosmic_to_obsseq.csh 20071002 ../cosmic yes no no
./cosmic_to_obsseq.csh 20071003 ../cosmic yes no no
./cosmic_to_obsseq.csh 20071004 ../cosmic yes no no
./cosmic_to_obsseq.csh 20071005 ../cosmic yes no no

# convert only.  assume all data already downloaded:
./cosmic_to_obsseq.csh 20071001 ../cosmic no yes no
./cosmic_to_obsseq.csh 20071002 ../cosmic no yes no
./cosmic_to_obsseq.csh 20071003 ../cosmic no yes no
./cosmic_to_obsseq.csh 20071004 ../cosmic no yes no
./cosmic_to_obsseq.csh 20071005 ../cosmic no yes no

# download and convert, not removing files:
./cosmic_to_obsseq.csh 20071001 ../cosmic yes yes no
./cosmic_to_obsseq.csh 20071002 ../cosmic yes yes no
./cosmic_to_obsseq.csh 20071003 ../cosmic yes yes no
./cosmic_to_obsseq.csh 20071004 ../cosmic yes yes no
./cosmic_to_obsseq.csh 20071005 ../cosmic yes yes no

# clean up only after verifying conversion worked:
./cosmic_to_obsseq.csh 20071001 ../cosmic no no yes
./cosmic_to_obsseq.csh 20071002 ../cosmic no no yes
./cosmic_to_obsseq.csh 20071003 ../cosmic no no yes
./cosmic_to_obsseq.csh 20071004 ../cosmic no no yes
./cosmic_to_obsseq.csh 20071005 ../cosmic no no yes

# download, convert, and clean up all in one go:
./cosmic_to_obsseq.csh 20071001 ../cosmic yes yes yes
./cosmic_to_obsseq.csh 20071002 ../cosmic yes yes yes
./cosmic_to_obsseq.csh 20071003 ../cosmic yes yes yes
./cosmic_to_obsseq.csh 20071004 ../cosmic yes yes yes
./cosmic_to_obsseq.csh 20071005 ../cosmic yes yes yes

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

