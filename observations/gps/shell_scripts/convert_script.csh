#!/bin/csh
#
# DART software - Copyright © 2004 - 2010 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# Main script:
# generate multiple days of gps observations
#
# calls the cosmic_to_obsseq script with 4 args:
#
#  - the date in YYYYMMDD format
#  - the working directory location
#  - whether to download the data automatically from the cosmic web site 
#     (downloading data requires signing up for a username/password to 
#     access the site, and then setting the username and password here 
#     before running this script.)  set to 'no' if the data has already been
#     downloaded separately before now.  set to 'both' to download both the
#     current day plus the next day (needed to get the first N hours of the
#     following day for assimilation windows centered on midnight).  set to
#     'next' if the current day's data is already downloaded.
#  - whether to delete the data automatically from the local disk after the
#     conversion is done.  values are 'no', 'curr', or 'both'.  'curr' deletes
#     the current day but leaves the following day.  'both' deletes both
#     the current and next day's data.
#

setenv cosmic_user xxx
setenv cosmic_pw   yyy

# assumes all data predownloaded, and will be deleted afterwards
# by hand.
./cosmic_to_obsseq.csh 20070930 ../cosmic no no
./cosmic_to_obsseq.csh 20071001 ../cosmic no no
./cosmic_to_obsseq.csh 20071002 ../cosmic no no
./cosmic_to_obsseq.csh 20071003 ../cosmic no no
./cosmic_to_obsseq.csh 20071004 ../cosmic no no
./cosmic_to_obsseq.csh 20071005 ../cosmic no no
./cosmic_to_obsseq.csh 20071006 ../cosmic no no
./cosmic_to_obsseq.csh 20071007 ../cosmic no no

## example of using both, curr and next for on-demand download and cleanup.
#./cosmic_to_obsseq.csh 20071001 ../cosmic both curr
#./cosmic_to_obsseq.csh 20071002 ../cosmic next curr
#./cosmic_to_obsseq.csh 20071003 ../cosmic next curr
#./cosmic_to_obsseq.csh 20071004 ../cosmic next both

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

