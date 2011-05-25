#!/bin/bash
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Adds a file name to all entries with path_names at the beginning -
# use this for adding new files into the DART mkmf framework.

for pathfile in `ls -1 path_names*`
do
  echo "$1" >> $pathfile
done

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

